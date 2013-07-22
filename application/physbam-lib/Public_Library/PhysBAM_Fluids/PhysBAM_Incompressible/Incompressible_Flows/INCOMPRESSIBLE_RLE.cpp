#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_SIMPLE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_TRANSFER_ITERATOR.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_LINEAR_PROFILE.h>
#include <PhysBAM_Tools/Grids_RLE_Advection/ADVECTION_FACE_VELOCITY_DONOR_CELL_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Boundaries/BOUNDARY_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/AVERAGING_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/FACE_LOOKUP_RLE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_RLE_GRID.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/EXTRAPOLATION_RLE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_RLE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_RLE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> INCOMPRESSIBLE_RLE<T_GRID>::
INCOMPRESSIBLE_RLE(const T_GRID& grid_input)
    :grid(grid_input),projection(grid,V),donor_cell_advection(*new ADVECTION_FACE_VELOCITY_DONOR_CELL_RLE<T>),mpi_grid(0),boundary_default(*new BOUNDARY_RLE<T_GRID,T>)
{
    Set_Custom_Boundary(boundary_default);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> INCOMPRESSIBLE_RLE<T_GRID>::
~INCOMPRESSIBLE_RLE()
{
    delete &donor_cell_advection;delete &boundary_default;
}
//#####################################################################
// Function Advance_One_Time_Step_Explicit_Part
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_RLE<T_GRID>::
Advance_One_Time_Step_Explicit_Part(const T dt,const T time,const int substeps,const bool implicit_viscosity,const ARRAY<T>* phi_ghost)
{
    assert(!projection.flame && !viscosity && !use_variable_viscosity && !implicit_viscosity);

    ARRAY<T> V_ghost(grid.number_of_faces,false);boundary->Fill_Ghost_Cells_Face(grid,V,V_ghost,time);
    ARRAY<T> V_save(V_ghost);

    // update convection
    LOG::Time("advecting velocity conservatively");
    donor_cell_advection.Euler_Step(grid,V,V_ghost,V_save,dt/substeps);
    for(int s=2;s<=substeps;s++){
        boundary->Apply_Boundary_Condition_Face(grid,V,time+s*dt/substeps);
        boundary->Fill_Ghost_Cells_Face(grid,V,V_ghost,time);
        donor_cell_advection.Euler_Step(grid,V,V_ghost,V_save,dt/substeps);}
    LOG::Time("advecting velocity nonconservatively");
    advection->Update_Advection_Equation_Face(grid,V,V_save,V_save,*boundary,dt,time);

    // update gravity and external forces
    LOG::Time("adding external forces");
    TV gravity_impulse=dt*gravity*downward_direction;
    T_GRID::template Face_Loop<Add_Impulse>(grid,V,gravity_impulse);
    if(use_force) for(int f=1;f<=grid.number_of_faces;f++)V(f)+=dt*force(f);

    boundary->Apply_Boundary_Condition_Face(grid,V,time+dt);
    grid.Put_Ghost_Faces((T)0,V);
}
//#####################################################################
// Function Add_Impulse
//#####################################################################
template<class T_GRID> template<class T_FACE> void INCOMPRESSIBLE_RLE<T_GRID>::
Add_Impulse::Apply(const T_GRID& grid,ARRAY<T>& V,const TV& impulse)
{
    assert(grid.long_run_faces_horizontal==2);
    for(T_FACE face(grid,0);face;face++){int f=face.Face();
        V(f)+=impulse[T_FACE::Axis()];
        if(face.Long()) V(f+1)+=impulse[T_FACE::Axis()];}
}
//#####################################################################
// Function Advance_One_Time_Step_Implicit_Part
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_RLE<T_GRID>::
Advance_One_Time_Step_Implicit_Part(const T dt,const T time,const bool implicit_viscosity,const ARRAY<T>* phi_ghost)
{
    assert(!implicit_viscosity);
    projection.Make_Divergence_Free(time,phi_ghost);
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR INCOMPRESSIBLE_RLE<T_GRID>::
CFL(const bool inviscid,const bool viscous_only) const
{
    assert((inviscid || !viscosity) && !viscous_only);
    ARRAY<TV> max_V_cell(grid.number_of_cells);
    T_GRID::template Face_Loop<Compute_Maximum_Cell_Velocity>(grid,V,max_V_cell);
    T max_V_norm=0;
    for(int c=1;c<=grid.number_of_cells;c++)max_V_norm=max(max_V_norm,max_V_cell(c).L1_Norm());
    T dt_convection=max_V_norm/grid.uniform_grid.dX.x;
    T dt_force=abs(gravity)*TV::Dot_Product(abs(downward_direction),grid.uniform_grid.one_over_dX);
    if(use_force){
        T max_force_norm=0;
        for(CELL_ITERATOR cell(grid,0);cell;cell++){
            T local_force_norm=0;
            for(int axis=1;axis<=T_GRID::dimension;axis++)local_force_norm+=maxabs(force(cell.First_Face_Index(axis)),force(cell.Second_Face_Index(axis)));
            max_force_norm=max(max_force_norm,local_force_norm);}
        dt_force+=max_force_norm/grid.uniform_grid.dX.x;}
    T dt_overall=(dt_convection+sqrt(sqr(dt_convection)+4*dt_force))/2;
    return 1/max(dt_overall,1/max_time_step);
}
//#####################################################################
// Function Compute_Maximum_Cell_Velocity
//#####################################################################
template<class T_GRID> template<class T_FACE> void INCOMPRESSIBLE_RLE<T_GRID>::
Compute_Maximum_Cell_Velocity::Apply(const T_GRID& grid,const ARRAY<T>& V,ARRAY<TV>& max_V_cell)
{
    int axis=T_FACE::Axis();
    for(T_FACE face(grid,0);face;face++){int f=face.Face(),c1=face.Cell(0),c2=face.Cell(1);
        T max_V=abs(V(f));if(face.Long()) max_V=max(max_V,abs(V(f+1)));
        max_V_cell(c1)[axis]=max(max_V_cell(c1)[axis],max_V);
        max_V_cell(c2)[axis]=max(max_V_cell(c2)[axis],max_V);}
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_RLE<T_GRID>::
Set_Dirichlet_Boundary_Conditions(const ARRAY<T>& phi,const T pressure)
{
    ARRAY<bool>& psi_D=projection.laplace.psi_D;
    for(CELL_ITERATOR cell(grid,0);cell;cell++){int c=cell.Cell();if(phi(c)>0){
        psi_D(c)=true;projection.p(c)=pressure;if(cell.Long()){psi_D(c+1)=true;projection.p(c+1)=pressure;}}}
    if(mpi_grid){
        mpi_grid->Exchange_Boundary_Cell_Data(psi_D,1,false);
        mpi_grid->Exchange_Boundary_Cell_Data(projection.p,1,false);}
}
//#####################################################################
// Function Extrapolate_Velocity_Across_Interface
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_RLE<T_GRID>::
Extrapolate_Velocity_Across_Interface(ARRAY<T>& phi_ghost,const bool enforce_divergence_free,const T band_width,const T_FACE_NEIGHBORS* face_neighbors_visible)
{
    boundary->Fill_Ghost_Cells_Face(grid,V,V,0); // TODO: use real time
    EXTRAPOLATION_RLE<T_GRID,T> extrapolate(grid,phi_ghost);extrapolate.Set_Band_Width(band_width+2);extrapolate.Set_Default_Value((T)0);
    // extrapolate velocity
    ARRAY<bool> fixed(projection.laplace.psi_N);T_GRID::template Face_Loop<Find_Fixed_Faces>(grid,phi_ghost,fixed);
    if(face_neighbors_visible) extrapolate.Set_Collision_Aware_Extrapolation(*face_neighbors_visible);
    extrapolate.Extrapolate_Faces(V,fixed);
    // make divergence free
    if(enforce_divergence_free) PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Find_Fixed_Faces
//#####################################################################
template<class T_GRID> template<class T_FACE> void INCOMPRESSIBLE_RLE<T_GRID>::
Find_Fixed_Faces::Apply(const T_GRID& grid,const ARRAY<T>& phi_ghost,ARRAY<bool>& fixed)
{
    for(T_FACE face(grid,grid.number_of_ghost_cells,false);face;face++)if(phi_ghost(face.Cell(0))<=0 || phi_ghost(face.Cell(1))<=0){
        int f=face.Face();fixed(f)=true;if(face.Long()) fixed(f+1)=true;}
}
//#####################################################################
// Function Transfer_Velocity
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_RLE<T_GRID>::
Transfer_Velocity(const T_GRID& new_grid)
{
    ARRAY<T> new_V(new_grid.number_of_faces);
    ARRAY<bool> new_valid_mask(new_grid.number_of_faces,false);ARRAYS_COMPUTATIONS::Fill(new_valid_mask,true);
    T_GRID::template Horizontal_Face_Loop<Transfer_Horizontal_Velocity>(grid,new_grid,V,new_V,valid_mask,new_valid_mask);
    Transfer_Vertical_Velocity(grid,new_grid,V,new_V,valid_mask,new_valid_mask);
    ARRAY<T>::Exchange_Arrays(V,new_V);ARRAY<bool>::Exchange_Arrays(valid_mask,new_valid_mask);
}
//#####################################################################
// Function Transfer_Horizontal_Velocity
//#####################################################################
template<class T_GRID> template<class T_FACE> void INCOMPRESSIBLE_RLE<T_GRID>::
Transfer_Horizontal_Velocity::Apply(const T_GRID& old_grid,const T_GRID& new_grid,const ARRAY<T>& V,ARRAY<T>& new_V,const ARRAY<bool>& valid_mask,ARRAY<bool>& new_valid_mask)
{
    if(old_grid.long_run_faces_horizontal==1)
        for(RLE_GRID_TRANSFER_ITERATOR<T_GRID,T_FACE> faces(old_grid,new_grid,old_grid.number_of_ghost_cells);faces;faces++){
            Transfer_Constant_Horizontal_Faces(faces,V,new_V);
            if(faces.source.Short() && faces.destination.Short()) new_valid_mask(faces.destination.Face())=valid_mask(faces.source.Face());}
    else
        for(RLE_GRID_TRANSFER_ITERATOR<T_GRID,T_FACE> faces(old_grid,new_grid,old_grid.number_of_ghost_cells);faces;faces++){
            Transfer_Linear_Horizontal_Faces(faces,V,new_V);
            if(faces.source.Short() && faces.destination.Short()) new_valid_mask(faces.destination.Face())=valid_mask(faces.source.Face());}
}
//#####################################################################
// Function Transfer_Vertical_Velocity
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_RLE<T_GRID>::
Transfer_Vertical_Velocity(const T_GRID& old_grid,const T_GRID& new_grid,const ARRAY<T>& V,ARRAY<T>& new_V,const ARRAY<bool>& valid_mask,ARRAY<bool>& new_valid_mask)
{
    for(RLE_GRID_TRANSFER_ITERATOR<T_GRID,FACE_Y_ITERATOR> faces(old_grid,new_grid,old_grid.number_of_ghost_cells);faces;faces++){
        Transfer_Vertical_Faces(faces,V,new_V);
        new_valid_mask(faces.destination.Face())=valid_mask(faces.source.Face());}
}
//#####################################################################
// Function Transfer_Velocity_Ghost
//#####################################################################
namespace{
struct Transfer_Velocity_Ghost_Helper{template<class T_FACE,class T_GRID,class T> static void Apply(const T_GRID& old_grid,const T_GRID& new_grid,const ARRAY<T>& V_ghost,ARRAY<T>& new_V_ghost)
{
    for(RLE_GRID_TRANSFER_ITERATOR<T_GRID,T_FACE> faces(old_grid,new_grid,old_grid.number_of_ghost_cells);faces;faces++)
        Transfer(faces,V_ghost,new_V_ghost);
}};}
template<class T_GRID> void INCOMPRESSIBLE_RLE<T_GRID>::
Transfer_Velocity_Ghost(const T_GRID& new_grid,ARRAY<T>& V_ghost) const
{
    ARRAY<T> new_V_ghost(new_grid.number_of_faces);
    T_GRID::template Face_Loop<Transfer_Velocity_Ghost_Helper>(grid,new_grid,V_ghost,new_V_ghost);
    ARRAY<T>::Exchange_Arrays(V_ghost,new_V_ghost);
}
//#####################################################################
template class INCOMPRESSIBLE_RLE<RLE_GRID_2D<float> >;
template class INCOMPRESSIBLE_RLE<RLE_GRID_3D<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INCOMPRESSIBLE_RLE<RLE_GRID_2D<double> >;
template class INCOMPRESSIBLE_RLE<RLE_GRID_3D<double> >;
#endif
#endif

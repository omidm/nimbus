//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM.h>
#include <PhysBAM_Dynamics/Parallel_Computation/REFINEMENT_THREADS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<T_GRID>::
PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM(const T_GRID& mac_grid,T_LEVELSET& levelset_input,const int scale,const T alpha_in,const bool use_surface_solve,const bool flame_input,const bool multiphase,const bool use_variable_beta,const bool use_poisson)
    :PROJECTION_REFINEMENT_UNIFORM<T_GRID>(mac_grid,scale,alpha_in,flame_input,multiphase,use_variable_beta,use_poisson),boundary(0),phi_boundary(0),levelset_projection(fine_grid),levelset(levelset_input),coarse_levelset(coarse_grid,coarse_phi),surface_solve(use_surface_solve),buffer(1)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<T_GRID>::
~PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM()
{
}
//#####################################################################
// Function Initialize_Grid
//#####################################################################
template<class T_GRID> void PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<T_GRID>::
Initialize_Grid(const T_GRID& mac_grid)
{
    BASE::Initialize_Grid(mac_grid);
    local_phi.Resize(local_grid.Domain_Indices(1));
    coarse_phi.Resize(coarse_grid.Domain_Indices(1));
    fast_local_projection.collidable_solver->Use_External_Level_Set(*new T_LEVELSET(local_grid,local_phi));    
    levelset_projection.elliptic_solver->Set_Relative_Tolerance((T)1e-7);
    levelset_projection.elliptic_solver->pcg.Set_Maximum_Iterations(400);
    levelset_projection.elliptic_solver->pcg.evolution_solver_type=krylov_solver_cg;
    levelset_projection.elliptic_solver->pcg.cg_restart_iterations=40;
    levelset_projection.elliptic_solver->pcg.Show_Results();
    levelset_projection.elliptic_solver->Solve_Neumann_Regions(false);
    levelset_projection.Initialize_Grid(fine_grid);
    levelset_projection.collidable_solver->Use_External_Level_Set(levelset);
    phi_ghost.Resize(mac_grid.Domain_Indices(max(1+coarse_scale,3)));
}
//#####################################################################
// Function Set_Coarse_Boundary_Conditions
//#####################################################################
template<class T_GRID> void PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<T_GRID>::
Set_Coarse_Boundary_Conditions(T_FACE_ARRAYS_SCALAR& coarse_face_velocities)
{
    int ghost=max(3,coarse_scale);
    phi_boundary->Fill_Ghost_Cells(fine_grid,levelset.phi,phi_ghost,0,0,ghost);
    for(int axis=1;axis<=T_GRID::dimension;axis++) for(int axis_side=1;axis_side<=2;axis_side++){
        int side=2*(axis-1)+axis_side;
        TV_INT interior_cell_offset=axis_side==1?TV_INT():-TV_INT::Axis_Vector(axis);
        TV_INT exterior_cell_offset=axis_side==1?-TV_INT::Axis_Vector(axis):TV_INT();
        TV_INT boundary_face_offset=axis_side==1?TV_INT::Axis_Vector(axis):-TV_INT::Axis_Vector(axis);
        if(solid_wall(axis)(axis_side)){
            for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_grid,1,T_GRID::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){
                TV_INT face=iterator.Face_Index()+boundary_face_offset;
                if(!Contains_Outside(face+interior_cell_offset,phi_ghost,0)){
                    if(coarse_face_velocities.Component(axis).Valid_Index(face)){
                        elliptic_solver->psi_N.Component(axis)(face)=true;coarse_face_velocities.Component(axis)(face)=0;}}
                else{TV_INT cell=face+exterior_cell_offset;
                    elliptic_solver->psi_D(cell)=true;p(cell)=0;}}}
        else
            for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_grid,1,T_GRID::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Face_Index()+interior_cell_offset;
                elliptic_solver->psi_D(cell)=true;p(cell)=0;}}
    Set_Beta_Face_For_Boundary_Conditions(coarse_face_velocities);
    if(false) if(POISSON_COLLIDABLE_UNIFORM<T_GRID>* poisson_collidable=dynamic_cast<POISSON_COLLIDABLE_UNIFORM<T_GRID>*>(poisson)){
        T_FACE_ARRAYS_SCALAR beta_face_new(coarse_grid);
        T_ARRAYS_SCALAR phi_ghost(coarse_grid.Domain_Indices(3),false);levelset.boundary->Fill_Ghost_Cells(coarse_grid,coarse_phi,phi_ghost,0,0,ghost);
        poisson_collidable->Find_Constant_beta(beta_face_new,phi_ghost);
        for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> index=iterator.Full_Index();poisson->beta_face(index)*=beta_face_new(index);}}
    for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_grid);iterator.Valid();iterator.Next()) if(Contains_Outside(iterator.Cell_Index(),levelset.phi,0)){
        elliptic_solver->psi_D(iterator.Cell_Index())=true;p(iterator.Cell_Index())=0;}
    if(coarse_mpi_grid){
        coarse_mpi_grid->Exchange_Boundary_Cell_Data(elliptic_solver->psi_D,1,false);
        coarse_mpi_grid->Exchange_Boundary_Cell_Data(p,1,false);}
}
//#####################################################################
// Function Set_Local_Boundary_Conditions
//#####################################################################
template<class T_GRID> bool PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<T_GRID>::
Set_Local_Boundary_Conditions(GRID<TV>& local_grid,PROJECTION_UNIFORM<GRID<TV> >& local_projection,TV_INT coarse_index)
{
    local_projection.elliptic_solver->psi_D.Fill(false);bool contains_surface=false;
    for(typename GRID<TV>::CELL_ITERATOR iterator(local_grid);iterator.Valid();iterator.Next()){
        TV_INT fine_cell=coarse_index*coarse_scale+iterator.Cell_Index()-coarse_scale*TV_INT::All_Ones_Vector();
        if(levelset.phi(fine_cell)>0){contains_surface=true;local_projection.elliptic_solver->psi_D(iterator.Cell_Index())=true;local_projection.p(iterator.Cell_Index())=0;}}
    return contains_surface;
}
//#####################################################################
// Function Set_Local_Phi_From_Fine_Phi
//#####################################################################
template<class T_GRID> void PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<T_GRID>::
Set_Local_Phi_From_Fine_Phi(GRID<TV>& local_mac_grid,ARRAY<T,TV_INT>& local_phi,const ARRAY<T,TV_INT>& fine_phi,TV_INT cell_index)
{
    for(typename GRID<TV>::CELL_ITERATOR iterator(local_mac_grid,1);iterator.Valid();iterator.Next()){
        TV_INT fine_index=(cell_index-TV_INT::All_Ones_Vector())*coarse_scale+iterator.Cell_Index();
        local_phi(iterator.Cell_Index())=fine_phi(fine_index);}
}
//#####################################################################
// Function Local_Projection_PCG
//#####################################################################
template<class T_GRID> void PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<T_GRID>::
Local_Projection_PCG(T_FACE_ARRAYS_SCALAR& fine_face_velocities,T_GRID& local_grid,T_FACE_ARRAYS_SCALAR& local_face_velocities,FAST_PROJECTION_DYNAMICS_UNIFORM<GRID<TV> >& local_projection,const T dt,const T time,TV_INT cell_index)
{
    if(surface_solve && Contains_Outside(cell_index,phi_ghost,buffer)) return;
    if(!surface_solve && !Contains_Inside(cell_index,phi_ghost,buffer)) return;
    Map_Fine_To_Local_Boundary_For_Cell(local_grid,local_face_velocities,fine_face_velocities,cell_index);
    Map_Fine_To_Local_Interior_For_Cell(local_grid,local_face_velocities,fine_face_velocities,cell_index,false);
    for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_grid);local_iterator.Valid();local_iterator.Next()){
        local_projection.p(local_iterator.Cell_Index())=0;}
    bool contains_solids=Map_Fine_To_Local_Boundaries_For_Cell(local_grid,local_projection.elliptic_solver->psi_N,cell_index);
    bool contains_surface=false;
    if(!surface_solve) contains_surface=Set_Local_Boundary_Conditions(local_grid,local_projection,cell_index);
    local_projection.elliptic_solver->Set_Neumann_Outer_Boundaries();
    local_projection.p*=dt;
    if(!surface_solve){
        Set_Local_Phi_From_Fine_Phi(local_grid,local_phi,levelset.phi,cell_index);
        local_projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method();}
    if(contains_solids || contains_surface) local_projection.Make_Divergence_Free(local_face_velocities,dt,time);
    else local_projection.Make_Divergence_Free_Fast(local_face_velocities,dt,time);
    local_projection.p/=dt;
    Map_Local_To_Fine_Interior_For_Cell(local_grid,local_face_velocities,fine_face_velocities,cell_index);
    if(!surface_solve) local_projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method(false);
}
//#####################################################################
// Function Contains_Outside
//#####################################################################
template<class T_GRID> bool PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<T_GRID>::
Contains_Outside(TV_INT cell_index,const ARRAY<T,TV_INT>& levelset_phi,int buffer)
{ 
    RANGE<TV_INT> domain=local_grid.Domain_Indices();domain.max_corner+=buffer*TV_INT::All_Ones_Vector();domain.min_corner-=buffer*TV_INT::All_Ones_Vector();
    for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_grid,domain);local_iterator.Valid();local_iterator.Next()){
        TV_INT fine_face=(cell_index-TV_INT::All_Ones_Vector())*coarse_scale+local_iterator.Cell_Index(); //Assume fine_grid is levelset_grid for now
        if(levelset_phi(fine_face)>0) return true;}
    return false;
}
//#####################################################################
// Function Contains_Inside
//#####################################################################
template<class T_GRID> bool PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<T_GRID>::
Contains_Inside(TV_INT cell_index,const ARRAY<T,TV_INT>& levelset_phi,int buffer)
{ 
    RANGE<TV_INT> domain=local_grid.Domain_Indices();domain.max_corner+=buffer*TV_INT::All_Ones_Vector();domain.min_corner-=buffer*TV_INT::All_Ones_Vector();
    for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_grid,domain);local_iterator.Valid();local_iterator.Next()){
        TV_INT fine_face=(cell_index-TV_INT::All_Ones_Vector())*coarse_scale+local_iterator.Cell_Index(); //Assume fine_grid is levelset_grid for now
        if(levelset_phi(fine_face)<=0) return true;}
    return false;
}
//#####################################################################
// Function Set_Coarse_Phi_From_Fine_Phi
//#####################################################################
template<class T_GRID> void PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<T_GRID>::
Set_Coarse_Phi_From_Fine_Phi(ARRAY<T,TV_INT>& coarse_phi,const ARRAY<T,TV_INT>& fine_phi)
{
    for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_grid,1);iterator.Valid();iterator.Next())
        coarse_phi(iterator.Cell_Index())=phi_interpolation.Clamped_To_Array(fine_grid,fine_phi,iterator.Location());
}
//#####################################################################
// Function Set_Levelset_Boundary_Conditions
//#####################################################################
template<class T_GRID> void PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<T_GRID>::
Set_Levelset_Boundary_Conditions(const GRID<TV>& levelset_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& levelset_velocities,const ARRAY<T,TV_INT>& levelset_phi,const T time)
{
    Map_Fine_To_Levelset_For_Constraints(levelset_velocities);
    for(int axis=1;axis<=TV::dimension;axis++) for(int axis_side=1;axis_side<=2;axis_side++){int side=2*(axis-1)+axis_side;
        TV_INT interior_cell_offset=axis_side==1?TV_INT():-TV_INT::Axis_Vector(axis);
        TV_INT exterior_cell_offset=axis_side==1?-TV_INT::Axis_Vector(axis):TV_INT();
        TV_INT boundary_face_offset=axis_side==1?TV_INT::Axis_Vector(axis):-TV_INT::Axis_Vector(axis);
        if(solid_wall(axis)(axis_side)){
            for(typename GRID<TV>::FACE_ITERATOR iterator(levelset_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){
                TV_INT face=iterator.Face_Index()+boundary_face_offset;
                if(levelset_phi(face+interior_cell_offset)<=0){
                    if(levelset_velocities.Component(axis).Valid_Index(face)){
                        levelset_projection.elliptic_solver->psi_N.Component(axis)(face)=true;levelset_velocities.Component(axis)(face)=0;}}
                else{TV_INT cell=face+exterior_cell_offset;levelset_projection.elliptic_solver->psi_D(cell)=true;levelset_projection.p(cell)=0;}}}
        else for(typename GRID<TV>::FACE_ITERATOR iterator(levelset_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Face_Index()+interior_cell_offset;
            levelset_projection.elliptic_solver->psi_D(cell)=true;levelset_projection.p(cell)=0;}}
    for(typename GRID<TV>::CELL_ITERATOR iterator(levelset_grid);iterator.Valid();iterator.Next()){
        if(levelset_phi(iterator.Cell_Index())>0){
            levelset_projection.elliptic_solver->psi_D(iterator.Cell_Index())=true;levelset_projection.p(iterator.Cell_Index())=0;}}
    for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_grid);iterator.Valid();iterator.Next()){ 
        if(!Contains_Outside(iterator.Cell_Index(),levelset_phi,buffer)) for(int axis=1;axis<=TV::dimension;axis++) for(int axis_side=1;axis_side<=2;axis_side++){
            TV_INT offset=(axis_side==1?-TV_INT::Axis_Vector(axis):TV_INT::Axis_Vector(axis));
            TV_INT face=iterator.Cell_Index()+(axis_side==1?TV_INT():TV_INT::Axis_Vector(axis));
            if(domain_boundary(axis)(axis_side) && ((iterator.Cell_Index()(axis)+offset(axis))<1 || (iterator.Cell_Index()(axis)+offset(axis))>coarse_grid.Counts()(axis))) continue;
            TV_INT adjacent_cell=iterator.Cell_Index()+offset;
            bool adjacent_outside=false;
            if(Contains_Outside(adjacent_cell,levelset_phi,buffer)) adjacent_outside=true;
            RANGE<TV_INT> domain;domain.min_corner=TV_INT::All_Ones_Vector();domain.max_corner=TV_INT::All_Ones_Vector()*coarse_scale;domain.max_corner(axis)=1;
            if(adjacent_outside) for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_grid,domain);local_iterator.Valid();local_iterator.Next()){
                TV_INT fine_face=coarse_scale*(face-TV_INT::All_Ones_Vector())+local_iterator.Cell_Index();
                levelset_projection.elliptic_solver->psi_N(FACE_INDEX<TV::dimension>(axis,fine_face))=true;}}}
}
//#####################################################################
// Function Map_Fine_To_Coarse
//#####################################################################
template<class T_GRID> void PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<T_GRID>::
Map_Fine_To_Levelset_For_Constraints(T_FACE_ARRAYS_SCALAR& face_velocities)
{
    for(typename GRID<TV>::FACE_ITERATOR iterator(fine_grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        if(fine_psi_N(face)) face_velocities(face)=face_velocities_save(face);}
}
//#####################################################################
// Function Map_Fine_To_Coarse
//#####################################################################
template<class T_GRID> void PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<T_GRID>::
Map_Fine_To_Coarse(T_FACE_ARRAYS_SCALAR& coarse_face_velocities,const T_FACE_ARRAYS_SCALAR& face_velocities)
{
    phi_boundary->Fill_Ghost_Cells(fine_grid,levelset.phi,phi_ghost,0,0,3);    
    Set_Coarse_Phi_From_Fine_Phi(coarse_phi,phi_ghost);
    collidable_solver->Use_External_Level_Set(coarse_levelset);
    BASE::Map_Fine_To_Coarse(coarse_face_velocities,face_velocities);
}
//#####################################################################
// Function Map_Coarse_To_Fine
//#####################################################################
template<class T_GRID> void PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<T_GRID>::
Map_Coarse_To_Fine(const T_FACE_ARRAYS_SCALAR& coarse_face_velocities,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    int ghost=max(buffer+coarse_scale,3);
    phi_boundary->Fill_Ghost_Cells(fine_grid,levelset.phi,phi_ghost,dt,time,ghost);    
    BASE::Map_Coarse_To_Fine(coarse_face_velocities,face_velocities,dt,time);
    collidable_solver->Use_External_Level_Set(levelset);
    if(surface_solve){
        for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_grid);iterator.Valid();iterator.Next()){
            for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_grid);local_iterator.Valid();local_iterator.Next()){
                TV_INT fine_index=(iterator.Cell_Index()-TV_INT::All_Ones_Vector())*coarse_scale+local_iterator.Cell_Index();
                levelset_projection.p(fine_index)=0;}}
        levelset_projection.elliptic_solver->psi_D.Fill(false);levelset_projection.elliptic_solver->psi_N=fine_psi_N;
        Set_Levelset_Boundary_Conditions(fine_grid,face_velocities,phi_ghost,time);
        if(fine_mpi_grid){
            levelset_projection.elliptic_solver->mpi_grid=fine_mpi_grid;
            fine_mpi_grid->Union_Common_Face_Data(levelset_projection.elliptic_solver->psi_N);
            fine_mpi_grid->Exchange_Boundary_Cell_Data(levelset_projection.elliptic_solver->psi_D,1,false);
            fine_mpi_grid->Exchange_Boundary_Cell_Data(levelset_projection.p,1,false);}
        levelset_projection.p*=dt;
        levelset_projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method();
        boundary->Apply_Boundary_Condition_Face(levelset_projection.p_grid,face_velocities,time+dt);        
        levelset_projection.Make_Divergence_Free(face_velocities,dt,time);
        levelset_projection.p/=dt;
        levelset_projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method(false);}
}
//#####################################################################
template class PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<GRID<VECTOR<float,1> > >;
template class PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<GRID<VECTOR<float,2> > >;
template class PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<GRID<VECTOR<double,1> > >;
template class PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<GRID<VECTOR<double,2> > >;
template class PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<GRID<VECTOR<double,3> > >;
#endif

//#####################################################################
// Copyright 2009-2011, Mridul Aanjaneya, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES
//#####################################################################
#include <PhysBAM_Dynamics/Coupled_Evolution/COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/BOUNDARY_OBJECT_EULER.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES<TV>::
COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES(T_GRID& grid_input,T_ARRAYS_DIMENSION_SCALAR& U_input,T_ARRAYS_BOOL& psi_input,
                                          MPI_UNIFORM_GRID<T_GRID>* mpi_grid_input):
    grid(grid_input),U(U_input),psi(psi_input),mpi_grid(mpi_grid_input),collision_bodies_affecting_fluid(0),object_boundary(new BOUNDARY_OBJECT_EULER<T_GRID>()),
    use_fast_marching(false),use_higher_order_solid_extrapolation(true),number_of_cells_to_extrapolate(7)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES<TV>::
~COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES()
{}
//#####################################################################
// Initialize_Solid_Fluid_Coupling
//#####################################################################
template<class TV> void COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES<TV>::
Initialize_Solid_Fluid_Coupling(GRID_BASED_COLLISION_GEOMETRY<T_GRID>* collision_bodies_affecting_fluid_input)
{
    phi_all_solids_negated.Resize(grid.Domain_Indices(number_of_cells_to_extrapolate));
    outside_fluid.Resize(grid.Domain_Indices(number_of_cells_to_extrapolate));
    collision_bodies_affecting_fluid=dynamic_cast<GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>*>(collision_bodies_affecting_fluid_input);
    collision_bodies_affecting_fluid->Initialize_Grids();
    collision_bodies_affecting_fluid->Rasterize_Objects();
}
//#####################################################################
// Update_Cut_Out_Grid
//#####################################################################
template<class TV> void COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES<TV>::
Update_Cut_Out_Grid()
{   
    Compute_Phi_Solids(0);
    for(CELL_ITERATOR iterator(grid,3);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        psi(cell_index)=!outside_fluid(cell_index);}
}
//#####################################################################
// Function Fill_Solid_Cells
//#####################################################################
template<class TV> void COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES<TV>::
Fill_Solid_Cells(const T dt,const T time)
{
    if(collision_bodies_affecting_fluid){
        int number_of_ghost_cells=mpi_grid?3:0;
        Compute_Phi_Solids(number_of_ghost_cells);
        Extrapolate_State_Into_Solids(phi_all_solids_negated,dt,time,number_of_ghost_cells,number_of_cells_to_extrapolate);}
}
//#####################################################################
// Function Compute_Phi_Solids
//#####################################################################
template<class TV> void COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES<TV>::
Compute_Phi_Solids(const int number_of_ghost_cells)
{
    ARRAY<TV_INT> seed_indices;
    phi_all_solids_negated.Fill(-FLT_MAX);
    outside_fluid.Fill(false);
    T_FACE_ARRAYS_BOOL kinematic_faces(grid.Domain_Indices(1));kinematic_faces.Fill(false);

    for(COLLISION_GEOMETRY_ID id(1);id<=collision_bodies_affecting_fluid->collision_geometry_collection.bodies.m;id++)
        if(collision_bodies_affecting_fluid->collision_geometry_collection.Is_Active(id)){
            COLLISION_GEOMETRY<TV>& collision_body=*(collision_bodies_affecting_fluid->collision_geometry_collection.bodies(id));
            T collision_thickness_over_two=(T).5*collision_bodies_affecting_fluid->collision_thickness;
            RANGE<TV_INT> bounding_box(grid.Clamp_To_Cell(collision_body.Axis_Aligned_Bounding_Box().Thickened(grid.dX.Max()*(T)2),number_of_cells_to_extrapolate));
            for(CELL_ITERATOR iterator(grid,bounding_box);iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();T phi_value=collision_body.Implicit_Geometry_Extended_Value(iterator.Location());
                if(collision_body.Inside(iterator.Location(),collision_thickness_over_two)) outside_fluid(cell_index)=true;
                phi_all_solids_negated(cell_index)=max(-phi_value,phi_all_solids_negated(cell_index));}}

    if(use_fast_marching){
        collision_bodies_affecting_fluid->outside_fluid=&outside_fluid;
        UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV> iterator_info(*collision_bodies_affecting_fluid);
        iterator_info.Initialize_Collision_Aware_Face_Iterator(outside_fluid,kinematic_faces,7,false);

        for(UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<TV> iterator(iterator_info);iterator.Valid();iterator.Next()){
            seed_indices.Append(iterator.First_Cell_Index());seed_indices.Append(iterator.Second_Cell_Index());}

        typename LEVELSET_POLICY<GRID<TV> >::LEVELSET levelset(grid,phi_all_solids_negated); // TODO(jontg): Make this a permanent member variable?
        FAST_MARCHING_METHOD_UNIFORM<GRID<TV> > fmm(levelset,number_of_cells_to_extrapolate);
        fmm.Fast_Marching_Method(phi_all_solids_negated,grid.dX.Max()*(T)number_of_cells_to_extrapolate,&seed_indices,true);}
}
//#####################################################################
// Function Extrapolate_State_Into_Solids
//#####################################################################
template<class TV> void COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES<TV>::
Extrapolate_State_Into_Solids(T_ARRAYS_SCALAR& phi_all_solids_negated,const T dt,const T time,const int number_of_ghost_cells,const int number_of_cells_to_extrapolate)
{
    T_ARRAYS_DIMENSION_SCALAR U_extrapolated(grid.Domain_Indices(number_of_ghost_cells));

    if(mpi_grid) boundary->Fill_Ghost_Cells(grid,U,U_extrapolated,(T)dt,(T)time,number_of_ghost_cells);
    else U_extrapolated=U;

    T_EXTRAPOLATION_SCALAR_DIMENSION extrapolate(grid,phi_all_solids_negated,U_extrapolated,number_of_ghost_cells);extrapolate.Set_Band_Width((T)number_of_cells_to_extrapolate);
    extrapolate.Extrapolate((T)0,false);
    T_ARRAYS_DIMENSION_SCALAR::Get(U,U_extrapolated);

    T band_width=number_of_cells_to_extrapolate*grid.dX.Max();
    T_LINEAR_INTERPOLATION_DIMENSION interpolation;
    T max_distance,object_velocity_normal_component;TV location,normal_direction,object_velocity,reflected_point;
    const RANGE<TV>& domain=RANGE<TV>::Intersect(grid.Ghost_Domain(number_of_ghost_cells),mpi_grid?mpi_grid->global_grid.Domain():grid.Domain());

    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();TV location=iterator.Location();
        if(outside_fluid(cell_index) && phi_all_solids_negated(cell_index)<band_width){
            max_distance=phi_all_solids_negated(cell_index)*(T)2;
            Get_Neumann_Data(location,max_distance,normal_direction,object_velocity_normal_component,reflected_point);
            COLLISION_GEOMETRY_ID body_id;
            if(use_higher_order_solid_extrapolation && domain.Inside(reflected_point,grid.dX.Max()*(T).5) && 
                    !collision_bodies_affecting_fluid->Implicit_Geometry_Lazy_Inside_Any_Body(reflected_point,body_id)){
                U(cell_index)=interpolation.Clamped_To_Array(grid,U_extrapolated,reflected_point);}
            object_boundary->Apply_Neumann_Boundary_Condition(U(cell_index),normal_direction,object_velocity_normal_component);}}
}
//#####################################################################
// Function Get_Neumann_Data 
//#####################################################################
template<class TV> void COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES<TV>::
Get_Neumann_Data(const TV& location,const T max_distance,TV& normal_direction,T& object_velocity_normal_component,
                 TV& reflected_point) const
{
    T distance;TV boundary_point;COLLISION_GEOMETRY_ID body_id;int simplex_id;
    boundary_point=collision_bodies_affecting_fluid->collision_geometry_collection.Closest_Boundary_Point(
            location,max_distance,distance,body_id,simplex_id);
    const COLLISION_GEOMETRY<TV>& collision_body=collision_bodies_affecting_fluid->collision_geometry_collection(body_id);
    reflected_point=location+((T)2*(boundary_point-location));
    normal_direction=boundary_point-location;normal_direction.Normalize();
    TV object_velocity=collision_body.Pointwise_Object_Velocity(simplex_id,boundary_point);
    object_velocity_normal_component=TV::Dot_Product(object_velocity,normal_direction);
}
//#####################################################################
template class COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES<VECTOR<float,1> >;
template class COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES<VECTOR<float,2> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES<VECTOR<double,1> >;
template class COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES<VECTOR<double,2> >;
#endif

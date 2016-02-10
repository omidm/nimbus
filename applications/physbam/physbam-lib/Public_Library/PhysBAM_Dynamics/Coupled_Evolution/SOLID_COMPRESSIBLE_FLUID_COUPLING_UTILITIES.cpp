//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Geometry/Basic_Geometry/POLYGON.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Computations/BOX_BOX_INTERSECTION_AREA.h>
#include <PhysBAM_Geometry/Basic_Geometry_Computations/BOX_POLYGON_INTERSECTION_AREA.h>
#include <PhysBAM_Geometry/Collisions_And_Grids/OBJECTS_IN_CELL.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/CUT_CELL_COMPUTATIONS.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/RIGID_GEOMETRY_RASTERIZATION_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/EULER_THINSHELL_ADVECTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
#include <PhysBAM_Dynamics/Forces_And_Torques/EULER_FLUID_FORCES.h>
#include <iomanip>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES(EULER_UNIFORM<T_GRID>& euler_input,MPI_UNIFORM_GRID<T_GRID>* mpi_grid_input):
    euler(euler_input),mpi_grid(mpi_grid_input),collision_bodies_affecting_fluid(0),thinshell(false),use_fast_marching(false),
    use_higher_order_solid_extrapolation(true),fluid_affects_solid(false),number_of_cells_to_extrapolate(7),
    solid_state(TV_DIMENSION()),euler_fluid_forces(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
~SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES()
{
    delete euler_fluid_forces;
}
//#####################################################################
// Initialize_Solid_Fluid_Coupling
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Initialize_Solid_Fluid_Coupling(GRID_BASED_COLLISION_GEOMETRY<T_GRID>* collision_bodies_affecting_fluid_input)
{
    phi_all_solids_negated.Resize(euler.grid.Domain_Indices(number_of_cells_to_extrapolate));
    outside_fluid.Resize(euler.grid.Domain_Indices(number_of_cells_to_extrapolate));
    uncovered_cells.Resize(euler.grid.Domain_Indices());
    collision_bodies_affecting_fluid=
        dynamic_cast<GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>*>(collision_bodies_affecting_fluid_input);
    collision_bodies_affecting_fluid->Initialize_Grids();
    collision_bodies_affecting_fluid->Update_Intersection_Acceleration_Structures(false); // NON-swept acceleration structures
    collision_bodies_affecting_fluid->Rasterize_Objects();
    if(fluid_affects_solid){
        if(!euler.timesplit){
            euler.compute_pressure_fluxes=true;
            solid_fluid_face_time_n.Resize(euler.grid,0);
            cells_inside_fluid_time_n.Resize(euler.grid.Domain_Indices(0));
            pressure_at_faces.Resize(euler.grid.Domain_Indices(0));}
        if(euler.timesplit && thinshell){
            U_n.Resize(euler.grid.Domain_Indices(3));
            accumulated_flux.Resize(euler.grid.Domain_Indices());
            near_interface.Resize(euler.grid.Domain_Indices(1));near_interface.Fill(false);
            cut_cells_n.Resize(euler.grid.Domain_Indices(1));cut_cells_n.Fill(0);
            cut_cells_np1.Resize(euler.grid.Domain_Indices(1));cut_cells_np1.Fill(0);
            cell_volumes_n.Resize(euler.grid.Domain_Indices(2));
            cell_volumes_np1.Resize(euler.grid.Domain_Indices(2));
            psi_n.Resize(euler.grid.Domain_Indices(1));
            psi_np1.Resize(euler.grid.Domain_Indices(1));

            uncovered_cells_n_p_half.Resize(euler.grid.Domain_Indices());
            cut_cells_n_p_half.Resize(euler.grid.Domain_Indices(1));cut_cells_n_p_half.Fill(0);
            cell_volumes_n_p_half.Resize(euler.grid.Domain_Indices(2));
            psi_n_p_half.Resize(euler.grid.Domain_Indices(1));}}
}
//#####################################################################
// Update_Cut_Out_Grid
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Update_Cut_Out_Grid()
{   
    if(euler.timesplit && !thinshell){
        Compute_Phi_Solids(0);
        for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            //euler.psi(cell_index)=!outside_fluid(cell_index);}
            euler.psi(cell_index)=phi_all_solids_negated(cell_index)<(T)2*euler.grid.dX.Max();}}
    else if(euler.timesplit && thinshell)
        for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            euler.psi(cell_index)=psi_np1(cell_index);}
    else{
        Compute_Phi_Solids(0);
        for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            euler.psi(cell_index)=!outside_fluid(cell_index);}}
}   
//#####################################################################
// Function Get_Neumann_Data 
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Get_Neumann_Data(const TV& location,const T max_distance,TV& normal_direction,T& object_velocity_normal_component,
                 TV& reflected_point) const
{
    assert(!euler.timesplit || !thinshell);
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
// Function Get_Neumann_Data 
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Get_Neumann_Data(const TV& location,const T max_distance,TV& object_velocity,TV& reflected_point) const
{
    assert(!euler.timesplit || !thinshell);
    T distance;TV boundary_point;COLLISION_GEOMETRY_ID body_id;int simplex_id;
    boundary_point=collision_bodies_affecting_fluid->collision_geometry_collection.Closest_Boundary_Point(
            location,max_distance,distance,body_id,simplex_id);
    const COLLISION_GEOMETRY<TV>& collision_body=collision_bodies_affecting_fluid->collision_geometry_collection(body_id);
    reflected_point=location+((T)2*(boundary_point-location));
    object_velocity=collision_body.Pointwise_Object_Velocity(simplex_id,boundary_point);
}
//#####################################################################
// Function Fill_Uncovered_Cells 
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Fill_Uncovered_Cells()
{
    if(!euler.timesplit || !thinshell){
        T_LINEAR_INTERPOLATION_DIMENSION interpolation;
        T max_distance,object_velocity_normal_component;TV location,normal_direction,object_velocity,reflected_point;
        for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            if(uncovered_cells(cell_index)){location=euler.grid.Center(cell_index);
                max_distance=phi_all_solids_negated(cell_index)*(T)2;
                if(euler.use_solid_velocity_in_ghost_cells){Get_Neumann_Data(
                                location,max_distance,object_velocity,reflected_point);normal_direction=object_velocity;}
                else Get_Neumann_Data(location,max_distance,normal_direction,object_velocity_normal_component,reflected_point);
                euler.U(cell_index)=interpolation.Clamped_To_Array(euler.grid,euler.U,reflected_point);
                euler.conservation->object_boundary->Apply_Neumann_Boundary_Condition(euler.U(cell_index),normal_direction,
                                                                                      object_velocity_normal_component);}}
        euler.Invalidate_Ghost_Cells();}
}
//#####################################################################
// Function Extrapolate_State_Into_Solids
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Extrapolate_State_Into_Solids(T_ARRAYS_SCALAR& phi_all_solids_negated,const int number_of_ghost_cells,
                              const int number_of_cells_to_extrapolate)
{
    assert(!euler.timesplit || !thinshell);
    T_ARRAYS_DIMENSION_SCALAR U_extrapolated(euler.grid.Domain_Indices(number_of_ghost_cells));

    if(mpi_grid){
        BOUNDARY_UNIFORM<GRID<TV>,TV_DIMENSION> boundary; // Constant extrapolation at non-mpi faces.
        BOUNDARY_MPI<GRID<TV>,TV_DIMENSION> mpi_boundary(mpi_grid,boundary);
        mpi_boundary.Fill_Ghost_Cells(euler.grid,euler.U,U_extrapolated,(T)0,(T)0,number_of_ghost_cells);}
    else U_extrapolated=euler.U;

    T_EXTRAPOLATION_SCALAR_DIMENSION extrapolate(euler.grid,phi_all_solids_negated,U_extrapolated,number_of_ghost_cells);
    extrapolate.Set_Band_Width((T)number_of_cells_to_extrapolate);
    extrapolate.Extrapolate((T)0,false);
    T_ARRAYS_DIMENSION_SCALAR::Get(euler.U,U_extrapolated);

    T band_width=number_of_cells_to_extrapolate*euler.grid.dX.Max();
    T_LINEAR_INTERPOLATION_DIMENSION interpolation;
    T max_distance,object_velocity_normal_component;TV location,normal_direction,object_velocity,reflected_point;
    const RANGE<TV>& domain=RANGE<TV>::Intersect(
            euler.grid.Ghost_Domain(number_of_ghost_cells),mpi_grid?mpi_grid->global_grid.Domain():euler.grid.Domain());

    for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();TV location=iterator.Location();
        if(outside_fluid(cell_index) && phi_all_solids_negated(cell_index)<band_width){
            max_distance=phi_all_solids_negated(cell_index)*(T)2;
            if(euler.use_solid_velocity_in_ghost_cells){
                Get_Neumann_Data(location,max_distance,object_velocity,reflected_point);normal_direction=object_velocity;}
            else Get_Neumann_Data(location,max_distance,normal_direction,object_velocity_normal_component,reflected_point);
            COLLISION_GEOMETRY_ID body_id;
            if(use_higher_order_solid_extrapolation && domain.Inside(reflected_point,euler.grid.dX.Max()*(T).5) && 
               !collision_bodies_affecting_fluid->Implicit_Geometry_Lazy_Inside_Any_Body(reflected_point,body_id)){
                euler.U(cell_index)=interpolation.Clamped_To_Array(euler.grid,U_extrapolated,reflected_point);}
            euler.conservation->object_boundary->Apply_Neumann_Boundary_Condition(euler.U(cell_index),normal_direction,
                                                                                  object_velocity_normal_component);}}
    euler.Invalidate_Ghost_Cells();
}
//#####################################################################
// Function Compute_Phi_Solids
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Compute_Phi_Solids(const int number_of_ghost_cells)
{
    LOG::SCOPE scope("Compute_Phi_Solids");
    ARRAY<TV_INT> seed_indices;
    phi_all_solids_negated.Fill(-FLT_MAX);
    outside_fluid.Fill(false);
    T_FACE_ARRAYS_BOOL kinematic_faces(euler.grid.Domain_Indices(1));kinematic_faces.Fill(false);

    for(COLLISION_GEOMETRY_ID id(1);id<=collision_bodies_affecting_fluid->collision_geometry_collection.bodies.m;id++)
        if(collision_bodies_affecting_fluid->collision_geometry_collection.Is_Active(id)){
            COLLISION_GEOMETRY<TV>& body=*(collision_bodies_affecting_fluid->collision_geometry_collection.bodies(id));
            T collision_thickness_over_two=(T).5*collision_bodies_affecting_fluid->collision_thickness;
            RANGE<TV_INT> bounding_box(euler.grid.Clamp_To_Cell(body.Axis_Aligned_Bounding_Box().Thickened(
                                    euler.grid.dX.Max()*(T)2),number_of_cells_to_extrapolate));
            for(CELL_ITERATOR iterator(euler.grid,bounding_box);iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();T phi_value=
                    body.Implicit_Geometry_Extended_Value(iterator.Location());
                if(body.Inside(iterator.Location(),collision_thickness_over_two)) outside_fluid(cell_index)=true;
                phi_all_solids_negated(cell_index)=max(-phi_value,phi_all_solids_negated(cell_index));}}

    if(use_fast_marching){
        collision_bodies_affecting_fluid->outside_fluid=&outside_fluid;
        UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV> iterator_info(*collision_bodies_affecting_fluid);
        iterator_info.Initialize_Collision_Aware_Face_Iterator(outside_fluid,kinematic_faces,7,false);

        for(UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<TV> iterator(iterator_info);iterator.Valid();iterator.Next()){
            seed_indices.Append(iterator.First_Cell_Index());seed_indices.Append(iterator.Second_Cell_Index());}

        typename LEVELSET_POLICY<GRID<TV> >::LEVELSET levelset(euler.grid,phi_all_solids_negated);
        FAST_MARCHING_METHOD_UNIFORM<GRID<TV> > fmm(levelset,number_of_cells_to_extrapolate);
        fmm.Fast_Marching_Method(phi_all_solids_negated,euler.grid.dX.Max()*(T)number_of_cells_to_extrapolate,
                                 &seed_indices,true);}
}
//#####################################################################
// Function Fill_Solid_Cells
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Fill_Solid_Cells(bool fill_pressure_only)
{
    if(euler.timesplit && thinshell) return;
    if(collision_bodies_affecting_fluid){
        if(thinshell) Fill_Uncovered_Cells();
        int number_of_ghost_cells=mpi_grid?3:0;
        Compute_Phi_Solids(number_of_ghost_cells);
        if(fill_pressure_only){LOG::filecout("Defunct code\n");exit(-1);}
        else Extrapolate_State_Into_Solids(phi_all_solids_negated,number_of_ghost_cells,number_of_cells_to_extrapolate);}
}
//#####################################################################
// Function Project_Fluid_Pressure_At_Neumann_Faces
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Project_Fluid_Pressure_At_Neumann_Faces(const T_ARRAYS_SCALAR& p_ghost,T_FACE_ARRAYS_SCALAR& p_face) const
{
    // Bp
    const RANGE<TV>& domain=euler.grid.domain;
    for(FACE_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
        TV_INT face_index=iterator.Face_Index(),first_cell_index=iterator.First_Cell_Index(),
               second_cell_index=iterator.Second_Cell_Index();
        if(euler.euler_projection.elliptic_solver->psi_N.Component(axis)(face_index)){
            int direction;
            // TODO: This doesn't work with thin-shells
            bool first_cell_inside_solid=euler.euler_projection.elliptic_solver->psi_D(iterator.First_Cell_Index()) ||
                !domain.Lazy_Inside(iterator.First_Cell_Center());
            bool second_cell_inside_solid=euler.euler_projection.elliptic_solver->psi_D(iterator.Second_Cell_Index()) ||
                !domain.Lazy_Inside(iterator.Second_Cell_Center());

            if(!first_cell_inside_solid && second_cell_inside_solid) direction=1; // solid on right
            else if(first_cell_inside_solid && !second_cell_inside_solid) direction=-1; // solid on left
            else if(!first_cell_inside_solid && !second_cell_inside_solid){
                direction=0;
                LOG::cerr<<"Warning: Neumann face has fluid on both sides, face_index="<<
                    face_index<<", axis="<<axis<<std::endl;}
            else direction=0; // solid on both sides
    
            if(direction!=0){
                if(direction==1) p_face.Component(axis)(face_index)=p_ghost(first_cell_index);
                else p_face.Component(axis)(face_index)=p_ghost(second_cell_index);}}}
}
//#####################################################################
// Function Apply_Isobaric_Fix
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Apply_Isobaric_Fix(const T dt,const T time)
{
    euler.Fill_Ghost_Cells(dt,time,3);
    for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(euler.psi(cell_index) && phi_all_solids_negated(cell_index)<0){
            bool encountered_neumann_face=false;TV_INT reference_point;
            for(int axis=1;axis<=T_GRID::dimension;axis++){
                TV_INT first_face_index=iterator.First_Face_Index(axis),second_face_index=iterator.Second_Face_Index(axis);
                if(euler.euler_projection.elliptic_solver->psi_N.Component(axis)(first_face_index) &&
                        euler.euler_projection.elliptic_solver->psi_N.Component(axis)(second_face_index)) continue;
                else if(euler.euler_projection.elliptic_solver->psi_N.Component(axis)(first_face_index)){
                    encountered_neumann_face=true;reference_point=iterator.Cell_Neighbor(2*axis);}
                else if(euler.euler_projection.elliptic_solver->psi_N.Component(axis)(second_face_index)){
                    encountered_neumann_face=true;reference_point=iterator.Cell_Neighbor(2*axis-1);}}
            if(encountered_neumann_face){
                std::stringstream ss;
                ss<<"ISOBARIC FIX: fixing cell "<<cell_index<<" with reference cell "<<reference_point<<std::endl;
                LOG::filecout(ss.str());
                T rho=euler.U(cell_index)(1);
                TV velocity=EULER<T_GRID>::Get_Velocity(euler.U,cell_index);
                T e=EULER<T_GRID>::e(euler.U,cell_index);
                T p_cell=euler.eos->p(rho,e);

                T rho_reference=euler.U_ghost(reference_point)(1);
                T e_reference=EULER<T_GRID>::e(euler.U_ghost,reference_point);
                T p_reference=euler.eos->p(rho_reference,e_reference);

                if(p_cell>p_reference) rho=rho_reference*sqrt(p_cell/p_reference); //isobaric fix
                else rho=rho_reference;
                //rho=rho_reference*pow(p_cell/p_reference,(T)(1/1.4)); //isobaric fix (constant entropy)
                e=euler.eos->e_From_p_And_rho(p_cell,rho);
                EULER<T_GRID>::Set_Euler_State_From_rho_velocity_And_internal_energy(euler.U,cell_index,rho,velocity,e);}}}
    euler.Invalidate_Ghost_Cells();
}
//#####################################################################
// Function Extract_Time_N_Data_For_Explicit_Fluid_Forces
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Extract_Time_N_Data_For_Explicit_Fluid_Forces()
{
    if(!fluid_affects_solid || euler.timesplit) return;
    T_ARRAYS_SCALAR p_approx(euler.grid.Domain_Indices(1));
    for(CELL_ITERATOR iterator(euler.grid,1);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        p_approx(cell_index)=euler.eos->p(euler.U_ghost(cell_index)(1),euler.e(euler.U_ghost,cell_index));}
    euler.euler_projection.Compute_Face_Pressure_From_Cell_Pressures(
            euler.grid,euler.U_ghost,euler.psi,pressure_at_faces,p_approx);

    solid_fluid_face_time_n.Fill(false);
    collision_bodies_affecting_fluid->Compute_Psi_N(solid_fluid_face_time_n,0);
    cells_inside_fluid_time_n=euler.psi;
}
//#####################################################################
// Function Snapshot_State
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Snapshot_State(const T_ARRAYS_DIMENSION_SCALAR& U_ghost)
{
    T_ARRAYS_DIMENSION_SCALAR::Put(U_ghost,U_n);
    for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next())
        if(!euler.psi(iterator.Cell_Index())) U_n(iterator.Cell_Index())=TV_DIMENSION();

    accumulated_flux.Fill(TV_DIMENSION());
    TV_DIMENSION material;
    for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next())
        material+=U_n(iterator.Cell_Index())*cell_volumes_np1(iterator.Cell_Index());
    std::stringstream ss;
    ss<<std::setprecision(16)<<"TOTAL ACCUMULATED MATERIAL = "<<material<<std::endl;
    LOG::filecout(ss.str());

    near_interface.Fill(false);
}
//#####################################################################
// Function Initialize_Collision_Data
//#####################################################################
namespace {
template<class T,int d>
void Adjust_Psi_And_Cut_Cells_For_Consistency(const GRID<VECTOR<T,d> >& grid,const int num_ghost_cells,
                                              const ARRAY<bool,VECTOR<int,d> >& psi,ARRAY<CUT_CELLS<T,d>*,
                                              VECTOR<int,d> >& cut_cells)
{
    typedef VECTOR<T,d> TV;typedef VECTOR<int,d> TV_INT;
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(!psi(cell_index) && cut_cells(cell_index) && cut_cells(cell_index)->dominant_element){
            for(int i=cut_cells(cell_index)->visibility(cut_cells(cell_index)->dominant_element).Size();i >= 1;--i)
                if(cut_cells(cell_index)->visibility(cut_cells(cell_index)->dominant_element)(i)==cell_index)
                    cut_cells(cell_index)->visibility(cut_cells(cell_index)->dominant_element).Remove_Index(i);
            cut_cells(cell_index)->dominant_element=0;}
        else if(psi(cell_index) && cut_cells(cell_index) && !cut_cells(cell_index)->dominant_element){
            std::stringstream ss;
            ss<<"WARNING: Adjusting inconsistency in cell "<<cell_index
                <<" which has psi = true, but no dominant element"<<std::endl;
            ARRAY<int> candidate_polygons;
            for(int i=1;i<=cut_cells(cell_index)->geometry.Size();++i)
                for(int j=1;j<=cut_cells(cell_index)->visibility(i).Size();++j){
                    if(psi(cut_cells(cell_index)->visibility(i)(j))) candidate_polygons.Append(i);}
            if(candidate_polygons.Size()==1){
                ss<<"WARNING: \tFound candidate correction; polygon "<<candidate_polygons(1)<<std::endl;
                cut_cells(cell_index)->visibility(candidate_polygons(1)).Append(cell_index);
                cut_cells(cell_index)->dominant_element=candidate_polygons(1);}
            else if(!candidate_polygons.Size()){
                ss<<"ERROR: \tCell "<<cell_index
                    <<" cannot locate ANY visible neighbors; don't know what to do!!!"<<std::endl;}
            else{
                for(int i=1;i<=candidate_polygons.Size();++i)
                    if(cut_cells(cell_index)->geometry(candidate_polygons(i)).Inside_Polygon(grid.Center(cell_index))){
                        ss<<"WARNING: \tFound candidate correction; polygon "<<candidate_polygons(i)<<std::endl;
                        cut_cells(cell_index)->visibility(candidate_polygons(i)).Append(cell_index);
                        cut_cells(cell_index)->dominant_element=candidate_polygons(i);break;}
                if(!cut_cells(cell_index)->dominant_element){ss<<"WARNING Cell "<<cell_index
                    <<"'s cell center does not lie in any cut cell polygons; picking based on psi neighbors!"<<std::endl;
                    for(int poly=1;!cut_cells(cell_index)->dominant_element && poly<=candidate_polygons.Size();++poly)
                        for(int neighbor=1;neighbor<=cut_cells(cell_index)->visibility(poly).Size();++neighbor)
                            if(psi(cut_cells(cell_index)->visibility(poly)(neighbor))){
                                cut_cells(cell_index)->visibility(poly).Append(cell_index);
                                cut_cells(cell_index)->dominant_element=poly;break;}}
                if(!cut_cells(cell_index)->dominant_element){ss<<"ERROR: Cell "<<cell_index
                    <<"'s cell center does not lie in any cut cell polygons!"<<std::endl;
                    CUT_CELLS<T,d>::Print_Debug_Information(cut_cells(cell_index)->geometry);}}
            LOG::filecout(ss.str());}}

    for(typename GRID<TV>::CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(cut_cells(cell_index)){
            for(int poly=cut_cells(cell_index)->geometry.Size();poly>=1;--poly){
                bool has_psi_neighbors=false;
                for(int neighbor=cut_cells(cell_index)->visibility(poly).Size();neighbor>=1;--neighbor){
                    if(psi(cut_cells(cell_index)->visibility(poly)(neighbor))) has_psi_neighbors=true;
                    else cut_cells(cell_index)->visibility(poly).Remove_Index(neighbor);}
                if(!has_psi_neighbors){
                    if(cut_cells(cell_index)->dominant_element > poly) --cut_cells(cell_index)->dominant_element;
                    else if(cut_cells(cell_index)->dominant_element == poly) cut_cells(cell_index)->dominant_element=0;
                    cut_cells(cell_index)->visibility.Remove_Index(poly);
                    cut_cells(cell_index)->visibility_nodes.Remove_Index(poly);
                    cut_cells(cell_index)->geometry.Remove_Index(poly);}}
            if(!cut_cells(cell_index)->geometry.Size()){delete cut_cells(cell_index);cut_cells(cell_index)=0;}}}
}

template<class T,int d>
void Compute_Psi(const GRID<VECTOR<T,d> >& grid,const int num_ghost_cells,
                 const COLLISION_GEOMETRY_COLLECTION<VECTOR<T,d> >& collision_geometry_collection,
                 ARRAY<bool,VECTOR<int,d> >& psi)
{
    LOG::SCOPE scope("Compute_Psi");
    typedef VECTOR<T,d> TV;typedef VECTOR<int,d> TV_INT;
    typedef UNIFORM_GRID_ITERATOR_CELL<TV> CELL_ITERATOR;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX T_SIMPLEX;
    psi.Fill(true);
    int dummy_index;
    for(COLLISION_GEOMETRY_ID i(1);i<=collision_geometry_collection.Size();i++)if(collision_geometry_collection.Is_Active(i)){
        T collision_thickness_over_two=(T).5*collision_geometry_collection.collision_body_thickness;
        if(DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV> const* object=
           dynamic_cast<DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV> const*>(&collision_geometry_collection(i))){
            if(!object->volume_object){
                BOX<TV>& box(*(object->object.bounding_box));
                TV_INT min_index=grid.Clamp_To_Cell(box.Minimum_Corner(),num_ghost_cells);
                TV_INT max_index=grid.Clamp_To_Cell(box.Maximum_Corner(),num_ghost_cells);
                int dummy_index;
                for(CELL_ITERATOR iterator(grid,RANGE<TV_INT>(min_index,max_index));iterator.Valid();iterator.Next())
                    if(object->object.Inside_Any_Simplex(iterator.Location(),dummy_index,collision_thickness_over_two)){
                        psi(iterator.Cell_Index())=false;}}
            else{
                for(int i=1;i<=object->volume_object->mesh.elements.m;i++){
                    const T_SIMPLEX& simplex=object->volume_object->Get_Element(i);
                    RANGE<TV> box(simplex.Bounding_Box());
                    TV_INT min_index=grid.Clamp_To_Cell(box.Minimum_Corner(),num_ghost_cells);
                    TV_INT max_index=grid.Clamp_To_Cell(box.Maximum_Corner(),num_ghost_cells);
                    for(CELL_ITERATOR iterator(grid,RANGE<TV_INT>(min_index,max_index));iterator.Valid();iterator.Next())
                        if(simplex.Inside(iterator.Location(),collision_thickness_over_two)){
                            psi(iterator.Cell_Index())=false;}}}}
        else if(RIGID_COLLISION_GEOMETRY_BASE<TV> const* object=
                dynamic_cast<RIGID_COLLISION_GEOMETRY_BASE<TV> const*>(&collision_geometry_collection(i))){
            const RIGID_BODY<TV>& rigid_body=dynamic_cast<const RIGID_BODY<TV>&>(object->rigid_geometry);
            if(!rigid_body.simplicial_object->mesh.incident_elements)
                rigid_body.simplicial_object->mesh.Initialize_Incident_Elements();
            if(!rigid_body.simplicial_object->mesh.adjacent_elements)
                rigid_body.simplicial_object->mesh.Initialize_Adjacent_Elements();
            RANGE<TV> box(object->Axis_Aligned_Bounding_Box());
            TV_INT min_index=grid.Clamp_To_Cell(box.Minimum_Corner(),num_ghost_cells);
            TV_INT max_index=grid.Clamp_To_Cell(box.Maximum_Corner(),num_ghost_cells);
            if(rigid_body.thin_shell){
                for(CELL_ITERATOR iterator(grid,RANGE<TV_INT>(min_index,max_index));iterator.Valid();iterator.Next()){
                    if(object->Inside_Any_Simplex(iterator.Location(),dummy_index)){psi(iterator.Cell_Index())=false;}}}
            else{
                for(CELL_ITERATOR iterator(grid,RANGE<TV_INT>(min_index,max_index));iterator.Valid();iterator.Next()){
                    if(object->Inside(iterator.Location(),collision_thickness_over_two)){psi(iterator.Cell_Index())=false;}}}}
        else PHYSBAM_FATAL_ERROR("Unrecognized collision body type");}
}

template<class T,int d>
void Compute_Cell_Volumes(const GRID<VECTOR<T,d> >& grid,const int num_ghost_cells,const ARRAY<bool,VECTOR<int,d> >& psi,
                          const ARRAY<CUT_CELLS<T,d>*,VECTOR<int,d> >& cut_cells,ARRAY<T,VECTOR<int,d> >& cell_volumes)
{
    typedef VECTOR<T,d> TV;typedef VECTOR<int,d> TV_INT;
    typedef UNIFORM_GRID_ITERATOR_CELL<TV> CELL_ITERATOR;
    cell_volumes.Fill(grid.Cell_Size());
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(cut_cells(cell_index))
            for(int poly=1;poly<=cut_cells(cell_index)->geometry.Size();++poly){
                POLYGON<TV>& polygon=cut_cells(cell_index)->geometry(poly);T volume=polygon.Area();
                if(poly==cut_cells(cell_index)->dominant_element){
                    cell_volumes(cell_index)+=volume-grid.Cell_Size();continue;}
                if(!cut_cells(cell_index)->visibility(poly).Size()){
                    std::stringstream ss;
                    ss<<"WARNING: No visible neighbors for cut cell "<<cell_index<<", poly "<<poly<<": "<<std::endl;
                    CUT_CELLS<T,TV::dimension>::Print_Debug_Information(cut_cells(cell_index)->geometry,
                                                                        iterator.Bounding_Box());
                    for(int i=1;i<=polygon.X.Size();++i)
                        ss<<polygon.X(i)<<", ";ss<<"discarding "<<volume<<" volume"<<std::endl;continue;LOG::filecout(ss.str());}
                volume /= cut_cells(cell_index)->visibility(poly).Size();
                for(int n=1;n<=cut_cells(cell_index)->visibility(poly).Size();++n){
                    cell_volumes(cut_cells(cell_index)->visibility(poly)(n)) += volume;}}}

    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(!psi(cell_index) || (cut_cells(cell_index) && !cut_cells(cell_index)->dominant_element)) cell_volumes(cell_index)=0;}
}

template<class T,int d>
void Compute_Swept_Cells(const GRID<VECTOR<T,d> >& grid,const int num_ghost_cells,const T dt,
                         const COLLISION_GEOMETRY_COLLECTION<VECTOR<T,d> >& collision_geometry_collection,
                         const ARRAY<bool,VECTOR<int,d> >& psi_n,const ARRAY<bool,VECTOR<int,d> >& psi_np1,
                         ARRAY<bool,VECTOR<int,d> >& swept_cells)
{
    LOG::SCOPE scope("Compute_Swept_Cells");
    typedef VECTOR<T,d> TV;typedef VECTOR<int,d> TV_INT;
    typedef UNIFORM_GRID_ITERATOR_CELL<TV> CELL_ITERATOR;
    swept_cells.Fill(false);
    for(COLLISION_GEOMETRY_ID id(1);id<=collision_geometry_collection.bodies.m;id++)
        if(collision_geometry_collection.Is_Active(id)){
            COLLISION_GEOMETRY<TV>& collision_body=*(collision_geometry_collection.bodies(id));
            RANGE<TV_INT> bounding_box(grid.Clamp_To_Cell(
                            collision_body.Axis_Aligned_Bounding_Box().Thickened(grid.dX.Max()*(T)2),num_ghost_cells));
            for(CELL_ITERATOR iterator(grid,bounding_box);iterator.Valid();iterator.Next()){
                const TV_INT cell_index=iterator.Cell_Index();
                if(collision_body.Any_Simplex_Crossover(iterator.Location(),iterator.Location(),dt))
                    swept_cells(cell_index)=true;
                else if(psi_n(cell_index) != psi_np1(cell_index)){
                    std::stringstream ss;
                    ss<<"WARNING: psi_n = "<<psi_n(cell_index)<<", psi_np1 = "<<psi_np1(cell_index)
                        <<" but "<<cell_index<<" is NOT a swept cell... fixing!"<<std::endl;swept_cells(cell_index)=true;LOG::filecout(ss.str());}}}
}
};
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Initialize_Collision_Data()
{
    CUT_CELL_COMPUTATIONS::Compute_Cut_Geometries(euler.grid,1,*collision_bodies_affecting_fluid,cut_cells_np1);
    Compute_Psi<T,TV::dimension>(euler.grid,1,collision_bodies_affecting_fluid->collision_geometry_collection,psi_np1);
    Adjust_Psi_And_Cut_Cells_For_Consistency(euler.grid,1,psi_np1,cut_cells_np1);
    Compute_Cell_Volumes<T,TV::dimension>(euler.grid,1,psi_np1,cut_cells_np1,cell_volumes_np1);
    uncovered_cells.Fill(false);
}
//#####################################################################
// Function Update_Np1_Collision_Data
//#####################################################################
namespace{
void Get_Neighbor_Offsets(ARRAY<VECTOR<int,1> >& neighbor_offsets)
{
    neighbor_offsets.Append(VECTOR<int,1>(-1));neighbor_offsets.Append(VECTOR<int,1>(0));
    neighbor_offsets.Append(VECTOR<int,1>(1));
}

void Get_Neighbor_Offsets(ARRAY<VECTOR<int,2> >& neighbor_offsets)
{
    neighbor_offsets.Append(VECTOR<int,2>(-1,-1));neighbor_offsets.Append(VECTOR<int,2>(-1, 0));
    neighbor_offsets.Append(VECTOR<int,2>(-1,1));neighbor_offsets.Append(VECTOR<int,2>( 0,-1));
    neighbor_offsets.Append(VECTOR<int,2>( 0, 0));neighbor_offsets.Append(VECTOR<int,2>( 0,1));
    neighbor_offsets.Append(VECTOR<int,2>( 1,-1));neighbor_offsets.Append(VECTOR<int,2>( 1, 0));
    neighbor_offsets.Append(VECTOR<int,2>( 1,1));
}

void Get_Neighbor_Offsets(ARRAY<VECTOR<int,3> >& neighbor_offsets)
{PHYSBAM_NOT_IMPLEMENTED("3-D cut cells not implemented!");}
};
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Update_Np1_Collision_Data(const T dt)
{
    LOG::SCOPE scope("Update_Np1_Collision_Data");
    collision_bodies_affecting_fluid->Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE);
    collision_bodies_affecting_fluid->Update_Intersection_Acceleration_Structures(
            true,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE,
            COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);
    collision_bodies_affecting_fluid->Rasterize_Objects();
    collision_bodies_affecting_fluid->Compute_Occupied_Blocks(false,(T)2*euler.grid.Minimum_Edge_Length(),5);

    cut_cells_n.Delete_Pointers_And_Clean_Memory();cut_cells_n.Resize(euler.grid.Domain_Indices(1));
    T_ARRAYS_CUT_CELLS::Put(cut_cells_np1,cut_cells_n);cut_cells_np1.Fill(0);
    CUT_CELL_COMPUTATIONS::Compute_Cut_Geometries(euler.grid,1,*collision_bodies_affecting_fluid,cut_cells_np1);

#if 0
    {LOG::SCOPE scope("Inside cells");LOG::cerr<<"Inside cells"<<std::endl;
    for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){const TV_INT& cell_index=iterator.Cell_Index();
        if(cut_cells_np1(cell_index)) for(int i=1;i<=cut_cells_np1(cell_index)->geometry.Size();++i)
            if(RIGID_COLLISION_GEOMETRY<TV> const* object=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV> const*>(
                            &(collision_bodies_affecting_fluid->collision_geometry_collection(COLLISION_GEOMETRY_ID(1))))){
                const RIGID_BODY<TV>& rigid_body=dynamic_cast<const RIGID_BODY<TV>&>(object->rigid_geometry);
                ARRAY<TV> candidate_nodes(cut_cells_np1(cell_index)->visibility_nodes(i));
                int dummy_index;
                for(int j=1;j<=candidate_nodes.Size();++j)
                    if(rigid_body.Simplex_Inside(candidate_nodes(j),
                                                 (T).5*collision_bodies_affecting_fluid->collision_thickness) &&
                       !object->Inside_Any_Simplex(candidate_nodes(j),dummy_index)){
                        LOG::cerr<<"GEOMETRY DEBUG INFO: "<<std::endl;
                        LOG::cerr<<"GEOMETRY DEBUG INFO: "<<candidate_nodes(j)<<std::endl;
                        CUT_CELLS<T,TV::dimension>::Print_Debug_Information(cut_cells_np1(cell_index)->geometry(i));break;}}}}
#endif

    T_ARRAYS_BOOL::Exchange_Arrays(psi_np1,psi_n);
    Compute_Psi<T,TV::dimension>(euler.grid,1,collision_bodies_affecting_fluid->collision_geometry_collection,psi_np1);
    Adjust_Psi_And_Cut_Cells_For_Consistency(euler.grid,1,psi_np1,cut_cells_np1);
    
    T_ARRAYS_SCALAR::Exchange_Arrays(cell_volumes_np1,cell_volumes_n);
    Compute_Cell_Volumes<T,TV::dimension>(euler.grid,1,psi_np1,cut_cells_np1,cell_volumes_np1);

    collision_bodies_affecting_fluid->Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);
    Compute_Swept_Cells<T,TV::dimension>(euler.grid,1,dt,collision_bodies_affecting_fluid->collision_geometry_collection,
                                         psi_n,psi_np1,uncovered_cells);

    ARRAY<TV_INT> neighbor_offsets;Get_Neighbor_Offsets(neighbor_offsets);
    CUT_CELL_COMPUTATIONS::Compute_Temporal_Neighbors<T,TV::dimension>(
            euler.grid,1,dt,neighbor_offsets,*collision_bodies_affecting_fluid,
            cut_cells_n,cut_cells_np1,psi_np1,visibility_np1);
}
//#####################################################################
// Function Compute_Intermediate_Solid_Position_Data
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Compute_Intermediate_Solid_Position_Data(const T dt)
{
    LOG::SCOPE scope("Compute_Intermediate_Solid_Position_Data");
    collision_bodies_affecting_fluid->Average_States(
            COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,
            COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_SAVED_NEW_STATE,
            COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE, (T).5);
    collision_bodies_affecting_fluid->Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE);
    collision_bodies_affecting_fluid->Update_Intersection_Acceleration_Structures(
            true,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE,
            COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);
    collision_bodies_affecting_fluid->Rasterize_Objects();
    collision_bodies_affecting_fluid->Compute_Occupied_Blocks(false,(T)2*euler.grid.Minimum_Edge_Length(),5);

    cut_cells_n_p_half.Delete_Pointers_And_Clean_Memory();
    cut_cells_n_p_half.Resize(euler.grid.Domain_Indices(1)); cut_cells_n_p_half.Fill(0);
    CUT_CELL_COMPUTATIONS::Compute_Cut_Geometries(euler.grid,1,*collision_bodies_affecting_fluid,cut_cells_n_p_half);

    Compute_Psi<T,TV::dimension>(euler.grid,1,collision_bodies_affecting_fluid->collision_geometry_collection,psi_n_p_half);
    Adjust_Psi_And_Cut_Cells_For_Consistency(euler.grid,1,psi_n_p_half,cut_cells_n_p_half);

    Compute_Cell_Volumes<T,TV::dimension>(euler.grid,1,psi_n_p_half,cut_cells_n_p_half,cell_volumes_n_p_half);

    collision_bodies_affecting_fluid->Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);
    Compute_Swept_Cells<T,TV::dimension>(euler.grid,1,dt,collision_bodies_affecting_fluid->collision_geometry_collection,
                                         psi_n,psi_n_p_half,uncovered_cells_n_p_half);

    ARRAY<TV_INT> neighbor_offsets;Get_Neighbor_Offsets(neighbor_offsets);
    CUT_CELL_COMPUTATIONS::Compute_Temporal_Neighbors<T,TV::dimension>(
            euler.grid,1,dt,neighbor_offsets,*collision_bodies_affecting_fluid,
            cut_cells_n,cut_cells_n_p_half,psi_n_p_half,visibility_n_p_half);

    collision_bodies_affecting_fluid->Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_SAVED_NEW_STATE);
    collision_bodies_affecting_fluid->Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE);
}
//#####################################################################
// Function Revert_Cells_Near_Interface
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Revert_Cells_Near_Interface(const int iteration_number)
{
    if(iteration_number==1){
        ARRAY<TV_INT> neighbor_offsets;Get_Neighbor_Offsets(neighbor_offsets);
        for(CELL_ITERATOR iterator(euler.grid,1);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            near_interface(cell_index) = (cut_cells_n(cell_index)!=0 || cut_cells_n_p_half(cell_index)!=0
                                          || cut_cells_np1(cell_index)!=0);}
        for(CELL_ITERATOR iterator(euler.grid,1);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            for(int i=1;i<=neighbor_offsets.Size();++i){TV_INT neighbor_index=cell_index+neighbor_offsets(i);
                if(!near_interface.Valid_Index(neighbor_index)) continue;
                if(cut_cells_n(cell_index))
                    for(int p=1;p<=cut_cells_n(cell_index)->geometry.Size() && !near_interface(neighbor_index);++p)
                        if(p != cut_cells_n(cell_index)->dominant_element
                           && cut_cells_n(cell_index)->visibility(p).Contains(neighbor_index))
                            near_interface(neighbor_index)=true;
                if(cut_cells_n_p_half(cell_index))
                    for(int p=1;p<=cut_cells_n_p_half(cell_index)->geometry.Size() && !near_interface(neighbor_index);++p)
                        if(p != cut_cells_n_p_half(cell_index)->dominant_element
                           && cut_cells_n_p_half(cell_index)->visibility(p).Contains(neighbor_index))
                        near_interface(neighbor_index)=true;
                if(cut_cells_np1(cell_index))
                    for(int p=1;p<=cut_cells_np1(cell_index)->geometry.Size() && !near_interface(neighbor_index);++p)
                        if(p != cut_cells_np1(cell_index)->dominant_element
                           && cut_cells_np1(cell_index)->visibility(p).Contains(neighbor_index))
                            near_interface(neighbor_index)=true;}}}
}
//#####################################################################
// Function Update_Cells_Near_Interface
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Update_Cells_Near_Interface(const T dt,const int rk_order,const int rk_substep)
{
    for(FACE_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){
        FACE_INDEX<TV::dimension> face_index=iterator.Full_Index();
        accumulated_flux(face_index) += euler.conservation->fluxes(face_index);}

    if(rk_order==2 && rk_substep==2) accumulated_flux *= (T).5;
    else if(rk_order==3){
        if(rk_substep==2) accumulated_flux *= (T).5;
        if(rk_substep==3) accumulated_flux *= (T)2./3;}

    const TV dt_over_dx=((rk_order==3 && rk_substep==2) ? dt/(T)2 : dt)*euler.grid.One_Over_DX();
    for(CELL_ITERATOR iterator(euler.grid,1);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(psi_n(cell_index) && near_interface(cell_index)) euler.U(cell_index)=U_n(cell_index);}

    TV_DIMENSION material_n=TV_DIMENSION();
    if(rk_order==rk_substep){
        for(CELL_ITERATOR iterator(euler.grid,1);iterator.Valid();iterator.Next()) if(psi_n(iterator.Cell_Index()))
            material_n+=U_n(iterator.Cell_Index())*cell_volumes_n(iterator.Cell_Index());}

#if 1
    for(CELL_ITERATOR iterator(euler.grid,1);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        T cell_size=euler.grid.Cell_Size();
        if(psi_n(cell_index) && !near_interface(cell_index)
           && (abs(cell_volumes_n(cell_index)-cell_size) > (T)10*std::numeric_limits<T>::epsilon()
               || abs(cell_volumes_np1(cell_index)-cell_size) > (T)10*std::numeric_limits<T>::epsilon())){
            std::stringstream ss;
            ss<<"ERROR: Cell "<<cell_index<<" should have cell volume "<<cell_size
                <<"; but instead has time n size "<<cell_volumes_n(cell_index)<<" and time np1 size "
              <<cell_volumes_np1(cell_index)<<std::endl;
            LOG::filecout(ss.str());}}

    {euler.Invalidate_Ghost_Cells();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After FSI Reset (part 0)",0,0);
    ARRAY<bool,FACE_INDEX<TV::dimension> > psi_n_copy(euler.euler_projection.elliptic_solver->psi_N);
    ARRAY<bool,TV_INT> psi_d_copy(euler.euler_projection.elliptic_solver->psi_D);
    euler.euler_projection.elliptic_solver->psi_N.Fill(false);
    euler.euler_projection.elliptic_solver->psi_D.Fill(false);
    ARRAY<TRIPLE<TV_INT,TV_INT,TV_DIMENSION> > hybrid_flux_data;
    Compute_Hybrid_Boundary_Fluxes(euler.grid,dt,near_interface,euler.psi,accumulated_flux,hybrid_flux_data);
    for(int i=1;i<=hybrid_flux_data.Size();++i){
        int axis=(hybrid_flux_data(i).y-hybrid_flux_data(i).x).Arg_Max();
        euler.euler_projection.elliptic_solver->psi_N.Component(axis)(hybrid_flux_data(i).y)=true;}
    for(CELL_ITERATOR iter(euler.grid);iter.Valid();iter.Next())
        euler.euler_projection.elliptic_solver->psi_D(iter.Cell_Index())=
            (euler.psi(iter.Cell_Index()) || psi_np1(iter.Cell_Index())) && near_interface(iter.Cell_Index());
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After FSI Reset",0,0);
    euler.euler_projection.elliptic_solver->psi_N=psi_n_copy;
    euler.euler_projection.elliptic_solver->psi_D=psi_d_copy;
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After FSI Reset (part 2)",0,0);
    }
#endif

    if(rk_order==3 && rk_substep==2)
        Advect_Near_Interface_Data(euler.grid,collision_bodies_affecting_fluid->collision_thickness,dt/(T)2,
                                   near_interface,uncovered_cells_n_p_half,accumulated_flux,cell_volumes_n_p_half,
                                   visibility_n_p_half,psi_n,cut_cells_n,U_n,psi_n_p_half,cut_cells_n_p_half,euler.U);
    else
        Advect_Near_Interface_Data(euler.grid,collision_bodies_affecting_fluid->collision_thickness,dt     ,
                                   near_interface,uncovered_cells,         accumulated_flux,cell_volumes_np1,
                                   visibility_np1,     psi_n,cut_cells_n,U_n,psi_np1,     cut_cells_np1,     euler.U);

#if 1
    TV_DIMENSION material_np1=TV_DIMENSION();

    TV_DIMENSION boundary_flux=TV_DIMENSION();
    if(rk_order==rk_substep){
        if(rk_order==rk_substep){
            for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
                if(psi_np1(cell_index)) material_np1+=
                    euler.U(cell_index)*((rk_order==3 && rk_substep==2) ?
                                         cell_volumes_n_p_half(cell_index) : cell_volumes_np1(cell_index));}}

        for(FACE_ITERATOR iterator(euler.grid,0,GRID<TV>::BOUNDARY_REGION);iterator.Valid();iterator.Next()){
            if(euler.psi.Valid_Index(iterator.First_Cell_Index())       && euler.psi(iterator.First_Cell_Index())
               && !near_interface(iterator.First_Cell_Index())){
                boundary_flux -= accumulated_flux(iterator.Full_Index())*dt_over_dx(iterator.Axis())*euler.grid.Cell_Size();}
            else if(euler.psi.Valid_Index(iterator.Second_Cell_Index()) && euler.psi(iterator.Second_Cell_Index())
                    && !near_interface(iterator.Second_Cell_Index())){
                boundary_flux += accumulated_flux(iterator.Full_Index())*dt_over_dx(iterator.Axis())*euler.grid.Cell_Size();}}
        std::stringstream ss;
        ss<<std::setprecision(16)<<"TOTAL ACCUMULATED BOUNDARY FLUX MATERIAL = "<<boundary_flux<<std::endl;LOG::filecout(ss.str());}

    if(rk_order==rk_substep &&
        (material_np1-material_n-boundary_flux).Magnitude_Squared() >= (T)10*std::numeric_limits<T>::epsilon()){
        std::stringstream ss;
        ss<<"ERROR: Conservation error in regular update; "<<rk_order<<", "<<rk_substep
            <<"\tdifference = "<<(material_np1-material_n-boundary_flux)<<std::endl
            <<"\t Time N   = "<<material_n<<std::endl<<"\t Time Np1 = "<<material_np1<<std::endl
          <<"\t Boundary Flux = "<<boundary_flux<<std::endl;LOG::filecout(ss.str());}
#endif

    euler.Invalidate_Ghost_Cells();
    if(rk_order==3 && rk_substep==2) accumulated_flux *= (T).5;

    PHYSBAM_DEBUG_WRITE_SUBSTEP("After FSI update",0,0);
}
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Compute_Post_Advected_Variables()
{
    for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(near_interface(cell_index)){
            std::stringstream ss;
            if(uncovered_cells(cell_index) && psi_np1(cell_index)){
                T one_over_c=euler.eos->one_over_c(euler.U(cell_index)(1),euler.e(euler.U(cell_index)));
#if 0
                ss<<"Computing 1/(rho c^2) for uncovered cell "<<cell_index<<" = "
                    <<one_over_c*one_over_c/(euler.U(cell_index)(1))<<std::endl;
#endif
                euler.euler_projection.one_over_rho_c_squared(cell_index) = one_over_c*one_over_c/(euler.U(cell_index)(1));}
#if 0
            if(euler.psi(cell_index)) ss<<"Updated Cell "<<cell_index<<" from "<<U_n(cell_index)<<" to "
                <<euler.U(cell_index)<<"\t\t 1/(rho c^2) = "<<euler.euler_projection.one_over_rho_c_squared(cell_index)
                <<std::endl;
#endif
            LOG::filecout(ss.str());
        }}
}
//#####################################################################
template class SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<VECTOR<float,1> >;
template class SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<VECTOR<float,2> >;
template class SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<VECTOR<double,1> >;
template class SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<VECTOR<double,2> >;
template class SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<VECTOR<double,3> >;
#endif

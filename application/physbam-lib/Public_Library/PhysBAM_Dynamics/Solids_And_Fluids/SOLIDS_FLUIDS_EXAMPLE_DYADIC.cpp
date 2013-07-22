#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_OCTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_QUADTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Collisions/GRID_BASED_COLLISION_GEOMETRY_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_DYADIC.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_OCTREE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_QUADTREE.h>
#include <PhysBAM_Dynamics/Interpolation/FACE_LOOKUP_FIRE_DYADIC.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_DYADIC.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::
SOLIDS_FLUIDS_EXAMPLE_DYADIC(const STREAM_TYPE stream_type,const typename FLUIDS_PARAMETERS<T_GRID>::TYPE type)
    :SOLIDS_FLUIDS_EXAMPLE<TV>(stream_type),fluids_parameters(type)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::
~SOLIDS_FLUIDS_EXAMPLE_DYADIC()
{}
//#####################################################################
// Function Add_Volumetric_Body_To_Fluid_Simulation
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::
Add_Volumetric_Body_To_Fluid_Simulation(RIGID_BODY<TV>& rigid_body,const bool affects_fluid,const bool affected_by_fluid)
{
    collision_body_affected_by_fluid.Append(affected_by_fluid);
    if(affected_by_fluid) PHYSBAM_FATAL_ERROR("TODO:not yet written for new coupling");
    if(affects_fluid){
        RIGID_COLLISION_GEOMETRY<TV>* collision_geometry=new RIGID_COLLISION_GEOMETRY<TV>(rigid_body);
        fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(collision_geometry,rigid_body.particle_index,true);}
}
//#####################################################################
// Function Add_Thin_Shell_To_Fluid_Simulation
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::
Add_Thin_Shell_To_Fluid_Simulation(RIGID_BODY<TV>& rigid_body,const bool affects_fluid,const bool affected_by_fluid)
{
    rigid_body.thin_shell=true;
    collision_body_affected_by_fluid.Append(affected_by_fluid);
    if(affected_by_fluid) PHYSBAM_FATAL_ERROR("TODO:not yet written for new coupling");
    if(affects_fluid){
        RIGID_COLLISION_GEOMETRY<TV>* collision_geometry=new RIGID_COLLISION_GEOMETRY<TV>(rigid_body);
        if(collision_geometry->collision_thickness<minimum_collision_thickness) collision_geometry->Set_Collision_Thickness(minimum_collision_thickness);
        fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(collision_geometry,rigid_body.particle_index,true);}
}
//#####################################################################
// Function Add_To_Fluid_Simulation
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::
Add_To_Fluid_Simulation(DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions,const bool affects_fluid,const bool affected_by_fluid)
{
    if(deformable_collisions.collision_thickness<minimum_collision_thickness) deformable_collisions.Set_Collision_Thickness(minimum_collision_thickness);
    collision_body_affected_by_fluid.Append(affected_by_fluid);
    if(affected_by_fluid) PHYSBAM_FATAL_ERROR("TODO:not yet written for new coupling");
    if(affects_fluid){
        fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(&deformable_collisions,0,false);
        deformable_collisions.Initialize_For_Thin_Shells_Fluid_Coupling();}
}
//#####################################################################
// Function Initialize_Solid_Fluid_Coupling
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::
Initialize_Solid_Fluid_Coupling()
{
    fluids_parameters.collision_bodies_affecting_fluid->Initialize_Grids();

#if 0
    bool initialize_valid_value_mask=!restart; // If we're restarting we'll be reading in the values anyway
    // TODO: update for new advections!
    fluids_parameters.semi_lagrangian_collidable_density.Initialize(*fluids_parameters.grid,initialize_valid_value_mask);
    fluids_parameters.semi_lagrangian_collidable_density.interpolation->Set_Default_Replacement_Value(fluids_parameters.density_container.ambient_density);
    fluids_parameters.semi_lagrangian_collidable_temperature.Initialize(*fluids_parameters.grid,initialize_valid_value_mask);
    fluids_parameters.semi_lagrangian_collidable_temperature.interpolation->Set_Default_Replacement_Value(fluids_parameters.temperature_container.ambient_temperature);
    fluids_parameters.semi_lagrangian_collidable_velocity.Initialize(*fluids_parameters.grid,initialize_valid_value_mask);
    if(fluids_parameters.water||fluids_parameters.fire) fluids_parameters.semi_lagrangian_collidable_phi.Initialize(*fluids_parameters.grid,initialize_valid_value_mask);
    if(fluids_parameters.water) fluids_parameters.semi_lagrangian_collidable_phi.interpolation->Set_Default_Replacement_Value((T)1e-5);
    if(fluids_parameters.fire){ // special treatment for fire because we advect different quantities using different velocity fields
        fluids_parameters.semi_lagrangian_collidable_phi.interpolation->Set_Default_Replacement_Value((T)1e-5);}
#endif
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::
Set_Dirichlet_Boundary_Conditions(const T time)
{
    if(fluids_parameters.water) fluids_parameters.incompressible.Set_Dirichlet_Boundary_Conditions(fluids_parameters.particle_levelset_evolution.phi);
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::
Get_Source_Velocities(const GEOMETRY& source,const TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity)
{
    for(FACE_ITERATOR iterator(*fluids_parameters.grid,fluids_parameters.grid->Map_All_Faces());iterator.Valid();iterator.Next()){
        int i=iterator.Face_Index();
        if(source.Lazy_Inside(world_to_source.Homogeneous_Times(iterator.Location()))){
            fluids_parameters.incompressible.projection.elliptic_solver->psi_N(i)=true;
            fluids_parameters.incompressible.projection.face_velocities(i)=constant_source_velocity[iterator.Axis()+1];}}
}
//#####################################################################
// Function Adjust_Phi_With_Source
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::
Adjust_Phi_With_Source(const GEOMETRY& source,const TRANSFORMATION_MATRIX& world_to_source)
{
    for(CELL_ITERATOR iterator(*fluids_parameters.grid,fluids_parameters.grid->number_of_ghost_cells);iterator.Valid();iterator.Next()){
        TV source_X=world_to_source.Homogeneous_Times(iterator.Location());
        if(source.Lazy_Inside(source_X))
            fluids_parameters.particle_levelset_evolution.phi(iterator.Cell_Index())=min(fluids_parameters.particle_levelset_evolution.phi(iterator.Cell_Index()),source.Signed_Distance(source_X));}
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::
Get_Source_Reseed_Mask(const GEOMETRY& source,const TRANSFORMATION_MATRIX& world_to_source,ARRAY<bool>*& cell_centered_mask,const bool reset_mask)
{
    ARRAY<CELL*>& cell_pointer_from_index=fluids_parameters.grid->Cell_Pointer_From_Index();
    T_GRID& grid=*fluids_parameters.grid;
    if(reset_mask){if(cell_centered_mask) delete cell_centered_mask;cell_centered_mask=new ARRAY<bool>(grid.number_of_cells);}
    T padding=3*grid.Minimum_Edge_Length();
    for(int i=1;i<=grid.number_of_cells;i++) if(!source.Outside(world_to_source.Homogeneous_Times(cell_pointer_from_index(i)->Center()),padding)) (*cell_centered_mask)(i)=true;
}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::
Adjust_Density_And_Temperature_With_Sources(const GEOMETRY& source,const TRANSFORMATION_MATRIX& world_to_source,const T source_density,const T source_temperature)
{
    T_GRID& grid=*fluids_parameters.grid;
    for(int i=1;i<=grid.number_of_nodes;i++)if(source.Lazy_Inside(world_to_source.Homogeneous_Times(grid.Node_Location(i)))){
        if(fluids_parameters.use_density) fluids_parameters.density_container.density(i)=source_density;
        if(fluids_parameters.use_temperature) fluids_parameters.temperature_container.temperature(i)=source_temperature;}
}
//#####################################################################
// Function Revalidate_Fluid_Scalars
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::
Revalidate_Fluid_Scalars()
{
    T_PARTICLE_LEVELSET_EVOLUTION& particle_levelset_evolution=fluids_parameters.particle_levelset_evolution;
    T_LEVELSET& levelset=particle_levelset_evolution.particle_levelset.levelset;
    T_LEVELSET_ADVECTION& levelset_advection=particle_levelset_evolution.levelset_advection;
    if(fluids_parameters.water || fluids_parameters.fire){
        assert(levelset_advection.nested_semi_lagrangian_collidable);
        levelset_advection.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(*fluids_parameters.grid,fluids_parameters.collidable_phi_replacement_value,levelset.phi);
    }
    if(fluids_parameters.use_density){
        assert(fluids_parameters.density_container.nested_semi_lagrangian_collidable);
        fluids_parameters.density_container.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(fluids_parameters.density_container.grid,
            fluids_parameters.density_container.ambient_density,fluids_parameters.density_container.density);}
    if(fluids_parameters.use_temperature){
        assert(fluids_parameters.temperature_container.nested_semi_lagrangian_collidable);
        fluids_parameters.temperature_container.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(fluids_parameters.temperature_container.grid,
            fluids_parameters.temperature_container.ambient_temperature,fluids_parameters.temperature_container.temperature);}
}
//#####################################################################
// Function Revalidate_Phi_After_Modify_Levelset
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::
Revalidate_Phi_After_Modify_Levelset()
{
    assert(fluids_parameters.water || fluids_parameters.fire);
    T_PARTICLE_LEVELSET_EVOLUTION& particle_levelset_evolution=fluids_parameters.particle_levelset_evolution;T_LEVELSET& levelset=particle_levelset_evolution.particle_levelset.levelset;
    T_LEVELSET_ADVECTION& levelset_advection=particle_levelset_evolution.levelset_advection;
    assert(levelset_advection.nested_semi_lagrangian_collidable);
    levelset_advection.nested_semi_lagrangian_collidable->cell_valid_points_current=levelset_advection.nested_semi_lagrangian_collidable->cell_valid_points_next;
    levelset_advection.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(*fluids_parameters.grid,fluids_parameters.collidable_phi_replacement_value,levelset.phi);
}
//#####################################################################
// Function Revalidate_Fluid_Velocity
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::
Revalidate_Fluid_Velocity()
{    
    assert(!fluids_parameters.incompressible.nested_nested_semi_lagrangian_fire_multiphase_collidable&&fluids_parameters.incompressible.nested_semi_lagrangian_collidable);
    if(fluids_parameters.incompressible.nested_semi_lagrangian_collidable) 
        fluids_parameters.incompressible.nested_semi_lagrangian_collidable->Average_To_Invalidated_Face(*fluids_parameters.grid,fluids_parameters.incompressible.projection.face_velocities);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::
Get_Object_Velocities(const T dt,const T time)
{
    PROJECTION_DYADIC<T_GRID>& projection=fluids_parameters.incompressible.projection;
    fluids_parameters.collision_bodies_affecting_fluid->Compute_Psi_N(projection.elliptic_solver->psi_N,&projection.face_velocities);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::
Get_Object_Velocities(PROJECTION_DYADIC<T_GRID>& projection,ARRAY<T> face_velocities,const T dt,const T time)
{
    fluids_parameters.collision_bodies_affecting_fluid->Compute_Psi_N(projection.elliptic_solver->psi_N,&face_velocities);
}
//#####################################################################
// Function Initialize_Occupied_Cell_For_Advection
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::
Initialize_Swept_Occupied_Blocks_For_Advection(const T dt,const T time,const bool include_removed_particle_velocities)
{
    T_PARTICLE_LEVELSET& particle_levelset=fluids_parameters.particle_levelset_evolution.particle_levelset;
    T maximum_fluid_speed=fluids_parameters.incompressible.projection.face_velocities.Maxabs(),maximum_particle_speed=0;
    if(fluids_parameters.fire){
        Get_Levelset_Velocity(*fluids_parameters.grid,fluids_parameters.particle_levelset_evolution.particle_levelset.levelset,fluids_parameters.particle_levelset_evolution.face_velocities,time);
        maximum_fluid_speed=max(maximum_fluid_speed,fluids_parameters.particle_levelset_evolution.face_velocities.Maxabs());}
    if(include_removed_particle_velocities){
        if(particle_levelset.use_removed_negative_particles) for(int i=1;i<=fluids_parameters.grid->number_of_cells;i++){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles=particle_levelset.removed_negative_particles(i);
            if(particles) maximum_particle_speed=max(maximum_particle_speed,particles->V.Maximum_Magnitude());}
        if(particle_levelset.use_removed_positive_particles) for(int i=1;i<=fluids_parameters.grid->number_of_cells;i++){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles=particle_levelset.removed_positive_particles(i);
                if(particles) maximum_particle_speed=max(maximum_particle_speed,particles->V.Maximum_Magnitude());}}
    LOG::cout<<"maximum_fluid_speed = "<<maximum_fluid_speed<<", maximum_particle_speed = "<<maximum_particle_speed<<std::endl;
    T max_particle_collision_distance=fluids_parameters.particle_levelset_evolution.particle_levelset.max_collision_distance_factor*fluids_parameters.grid->Minimum_Edge_Length();
    fluids_parameters.collision_bodies_affecting_fluid->Compute_Occupied_Blocks(true,fluids_parameters.cfl*dt*max(maximum_fluid_speed,maximum_particle_speed)+
                                                                              2*max_particle_collision_distance+(T).5*fluids_parameters.grid->Minimum_Edge_Length(),10);
//    fluids_parameters.collision_bodies_affecting_fluid->Compute_Occupied_Cells(fluids_parameters.grid,occupied_cell_for_advection,true,
//    2*dt*max(maximum_fluid_speed,maximum_particle_speed)+2*max_collision_distance+fluids_parameters.grid.Minimum_Edge_Length(),10);
}
//#####################################################################
// Function Read_Output_Files_Fluids
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::
Read_Output_Files_Fluids(const int frame)
{
    fluids_parameters.Read_Output_Files(stream_type,output_directory,frame);
/*
    if(fluids_parameters.smoke||fluids_parameters.fire||fluids_parameters.water){
        if(fluids_parameters.use_solid_fluid_coupling){std::string filename;
            if(fluids_parameters.smoke && fluids_parameters.use_density && fluids_parameters.semi_lagrangian_collidable_density.use_valid_mask){
                filename=STRING_UTILITIES::string_sprintf("%s/%d/density_valid_mask",output_directory,frame);
                if(FILE_UTILITIES::File_Exists(filename)){
                    LOG::cout << "Reading " << filename << std::endl;
                    FILE_UTILITIES::Read_From_File(stream_type,filename,fluids_parameters.semi_lagrangian_collidable_density.valid_points_current);}}
            if(fluids_parameters.smoke && fluids_parameters.use_temperature && fluids_parameters.semi_lagrangian_collidable_temperature.use_valid_mask){
                filename=STRING_UTILITIES::string_sprintf("%s/%d/temperature_valid_mask",output_directory,frame);
                if(FILE_UTILITIES::File_Exists(filename)){
                    LOG::cout << "Reading " << filename << std::endl;
                    FILE_UTILITIES::Read_From_File(stream_type,filename,fluids_parameters.semi_lagrangian_collidable_temperature.valid_points_current);}}
            if(fluids_parameters.semi_lagrangian_collidable_velocity.use_valid_mask){
                filename=STRING_UTILITIES::string_sprintf("%s/%d/velocity_valid_mask",output_directory,frame);
                if(FILE_UTILITIES::File_Exists(filename)){
                    LOG::cout << "Reading " << filename << std::endl;
                    FILE_UTILITIES::Read_From_File(stream_type,filename,fluids_parameters.semi_lagrangian_collidable_velocity.valid_points_current);}}
            if((fluids_parameters.water || fluids_parameters.fire) && fluids_parameters.semi_lagrangian_collidable_phi.use_valid_mask){
                filename=STRING_UTILITIES::string_sprintf("%s/%d/phi_valid_mask",output_directory,frame);
                if(FILE_UTILITIES::File_Exists(filename)){
                    LOG::cout << "Reading " << filename << std::endl;
                    FILE_UTILITIES::Read_From_File(stream_type,filename,fluids_parameters.semi_lagrangian_collidable_phi.valid_points_current);}}}}*/
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::
Write_Output_Files(const int frame) const
{
    FILE_UTILITIES::Create_Directory(output_directory);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame),prefix=output_directory+"/"+T_GRID::name;
    FILE_UTILITIES::Create_Directory(output_directory+"/"+f);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    Write_Frame_Title(frame);
    solids_parameters.Write_Output_Files(stream_type,output_directory,first_frame,frame);
    fluids_parameters.Write_Output_Files(stream_type,output_directory,first_frame,frame);
    if(fluids_parameters.smoke||fluids_parameters.fire||fluids_parameters.water){
        if(fluids_parameters.solid_affects_fluid && fluids_parameters.fluid_affects_solid){
            if(fluids_parameters.write_debug_data){
                /*
                FILE_UTILITIES::Write_To_File(stream_type,prefix+"_thin_shells_grid_visibility."+f,fluids_parameters.collision_bodies_affecting_fluid->node_neighbors_visible,
                    fluids_parameters.collision_bodies_affecting_fluid->face_corners_visible_from_face_center);
                FILE_UTILITIES::Write_To_File(stream_type,prefix+"_node_neighbors_visible."+f,fluids_parameters.collision_bodies_affecting_fluid->node_neighbors_visible);
                if(fluids_parameters.smoke && fluids_parameters.use_density && fluids_parameters.semi_lagrangian_collidable_density.use_valid_mask)
                    FILE_UTILITIES::Write_To_File(stream_type,prefix+"_density_valid_mask."+f,fluids_parameters.semi_lagrangian_collidable_density.valid_points_current);
                if(fluids_parameters.smoke && fluids_parameters.use_temperature && fluids_parameters.semi_lagrangian_collidable_temperature.use_valid_mask)
                    FILE_UTILITIES::Write_To_File(stream_type,prefix+"_temperature_valid_mask."+f,fluids_parameters.semi_lagrangian_collidable_temperature.valid_points_current);
                if(fluids_parameters.semi_lagrangian_collidable_velocity.use_valid_mask)
                    FILE_UTILITIES::Write_To_File(stream_type,prefix+"_velocity_valid_mask."+f,fluids_parameters.semi_lagrangian_collidable_velocity.valid_points_current);
                if((fluids_parameters.water || fluids_parameters.fire) && fluids_parameters.semi_lagrangian_collidable_phi.use_valid_mask)
                    FILE_UTILITIES::Write_To_File(stream_type,prefix+"_phi_valid_mask."+f,fluids_parameters.semi_lagrangian_collidable_phi.valid_points_current);
                if((fluids_parameters.water || fluids_parameters.fire) && fluids_parameters.semi_lagrangian_collidable_phi.use_valid_mask)
                    FILE_UTILITIES::Write_To_File(stream_type,prefix+"_phi_valid_next_mask."+f,fluids_parameters.semi_lagrangian_collidable_phi.valid_points_next);*/}}}
}
//#####################################################################
#define COMMA ,
#define P(...) __VA_ARGS__
#define INSTANTIATION_HELPER_T_d_SHAPE(T,d,T_GRID,SHAPE) \
    template void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::Get_Source_Velocities(const SHAPE&,const MATRIX<T,d+1>&,const VECTOR<T,d>&); \
    template void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::Adjust_Phi_With_Source(const SHAPE&,const MATRIX<T,d+1>&); \
    template void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::Get_Source_Reseed_Mask(const SHAPE& source,const MATRIX<T,d+1>&,REBIND<GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR,bool>::TYPE*&,const bool); \
    template void SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>::Adjust_Density_And_Temperature_With_Sources(const SHAPE&,const MATRIX<T,d+1>&,const T,const T);
#define INSTANTIATION_HELPER_SHAPE(SHAPE) INSTANTIATION_HELPER_T_d_SHAPE(P(SHAPE)::VECTOR_T::SCALAR,P(SHAPE)::VECTOR_T::m,DYADIC_GRID_POLICY<P(SHAPE)::VECTOR_T>::DYADIC_GRID,P(SHAPE))
#define INSTANTIATION_HELPER(T) \
    template class SOLIDS_FLUIDS_EXAMPLE_DYADIC<OCTREE_GRID<T> >; \
    template class SOLIDS_FLUIDS_EXAMPLE_DYADIC<QUADTREE_GRID<T> >; \
    INSTANTIATION_HELPER_SHAPE(P(BOX<VECTOR<T COMMA 2> >)) \
    INSTANTIATION_HELPER_SHAPE(P(BOX<VECTOR<T COMMA 3> >)) \
    INSTANTIATION_HELPER_SHAPE(P(SPHERE<VECTOR<T,2> >)) \
    INSTANTIATION_HELPER_SHAPE(P(SPHERE<VECTOR<T,3> >)) \
    INSTANTIATION_HELPER_SHAPE(CYLINDER<T>)
INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif
#endif

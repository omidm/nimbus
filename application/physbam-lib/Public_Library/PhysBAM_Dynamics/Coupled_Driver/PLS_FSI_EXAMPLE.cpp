//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_MATRIX_1X1.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_OCTREE.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_QUADTREE.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/LEVELSET_RLE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_1D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_3D.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Uniform/READ_WRITE_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Fluids/Coupled_Evolution/COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/FLUID_COLLISION_BODY_INACCURATE_UNION.h>
#include <PhysBAM_Dynamics/Coupled_Driver/PLS_FSI_EXAMPLE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
#include <PhysBAM_Dynamics/Forces_And_Torques/EULER_FLUID_FORCES.h>
#include <PhysBAM_Dynamics/Geometry/GENERAL_GEOMETRY_FORWARD.h>
#include <PhysBAM_Dynamics/Incompressible_Flows/INCOMPRESSIBLE_MULTIPHASE_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Particles/PARTICLES_FORWARD.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_PARAMETERS.h>
#include <stdexcept>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV_PFE> PLS_FSI_EXAMPLE<TV_PFE>::
PLS_FSI_EXAMPLE(const STREAM_TYPE stream_type,const int number_of_regions)
    :BASE((Initialize_Particles(),stream_type)),solids_parameters(*new SOLIDS_PARAMETERS<TV_PFE>),solids_fluids_parameters(*new SOLIDS_FLUIDS_PARAMETERS<TV_PFE>(this)),
    solid_body_collection(*new SOLID_BODY_COLLECTION<TV_PFE>(this,0)),solids_evolution(new NEWMARK_EVOLUTION<TV_PFE>(solids_parameters,solid_body_collection)),
    fluids_parameters(number_of_regions,fluids_parameters.WATER,incompressible_fluid_container),resolution(0),convection_order(1),use_pls_evolution_for_structure(false),
    two_phase(false)
{
    Initialize_Read_Write_General_Structures();
    Set_Minimum_Collision_Thickness();
    Set_Write_Substeps_Level(-1);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV_PFE> PLS_FSI_EXAMPLE<TV_PFE>::
~PLS_FSI_EXAMPLE()
{
    fluids_parameters.projection=0;
    delete solids_evolution;
    delete &solid_body_collection;
    delete &solids_parameters;
    delete &solids_fluids_parameters;
}
//#####################################################################
// Function Register_Options
//#####################################################################
template<class TV_PFE> void PLS_FSI_EXAMPLE<TV_PFE>::
Register_Options()
{
    if(!parse_args) return;
    BASE::Register_Options();
    parse_args->Add_String_Argument("-params","","parameter file");
    parse_args->Add_Double_Argument("-solidscfl",.9,"solids CFL");
    parse_args->Add_Option_Argument("-solidscg","Use CG for time integration");
    parse_args->Add_Option_Argument("-solidscr","Use CONJUGATE_RESIDUAL for time integration");
    parse_args->Add_Option_Argument("-solidssymmqmr","Use SYMMQMR for time integration");
    parse_args->Add_Double_Argument("-rigidcfl",.5,"rigid CFL");
    parse_args->Add_Option_Argument("-skip_debug_data","turn off file io for debug data");
    parse_args->Add_Integer_Argument("-resolution",1);
}
//#####################################################################
// Function Parse_Options
//#####################################################################
template<class TV_PFE> void PLS_FSI_EXAMPLE<TV_PFE>::
Parse_Options()
{
    BASE::Parse_Options();
    if(parse_args->Is_Value_Set("-skip_debug_data")) fluids_parameters.write_debug_data=false;
    resolution=parse_args->Get_Integer_Value("-resolution");
}
//#####################################################################
// Function Add_Volumetric_Body_To_Fluid_Simulation
//#####################################################################
template<class TV_PFE> void PLS_FSI_EXAMPLE<TV_PFE>::
Add_Volumetric_Body_To_Fluid_Simulation(RIGID_BODY<TV_PFE>& rigid_body,bool add_collision,bool add_coupling)
{
    RIGID_COLLISION_GEOMETRY<TV_PFE>* collision_geometry=new RIGID_COLLISION_GEOMETRY<TV_PFE>(rigid_body);
    if(add_collision) fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(collision_geometry,rigid_body.particle_index,true);
    if(add_coupling)
        if(SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV_PFE>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV_PFE>*>(solids_evolution))
            coupled_evolution->iterator_info.coupling_bodies.Append(collision_geometry);
    if(rigid_body.simplicial_object) rigid_body.simplicial_object->Initialize_Hierarchy();
}
//#####################################################################
// Function Add_Thin_Shell_To_Fluid_Simulation
//#####################################################################
template<class TV_PFE> void PLS_FSI_EXAMPLE<TV_PFE>::
Add_Thin_Shell_To_Fluid_Simulation(RIGID_BODY<TV_PFE>& rigid_body,bool add_collision,bool add_coupling)
{
    RIGID_COLLISION_GEOMETRY<TV_PFE>* collision_geometry=new RIGID_COLLISION_GEOMETRY<TV_PFE>(rigid_body);
    if(add_collision) fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(collision_geometry,rigid_body.particle_index,true);
    if(add_coupling)
        if(SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV_PFE>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV_PFE>*>(solids_evolution))
            coupled_evolution->iterator_info.coupling_bodies.Append(collision_geometry);
    rigid_body.thin_shell=true;
    if(collision_geometry->collision_thickness<minimum_collision_thickness) collision_geometry->Set_Collision_Thickness(minimum_collision_thickness);
    if(rigid_body.simplicial_object) rigid_body.simplicial_object->Initialize_Hierarchy();
}
//#####################################################################
// Function Add_To_Fluid_Simulation
//#####################################################################
template<class TV_PFE> void PLS_FSI_EXAMPLE<TV_PFE>::
Add_To_Fluid_Simulation(DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV_PFE>& deformable_collisions,bool add_collision,bool add_coupling)
{
    if(add_collision) fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(&deformable_collisions,0,false);
    if(add_coupling)
        if(SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV_PFE>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV_PFE>*>(solids_evolution))
            coupled_evolution->iterator_info.coupling_bodies.Append(&deformable_collisions);
    if(deformable_collisions.collision_thickness<minimum_collision_thickness) deformable_collisions.Set_Collision_Thickness(minimum_collision_thickness);
    deformable_collisions.Initialize_For_Thin_Shells_Fluid_Coupling();
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
template<class TV_PFE> template<class GEOMETRY> void PLS_FSI_EXAMPLE<TV_PFE>::
Get_Source_Velocities(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const TV_PFE& constant_source_velocity)
{
    for(UNIFORM_GRID_ITERATOR_FACE<TV_PFE> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) if(source.Lazy_Inside(world_to_source.Homogeneous_Times(iterator.Location()))){
        int axis=iterator.Axis();fluids_parameters.incompressible->projection.elliptic_solver->psi_N.Component(axis)(iterator.Face_Index())=true;
        incompressible_fluid_container.face_velocities.Component(axis)(iterator.Face_Index())=constant_source_velocity[axis];}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
template<class TV_PFE> template<class GEOMETRY> void PLS_FSI_EXAMPLE<TV_PFE>::
Get_Source_Velocities(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const TV_PFE& constant_source_velocity,const ARRAY<bool,FACE_INDEX<TV_PFE::m> >& invalid_mask)
{
    for(UNIFORM_GRID_ITERATOR_FACE<TV_PFE> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){const int axis=iterator.Axis();const TV_INT face_index=iterator.Face_Index();
        if(!invalid_mask(axis,face_index) && source.Lazy_Inside(world_to_source.Homogeneous_Times(iterator.Location()))){
            fluids_parameters.incompressible->projection.elliptic_solver->psi_N.Component(axis)(face_index)=true;
            incompressible_fluid_container.face_velocities.Component(axis)(face_index)=constant_source_velocity[axis];}}
}
//#####################################################################
// Function Adjust_Phi_With_Source
//#####################################################################
template<class TV_PFE> template<class GEOMETRY> void PLS_FSI_EXAMPLE<TV_PFE>::
Adjust_Phi_With_Source(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source)
{
    for(UNIFORM_GRID_ITERATOR_CELL<TV_PFE> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){TV_PFE source_X=world_to_source.Homogeneous_Times(iterator.Location());
        if(source.Lazy_Inside(source_X)) 
            fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index())=min(fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index()),
                                                                                          source.Signed_Distance(source_X));}
}
//#####################################################################
// Function Adjust_Phi_With_Source
//#####################################################################
template<class TV_PFE> template<class GEOMETRY> void PLS_FSI_EXAMPLE<TV_PFE>::
Adjust_Phi_With_Source(const GEOMETRY& source,const int region,const T_TRANSFORMATION_MATRIX& world_to_source)
{
    T bandwidth=3*fluids_parameters.grid->Minimum_Edge_Length();
    ARRAY<ARRAY<T,TV_INT> >& phis=fluids_parameters.particle_levelset_evolution_multiple->phis;
    for(UNIFORM_GRID_ITERATOR_CELL<TV_PFE> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
        TV_PFE source_X=world_to_source.Homogeneous_Times(iterator.Location());
        if(source.Inside(source_X,-bandwidth)){
            T source_signed_distance=source.Signed_Distance(source_X);
            for(int i=1;i<=fluids_parameters.number_of_regions;i++){
                if(i==region) phis(i)(iterator.Cell_Index())=min(phis(i)(iterator.Cell_Index()),source_signed_distance);
                else phis(i)(iterator.Cell_Index())=max(phis(i)(iterator.Cell_Index()),-source_signed_distance);}}}
}
//#####################################################################
// Function Revalidate_Fluid_Scalars
//#####################################################################
template<class TV_PFE> void PLS_FSI_EXAMPLE<TV_PFE>::
Revalidate_Fluid_Scalars()
{
    T_FAST_LEVELSET& levelset=fluids_parameters.particle_levelset_evolution->Levelset(1);
    T_FAST_LEVELSET_ADVECTION& levelset_advection=fluids_parameters.particle_levelset_evolution->Levelset_Advection(1);
    if(levelset_advection.nested_semi_lagrangian_collidable)
        levelset_advection.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(*fluids_parameters.grid,fluids_parameters.collidable_phi_replacement_value,levelset.phi);
}
//#####################################################################
// Function Revalidate_Phi_After_Modify_Levelset
//#####################################################################
template<class TV_PFE> void PLS_FSI_EXAMPLE<TV_PFE>::
Revalidate_Phi_After_Modify_Levelset()
{
    T_FAST_LEVELSET& levelset=fluids_parameters.particle_levelset_evolution->Levelset(1);
    T_FAST_LEVELSET_ADVECTION& levelset_advection=fluids_parameters.particle_levelset_evolution->Levelset_Advection(1);
    if(levelset_advection.nested_semi_lagrangian_collidable){
        levelset_advection.nested_semi_lagrangian_collidable->cell_valid_points_current=levelset_advection.nested_semi_lagrangian_collidable->cell_valid_points_next;
        levelset_advection.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(*fluids_parameters.grid,fluids_parameters.collidable_phi_replacement_value,levelset.phi);}
}
//#####################################################################
// Function Revalidate_Fluid_Velocity
//#####################################################################
template<class TV_PFE> void PLS_FSI_EXAMPLE<TV_PFE>::
Revalidate_Fluid_Velocity(ARRAY<T,FACE_INDEX<TV_PFE::m> >& face_velocities)
{
    if(fluids_parameters.incompressible->nested_semi_lagrangian_collidable) 
        fluids_parameters.incompressible->nested_semi_lagrangian_collidable->Average_To_Invalidated_Face(*fluids_parameters.grid,face_velocities);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
template<class TV_PFE> void PLS_FSI_EXAMPLE<TV_PFE>::
Get_Object_Velocities(LAPLACE_UNIFORM<GRID<TV_PFE> >* elliptic_solver,ARRAY<T,FACE_INDEX<TV_PFE::m> >& face_velocities,const T dt,const T time)
{
    if(fluids_parameters.solid_affects_fluid && (fluids_parameters.fluid_affects_solid || fluids_parameters.use_slip)){
        if(!fluids_parameters.use_slip){
           SOLID_FLUID_COUPLED_EVOLUTION<TV_PFE>& coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION<TV_PFE>&>(*solids_evolution);
           coupled_evolution.Apply_Solid_Boundary_Conditions(time,false,face_velocities);}
        else{
            SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV_PFE>& coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV_PFE>&>(*solids_evolution);
            coupled_evolution.Get_Coupled_Faces_And_Interpolated_Solid_Velocities(time,elliptic_solver->psi_N,face_velocities);}}
    else fluids_parameters.collision_bodies_affecting_fluid->Compute_Psi_N(elliptic_solver->psi_N,&face_velocities);
}
//#####################################################################
// Function Get_Levelset_Velocity
//#####################################################################
template<class TV_PFE> void PLS_FSI_EXAMPLE<TV_PFE>::
Get_Levelset_Velocity(const GRID<TV_PFE>& grid,LEVELSET_MULTIPLE_UNIFORM<GRID<TV_PFE> >& levelset_multiple,ARRAY<T,FACE_INDEX<TV_PFE::m> >& V_levelset,const T time) const
{
    ARRAY<T,FACE_INDEX<TV_PFE::m> >::Put(incompressible_fluid_container.face_velocities,V_levelset);
}
//#####################################################################
// Function Initialize_Swept_Occupied_Blocks_For_Advection
//#####################################################################
template<class TV_PFE> void PLS_FSI_EXAMPLE<TV_PFE>::
Initialize_Swept_Occupied_Blocks_For_Advection(const T dt,const T time,const ARRAY<T,FACE_INDEX<TV_PFE::m> >& face_velocities)
{
    GRID<TV_PFE>& grid=*fluids_parameters.grid;
    T maximum_fluid_speed=face_velocities.Maxabs().Max(),maximum_particle_speed=0;
    PARTICLE_LEVELSET_UNIFORM<GRID<TV_PFE> >& particle_levelset=fluids_parameters.particle_levelset_evolution->Particle_Levelset(1);
    if(particle_levelset.use_removed_negative_particles) for(UNIFORM_GRID_ITERATOR_CELL<TV_PFE> iterator(particle_levelset.levelset.grid);iterator.Valid();iterator.Next()){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV_PFE>* particles=particle_levelset.removed_negative_particles(iterator.Cell_Index());
            if(particles) maximum_particle_speed=max(maximum_particle_speed,ARRAYS_COMPUTATIONS::Maximum_Magnitude(particles->V));}
    if(particle_levelset.use_removed_positive_particles) for(UNIFORM_GRID_ITERATOR_CELL<TV_PFE> iterator(particle_levelset.levelset.grid);iterator.Valid();iterator.Next()){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV_PFE>* particles=particle_levelset.removed_positive_particles(iterator.Cell_Index());
            if(particles) maximum_particle_speed=max(maximum_particle_speed,ARRAYS_COMPUTATIONS::Maximum_Magnitude(particles->V));}
    T max_particle_collision_distance=0;
    max_particle_collision_distance=max(max_particle_collision_distance,fluids_parameters.particle_levelset_evolution->Particle_Levelset(1).max_collision_distance_factor*grid.dX.Max());
    fluids_parameters.collision_bodies_affecting_fluid->Compute_Occupied_Blocks(true,
        dt*max(maximum_fluid_speed,maximum_particle_speed)+2*max_particle_collision_distance+(T).5*fluids_parameters.p_grid.dX.Max(),10);
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
template<class TV_PFE> void PLS_FSI_EXAMPLE<TV_PFE>::
Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<TV_PFE>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time)
{}
//#####################################################################
// Function Read_Output_Files_Fluids
//#####################################################################
template<class TV_PFE> void PLS_FSI_EXAMPLE<TV_PFE>::
Read_Output_Files_Fluids(const int frame)
{
    fluids_parameters.Read_Output_Files(stream_type,output_directory,frame);
    incompressible_fluid_container.Read_Output_Files(stream_type,output_directory,frame);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV_PFE> void PLS_FSI_EXAMPLE<TV_PFE>::
Write_Output_Files(const int frame) const
{
    FILE_UTILITIES::Create_Directory(output_directory);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Create_Directory(output_directory+"/"+f);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    Write_Frame_Title(frame);
    solid_body_collection.Write(stream_type,output_directory,frame,first_frame,solids_parameters.write_static_variables_every_frame,
        solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies,solids_parameters.write_deformable_body,solids_parameters.write_from_every_process,
        solids_parameters.triangle_collision_parameters.output_interaction_pairs);
    if(NEWMARK_EVOLUTION<TV_PFE>* newmark=dynamic_cast<NEWMARK_EVOLUTION<TV_PFE>*>(solids_evolution))
        newmark->Write_Position_Update_Projection_Data(stream_type,output_directory+"/"+f+"/");
    
    fluids_parameters.Write_Output_Files(stream_type,output_directory,first_frame,frame);
    incompressible_fluid_container.Write_Output_Files(stream_type,output_directory,frame);
}
//#####################################################################
// Function Log_Parameters 
//#####################################################################
template<class TV_PFE> void PLS_FSI_EXAMPLE<TV_PFE>::
Log_Parameters() const
{       
    LOG::SCOPE scope("PLS_FSI_EXAMPLE parameters");
    BASE::Log_Parameters();
    std::stringstream ss;
    ss<<"minimum_collision_thickness="<<minimum_collision_thickness<<std::endl;
    LOG::filecout(ss.str());
    fluids_parameters.Log_Parameters();
}
//#####################################################################
// Function Read_Output_Files_Solids
//#####################################################################
template<class TV_PFE> void PLS_FSI_EXAMPLE<TV_PFE>::
Read_Output_Files_Solids(const int frame)
{
    solid_body_collection.Read(stream_type,output_directory,frame,frame,solids_parameters.write_static_variables_every_frame,solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies,
        solids_parameters.write_deformable_body,solids_parameters.write_from_every_process);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    //if(NEWMARK_EVOLUTION<TV_PFE>* newmark=dynamic_cast<NEWMARK_EVOLUTION<TV_PFE>*>(solids_evolution))
    //    newmark->Read_Position_Update_Projection_Data(stream_type,output_directory+"/"+f+"/");
}
//#####################################################################
// Function Parse_Late_Options
//#####################################################################
template<class TV_PFE> void PLS_FSI_EXAMPLE<TV_PFE>::
Parse_Late_Options()
{
    if(!parse_args) return;
    BASE::Parse_Late_Options();
    if(parse_args->Is_Value_Set("-solidscfl")) solids_parameters.cfl=(T)parse_args->Get_Double_Value("-solidscfl");
    if(parse_args->Is_Value_Set("-rigidcfl")) solids_parameters.rigid_body_evolution_parameters.rigid_cfl=(T)parse_args->Get_Double_Value("-rigidcfl");
    if(parse_args->Is_Value_Set("-solidscg")) solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cg;
    if(parse_args->Is_Value_Set("-solidscr")) solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cr;
    if(parse_args->Is_Value_Set("-solidssymmqmr")) solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_symmqmr;
}
//#####################################################################
// Function Adjust_Particle_For_Domain_Boundaries
//#####################################################################
template<class TV_PFE> void PLS_FSI_EXAMPLE<TV_PFE>::
Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV_PFE>& particles,const int index,TV_PFE& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time)
{
    // remove this for speed - don't call the function for the other particles
    if(particle_type==PARTICLE_LEVELSET_POSITIVE || particle_type==PARTICLE_LEVELSET_REMOVED_POSITIVE) return;

    TV_PFE& X=particles.X(index);TV_PFE X_new=X+dt*V;
    T max_collision_distance=fluids_parameters.particle_levelset_evolution->particle_levelset.Particle_Collision_Distance(particles.quantized_collision_distance(index));
    T min_collision_distance=fluids_parameters.particle_levelset_evolution->particle_levelset.min_collision_distance_factor*max_collision_distance;
    TV_PFE min_corner=fluids_parameters.grid->domain.Minimum_Corner(),max_corner=fluids_parameters.grid->domain.Maximum_Corner();
    for(int axis=1;axis<=TV_PFE::m;axis++){
        if(fluids_parameters.domain_walls[axis][1] && X_new[axis]<min_corner[axis]+max_collision_distance){
            T collision_distance=X[axis]-min_corner[axis];
            if(collision_distance>max_collision_distance)collision_distance=X_new[axis]-min_corner[axis];
            collision_distance=max(min_collision_distance,collision_distance);
            X_new[axis]+=max((T)0,min_corner[axis]-X_new[axis]+collision_distance);
            V[axis]=max((T)0,V[axis]);X=X_new-dt*V;}
        if(fluids_parameters.domain_walls[axis][2] && X_new[axis]>max_corner[axis]-max_collision_distance){
            T collision_distance=max_corner[axis]-X[axis];
            if(collision_distance>max_collision_distance) collision_distance=max_corner[axis]-X_new[axis];
            collision_distance=max(min_collision_distance,collision_distance);
            X_new[axis]-=max((T)0,X_new[axis]-max_corner[axis]+collision_distance);
            V[axis]=min((T)0,V[axis]);X=X_new-dt*V;}}
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
template<class TV_PFE> void PLS_FSI_EXAMPLE<TV_PFE>::
Update_Fluid_Parameters(const T dt,const T time)
{
    fluids_parameters.incompressible->Set_Gravity(fluids_parameters.gravity,fluids_parameters.gravity_direction);
    fluids_parameters.incompressible->Set_Body_Force(fluids_parameters.use_body_force);
    fluids_parameters.incompressible->projection.Use_Non_Zero_Divergence(fluids_parameters.use_non_zero_divergence);
    fluids_parameters.incompressible->projection.elliptic_solver->Solve_Neumann_Regions(fluids_parameters.solve_neumann_regions);
    fluids_parameters.incompressible->projection.elliptic_solver->solve_single_cell_neumann_regions=fluids_parameters.solve_single_cell_neumann_regions;
    fluids_parameters.incompressible->Use_Explicit_Part_Of_Implicit_Viscosity(fluids_parameters.use_explicit_part_of_implicit_viscosity);
    if(fluids_parameters.implicit_viscosity && fluids_parameters.implicit_viscosity_iterations)
        fluids_parameters.incompressible->Set_Maximum_Implicit_Viscosity_Iterations(fluids_parameters.implicit_viscosity_iterations);
    fluids_parameters.incompressible->Use_Variable_Vorticity_Confinement(fluids_parameters.use_variable_vorticity_confinement);
    fluids_parameters.incompressible->Set_Surface_Tension(fluids_parameters.surface_tension);
    fluids_parameters.incompressible->Set_Variable_Surface_Tension(fluids_parameters.variable_surface_tension);
    fluids_parameters.incompressible->Set_Viscosity(fluids_parameters.viscosity);
    fluids_parameters.incompressible->Set_Variable_Viscosity(fluids_parameters.variable_viscosity);
    fluids_parameters.incompressible->projection.Set_Density(fluids_parameters.density);
}
//#####################################################################
template class PLS_FSI_EXAMPLE<VECTOR<float,1> >;
template class PLS_FSI_EXAMPLE<VECTOR<float,2> >;
template class PLS_FSI_EXAMPLE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PLS_FSI_EXAMPLE<VECTOR<double,1> >;
template class PLS_FSI_EXAMPLE<VECTOR<double,2> >;
template class PLS_FSI_EXAMPLE<VECTOR<double,3> >;
#endif

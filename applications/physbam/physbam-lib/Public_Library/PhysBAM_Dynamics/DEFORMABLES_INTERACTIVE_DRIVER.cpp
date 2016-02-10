//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLES_INTERACTIVE_DRIVER
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLES_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformables_Evolution/DEFORMABLES_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Dynamics/DEFORMABLES_INTERACTIVE_DRIVER.h>
#ifdef WIN32
    #include <windows.h>
#endif
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLES_INTERACTIVE_DRIVER<TV>::
DEFORMABLES_INTERACTIVE_DRIVER(DEFORMABLES_EXAMPLE<TV>& example_input)
    :BASE(example_input),example(example_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DEFORMABLES_INTERACTIVE_DRIVER<TV>::
~DEFORMABLES_INTERACTIVE_DRIVER()
{}
//#####################################################################
// Function Advance_To_Target_Time
//#####################################################################
template<class TV> void DEFORMABLES_INTERACTIVE_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Advance_To_Target_Time Start",2,2);
    DEFORMABLES_PARAMETERS<TV>& deformables_parameters=example.deformables_parameters;
    TRIANGLE_COLLISION_PARAMETERS<TV>& triangle_collision_parameters=deformables_parameters.triangle_collision_parameters;
    DEFORMABLE_OBJECT_COLLISION_PARAMETERS<TV>& deformable_object_collision_parameters=deformables_parameters.deformable_object_collision_parameters;
    DEFORMABLES_EVOLUTION_CALLBACKS<TV>* deformables_evolution_callbacks=example.deformables_evolution->deformables_evolution_callbacks; //TODO: Fix this    

    triangle_collision_parameters.steps_since_self_collision_free=0;
    bool done=false;for(int substep=1;!done;substep++){
        std::stringstream ss;ss<<"substep "<<substep;
        LOG::SCOPE scope("SUBSTEP",ss.str());

        //setup
        if(triangle_collision_parameters.perform_self_collision && triangle_collision_parameters.temporary_enable_collisions){
            deformables_evolution_callbacks->Self_Collisions_Begin_Callback(time,substep);
            triangle_collision_parameters.repulsion_pair_update_count=0;
            example.deformable_body_collection.triangle_repulsions_and_collisions_geometry.Save_Self_Collision_Free_State();
            if((triangle_collision_parameters.topological_hierarchy_build_count++)%triangle_collision_parameters.topological_hierarchy_build_frequency==0){
                example.deformable_body_collection.triangle_repulsions_and_collisions_geometry.Build_Topological_Structure_Of_Hierarchies();}
            triangle_collision_parameters.self_collision_free_time=time;}
        deformables_evolution_callbacks->Preprocess_Solids_Substep(time,substep);
        if(deformable_object_collision_parameters.use_spatial_partition_for_levelset_collision_objects) // TODO - ANDY - why is this needed; move this to the right places inside solids evolution 
            example.collision_body_list.Update_Spatial_Partition(deformable_object_collision_parameters.spatial_partition_voxel_size_heuristic,
                deformable_object_collision_parameters.spatial_partition_number_of_cells,deformable_object_collision_parameters.spatial_partition_voxel_size_scale_factor);

        // solids dt
        T solids_dt=FLT_MAX;
        //solids_dt=min(solids_dt,example.deformable_body_collection.CFL(deformables_parameters.verbose_dt));
        //solids_dt=min(solids_dt,deformables_evolution_callbacks->Constraints_CFL());
        deformables_evolution_callbacks->Limit_Solids_Dt(solids_dt,time);
        if(example.deformable_body_collection.mpi_solids)
            solids_dt=example.deformable_body_collection.mpi_solids->Reduce_Min_Global(solids_dt);

        if(example.fixed_dt) solids_dt=example.fixed_dt;
        T dt=solids_dt;
        EXAMPLE<TV>::Clamp_Time_Step_With_Target_Time(time,target_time,dt,done,deformables_parameters.min_dt);

        //step
        example.Preprocess_Substep(dt,time);

        PHYSBAM_DEBUG_WRITE_SUBSTEP("Advance_To_Target_Time 1",2,2);
        if((triangle_collision_parameters.repulsion_pair_update_count++)%triangle_collision_parameters.repulsion_pair_update_frequency==0){
            example.deformable_body_collection.triangle_repulsions.Update_Faces_And_Hierarchies_With_Collision_Free_Positions(&example.deformable_body_collection.particles.X);
            example.deformable_body_collection.triangle_repulsions.Compute_Interaction_Pairs(example.deformable_body_collection.particles.X);}
        PHYSBAM_DEBUG_WRITE_SUBSTEP("Advance_To_Target_Time 2",2,2);
        example.deformables_evolution->kinematic_evolution.Set_External_Positions(example.rigid_geometry_collection.particles.X,example.rigid_geometry_collection.particles.rotation,time);
        example.deformables_evolution->kinematic_evolution.Set_External_Velocities(example.rigid_geometry_collection.particles.V,example.rigid_geometry_collection.particles.angular_velocity,time,time);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("Advance_To_Target_Time 3",2,2);
        example.deformables_evolution->Advance_One_Time_Step_Position(dt,time);

        if(triangle_collision_parameters.perform_self_collision && triangle_collision_parameters.temporary_enable_collisions){
            LOG::SCOPE scope("adjust velocity for self repulsion and self collisions");
            int repulsions,collisions_found;
            example.deformables_evolution->Adjust_Velocity_For_Self_Repulsion_And_Self_Collisions(dt,time,repulsions,collisions_found,false);
            triangle_collision_parameters.steps_since_self_collision_free=0;}

        example.deformables_evolution->Advance_One_Time_Step_Velocity(dt,time);
    
        example.deformables_evolution->time+=dt;
        deformables_evolution_callbacks->Postprocess_Solids_Substep(example.deformables_evolution->time,substep);
        deformables_evolution_callbacks->Apply_Constraints(dt,example.deformables_evolution->time);

        example.Postprocess_Substep(dt,time);
        time+=dt;}
}
//#####################################################################
// Function Execute_Main_Program
//#####################################################################
template<class TV> void DEFORMABLES_INTERACTIVE_DRIVER<TV>::
Execute_Main_Program()
{
    {LOG::SCOPE scope("INITIALIZING","Initializing");
    Initialize();
    example.Initialize_Bodies();
    Update_Structures();
    example.Post_Initialization();
    example.Log_Parameters();
    if(!example.restart) Write_Output_Files(example.first_frame);}
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void DEFORMABLES_INTERACTIVE_DRIVER<TV>::
Initialize()
{
    if(example.restart){current_frame=example.restart_frame;Read_Time(current_frame);}else current_frame=example.first_frame;
    output_number=current_frame;
    time=example.Time_At_Frame(current_frame);

    example.deformables_evolution->Set_Deformables_Evolution_Callbacks(example);
    example.deformables_evolution->time=time;
    since_last_frame=std::clock();
}
//#####################################################################
// Function Update_Structures
//#####################################################################
template<class TV> void DEFORMABLES_INTERACTIVE_DRIVER<TV>::
Update_Structures()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=example.deformable_body_collection;
    DEFORMABLE_OBJECT_COLLISION_PARAMETERS<TV>& deformable_object_collision_parameters=example.deformables_parameters.deformable_object_collision_parameters;
    
    if(example.restart){
        LOG::SCOPE scope("reading solids data");
        example.Read_Output_Files_Solids(example.restart_frame);
        example.deformables_evolution->time=time=example.Time_At_Frame(example.restart_frame);}
    
    deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Positions();
    deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities();

    deformable_body_collection.Update_Simulated_Particles();
    example.rigid_geometry_collection.Update_Kinematic_Particles();

    deformable_body_collection.collisions.Initialize_Object_Collisions(deformable_object_collision_parameters.collide_with_interior,deformable_object_collision_parameters.collision_tolerance,
        deformable_object_collision_parameters.use_spatial_partition_for_levelset_collision_objects,deformable_object_collision_parameters.disable_multiple_levelset_collisions,
        deformable_object_collision_parameters.maximum_levelset_collision_projection_velocity);
    deformable_body_collection.collisions.collision_body_list.Update_Spatial_Partition(deformable_object_collision_parameters.spatial_partition_voxel_size_heuristic,
        deformable_object_collision_parameters.spatial_partition_number_of_cells,deformable_object_collision_parameters.spatial_partition_voxel_size_scale_factor);

    deformable_body_collection.Set_CFL_Number(example.deformables_parameters.cfl);
}
//#####################################################################
// Function Update_Positions_From_Interactive
//#####################################################################
template<class TV> void DEFORMABLES_INTERACTIVE_DRIVER<TV>::Update_Positions_From_Interactive()
{ 
    if(visualizer->current_selection && visualizer->current_selection->type==OPENGL_SELECTION::COMPONENT_RIGID_BODIES_3D){
        OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T> *real_selection=dynamic_cast<OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T>*>(visualizer->current_selection);
        example.rigid_geometry_collection.particles.X(real_selection->body_id)[1]=visualizer->opengl_world.object_location(1);
        example.rigid_geometry_collection.particles.X(real_selection->body_id)[2]=visualizer->opengl_world.object_location(2);
        example.rigid_geometry_collection.particles.X(real_selection->body_id)[3]=visualizer->opengl_world.object_location(3);
        rigid_component->Reinitialize(true,false);}
}
//#####################################################################
// Function Simulate_To_Frame
//#####################################################################
template<class TV> void DEFORMABLES_INTERACTIVE_DRIVER<TV>::
Simulate_To_Frame(const int frame_input)
{
    while(frame_input<0 || current_frame<frame_input){
        std::stringstream ss;ss<<"Frame "<<current_frame+1;
        LOG::SCOPE scope("FRAME",ss.str());
        Preprocess_Frame(current_frame+1);
        if(visualizer->opengl_world.drag_current_selection) Update_Positions_From_Interactive();
        example.deformables_evolution->kinematic_evolution.Get_Current_Kinematic_Keyframes(example.Time_At_Frame(current_frame+1)-time,time);
        Advance_To_Target_Time(example.Time_At_Frame(current_frame+1));
        if(visualizer){
            now=std::clock();
            if((now-since_last_frame)<(1./example.frame_rate)*CLOCKS_PER_SEC){
#ifdef WIN32
                Sleep((DWORD)(1000*(1./example.frame_rate)-(now-since_last_frame)/CLOCKS_PER_SEC));
#else
                sleep((1./example.frame_rate)-(now-since_last_frame)/CLOCKS_PER_SEC);
#endif 
                now=std::clock();}
            since_last_frame=now;
            visualizer->Set_Frame(current_frame);}
        Postprocess_Frame(++current_frame);
        if(example.write_output_files && example.write_substeps_level==-1) Write_Output_Files(current_frame);
        else if(example.write_substeps_level!=-1) Write_Substep(STRING_UTILITIES::string_sprintf("END Frame %d",current_frame),0,example.write_substeps_level);
        LOG::cout<<"TIME = "<<time<<std::endl;}
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void DEFORMABLES_INTERACTIVE_DRIVER<TV>::
Write_Output_Files(const int frame)
{
    FILE_UTILITIES::Create_Directory(example.output_directory);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Create_Directory(example.output_directory+"/"+f);
    FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
    Write_First_Frame(frame);
    example.Write_Output_Files(frame);
    Write_Time(frame);
    Write_Last_Frame(frame);
}
//#####################################################################
template class DEFORMABLES_INTERACTIVE_DRIVER<VECTOR<float,1> >;
template class DEFORMABLES_INTERACTIVE_DRIVER<VECTOR<float,2> >;
template class DEFORMABLES_INTERACTIVE_DRIVER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DEFORMABLES_INTERACTIVE_DRIVER<VECTOR<double,1> >;
template class DEFORMABLES_INTERACTIVE_DRIVER<VECTOR<double,2> >;
template class DEFORMABLES_INTERACTIVE_DRIVER<VECTOR<double,3> >;
#endif

//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_DRIVER
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/SOLIDS_DRIVER.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLIDS_DRIVER<TV>::
SOLIDS_DRIVER(SOLIDS_EXAMPLE<TV>& example_input)
    :BASE(example_input),example(example_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLIDS_DRIVER<TV>::
~SOLIDS_DRIVER()
{}
//#####################################################################
// Function Execute_Main_Program
//#####################################################################
template<class TV> void SOLIDS_DRIVER<TV>::
Execute_Main_Program()
{
    {LOG::SCOPE scope("INITIALIZING","Initializing");
    Initialize();
    example.Post_Initialization();
    example.Log_Parameters();
    if(!example.restart) Write_Output_Files(example.first_frame);}
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Function Simulate_To_Frame
//#####################################################################
template<class TV> void SOLIDS_DRIVER<TV>::
Simulate_To_Frame(const int frame_input)
{
    while(current_frame<frame_input){
        LOG::SCOPE scope("FRAME",STRING_UTILITIES::string_sprintf("Frame %d",current_frame+1));
        Preprocess_Frame(current_frame+1);
        example.solids_evolution->kinematic_evolution.Get_Current_Kinematic_Keyframes(example.Time_At_Frame(current_frame+1)-time,time);
        Advance_To_Target_Time(example.Time_At_Frame(current_frame+1));
        Postprocess_Frame(++current_frame);
        if(example.write_output_files && example.write_substeps_level==-1) Write_Output_Files(current_frame);
        else if(example.write_substeps_level!=-1) Write_Substep(STRING_UTILITIES::string_sprintf("END Frame %d",current_frame),0,example.write_substeps_level);
        {std::stringstream ss;ss<<"TIME = "<<time<<std::endl;LOG::filecout(ss.str());}}
}
//#####################################################################
// Function Read_Time
//#####################################################################
template<class TV> void SOLIDS_DRIVER<TV>::
Read_Time(const int frame)
{
    time=example.Time_At_Frame(frame);
    std::string filename=STRING_UTILITIES::string_sprintf("%s/%d/time",example.output_directory.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(filename)){
        T corrected_time;
        FILE_UTILITIES::Read_From_File(example.stream_type,filename,corrected_time);
        if(abs(time-corrected_time)>(T)1e-4*abs(time)){ // only adjust time if significantly different from default in order to get deterministic restarts
            time=corrected_time;
            // adjust initial time so that Simulate_To_Frame() returns correct time (essential when writing substeps)
            example.initial_time=time-(frame-example.first_frame)/example.frame_rate;}}
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void SOLIDS_DRIVER<TV>::
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
// Function Initialize
//#####################################################################
template<class TV> void SOLIDS_DRIVER<TV>::
Initialize()
{
    if(example.auto_restart){
        std::string last_frame_file=example.output_directory+"/common/last_frame";
        int last_frame;FILE_UTILITIES::Read_From_Text_File(last_frame_file,last_frame);
        example.restart=true;example.restart_frame=last_frame;
        {std::stringstream ss;ss<<"Auto Restart from frame "<<last_frame<<" (from file "<<last_frame_file<<")"<<std::endl;LOG::filecout(ss.str());}}
    if(example.restart){current_frame=example.restart_frame;Read_Time(current_frame);}else current_frame=example.first_frame;
    output_number=current_frame;
    time=example.Time_At_Frame(current_frame);

    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    solids_evolution.Set_Solids_Evolution_Callbacks(example);
    example.Initialize_Bodies();

    example.Parse_Late_Options();
    solids_evolution.time=time;

    if(example.restart){
        LOG::SCOPE scope("reading solids data");
        example.Read_Output_Files_Solids(example.restart_frame);
        solids_evolution.time=time=example.Time_At_Frame(example.restart_frame);}

    solids_evolution.Initialize_Deformable_Objects(example.frame_rate,example.restart);

    solids_evolution.Initialize_Rigid_Bodies(example.frame_rate,example.restart);
}
//#####################################################################
// Function Advance_To_Target_Time
//#####################################################################
template<class TV> void SOLIDS_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    example.solids_parameters.triangle_collision_parameters.steps_since_self_collision_free=0;
    bool done=false;for(int substep=1;!done;substep++){
        LOG::SCOPE scope("SUBSTEP",STRING_UTILITIES::string_sprintf("substep %d",substep));
        
        Setup_Solids(time,substep);
        T dt=Compute_Dt(time,target_time,done);
        example.Preprocess_Substep(dt,time);
        Write_Substep("solid position update",substep,1);
        Solid_Position_Update(dt,substep);
        Write_Substep("solid velocity update",substep,1);
        Solid_Velocity_Update(dt,substep,done);
        example.Postprocess_Substep(dt,time);
        time+=dt;
        Write_Substep(STRING_UTILITIES::string_sprintf("END Substep %d",substep),substep,0);}
}
//#####################################################################
// Function Setup_Solids
//#####################################################################
template<class TV> void SOLIDS_DRIVER<TV>::
Setup_Solids(const T time,const int substep)
{
    SOLIDS_PARAMETERS<TV>& solids_parameters=example.solids_parameters;
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callbacks=solids_evolution.solids_evolution_callbacks;

    if(solids_parameters.triangle_collision_parameters.perform_self_collision && solids_parameters.triangle_collision_parameters.temporary_enable_collisions){
        solids_evolution_callbacks->Self_Collisions_Begin_Callback(time,substep);
        solids_parameters.triangle_collision_parameters.repulsion_pair_update_count=0;
        example.solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.Save_Self_Collision_Free_State();
        if((solids_parameters.triangle_collision_parameters.topological_hierarchy_build_count++)%solids_parameters.triangle_collision_parameters.topological_hierarchy_build_frequency==0){
            LOG::SCOPE scope("hierarchybuild","Building Hierarchy Topology");
            example.solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.Build_Topological_Structure_Of_Hierarchies();}
        solids_parameters.triangle_collision_parameters.self_collision_free_time=time;}

    if(solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies) example.solid_body_collection.rigid_body_collection.Reset_Impulse_Accumulators();
    solids_evolution_callbacks->Preprocess_Solids_Substep(time,substep);
    if(solids_parameters.deformable_object_collision_parameters.use_spatial_partition_for_levelset_collision_objects) // TODO - ANDY - why is this needed??? TODO: move this to the right places inside solids evolution 
        example.solid_body_collection.collision_body_list.Update_Spatial_Partition(solids_parameters.deformable_object_collision_parameters.spatial_partition_voxel_size_heuristic,
            solids_parameters.deformable_object_collision_parameters.spatial_partition_number_of_cells,solids_parameters.deformable_object_collision_parameters.spatial_partition_voxel_size_scale_factor);
    example.solid_body_collection.Update_Time_Varying_Material_Properties(time);
}
//#####################################################################
// Function Compute_Dt
//#####################################################################
template<class TV> typename TV::SCALAR SOLIDS_DRIVER<TV>::
Compute_Dt(const T time,const T target_time,bool& done)
{
    SOLIDS_PARAMETERS<TV>& solids_parameters=example.solids_parameters;
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callbacks=solids_evolution.solids_evolution_callbacks;

    // solids dt
    T solids_dt=FLT_MAX;
    if(!example.fixed_dt){
        if(solids_evolution.Use_CFL()) solids_dt=min(solids_dt,example.solid_body_collection.CFL(solids_parameters.verbose_dt));
        solids_dt=min(solids_dt,solids_evolution_callbacks->Constraints_CFL());
        if(solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies)
            solids_dt=min(solids_dt,example.solid_body_collection.rigid_body_collection.CFL_Rigid(solids_parameters.rigid_body_evolution_parameters,solids_parameters.verbose_dt));
        solids_evolution_callbacks->Limit_Solids_Dt(solids_dt,time);
        if(example.solid_body_collection.deformable_body_collection.mpi_solids)
            solids_dt=example.solid_body_collection.deformable_body_collection.mpi_solids->Reduce_Min_Global(solids_dt);}
    else solids_dt=example.fixed_dt;

    {std::stringstream ss;ss<<"dt = solids_dt = "<<solids_dt<<std::endl;LOG::filecout(ss.str());}
    if(example.abort_when_dt_below && solids_dt<example.abort_when_dt_below) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("dt too small (%g < %g)",solids_dt,example.abort_when_dt_below));
    done=false;
    EXAMPLE<TV>::Clamp_Time_Step_With_Target_Time(time,target_time,solids_dt,done,solids_parameters.min_dt);
    return solids_dt;
}
//#####################################################################
// Function Solid_Position_Update
//#####################################################################
template<class TV> void SOLIDS_DRIVER<TV>::
Solid_Position_Update(const T dt,const int substep)
{
     LOG::SCOPE scope("solids position update");

     SOLIDS_PARAMETERS<TV>& solids_parameters=example.solids_parameters;
     SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
     DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=example.solid_body_collection.deformable_body_collection;

     if((solids_parameters.triangle_collision_parameters.repulsion_pair_update_count++)%solids_parameters.triangle_collision_parameters.repulsion_pair_update_frequency==0){
         example.solid_body_collection.deformable_body_collection.triangle_repulsions.Update_Faces_And_Hierarchies_With_Collision_Free_Positions(&deformable_body_collection.particles.X);
         example.solid_body_collection.deformable_body_collection.triangle_repulsions.Compute_Interaction_Pairs(deformable_body_collection.particles.X);}
     solids_evolution.kinematic_evolution.Set_External_Positions(example.solid_body_collection.rigid_body_collection.rigid_body_particle.X,example.solid_body_collection.rigid_body_collection.rigid_body_particle.rotation,time);
     solids_evolution.kinematic_evolution.Set_External_Velocities(example.solid_body_collection.rigid_body_collection.rigid_body_particle.V,example.solid_body_collection.rigid_body_collection.rigid_body_particle.angular_velocity,time,time);
     example.solid_body_collection.rigid_body_collection.Update_Angular_Momentum();
     solids_evolution.Advance_One_Time_Step_Position(dt,time,true);

     if(solids_parameters.triangle_collision_parameters.perform_self_collision && solids_parameters.triangle_collision_parameters.temporary_enable_collisions){
         LOG::SCOPE scope("adjust velocity for self repulsion and self collisions");
         int repulsions,collisions_found;
         solids_evolution.Adjust_Velocity_For_Self_Repulsion_And_Self_Collisions(dt,time,repulsions,collisions_found,false);
         solids_parameters.triangle_collision_parameters.steps_since_self_collision_free=0;}

     Write_Substep("solid position updated",0,1);
}
//#####################################################################
// Function Solid_Velocity_Update
//#####################################################################
template<class TV> void SOLIDS_DRIVER<TV>::
Solid_Velocity_Update(const T dt,const int substep,const bool done)
{
    LOG::SCOPE scope("solids velocity update");

    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callbacks=solids_evolution.solids_evolution_callbacks;

    solids_evolution.Advance_One_Time_Step_Velocity(dt,time,true);
    solids_evolution.time+=dt;
    solids_evolution_callbacks->Postprocess_Solids_Substep(solids_evolution.time,substep);
    solids_evolution_callbacks->Apply_Constraints(dt,solids_evolution.time);
}
//#####################################################################
template class SOLIDS_DRIVER<VECTOR<float,1> >;
template class SOLIDS_DRIVER<VECTOR<float,2> >;
template class SOLIDS_DRIVER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SOLIDS_DRIVER<VECTOR<double,1> >;
template class SOLIDS_DRIVER<VECTOR<double,2> >;
template class SOLIDS_DRIVER<VECTOR<double,3> >;
#endif

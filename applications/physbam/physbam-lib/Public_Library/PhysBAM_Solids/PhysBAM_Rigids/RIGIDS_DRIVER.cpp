//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_DRIVER
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_IMPULSE_ACCUMULATOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/RIGIDS_DRIVER.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigids_Evolution/RIGIDS_EVOLUTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGIDS_DRIVER<TV>::
RIGIDS_DRIVER(RIGIDS_EXAMPLE<TV>& example_input)
    :BASE(example_input),example(example_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGIDS_DRIVER<TV>::
~RIGIDS_DRIVER()
{}
//#####################################################################
// Function Rigid_Cluster_Fracture
//#####################################################################
template<class TV> void RIGIDS_DRIVER<TV>::
Rigid_Cluster_Fracture(const T dt_full_advance,const T dt_cfl,const int substep)
{
    RIGIDS_EVOLUTION<TV>& rigids_evolution=*example.rigids_evolution;
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=example.rigid_body_collection.rigid_body_cluster_bindings;
    ARRAY<int> active_clusters;

    if(rigid_bindings.callbacks && ((substep-1)%example.rigids_parameters.rigid_cluster_fracture_frequency)==0 && rigid_bindings.Size()){
        T dt=min(dt_cfl*example.rigids_parameters.rigid_cluster_fracture_frequency,dt_full_advance);
        rigid_bindings.Deactivate_And_Return_Clusters(active_clusters,0);
        example.rigid_body_collection.Update_Simulated_Particles();

        rigid_bindings.callbacks->Pre_Advance_Unclustered(dt,time);
        example.rigids_evolution->kinematic_evolution.Set_External_Positions(example.rigid_body_collection.rigid_body_particle.X,
            example.rigid_body_collection.rigid_body_particle.rotation,time);
        example.rigids_evolution->kinematic_evolution.Set_External_Velocities(example.rigid_body_collection.rigid_body_particle.V,example.rigid_body_collection.rigid_body_particle.angular_velocity,time,time);
        example.rigid_body_collection.Update_Angular_Momentum();
        rigids_evolution.Advance_One_Time_Step_Position(dt,time,true);
        rigid_bindings.callbacks->Post_Advance_Unclustered(dt,time);
        rigid_bindings.callbacks->Compute_New_Clusters_Based_On_Unclustered_Strain();

        example.rigids_evolution->Restore_Position(rigids_evolution.rigid_X_save,rigids_evolution.rigid_rotation_save);
        rigid_bindings.Reactivate_Bindings(active_clusters);

        rigid_bindings.callbacks->Create_New_Clusters();
        example.rigid_body_collection.Update_Simulated_Particles();

        if(example.rigids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies) example.rigid_body_collection.Reset_Impulse_Accumulators();
        rigids_evolution.rigids_evolution_callbacks->Preprocess_Solids_Substep(time,substep);
        if(example.rigids_parameters.use_spatial_partition_for_levelset_collision_objects) // TODO - ANDY - why is this needed??? TODO: move this to the right places inside solids evolution 
            example.collision_body_list.Update_Spatial_Partition(example.rigids_parameters.spatial_partition_voxel_size_heuristic,
                example.rigids_parameters.spatial_partition_number_of_cells,example.rigids_parameters.spatial_partition_voxel_size_scale_factor);
    }
}
//#####################################################################
// Function Advance_To_Target_Time
//#####################################################################
template<class TV> void RIGIDS_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    RIGIDS_PARAMETERS<TV>& rigids_parameters=example.rigids_parameters;
    RIGIDS_EVOLUTION<TV>& rigids_evolution=*example.rigids_evolution;
    RIGIDS_EVOLUTION_CALLBACKS<TV>* rigids_evolution_callbacks=rigids_evolution.rigids_evolution_callbacks;
    
    bool done=false;for(int substep=1;!done;substep++){
        LOG::SCOPE scope("SUBSTEP",STRING_UTILITIES::string_sprintf("substep %d",substep));
        
        if(rigids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies)
            example.rigid_body_collection.Reset_Impulse_Accumulators();
        
        rigids_evolution_callbacks->Preprocess_Solids_Substep(time,substep);
        if(example.mpi_rigids || rigids_parameters.use_spatial_partition_for_levelset_collision_objects) {// TODO - ANDY - why is this needed??? TODO: move this to the right places
            // inside solids evolution 
            example.collision_body_list.Update_Spatial_Partition(rigids_parameters.spatial_partition_voxel_size_heuristic,
                rigids_parameters.spatial_partition_number_of_cells,rigids_parameters.spatial_partition_voxel_size_scale_factor);}
    
        // solids dt
        T solids_dt=FLT_MAX;
        if(!example.fixed_dt){
            if(rigids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies)
                solids_dt=min(solids_dt,example.rigid_body_collection.CFL_Rigid(rigids_parameters.rigid_body_evolution_parameters,rigids_parameters.verbose_dt));
            rigids_evolution_callbacks->Limit_Solids_Dt(solids_dt,time);}
        else solids_dt=example.fixed_dt;
        T dt=solids_dt;

        if(example.mpi_rigids) example.mpi_rigids->Synchronize_Dt(dt);
        {std::stringstream ss;ss<<"dt = solids_dt = "<<dt<<std::endl;LOG::filecout(ss.str());}

        done=false;
        EXAMPLE<TV>::Clamp_Time_Step_With_Target_Time(time,target_time,dt,done,rigids_parameters.min_dt);

        Rigid_Cluster_Fracture(target_time-time,dt,substep);
            
        example.Preprocess_Substep(dt,time);
        rigids_evolution.kinematic_evolution.Set_External_Positions(example.rigid_body_collection.rigid_body_particle.X,
            example.rigid_body_collection.rigid_body_particle.rotation,time);
        rigids_evolution.kinematic_evolution.Set_External_Velocities(example.rigid_body_collection.rigid_body_particle.V,example.rigid_body_collection.rigid_body_particle.angular_velocity,time,time);
        example.rigid_body_collection.Update_Angular_Momentum();

        rigids_evolution.Advance_One_Time_Step_Position(dt,time,true);

        rigids_evolution.Advance_One_Time_Step_Velocity(dt,time,true);
        
        rigids_evolution.time+=dt;
        rigids_evolution_callbacks->Postprocess_Solids_Substep(rigids_evolution.time,substep);
        rigids_evolution_callbacks->Apply_Constraints(dt,rigids_evolution.time);

        example.Postprocess_Substep(dt,time);

        time+=dt;
    }
}
//#####################################################################
// Function Execute_Main_Program
//#####################################################################
template<class TV> void RIGIDS_DRIVER<TV>::
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
// Function Initialize
//#####################################################################
template<class TV> void RIGIDS_DRIVER<TV>::
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

    RIGIDS_PARAMETERS<TV>& rigids_parameters=example.rigids_parameters;
    RIGIDS_EVOLUTION<TV>& rigids_evolution=*example.rigids_evolution;
    rigids_evolution.Set_Rigids_Evolution_Callbacks(example);
    rigids_evolution.mpi_rigids=example.mpi_rigids;

    example.Initialize_Bodies();
    example.Parse_Late_Options();
    rigids_evolution.time=time;
    if(example.mpi_rigids) for(int j=1;j<=example.rigid_body_collection.rigid_body_particle.array_collection->Size();j++) if(example.rigid_body_collection.rigid_geometry_collection.Is_Active(j))
        example.rigid_body_collection.Rigid_Body(j).impulse_accumulator=new RIGID_BODY_IMPULSE_ACCUMULATOR<TV,TV::dimension-1>(example.rigid_body_collection.Rigid_Body(j));

    if(example.restart){
        LOG::SCOPE scope("reading solids data");
        example.Read_Output_Files_Solids(example.restart_frame);
        rigids_evolution.time=time=example.Time_At_Frame(example.restart_frame);}

    example.rigid_body_collection.rigid_geometry_collection.collision_body_list->Update_Spatial_Partition(rigids_parameters.spatial_partition_voxel_size_heuristic,
        rigids_parameters.spatial_partition_number_of_cells,rigids_parameters.spatial_partition_voxel_size_scale_factor);
    if(!example.restart && example.mpi_rigids)
        example.mpi_rigids->Initialize(example.rigid_body_collection,example.processes_per_dimension,
            *example.rigid_body_collection.rigid_geometry_collection.collision_body_list->spatial_partition);
    example.rigid_body_collection.Update_Simulated_Particles();
    rigids_evolution.Initialize_Rigid_Bodies(example.frame_rate,example.restart);
}
//#####################################################################
// Function Simulate_To_Frame
//#####################################################################
template<class TV> void RIGIDS_DRIVER<TV>::
Simulate_To_Frame(const int frame_input)
{
    while(current_frame<example.last_frame){
        LOG::SCOPE scope("FRAME",STRING_UTILITIES::string_sprintf("Frame %d",current_frame+1));
        Preprocess_Frame(current_frame+1);
        example.rigids_evolution->kinematic_evolution.Get_Current_Kinematic_Keyframes(example.Time_At_Frame(current_frame+1)-time,time);
        Advance_To_Target_Time(example.Time_At_Frame(current_frame+1));
        Postprocess_Frame(++current_frame);
        if(example.write_output_files && example.write_substeps_level==-1) Write_Output_Files(current_frame);
        else if(example.write_substeps_level!=-1) Write_Substep(STRING_UTILITIES::string_sprintf("END Frame %d",current_frame),0,example.write_substeps_level);
        {std::stringstream ss;ss<<"TIME = "<<time<<std::endl;LOG::filecout(ss.str());}
    }
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void RIGIDS_DRIVER<TV>::
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
template class RIGIDS_DRIVER<VECTOR<float,1> >;
template class RIGIDS_DRIVER<VECTOR<float,2> >;
template class RIGIDS_DRIVER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGIDS_DRIVER<VECTOR<double,1> >;
template class RIGIDS_DRIVER<VECTOR<double,2> >;
template class RIGIDS_DRIVER<VECTOR<double,3> >;
#endif

//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BW_DRIVER
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/BW_DRIVER.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/BW_BACKWARD_EULER_SYSTEM.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BW_DRIVER<TV>::
BW_DRIVER(SOLIDS_EXAMPLE<TV>& example_input)
    :BASE(example_input),example(example_input),bw_collisions(example.solid_body_collection)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BW_DRIVER<TV>::
~BW_DRIVER()
{}
//#####################################################################
// Function Execute_Main_Program
//#####################################################################
template<class TV> void BW_DRIVER<TV>::
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
template<class TV> void BW_DRIVER<TV>::
Simulate_To_Frame(const int frame_input)
{
    while(current_frame<frame_input){
        LOG::SCOPE scope("FRAME","Frame %d",current_frame+1);
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
template<class TV> void BW_DRIVER<TV>::
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
template<class TV> void BW_DRIVER<TV>::
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
template<class TV> void BW_DRIVER<TV>::
Initialize()
{
    current_frame=example.first_frame;
    output_number=current_frame;
    time=example.Time_At_Frame(current_frame);

    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    solids_evolution.Set_Solids_Evolution_Callbacks(example);
    example.Initialize_Bodies();

    example.Parse_Late_Options();
    solids_evolution.time=time;

    solids_evolution.Initialize_Deformable_Objects(example.frame_rate,example.restart);

    solids_evolution.Initialize_Rigid_Bodies(example.frame_rate,example.restart);
}
//#####################################################################
// Function Advance_To_Target_Time
//#####################################################################
template<class TV> void BW_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;for(int substep=1;!done;substep++){
        LOG::SCOPE scope("SUBSTEP","substep %d",substep);
        
        if(!example.fixed_dt) PHYSBAM_FATAL_ERROR("Currently only works if using a fixed time step");
        T dt=example.fixed_dt;
        EXAMPLE<TV>::Clamp_Time_Step_With_Target_Time(time,target_time,dt,done,dt);
        example.Preprocess_Substep(dt,time);
        Advance_Substep(dt,time);
        example.Postprocess_Substep(dt,time);
        time+=dt;
        Write_Substep(STRING_UTILITIES::string_sprintf("END Substep %d",substep),substep,0);}
}
//#####################################################################
// Function Advance_Substep
//#####################################################################
template<class TV> void BW_DRIVER<TV>::
Advance_Substep(const T dt,const T time)
{
    Diagnostics(dt,time,1,"begin step");
    SOLID_BODY_COLLECTION<TV>& solid_body_collection=example.solid_body_collection;
    BW_BACKWARD_EULER_SYSTEM<TV> system(solid_body_collection,bw_collisions,dt,time);
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particle;

    bw_collisions.Detect_Cloth_Body_Contact();

    // Form the right hand side
    B_full.Resize(particles.array_collection->Size(),false,false);
    rigid_B_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
    GENERALIZED_VELOCITY<TV> B(B_full,rigid_B_full,solid_body_collection);
    ARRAYS_COMPUTATIONS::Fill(B_full,TV());ARRAYS_COMPUTATIONS::Fill(rigid_B_full,TWIST<TV>());
    // Update Position Based State
    solid_body_collection.Update_Position_Based_State(time,true);
    // Add implicit velocity dependent forces
    ARRAY<TWIST<TV> > twist(rigid_body_particles.V.Size());
    for(int i=1;i<=twist.Size();i++) twist(i)=TWIST<TV>(rigid_body_particles.V(i),rigid_body_particles.angular_velocity(i));
    GENERALIZED_VELOCITY<TV> V_not(particles.V,twist,solid_body_collection);
    solid_body_collection.Implicit_Velocity_Independent_Forces(V_not.V.array,V_not.rigid_V.array,B.V.array,B.rigid_V.array,dt,time);
    // Add velocity independent forces
    solid_body_collection.Add_Velocity_Independent_Forces(B_full,rigid_B_full,time);
    B_full*=(dt*solid_body_collection.deformable_body_collection.particles.one_over_mass);

    // TODO make sure the V we pass is delta V
    F_full.Resize(particles.array_collection->Size(),false,false);rigid_F_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
    R_full.Resize(particles.array_collection->Size(),false,false);rigid_R_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
    S_full.Resize(particles.array_collection->Size(),false,false);rigid_S_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
    V_full.Resize(particles.array_collection->Size(),false,false);rigid_V_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
    GENERALIZED_VELOCITY<TV> V(V_full,rigid_V_full,solid_body_collection),F(F_full,rigid_F_full,solid_body_collection),R(R_full,rigid_R_full,solid_body_collection),
        S(S_full,rigid_S_full,solid_body_collection),AR(AR_full,rigid_AR_full,solid_body_collection);

    static int solve_id=0;solve_id++;
    if(example.solids_parameters.implicit_solve_parameters.print_matrix){
        {std::stringstream ss;ss << "solve id " << solve_id << std::endl;LOG::filecout(ss.str());}
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("M-%i.txt",solve_id).c_str()).Write("M",system,S,R);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("b-%i.txt",solve_id).c_str()).Write("b",B);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("P-%i.txt",solve_id).c_str()).Write_Projection("P",system,S);}
    if(example.solids_parameters.implicit_solve_parameters.test_system) system.Test_System(S,R,F);

    static CONJUGATE_GRADIENT<T> cg;
    KRYLOV_SOLVER<T>* solver=&cg;
    solver->print_residuals=solid_body_collection.print_residuals;
    solver->Solve(system,V,B,F,S,R,AR,example.solids_parameters.implicit_solve_parameters.cg_tolerance,1,example.solids_parameters.implicit_solve_parameters.cg_iterations);

    if(example.solids_parameters.implicit_solve_parameters.print_matrix) OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("x-%i.txt",solve_id).c_str()).Write("x",V);

    bw_collisions.Remove_Separating_Cloth_Body_Contacts(system,R,B,V,F);

    particles.V+=V_full; // TODO do rigid as well
    particles.X+=dt*particles.V; // TODO do rigid as well
    for(int i=1;i<=twist.Size();i++){rigid_body_particles.V(i)=twist(i).linear;rigid_body_particles.angular_velocity(i)=twist(i).angular;}
    Diagnostics(dt,time+dt,100,"finish step");
}
//#####################################################################
// Function Diagnostics
//#####################################################################
template<class TV> void BW_DRIVER<TV>::
Diagnostics(const T dt,const T time,int step,const char* description)
{
    example.solid_body_collection.Print_Energy(time,step);
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Finished step %i (%s).  dt=%f time=%f",step,description,dt,time),2,3);
}
//#####################################################################
template class BW_DRIVER<VECTOR<float,1> >;
template class BW_DRIVER<VECTOR<float,2> >;
template class BW_DRIVER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BW_DRIVER<VECTOR<double,1> >;
template class BW_DRIVER<VECTOR<double,2> >;
template class BW_DRIVER<VECTOR<double,3> >;
#endif

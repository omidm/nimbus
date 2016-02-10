//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_DRIVER
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLIDS_FLUIDS_DRIVER<TV>::
SOLIDS_FLUIDS_DRIVER(SOLIDS_FLUIDS_EXAMPLE<TV>& example_input)
    :BASE(example_input),example(example_input),project_at_frame_boundaries(true)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLIDS_FLUIDS_DRIVER<TV>::
~SOLIDS_FLUIDS_DRIVER()
{}
//#####################################################################
// Function Execute_Main_Program
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_DRIVER<TV>::
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
template<class TV> void SOLIDS_FLUIDS_DRIVER<TV>::
Initialize()
{
    if(example.auto_restart){
        std::string last_frame_file=example.output_directory+"/common/last_frame";
        int last_frame;FILE_UTILITIES::Read_From_Text_File(last_frame_file,last_frame);
        example.restart=true;example.restart_frame=last_frame;
        std::stringstream ss;ss<<"Auto Restart from frame "<<last_frame<<" (from file "<<last_frame_file<<")"<<std::endl;LOG::filecout(ss.str());}
    if(example.restart){current_frame=example.restart_frame;Read_Time(current_frame);}else current_frame=example.first_frame;
    output_number=current_frame;
    time=example.Time_At_Frame(current_frame);
    example.solid_body_collection.deformable_body_collection.mpi_solids=example.solid_body_collection.deformable_body_collection.mpi_solids;
}
//#####################################################################
// Function Simulate_To_Frame
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_DRIVER<TV>::
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
        std::stringstream ss;ss<<"TIME = "<<time<<std::endl;LOG::filecout(ss.str());}
}
//#####################################################################
template class SOLIDS_FLUIDS_DRIVER<VECTOR<float,1> >;
template class SOLIDS_FLUIDS_DRIVER<VECTOR<float,2> >;
template class SOLIDS_FLUIDS_DRIVER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SOLIDS_FLUIDS_DRIVER<VECTOR<double,1> >;
template class SOLIDS_FLUIDS_DRIVER<VECTOR<double,2> >;
template class SOLIDS_FLUIDS_DRIVER<VECTOR<double,3> >;
#endif

//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_DRIVER
//#####################################################################
#ifndef __SOLIDS_FLUIDS_DRIVER__
#define __SOLIDS_FLUIDS_DRIVER__    

#include <PhysBAM_Tools/Ordinary_Differential_Equations/DRIVER.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE.h>
namespace PhysBAM{

template<class TV>
class SOLIDS_FLUIDS_DRIVER:public DRIVER<TV>,public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef DRIVER<TV> BASE;
public:
    using BASE::output_number;using BASE::time;using BASE::Write_Output_Files;using BASE::Read_Time;

    SOLIDS_FLUIDS_EXAMPLE<TV>& example;
    bool project_at_frame_boundaries;
    int current_frame;
    T next_dt; // for fluid time stepping
    bool next_done; // for fluid time stepping

    SOLIDS_FLUIDS_DRIVER(SOLIDS_FLUIDS_EXAMPLE<TV>& example_input);
    virtual ~SOLIDS_FLUIDS_DRIVER();

    virtual void Preprocess_Frame(const int frame)
    {if(example.substeps_delay_frame==frame){example.Set_Write_Substeps_Level(example.substeps_delay_level);output_number=frame-1;}
    example.Preprocess_Frame(frame);}
    
    virtual void Postprocess_Frame(const int frame)
    {example.Postprocess_Frame(frame);}

//#####################################################################
    void Advance_To_Target_Time(const T target_time){};
    void Execute_Main_Program() PHYSBAM_OVERRIDE;
    void Initialize() PHYSBAM_OVERRIDE;
    void Simulate_To_Frame(const int frame_input) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

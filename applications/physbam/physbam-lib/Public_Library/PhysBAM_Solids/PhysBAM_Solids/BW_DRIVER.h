//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BW_DRIVER
//#####################################################################
#ifndef __BW_DRIVER__
#define __BW_DRIVER__    

#include <PhysBAM_Tools/Ordinary_Differential_Equations/DRIVER.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/BW_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/SOLIDS_EXAMPLE.h>

namespace PhysBAM{

template<class TV>
class BW_DRIVER:public DRIVER<TV>,public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef DRIVER<TV> BASE;
public:
    using BASE::output_number;using BASE::time;using BASE::Write_Time;using BASE::Write_First_Frame;using BASE::Write_Last_Frame;using BASE::Write_Substep;

    SOLIDS_EXAMPLE<TV>& example;
    int current_frame;
    ARRAY<TV> V_full,F_full,R_full,B_full,S_full,AR_full;
    ARRAY<TWIST<TV> > rigid_V_full,rigid_F_full,rigid_R_full,rigid_S_full,rigid_B_full,rigid_AR_full;
    BW_COLLISIONS<TV> bw_collisions;

    BW_DRIVER(SOLIDS_EXAMPLE<TV>& example_input);
    virtual ~BW_DRIVER();

    virtual void Preprocess_Frame(const int frame)
    {if(example.substeps_delay_frame==frame){example.Set_Write_Substeps_Level(example.substeps_delay_level);output_number=frame-1;}
    example.Preprocess_Frame(frame);}

    virtual void Postprocess_Frame(const int frame)
    {example.Postprocess_Frame(frame);}

//#####################################################################
    void Advance_To_Target_Time(const T target_time) PHYSBAM_OVERRIDE;
    void Execute_Main_Program() PHYSBAM_OVERRIDE;
    void Initialize() PHYSBAM_OVERRIDE;
    void Simulate_To_Frame(const int frame_input) PHYSBAM_OVERRIDE;
    void Write_Output_Files(const int frame) PHYSBAM_OVERRIDE;
    virtual void Read_Time(const int frame);
    void Advance_Substep(const T dt,const T time);
    void Diagnostics(const T dt,const T time,int step,const char* description);
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLES_DRIVER
//#####################################################################
#ifndef __DEFORMABLES_DRIVER__
#define __DEFORMABLES_DRIVER__    

#include <PhysBAM_Tools/Ordinary_Differential_Equations/DRIVER.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/DEFORMABLES_EXAMPLE.h>
namespace PhysBAM{

template<class TV>
class DEFORMABLES_DRIVER:public DRIVER<TV>,public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef DRIVER<TV> BASE;
public:
    using BASE::output_number;using BASE::time;using BASE::Write_First_Frame;using BASE::Write_Last_Frame;using BASE::Write_Time;using BASE::Read_Time;

    DEFORMABLES_EXAMPLE<TV>& example;
    int current_frame;

    DEFORMABLES_DRIVER(DEFORMABLES_EXAMPLE<TV>& example_input);
    virtual ~DEFORMABLES_DRIVER();

    virtual void Preprocess_Frame(const int frame)
    {if(example.substeps_delay_frame==frame){example.Set_Write_Substeps_Level(example.substeps_delay_level);output_number=frame-1;}
    example.Preprocess_Frame(frame);}
    
    virtual void Postprocess_Frame(const int frame)
    {example.Postprocess_Frame(frame);}

//#####################################################################
    void Advance_To_Target_Time(const T target_time) PHYSBAM_OVERRIDE;
    void Execute_Main_Program() PHYSBAM_OVERRIDE;
    void Initialize() PHYSBAM_OVERRIDE;
    void Update_Structures();
    void Simulate_To_Frame(const int frame_input) PHYSBAM_OVERRIDE;
    void Write_Output_Files(const int frame) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_INTERACTIVE_DRIVER
//#####################################################################
#ifndef __RIGIDS_INTERACTIVE_DRIVER__
#define __RIGIDS_INTERACTIVE_DRIVER__    

#include <PhysBAM_Tools/Ordinary_Differential_Equations/DRIVER.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/RIGIDS_EXAMPLE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/ANIMATED_VISUALIZATION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Rigids_Components/OPENGL_COMPONENT_POLICY.h>
#include <time.h>
namespace PhysBAM{

template<class TV>
class RIGIDS_INTERACTIVE_DRIVER:public DRIVER<TV>,public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef DRIVER<TV> BASE;
    typedef typename OPENGL_COMPONENT_POLICY<TV>::RIGID_BODY_COLLECTION T_COMPONENT;
public:
    using BASE::output_number;using BASE::time;using BASE::Write_First_Frame;using BASE::Write_Last_Frame;using BASE::Write_Time;using BASE::Read_Time;

    RIGIDS_EXAMPLE<TV>& example;
    int current_frame;
    ANIMATED_VISUALIZATION* visualizer;
    T_COMPONENT* rigid_component;
    std::clock_t since_last_frame,now;

    RIGIDS_INTERACTIVE_DRIVER(RIGIDS_EXAMPLE<TV>& example_input);
    virtual ~RIGIDS_INTERACTIVE_DRIVER();

    virtual void Preprocess_Frame(const int frame)
    {if(example.substeps_delay_frame==frame){example.Set_Write_Substeps_Level(example.substeps_delay_level);output_number=frame-1;}
    example.Preprocess_Frame(frame);}
    
    virtual void Postprocess_Frame(const int frame)
    {example.Postprocess_Frame(frame);}

//#####################################################################
    void Rigid_Cluster_Fracture(const T dt_full_advance,const T dt_cfl,const int substep);
    void Advance_To_Target_Time(const T target_time) PHYSBAM_OVERRIDE;
    void Execute_Main_Program() PHYSBAM_OVERRIDE;
    void Initialize() PHYSBAM_OVERRIDE;
    void Simulate_To_Frame(const int frame_input) PHYSBAM_OVERRIDE;
    void Write_Output_Files(const int frame) PHYSBAM_OVERRIDE;
    void Update_Positions_From_Interactive();
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PLS_INTERACTIVE_DRIVER__
#define __PLS_INTERACTIVE_DRIVER__
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/ANIMATED_VISUALIZATION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Rigids_Components/OPENGL_COMPONENT_POLICY.h>
#ifndef WIN32
    #define USE_UNIX_TIME 1
    #include <sys/time.h>    
#endif
namespace PhysBAM{


template<class TV> class PLS_EXAMPLE;

template<class TV>
class PLS_INTERACTIVE_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;
    typedef typename ADVECTION_POLICY<GRID<TV> >::ADVECTION_SEMI_LAGRANGIAN_SCALAR T_ADVECTION_SEMI_LAGRANGIAN_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename OPENGL_COMPONENT_POLICY<TV>::RIGID_GEOMETRY_COLLECTION T_RIGIDS_COMPONENT;    
    typedef typename OPENGL_COMPONENT_POLICY<TV>::LEVELSET T_LEVELSET_COMPONENT;    
    typedef typename OPENGL_COMPONENT_POLICY<TV>::MAC_VELOCITY_FIELD T_MAC_COMPONENT;

protected:
    int current_frame;
    T time;
    int output_number;

    PLS_EXAMPLE<TV>& example;
public:
    ANIMATED_VISUALIZATION* visualizer;
    T_RIGIDS_COMPONENT* rigid_component;
    T_MAC_COMPONENT* mac_component;
    T_LEVELSET_COMPONENT* levelset_component;
#if USE_UNIX_TIME
    timeval since_last_frame,now;
#endif

    PLS_INTERACTIVE_DRIVER(PLS_EXAMPLE<TV>& example);
    virtual ~PLS_INTERACTIVE_DRIVER();
    
    void Update_Positions_From_Interactive();
    void Execute_Main_Program();
    void Initialize();
    void Advance_To_Target_Time(const T target_time);
    void Simulate_To_Frame(const int frame_input);
    void Write_Output_Files(const int frame);
    void Write_Substep(const std::string& title,const int substep,const int level=0);

//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __INCOMPRESSIBLE_INTERACTIVE_DRIVER__
#define __INCOMPRESSIBLE_INTERACTIVE_DRIVER__
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_UNIFORM_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Solids_Geometry_Evolution/KINEMATIC_EVOLUTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/ANIMATED_VISUALIZATION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Rigids_Components/OPENGL_COMPONENT_POLICY.h>
#include <time.h>
namespace PhysBAM{

template<class TV> class INCOMPRESSIBLE_EXAMPLE;

template<class TV>
class INCOMPRESSIBLE_INTERACTIVE_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename OPENGL_COMPONENT_POLICY<TV>::RIGID_GEOMETRY_COLLECTION T_RIGIDS_COMPONENT;    
    typedef typename OPENGL_COMPONENT_POLICY<TV>::SCALAR_FIELD T_SCALAR_COMPONENT;    
    typedef typename OPENGL_COMPONENT_POLICY<TV>::MAC_VELOCITY_FIELD T_MAC_COMPONENT;    

protected:
    int current_frame;
    T time;
    int output_number;

    INCOMPRESSIBLE_EXAMPLE<TV>& example;
    KINEMATIC_EVOLUTION<TV> kinematic_evolution;
public:
    ANIMATED_VISUALIZATION* visualizer;
    T_RIGIDS_COMPONENT* rigid_component;
    T_MAC_COMPONENT* mac_component;
    T_SCALAR_COMPONENT* scalar_component;
    std::clock_t since_last_frame,now;

    INCOMPRESSIBLE_INTERACTIVE_DRIVER(INCOMPRESSIBLE_EXAMPLE<TV>& example);
    virtual ~INCOMPRESSIBLE_INTERACTIVE_DRIVER();

    void Update_Positions_From_Interactive();
    void Scalar_Advance(const T dt,const T time);
    void Convect(const T dt,const T time);
    void Add_Forces(const T dt,const T time);
    void Project(const T dt,const T time);
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

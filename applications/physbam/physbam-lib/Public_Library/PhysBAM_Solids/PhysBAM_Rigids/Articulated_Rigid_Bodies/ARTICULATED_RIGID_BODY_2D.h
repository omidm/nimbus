//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Craig Schroeder, Tamar Shinar, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARTICULATED_RIGID_BODY_2D
//#####################################################################
#ifndef __ARTICULATED_RIGID_BODY_2D__
#define __ARTICULATED_RIGID_BODY_2D__

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY.h>
namespace PhysBAM{

template<class TV> class JOINT_FUNCTION;
template<class T_input>
class ARTICULATED_RIGID_BODY<VECTOR<T_input,2> >:public ARTICULATED_RIGID_BODY_BASE<VECTOR<T_input,2> >
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;
public:
    typedef ARTICULATED_RIGID_BODY_BASE<TV> BASE;
    using BASE::rigid_body_collection;using BASE::joint_mesh;using BASE::muscle_list;using BASE::muscle_activations;using BASE::iterative_tolerance;
    using BASE::max_iterations;using BASE::poststabilization_iterations;using BASE::use_epsilon_scale;using BASE::line_search_interval_tolerance;using BASE::max_line_search_iterations;
    using BASE::use_pd_actuators;using BASE::use_muscle_actuators;using BASE::actuation_iterations;using BASE::Update_Joint_Frames;using BASE::Parent_Id;using BASE::Child_Id;
    using BASE::Parent;using BASE::Child;using BASE::Apply_Poststabilization;using BASE::Apply_Poststabilization_To_Joint;

    ARRAY<ARRAY<PAIR<int,T> > > muscles_crossing_joints; // muscle id and moment arm

    ARTICULATED_RIGID_BODY(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input);

    virtual ~ARTICULATED_RIGID_BODY();

//#####################################################################
    void Compute_Position_Based_State(const T dt,const T time);
    JOINT_FUNCTION<TV>* Create_Joint_Function(const JOINT_ID joint_id);
    void Post_Stabilization_With_Actuation(const JOINT_ID joint_id);
    VECTOR<T,1> Compute_Target_PD_Angular_Impulse(const JOINT_ID joint_id);
    void Solve_Velocities_for_PD(const T time,const T dt,bool test_system,bool print_matrix) PHYSBAM_OVERRIDE;
    void Output_Articulation_Points(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2004-2007, Kevin Der, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARTICULATED_RIGID_BODY_3D
//#####################################################################
#ifndef __ARTICULATED_RIGID_BODY_3D__
#define __ARTICULATED_RIGID_BODY_3D__

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY.h>
namespace PhysBAM{

template<class TV> class JOINT_FUNCTION;
template<class T_input>
class ARTICULATED_RIGID_BODY<VECTOR<T_input,3> >:public ARTICULATED_RIGID_BODY_BASE<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef ARTICULATED_RIGID_BODY_BASE<TV> BASE;
    using BASE::rigid_body_collection;using BASE::joint_mesh;using BASE::muscle_list;using BASE::muscle_activations;using BASE::iterative_tolerance;
    using BASE::max_iterations;using BASE::poststabilization_iterations;using BASE::use_epsilon_scale;using BASE::line_search_interval_tolerance;using BASE::max_line_search_iterations;
    using BASE::use_pd_actuators;using BASE::use_muscle_actuators;using BASE::actuation_iterations;using BASE::enforce_nonnegative_activations;using BASE::clamp_negative_activations;
    using BASE::activation_optimization_iterations;using BASE::global_post_stabilization;using BASE::verbose;using BASE::use_angular_damping;using BASE::Update_Joint_Frames;
    using BASE::Parent_Id;using BASE::Child_Id;using BASE::Parent;using BASE::Child;using BASE::Apply_Poststabilization;using BASE::zero_row_tolerance;

private:
    VECTOR_ND<T> last_muscle_actuations;
    MATRIX_MXN<T> global_post_stabilization_matrix_11,global_post_stabilization_matrix_12,global_post_stabilization_matrix_21,global_post_stabilization_matrix_22;
    ARRAY<int> joint_constrained_dimensions,joint_muscle_control_dimensions;
    ARRAY<int> joint_offset_in_post_stabilization_matrix,joint_offset_in_muscle_control_matrix;
    ARRAY<MATRIX_MXN<T> > joint_angular_constraint_matrix,joint_angular_muscle_control_matrix;
    ARRAY<TV> precomputed_desired_damped_angular_velocities;
public:

    ARTICULATED_RIGID_BODY(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input);

    ~ARTICULATED_RIGID_BODY();

//#####################################################################
    void Compute_Position_Based_State(const T dt,const T time);
    static void Solve_Minimum_Norm_Solution_For_Linear_Constraints(MATRIX_MXN<T>& A,const VECTOR_ND<T>& b,VECTOR_ND<T>& x,const T zero_row_tolerance,const bool verbose=false);
    void Solve_For_Muscle_Control(MATRIX_MXN<T>& A,const VECTOR_ND<T>& b,VECTOR_ND<T>& x,const T dt);
    void Solve_Velocities_for_PD(const T time,const T dt,bool test_system,bool print_matrix) PHYSBAM_OVERRIDE;
    JOINT_FUNCTION<TV>* Create_Joint_Function(const JOINT_ID joint_id);
    void Compute_Desired_PD_Velocity(const T dt,const T time) PHYSBAM_OVERRIDE;
    void Output_Articulation_Points(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

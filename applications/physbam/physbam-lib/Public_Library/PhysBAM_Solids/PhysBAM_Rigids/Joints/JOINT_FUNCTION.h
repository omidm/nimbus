//#####################################################################
// Copyright 2005-2009, Ron Fedkiw, Jared Go, Eran Juandelam, Michael Lentine, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class JOINT_FUNCTION
//#####################################################################
#ifndef _JOINT_FUNCTION__
#define _JOINT_FUNCTION__

#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_1D.h>
#include <cassert>
#include <cmath>
namespace PhysBAM{

using ::std::exp;
using ::std::sqrt;

template<class TV> class JOINT;
template<class TV> class RIGID_BODY;
template<class T,class TV> class INTERPOLATION_CURVE;
template<class TV>
class JOINT_FUNCTION
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;
public:
    T k_p,k_v;
    T omega;
    bool muscle_control;
    JOINT<TV>& joint;
    RIGID_BODY<TV> *parent,*child;
    ROTATION<TV> target_angle;
    T_SPIN target_angular_velocity;
    // Our pd targets at the end of this step
    T_SPIN desired_angular_velocity;
    INTERPOLATION_CURVE<T,ROTATION<TV> >* interpolation_curve;
    bool target_position;
    bool active;

    JOINT_FUNCTION(JOINT<TV>& joint_input,RIGID_BODY<TV>& parent_input,RIGID_BODY<TV>& child_input)
        :muscle_control(false),joint(joint_input),parent(&parent_input),child(&child_input),interpolation_curve(0),target_position(true),active(true)
    {
        Set_k_p(1);
    }

    ~JOINT_FUNCTION()
    {}

    void Set_k_p(const T k_p_input)
    {k_p=k_p_input;omega=-sqrt(k_p);k_v=-2*omega;}

    void Set_Target_Angle(const ROTATION<TV>& target_angle_input)
    {target_position=true;target_angle=target_angle_input;target_angular_velocity=T_SPIN();}

    void Set_Target_Angular_Velocity(const T_SPIN& target_angular_velocity_input)
    {target_position=false;target_angular_velocity=target_angular_velocity_input;}

    VECTOR<T,1> PD_Error(const VECTOR<double,1> angle_error,const VECTOR<double,1> angular_velocity_error,const T dt)
    {return VECTOR<T,1>(angle_error+(angular_velocity_error-angle_error*(double)omega)*(double)dt)*exp(omega*dt);} // use doubles for higher precision

//#####################################################################
    ROTATION<TV> Angle();
    T_SPIN Angular_Velocity() const;
    ROTATION<TV> Target_Angle(const T time) const;
    T_SPIN Target_Angular_Velocity(const T time) const;
    void Compute_Desired_PD_Velocity(const T dt,const T time);
    void Prismatic_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& constrained_matrix,MATRIX_MXN<T>* unconstrained_matrix=0) const;
    void Angular_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& angular_constraint_matrix,MATRIX_MXN<T>* angular_unconstrained_matrix=0) const;
//#####################################################################
};
}
#endif

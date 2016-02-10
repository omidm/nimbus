//#####################################################################
// Copyright 2005-2007, Eran Guendelman, Michael Lentine, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Math_Tools/wrap.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
using namespace PhysBAM;
//#####################################################################
// Function Compute_Desired_PD_Velocity_Helper
//#####################################################################
template<class T> void
Compute_Desired_PD_Velocity_Helper(const T dt,const T time,JOINT_FUNCTION<VECTOR<T,2> >& jf)
{
    typedef VECTOR<T,1> T_SPIN;
    if(!jf.target_position){T_SPIN target=jf.Target_Angular_Velocity(time);jf.desired_angular_velocity=(jf.Angular_Velocity()-target)*exp(-jf.k_v*dt)+target;return;}
    T_SPIN target_minus_current_angle(wrap(jf.Target_Angle(time).Angle()-jf.Angle().Angle(),(T)-pi,(T)pi));
    T_SPIN target_minus_current_angular_velocity=jf.Target_Angular_Velocity(time)-jf.Angular_Velocity();
    T_SPIN desired_new_error=jf.PD_Error(VECTOR<double,1>(target_minus_current_angle),VECTOR<double,1>(target_minus_current_angular_velocity),dt);
    T_SPIN next_target_minus_current_angle(wrap(jf.Target_Angle(time+dt).Angle()-jf.Angle().Angle(),(T)-pi,(T)pi));
    jf.desired_angular_velocity=(next_target_minus_current_angle-desired_new_error)/dt;
}
//#####################################################################
// Function Compute_Desired_PD_Velocity_Helper
//#####################################################################
template<class T> void
Compute_Desired_PD_Velocity_Helper(const T dt,const T time,JOINT_FUNCTION<VECTOR<T,3> >& jf)
{
    typedef VECTOR<T,3> T_SPIN;typedef VECTOR<T,3> TV;
    ROTATION<TV> Fp_wj=jf.parent->Rotation()*jf.joint.F_pj().r;
    if(!jf.target_position){T_SPIN target=Fp_wj.Rotate(jf.Target_Angular_Velocity(time));jf.desired_angular_velocity=(jf.Angular_Velocity()-target)*exp(-jf.k_v*dt)+target;return;}
    T_SPIN rotation_axis_to_current_target=Fp_wj.Rotate((jf.Target_Angle(time)*jf.Angle().Inverse()).Rotation_Vector());
    T_SPIN rotation_axis_to_next_target=Fp_wj.Rotate((jf.Target_Angle(time+dt)*jf.Angle().Inverse()).Rotation_Vector());
    T target_minus_current_angle=wrap(rotation_axis_to_current_target.Normalize(),(T)-pi,(T)pi);
    T next_target_minus_current_angle=wrap(rotation_axis_to_next_target.Normalize(),(T)-pi,(T)pi);
    // target angular velocity in parent space; angular velocity already in world space
    T target_minus_current_angular_velocity=T_SPIN::Dot_Product(Fp_wj.Rotate(jf.Target_Angular_Velocity(time))-jf.Angular_Velocity(),rotation_axis_to_current_target);
    T desired_new_error=jf.PD_Error(VECTOR<double,1>(target_minus_current_angle),VECTOR<double,1>(target_minus_current_angular_velocity),dt).x;
    T_SPIN extra_angular_velocity=jf.Angular_Velocity().Projected_Orthogonal_To_Unit_Direction(rotation_axis_to_next_target);
    jf.desired_angular_velocity=((next_target_minus_current_angle-desired_new_error)/dt)*rotation_axis_to_next_target+extra_angular_velocity*exp(-jf.k_v*dt);
}
//#####################################################################
// Function Compute_Desired_PD_Velocity
//#####################################################################
template<class TV> void JOINT_FUNCTION<TV>::
Compute_Desired_PD_Velocity(const T dt,const T time)
{
    Compute_Desired_PD_Velocity_Helper(dt,time,*this);
}
//#####################################################################
// Function Prismatic_Constraint_Matrix
//#####################################################################
template<class TV> void JOINT_FUNCTION<TV>::
Prismatic_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& constrained_matrix,MATRIX_MXN<T>* unconstrained_matrix) const
{
}
//#####################################################################
// Function Angular_Constraint_Matrix
//#####################################################################
template<class TV> void JOINT_FUNCTION<TV>::
Angular_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& angular_constraint_matrix,MATRIX_MXN<T>* angular_unconstrained_matrix) const
{
    angular_constraint_matrix=MATRIX<T,T_SPIN::dimension>::Identity_Matrix();
    if(angular_unconstrained_matrix) angular_unconstrained_matrix->Resize(T_SPIN::dimension,0);
}
template<class TV> ROTATION<TV> JOINT_FUNCTION<TV>::
Angle()
{
    joint.Set_Joint_Frame(joint.Compute_Current_Joint_Frame(*parent,*child));
    return joint.J.r;
}
template<class TV> typename TV::SPIN JOINT_FUNCTION<TV>::
Angular_Velocity() const
{
    return RIGID_BODY<TV>::Relative_Angular_Velocity(*child,*parent); // child velocity w.r.t. parent velocity!
}
template<class TV> ROTATION<TV> JOINT_FUNCTION<TV>::
Target_Angle(const T time) const
{
    if(interpolation_curve) return interpolation_curve->Value(time);
    return target_angle;
}
template<class TV> typename TV::SPIN JOINT_FUNCTION<TV>::
Target_Angular_Velocity(const T time) const
{
    if(interpolation_curve) return interpolation_curve->Derivative(time);
    return target_angular_velocity;
}
//#####################################################################
template class JOINT_FUNCTION<VECTOR<float,2> >;
template class JOINT_FUNCTION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class JOINT_FUNCTION<VECTOR<double,2> >;
template class JOINT_FUNCTION<VECTOR<double,3> >;
#endif

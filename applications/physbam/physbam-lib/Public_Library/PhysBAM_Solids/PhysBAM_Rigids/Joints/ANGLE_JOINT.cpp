//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Eran Guendelman, Michael Lentine, Craig Schroeder, Tamar Shinar, Jonathan Su, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/ANGLE_JOINT.h>
using namespace PhysBAM;
//#####################################################################
// Destructor
//#####################################################################
template<class TV> ANGLE_JOINT<TV>::
~ANGLE_JOINT()
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_input> ANGLE_JOINT<VECTOR<T_input,1> >::
~ANGLE_JOINT()
{}
//#####################################################################
// Function Has_Angular_Constraint
//#####################################################################
template<class TV> bool ANGLE_JOINT<TV>::
Has_Angular_Constraint() const
{
    return TV::dimension==3 || constrain_angle;
}
//#####################################################################
// Function Constrain_Prismatically
//#####################################################################
template<class TV> void ANGLE_JOINT<TV>::
Constrain_Prismatically(TV& translation) const
{
    translation=TV();
}
//#####################################################################
// Function Constrain_Angles
//#####################################################################
template<class TV> void ANGLE_JOINT<TV>::
Constrain_Angles(T_SPIN& angles) const
{
    T_SPIN v;
    if(constrain_angle) v.x=clamp(angles.x,angle_min,angle_max);
    else v.x=angles.x;
    angles=v;
}
//#####################################################################
// Function Angular_Constraint_Matrix
//#####################################################################
template<class TV> void ANGLE_JOINT<TV>::
Angular_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& angular_constraint_matrix,MATRIX_MXN<T>* angular_unconstrained_matrix) const
{
    VECTOR<bool,T_SPIN::dimension> constrain(VECTOR<bool,T_SPIN::dimension>::All_Ones_Vector());
    constrain.x=constrain_angle && angle_min>=angle_max;
    Constraint_Matrix_Helper(parent_frame.r*F_pj().r,angular_constraint_matrix,angular_unconstrained_matrix,constrain);
}
//#####################################################################
// Function Constrain_Relative_Linear_Velocity
//#####################################################################
template<class TV> void ANGLE_JOINT<TV>::
Constrain_Relative_Linear_Velocity(const FRAME<TV>& parent_frame,TV& relative_linear_velocity) const
{
    relative_linear_velocity=TV();
}
//#####################################################################
// Function Update_State_From_Joint_Frame
//#####################################################################
template<class TV> void ANGLE_JOINT<TV>::
Update_State_From_Joint_Frame(const bool enforce_constraints)
{
    T_SPIN angles=J.r.Euler_Angles(); 
    if(enforce_constraints){Constrain_Angles(angles);
        Constrain_Prismatically(J.t);}
    J.r=ROTATION<TV>::From_Euler_Angles(angles);
    J_inverse=J.Inverse();
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d,s) \
    template ANGLE_JOINT<VECTOR<T,d> >::~ANGLE_JOINT(); \
    template bool ANGLE_JOINT<VECTOR<T,d> >::Has_Angular_Constraint() const; \
    template void ANGLE_JOINT<VECTOR<T,d> >::Constrain_Prismatically(VECTOR<T,d>& translation) const; \
    template void ANGLE_JOINT<VECTOR<T,d> >::Constrain_Angles(VECTOR<T,s>& angles) const; \
    template void ANGLE_JOINT<VECTOR<T,d> >::Angular_Constraint_Matrix(const FRAME<VECTOR<T,d> >& parent_frame,MATRIX_MXN<T>& angular_constraint_matrix, \
        MATRIX_MXN<T>* angular_unconstrained_matrix) const; \
    template void ANGLE_JOINT<VECTOR<T,d> >::Constrain_Relative_Linear_Velocity(const FRAME<VECTOR<T,d> >& parent_frame,VECTOR<T,d>& relative_linear_velocity) const; \
    template void ANGLE_JOINT<VECTOR<T,d> >::Update_State_From_Joint_Frame(const bool enforce_constraints);

INSTANTIATION_HELPER(float,2,1)
INSTANTIATION_HELPER(float,3,3)
template class ANGLE_JOINT<VECTOR<float,1> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,2,1)
INSTANTIATION_HELPER(double,3,3)
template class ANGLE_JOINT<VECTOR<double,1> >;
#endif

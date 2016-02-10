//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Michael Lentine, Craig Schroeder, Tamar Shinar, Jonathan Su, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
using namespace PhysBAM;
//#####################################################################
// Destructor
//#####################################################################
template<class TV> POINT_JOINT<TV>::
~POINT_JOINT()
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_input> POINT_JOINT<VECTOR<T_input,1> >::
~POINT_JOINT()
{}
//#####################################################################
// Function Constrain_Relative_Linear_Velocity
//#####################################################################
template<class TV> void POINT_JOINT<TV>::
Constrain_Relative_Linear_Velocity(const FRAME<TV>& parent_frame,TV& relative_linear_velocity) const
{
    relative_linear_velocity=TV();
}
//#####################################################################
// Function Prismatic_Component_Translation
//#####################################################################
template<class TV> TV POINT_JOINT<TV>::
Prismatic_Component_Translation() const
{
    return prismatic_translation;
}
//#####################################################################
// Function Constrain_Prismatically
//#####################################################################
template<class TV> void POINT_JOINT<TV>::
Constrain_Prismatically(TV& translation) const
{
    translation=prismatic_translation;
}
//#####################################################################
// Function Has_Angular_Constraint
//#####################################################################
template<class TV> bool POINT_JOINT<TV>::
Has_Angular_Constraint() const
{
    return constrain.Contains(true);
}
//#####################################################################
// Function Constrain_Angles
//#####################################################################
template<class TV> void POINT_JOINT<TV>::
Constrain_Angles(T_SPIN& angles_input) const
{
    for(int i=1;i<=T_SPIN::dimension;i++) angles_input(i)=wrap(angles_input(i),(T)-pi,(T)pi);
    angles_input=rotation_limits.Clamp(angles_input);
}
//#####################################################################
// Function Update_State_From_Joint_Frame
//#####################################################################
template<class TV> void POINT_JOINT<TV>::
Update_State_From_Joint_Frame(const bool enforce_constraints)
{
    if(enforce_constraints){
        T_SPIN angles=J.r.Euler_Angles();
        Constrain_Angles(angles);
        J.r=ROTATION<TV>::From_Euler_Angles(angles);
        Constrain_Prismatically(J.t);}
    J_inverse=J.Inverse();
}
//#####################################################################
// Function Angular_Constraint_Matrix
//#####################################################################
template<class TV> void POINT_JOINT<TV>::
Angular_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& constrained_matrix,MATRIX_MXN<T>* unconstrained_matrix) const
{
    Constraint_Matrix_Helper(parent_frame.r*F_pj().r,constrained_matrix,unconstrained_matrix,T_SPIN::Componentwise_Greater_Equal(rotation_limits.min_corner,
        rotation_limits.max_corner));
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d,s) \
    template POINT_JOINT<VECTOR<T,d> >::~POINT_JOINT(); \
    template void POINT_JOINT<VECTOR<T,d> >::Constrain_Relative_Linear_Velocity(const FRAME<VECTOR<T,d> >& parent_frame,VECTOR<T,d>& relative_linear_velocity) const; \
    template VECTOR<T,d>  POINT_JOINT<VECTOR<T,d> >::Prismatic_Component_Translation() const; \
    template void POINT_JOINT<VECTOR<T,d> >::Constrain_Prismatically(VECTOR<T,d>& translation) const; \
    template bool POINT_JOINT<VECTOR<T,d> >::Has_Angular_Constraint() const; \
    template void POINT_JOINT<VECTOR<T,d> >::Constrain_Angles(VECTOR<T,s>& angles_input) const; \
    template void POINT_JOINT<VECTOR<T,d> >::Update_State_From_Joint_Frame(const bool enforce_constraints); \
    template void POINT_JOINT<VECTOR<T,d> >::Angular_Constraint_Matrix(const FRAME<VECTOR<T,d> >& parent_frame,MATRIX_MXN<T>& constrained_matrix, \
        MATRIX_MXN<T>* unconstrained_matrix) const;

INSTANTIATION_HELPER(float,2,1)
INSTANTIATION_HELPER(float,3,3)
template class POINT_JOINT<VECTOR<float,1> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,2,1)
INSTANTIATION_HELPER(double,3,3)
template class POINT_JOINT<VECTOR<double,1> >;
#endif

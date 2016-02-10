//#####################################################################
// Copyright 2007, Craig Schroeder, Tamar Shinar, Jonathan Su
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/PRISMATIC_TWIST_JOINT.h>
using namespace PhysBAM;
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PRISMATIC_TWIST_JOINT<TV>::
~PRISMATIC_TWIST_JOINT()
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_input> PRISMATIC_TWIST_JOINT<VECTOR<T_input,1> >::
~PRISMATIC_TWIST_JOINT()
{}
template<class TV> bool PRISMATIC_TWIST_JOINT<TV>::
Has_Prismatic_Constraint() const
{
    return constrain.Number_True()>0;
}
template<class TV> void PRISMATIC_TWIST_JOINT<TV>::
Constrain_Prismatically(TV& translation) const
{
    for(int i=1;i<=TV::dimension;i++) if(constrain(i)) translation(i)=clamp(translation(i),prismatic_min(i),prismatic_max(i));
}
template<class TV> void PRISMATIC_TWIST_JOINT<TV>::
Constrain_Relative_Linear_Velocity(const FRAME<TV>& parent_frame,TV& relative_linear_velocity) const
{
    ROTATION<TV> joint_orientation=parent_frame.r*F_pj().r;
    TV_BOOL equality_constraint=Equality_Constraint();
    for(int i=1;i<=TV::dimension;i++) if(equality_constraint(i)){
        TV u=joint_orientation.Rotated_Axis(i);
        relative_linear_velocity-=TV::Dot_Product(relative_linear_velocity,u)*u;}
}
template<class TV> TV PRISMATIC_TWIST_JOINT<TV>::
Prismatic_Component_Translation() const
{
    TV current_translation=J.t;
    Constrain_Prismatically(current_translation);
    return current_translation;
}
template<class TV> void PRISMATIC_TWIST_JOINT<TV>::
Prismatic_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& constrained_matrix,MATRIX_MXN<T>* unconstrained_matrix) const
{
    Constraint_Matrix_Helper(parent_frame.r*F_pj().r,constrained_matrix,unconstrained_matrix,Equality_Constraint());
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d,s) \
    template PRISMATIC_TWIST_JOINT<VECTOR<T,d> >::~PRISMATIC_TWIST_JOINT(); \
    template bool PRISMATIC_TWIST_JOINT<VECTOR<T,d> >::Has_Prismatic_Constraint() const; \
    template void PRISMATIC_TWIST_JOINT<VECTOR<T,d> >::Constrain_Prismatically(VECTOR<T,d>& translation) const; \
    template void PRISMATIC_TWIST_JOINT<VECTOR<T,d> >::Prismatic_Constraint_Matrix(const FRAME<VECTOR<T,d> >& parent_frame,MATRIX_MXN<T>& constrained_matrix, \
        MATRIX_MXN<T>* unconstrained_matrix) const; \
    template void PRISMATIC_TWIST_JOINT<VECTOR<T,d> >::Constrain_Relative_Linear_Velocity(const FRAME<VECTOR<T,d> >& parent_frame,VECTOR<T,d>& relative_linear_velocity) const; \
    template VECTOR<T,d> PRISMATIC_TWIST_JOINT<VECTOR<T,d> >::Prismatic_Component_Translation() const;

INSTANTIATION_HELPER(float,2,1)
INSTANTIATION_HELPER(float,3,3)
template class PRISMATIC_TWIST_JOINT<VECTOR<float,1> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,2,1)
INSTANTIATION_HELPER(double,3,3)
template class PRISMATIC_TWIST_JOINT<VECTOR<double,1> >;
#endif

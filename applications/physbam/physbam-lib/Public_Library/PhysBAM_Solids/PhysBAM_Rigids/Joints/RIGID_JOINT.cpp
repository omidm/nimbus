//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Eran Guendelman, Michael Lentine, Craig Schroeder, Tamar Shinar, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Matrices/MATRIX_0X0.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/RIGID_JOINT.h>
using namespace PhysBAM;
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_JOINT<TV>::
~RIGID_JOINT()
{}
//#####################################################################
// Function Has_Angular_Constraint
//#####################################################################
template<class TV> bool RIGID_JOINT<TV>::
Has_Angular_Constraint() const
{
    return true;
}
//#####################################################################
// Function Constrain_Prismatically
//#####################################################################
template<class TV> void RIGID_JOINT<TV>::
Constrain_Prismatically(TV& translation) const
{
    translation=prismatic_translation;
}
//#####################################################################
// Function Constrain_Relative_Angular_Velocity
//#####################################################################
template<class TV> void RIGID_JOINT<TV>::
Constrain_Relative_Angular_Velocity(const FRAME<TV>& parent_frame,T_SPIN& relative_angular_velocity) const
{
    relative_angular_velocity=T_SPIN();
}
//#####################################################################
// Function Constrain_Relative_Linear_Velocity
//#####################################################################
template<class TV> void RIGID_JOINT<TV>::
Constrain_Relative_Linear_Velocity(const FRAME<TV>& parent_frame,TV& relative_linear_velocity) const
{
    relative_linear_velocity=TV();
}
//#####################################################################
// Function Prismatic_Component_Translation
//#####################################################################
template<class TV> TV RIGID_JOINT<TV>::
Prismatic_Component_Translation() const
{
    return prismatic_translation;
}
//#####################################################################
// Function Constrain_Angles
//#####################################################################
template<class TV> void RIGID_JOINT<TV>::
Constrain_Angles(T_SPIN& angles) const
{
    angles=T_SPIN();
}
//#####################################################################
// Function Update_State_From_Joint_Frame
//#####################################################################
template<class TV> void RIGID_JOINT<TV>::
Update_State_From_Joint_Frame(const bool enforce_constraints)
{
    if(TV::m==2){
        if(enforce_constraints){
            T_SPIN rotation_angle=J.r.Euler_Angles();
            Constrain_Angles(rotation_angle);
            J.r=ROTATION<TV>::From_Rotation_Vector(rotation_angle);
            Constrain_Prismatically(J.t);}
        J_inverse=J.Inverse();}
}
//#####################################################################
// Function Angular_Constraint_Matrix
//#####################################################################
template<class TV> void RIGID_JOINT<TV>::
Angular_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& angular_constraint_matrix,MATRIX_MXN<T>* angular_unconstrained_matrix) const
{
    angular_constraint_matrix=MATRIX<T,T_SPIN::dimension>::Identity_Matrix();
    if(angular_unconstrained_matrix) angular_unconstrained_matrix->Resize(T_SPIN::dimension,0);
}
//#####################################################################

template class RIGID_JOINT<VECTOR<float,1> >;
template class RIGID_JOINT<VECTOR<float,2> >;
template class RIGID_JOINT<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_JOINT<VECTOR<double,1> >;
template class RIGID_JOINT<VECTOR<double,2> >;
template class RIGID_JOINT<VECTOR<double,3> >;
#endif

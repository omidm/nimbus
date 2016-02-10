//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Eran Guendelman, Michael Lentine, Craig Schroeder, Tamar Shinar, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> JOINT<TV>::
JOINT()
    :joint_function(0),global_post_stabilization(true),primary_point_of_bend_joint(false),secondary_point_of_bend_joint(false),angular_damping(0),impulse_accumulator(0)
{
    for(int i=1;i<=dof;i++) control_dof(i)=false;
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> JOINT<TV>::
~JOINT()
{}
//#####################################################################
// Function Joint_Error
//#####################################################################
template<class TV> FRAME<TV> JOINT<TV>::
Joint_Error() const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Compute_Current_Joint_Frame
//#####################################################################
template<class TV> FRAME<TV> JOINT<TV>::
Compute_Current_Joint_Frame(const RIGID_BODY<TV>& parent,const RIGID_BODY<TV>& child) const
{
    return frame_jp*parent.Frame().Inverse()*child.Frame()*frame_cj;
}
//#####################################################################
// Function F_pc
//#####################################################################
template<class TV> FRAME<TV> JOINT<TV>::
F_pc() const
{
    return FRAME<TV>(frame_pj*J*frame_jc);
}
//#####################################################################
// Function F_cp
//#####################################################################
template<class TV> FRAME<TV> JOINT<TV>::
F_cp() const
{
    return FRAME<TV>(frame_cj*J_inverse*frame_jp);
}
//#####################################################################
// Function F_cj
//#####################################################################
template<class TV> FRAME<TV> JOINT<TV>::
F_cj() const
{
    return frame_cj;
}
//#####################################################################
// Function F_jc
//#####################################################################
template<class TV> FRAME<TV> JOINT<TV>::
F_jc() const
{
    return frame_jc;
}
//#####################################################################
// Function F_pj
//#####################################################################
template<class TV> FRAME<TV> JOINT<TV>::
F_pj() const
{
    return frame_pj;
}
//#####################################################################
// Function F_jp
//#####################################################################
template<class TV> FRAME<TV> JOINT<TV>::
F_jp() const
{
    return frame_jp;
}
//#####################################################################
// Function Set_Parent_To_Joint_Frame
//#####################################################################
template<class TV> void JOINT<TV>::
Set_Parent_To_Joint_Frame(const FRAME<TV>& F_jp_new)
{
    frame_jp=F_jp_new;frame_pj=frame_jp.Inverse();
}
//#####################################################################
// Function Set_Joint_To_Parent_Frame
//#####################################################################
template<class TV> void JOINT<TV>::
Set_Joint_To_Parent_Frame(const FRAME<TV>& F_pj_new)
{
    frame_pj=F_pj_new;frame_jp=frame_pj.Inverse();
}
//#####################################################################
// Function Set_Child_To_Joint_Frame
//#####################################################################
template<class TV> void JOINT<TV>::
Set_Child_To_Joint_Frame(const  FRAME<TV>& F_jc_new)
{
    frame_jc=F_jc_new;frame_cj=frame_jc.Inverse();
}
//#####################################################################
// Function Set_Joint_To_Child_Frame
//#####################################################################
template<class TV> void JOINT<TV>::
Set_Joint_To_Child_Frame(const FRAME<TV>& F_cj_new)
{
    frame_cj=F_cj_new;frame_jc=frame_cj.Inverse();
}
//#####################################################################
// Function Has_Angular_Constraint
//#####################################################################
template<class TV> bool JOINT<TV>::
Has_Angular_Constraint() const
{
    return false;
}
//#####################################################################
// Function Has_Prismatic_Constraint
//#####################################################################
template<class TV> bool JOINT<TV>::
Has_Prismatic_Constraint() const
{
    return true;
}
//#####################################################################
// Function Prismatic_Component_Translation
//#####################################################################
template<class TV> TV JOINT<TV>::
Prismatic_Component_Translation() const
{
    return TV();
}
//#####################################################################
// Function Update_State_From_Joint_Frame
//#####################################################################
template<class TV> void JOINT<TV>::
Update_State_From_Joint_Frame(const bool enforce_constraints)
{}
//#####################################################################
// Function Constrain_Relative_Angular_Velocity
//#####################################################################
template<class TV> void JOINT<TV>::
Constrain_Relative_Angular_Velocity(const FRAME<TV>& parent_frame,T_SPIN& relative_angular_velocity) const
{}
//#####################################################################
// Function Constrain_Relative_Linear_Velocity
//#####################################################################
template<class TV> void JOINT<TV>::
Constrain_Relative_Linear_Velocity(const FRAME<TV>& parent_frame,TV& relative_linear_velocity) const
{}
//#####################################################################
// Function Constrain_Prismatically
//#####################################################################
template<class TV> void JOINT<TV>::
Constrain_Prismatically(TV& translation) const
{}
//#####################################################################
// Function Constrain_Angles
//#####################################################################
template<class TV> void JOINT<TV>::
Constrain_Angles(T_SPIN& angles) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Prismatic_Constraint_Matrix
//#####################################################################
template<class TV> void JOINT<TV>::
Prismatic_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& constrained_matrix,MATRIX_MXN<T>* unconstrained_matrix) const
{
    constrained_matrix=MATRIX<T,TV::dimension>::Identity_Matrix();if(unconstrained_matrix) unconstrained_matrix->Resize(TV::dimension,0);
}
//#####################################################################
// Function Angular_Constraint_Matrix
//#####################################################################
template<class TV> void JOINT<TV>::
Angular_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& angular_constraint_matrix,MATRIX_MXN<T>* angular_unconstrained_matrix) const
{
    if(TV::dimension==3) PHYSBAM_FATAL_ERROR(); // must be overridden for 3D
}
//#####################################################################
// Function Constraint_Matrix_Helper
//#####################################################################
template<class TV> template<int d> void JOINT<TV>::
Constraint_Matrix_Helper(const ROTATION<TV>& orientation,MATRIX_MXN<T>& constrained_matrix,MATRIX_MXN<T>* unconstrained_matrix,const VECTOR<bool,d>& constrain) const
{
    MATRIX<T,d> R(orientation.Rotation_Matrix());
    int constrained_axes=constrain.Number_True();
    constrained_matrix.Resize(d,constrained_axes);
    for(int i=1,k=1;i<=d;i++) if(constrain(i)) constrained_matrix.Set_Column(k++,R.Column(i));
    if(unconstrained_matrix){
        unconstrained_matrix->Resize(d,d-constrained_axes);
        for(int i=1,k=1;i<=d;i++) if(!constrain(i)) unconstrained_matrix->Set_Column(k++,R.Column(i));}
}
//#####################################################################
// Function Constraint_Matrix_Helper
//#####################################################################
template<class TV> void JOINT<TV>::
Constraint_Matrix_Helper(const ROTATION<TV>& orientation,MATRIX_MXN<T>& constrained_matrix,MATRIX_MXN<T>* unconstrained_matrix,const VECTOR<bool,1>& constrain) const
{
    if(constrain.x){
        constrained_matrix=MATRIX<T,1>::Identity_Matrix();
        if(unconstrained_matrix) unconstrained_matrix->Resize(1,0);}
    else{
        constrained_matrix.Resize(1,0);
        if(unconstrained_matrix) *unconstrained_matrix=MATRIX<T,1>::Identity_Matrix();}
}
//#####################################################################
// Function Location
//#####################################################################
template<class TV> TV JOINT<TV>::
Location(const RIGID_BODY<TV>& parent,const RIGID_BODY<TV>& child) const
{
    return (T).5*(parent.World_Space_Point(frame_pj*Prismatic_Component_Translation())+child.World_Space_Point(frame_cj.t));
}
//#####################################################################
// Function Prismatic_Projection_Matrix
//#####################################################################
template<class TV> MATRIX<typename TV::SCALAR,TV::m> JOINT<TV>::
Prismatic_Projection_Matrix(const FRAME<TV>& parent_frame) const
{
    MATRIX_MXN<T> M;
    Prismatic_Constraint_Matrix(parent_frame,M);
    return MATRIX<T,TV::m>(M.Times_Transpose(M));
}
//#####################################################################
// Function Angular_Projection_Matrix
//#####################################################################
template<class TV> MATRIX<typename TV::SCALAR,JOINT<TV>::T_SPIN::m> JOINT<TV>::
Angular_Projection_Matrix(const FRAME<TV>& parent_frame) const
{
    MATRIX_MXN<T> M;
    Angular_Constraint_Matrix(parent_frame,M);
    return MATRIX<T,T_SPIN::m>(M.Times_Transpose(M));
}
//#####################################################################
template class JOINT<VECTOR<float,1> >;
template class JOINT<VECTOR<float,2> >;
template class JOINT<VECTOR<float,3> >;
template void JOINT<VECTOR<float,2> >::Constraint_Matrix_Helper<2>(ROTATION<VECTOR<float,2> > const&,MATRIX_MXN<float>&,MATRIX_MXN<float>*,VECTOR<bool,2> const&) const;
template void JOINT<VECTOR<float,3> >::Constraint_Matrix_Helper<3>(ROTATION<VECTOR<float,3> > const&,MATRIX_MXN<float>&,MATRIX_MXN<float>*,VECTOR<bool,3> const&) const;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class JOINT<VECTOR<double,1> >;
template class JOINT<VECTOR<double,2> >;
template class JOINT<VECTOR<double,3> >;
template void JOINT<VECTOR<double,2> >::Constraint_Matrix_Helper<2>(ROTATION<VECTOR<double,2> > const&,MATRIX_MXN<double>&,MATRIX_MXN<double>*,VECTOR<bool,2> const&) const;
template void JOINT<VECTOR<double,3> >::Constraint_Matrix_Helper<3>(ROTATION<VECTOR<double,3> > const&,MATRIX_MXN<double>&,MATRIX_MXN<double>*,VECTOR<bool,3> const&) const;
#endif

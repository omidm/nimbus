//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Eran Guendelman, Michael Lentine, Craig Schroeder, Tamar Shinar, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class JOINT
//#####################################################################
#ifndef __JOINT__
#define __JOINT__

#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_ID.h>
namespace PhysBAM{

template<class TV> class ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR;
template<class TV> class JOINT_FUNCTION;
template<class TV> class RIGID_BODY;
template<class TV>
class JOINT:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    typedef typename TV::SPIN T_SPIN;
    enum{dof=T_SPIN::dimension};

    FRAME<TV> frame_pj,frame_jp; // joint to parent and parent to joint, respectively, chosen so that the rest state is acheived with J=I
    FRAME<TV> frame_cj,frame_jc; // joint to child and child to joint, respectively, chosen to align the child's twist axis with the x-axis, and the bending axis with the z-axis
    FRAME<TV> J,J_inverse; // J=I is the rest state -- seem to only be used for kinematics
    FRAME<TV> internal_frame;
    JOINT_ID id_number;
    std::string name; // not saved to file - for convenience and debugging
    JOINT_FUNCTION<TV>* joint_function;
    bool global_post_stabilization;
    bool primary_point_of_bend_joint,secondary_point_of_bend_joint;
    T angular_damping;
    TWIST<TV> current_impulse;
    ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* impulse_accumulator;
    VECTOR<bool,dof> control_dof;

    JOINT();

    virtual ~JOINT();

    void Set_Joint_Frame(const FRAME<TV>& J_new,const bool enforce_constraints=true)
    {J=J_new;Update_State_From_Joint_Frame(enforce_constraints);}

    void Set_Joint_Inverse_Frame(const FRAME<TV>& J_inverse_new)
    {J_inverse=J_inverse_new;J=J_inverse.Inverse();Update_State_From_Joint_Frame();}

    void Set_Name(const std::string& name_input)
    {name=name_input;}

    void Constrain_Relative_Twist(const FRAME<TV>& parent_frame,TWIST<TV>& relative_twist) const
    {Constrain_Relative_Angular_Velocity(parent_frame,relative_twist.angular);Constrain_Relative_Linear_Velocity(parent_frame,relative_twist.linear);}

    void Set_Id_Number(const JOINT_ID id)
    {id_number=id;}

    void Set_Joint_Function(JOINT_FUNCTION<TV>* joint_function_input)
    {joint_function=joint_function_input;}

//#####################################################################
    FRAME<TV> Joint_Error() const;
    virtual VECTOR<bool,T_SPIN::dimension> Angular_Constraints() const{PHYSBAM_FATAL_ERROR();};
    TV Location(const RIGID_BODY<TV>& parent,const RIGID_BODY<TV>& child) const;
    FRAME<TV> Compute_Current_Joint_Frame(const RIGID_BODY<TV>& parent,const RIGID_BODY<TV>& child) const;
    FRAME<TV> F_pc() const;
    FRAME<TV> F_cp() const;
    FRAME<TV> F_cj() const;
    FRAME<TV> F_jc() const;
    FRAME<TV> F_pj() const;
    FRAME<TV> F_jp() const;
    virtual void Set_Parent_To_Joint_Frame(const FRAME<TV>& F_jp_new);
    virtual void Set_Joint_To_Parent_Frame(const FRAME<TV>& F_pj_new);
    virtual void Set_Child_To_Joint_Frame(const  FRAME<TV>& F_jc_new);
    virtual void Set_Joint_To_Child_Frame(const FRAME<TV>& F_cj_new);
    virtual bool Has_Angular_Constraint() const;
    virtual bool Has_Prismatic_Constraint() const;
    virtual TV Prismatic_Component_Translation() const;
    virtual void Update_State_From_Joint_Frame(const bool enforce_constraints=true);
    virtual void Constrain_Relative_Angular_Velocity(const FRAME<TV>& parent_frame,T_SPIN& relative_angular_velocity) const;
    virtual void Constrain_Relative_Linear_Velocity(const FRAME<TV>& parent_frame,TV& relative_linear_velocity) const;
    virtual void Constrain_Prismatically(TV& translation) const;
    virtual void Constrain_Angles(T_SPIN& angles) const;
    virtual void Prismatic_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& constrained_matrix,MATRIX_MXN<T>* unconstrained_matrix=0) const;
    virtual void Angular_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& angular_constraint_matrix,MATRIX_MXN<T>* angular_unconstrained_matrix=0) const;
    MATRIX<T,TV::m> Prismatic_Projection_Matrix(const FRAME<TV>& parent_frame) const;
    MATRIX<T,T_SPIN::m> Angular_Projection_Matrix(const FRAME<TV>& parent_frame) const;
protected:
    template<int d>
    void Constraint_Matrix_Helper(const ROTATION<TV>& orientation,MATRIX_MXN<T>& constrained_matrix,MATRIX_MXN<T>* unconstrained_matrix,const VECTOR<bool,d>& constrain) const;
    void Constraint_Matrix_Helper(const ROTATION<TV>& orientation,MATRIX_MXN<T>& constrained_matrix,MATRIX_MXN<T>* unconstrained_matrix,const VECTOR<bool,1>& constrain) const;
//#####################################################################
};
}
#include <PhysBAM_Solids/PhysBAM_Rigids/Read_Write/Joints/READ_WRITE_JOINT.h>
#endif

//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Michael Lentine, Craig Schroeder, Tamar Shinar, Jonathan Su, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POINT_JOINT
//#####################################################################
#ifndef __POINT_JOINT__
#define __POINT_JOINT__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT.h>
namespace PhysBAM{

template<class TV>
class POINT_JOINT:public JOINT<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    enum{dof=T_SPIN::dimension};
    typedef VECTOR<bool,dof> BOOL_VECTOR;
public:
    using JOINT<TV>::J;using JOINT<TV>::J_inverse;using JOINT<TV>::F_pj;

    bool prismatic_component;
    TV prismatic_translation;
    RANGE<VECTOR<T,dof> > rotation_limits;
    BOOL_VECTOR constrain;

    POINT_JOINT()
        :prismatic_component(false),rotation_limits(rotation_limits.Full_Box())
    {}

    virtual ~POINT_JOINT();

    void Set_Prismatic_Component_Translation(const TV& translation)
    {prismatic_translation=translation;prismatic_component=true;}

    void Use_Rotation_Constraint(const T min_input=-pi,const T max_input=pi)
    {STATIC_ASSERT(dof==1);constrain.x=true;rotation_limits.min_corner.x=min_input;rotation_limits.max_corner.x=max_input;}

    void Use_Twist_Constraint(const T min_input=-pi,const T max_input=pi)
    {STATIC_ASSERT(dof==3);constrain.x=true;rotation_limits.min_corner.x=min_input;rotation_limits.max_corner.x=max_input;}

    void Use_Phi_Constraint(const T min_input=-pi*.5,const T max_input=pi*.5)
    {STATIC_ASSERT(dof==3);constrain.y=true;rotation_limits.min_corner.y=min_input;rotation_limits.max_corner.y=max_input;}

    void Use_Theta_Constraint(const T min_input=-pi,const T max_input=pi)
    {STATIC_ASSERT(dof==3);constrain.z=true;rotation_limits.min_corner.z=min_input;rotation_limits.max_corner.z=max_input;}

    void Remove_Rotation_Constraint()
    {STATIC_ASSERT(dof==1);constrain.x=false;rotation_limits.min_corner.x=-FLT_MAX;rotation_limits.max_corner.x=FLT_MAX;}

    void Remove_Twist_Constraint()
    {STATIC_ASSERT(dof==3);constrain.x=false;rotation_limits.min_corner.x=-FLT_MAX;rotation_limits.max_corner.x=FLT_MAX;}

    void Remove_Phi_Constraint()
    {STATIC_ASSERT(dof==3);constrain.y=false;rotation_limits.min_corner.y=-FLT_MAX;rotation_limits.max_corner.y=FLT_MAX;}

    void Remove_Theta_Constraint()
    {STATIC_ASSERT(dof==3);constrain.z=false;rotation_limits.min_corner.z=-FLT_MAX;rotation_limits.max_corner.z=FLT_MAX;}

    VECTOR<bool,T_SPIN::dimension> Angular_Constraints() const
    {return T_SPIN::Componentwise_Greater_Equal(rotation_limits.min_corner,rotation_limits.max_corner);}

//#####################################################################
    void Constrain_Relative_Linear_Velocity(const FRAME<TV>& parent_frame,TV& relative_linear_velocity) const PHYSBAM_OVERRIDE;
    TV Prismatic_Component_Translation() const PHYSBAM_OVERRIDE;
    void Constrain_Prismatically(TV& translation) const PHYSBAM_OVERRIDE;
    bool Has_Angular_Constraint() const PHYSBAM_OVERRIDE;
    void Constrain_Angles(T_SPIN& angles_input) const PHYSBAM_OVERRIDE;
    void Update_State_From_Joint_Frame(const bool enforce_constraints=true) PHYSBAM_OVERRIDE;
    void Angular_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& constrained_matrix,MATRIX_MXN<T>* unconstrained_matrix=0) const PHYSBAM_OVERRIDE;
//#####################################################################
};

template<class T_input>
class POINT_JOINT<VECTOR<T_input,1> >:public JOINT<VECTOR<T_input,1> >
{
    typedef T_input T;
public:
    POINT_JOINT()
    {
        PHYSBAM_FATAL_ERROR();
    }

    virtual ~POINT_JOINT();
};
}
#endif

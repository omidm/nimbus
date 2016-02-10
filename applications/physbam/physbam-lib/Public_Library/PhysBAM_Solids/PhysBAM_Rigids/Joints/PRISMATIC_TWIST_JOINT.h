//#####################################################################
// Copyright 2007, Craig Schroeder, Tamar Shinar, Jonathan Su
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PRISMATIC_TWIST_JOINT
//#####################################################################
#ifndef __PRISMATIC_TWIST_JOINT__
#define __PRISMATIC_TWIST_JOINT__

#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/ANGLE_JOINT.h>
namespace PhysBAM{

template<class TV>
class PRISMATIC_TWIST_JOINT:public ANGLE_JOINT<TV>
{
    typedef typename TV::SCALAR T;typedef ANGLE_JOINT<TV> BASE;typedef typename TV::SPIN T_SPIN;
    using BASE::J;using BASE::F_pj;typedef typename TV::template REBIND<bool>::TYPE TV_BOOL;

public:
    TV_BOOL constrain;
    TV prismatic_min,prismatic_max;

    PRISMATIC_TWIST_JOINT()
    {constrain.x=true;}

    virtual ~PRISMATIC_TWIST_JOINT();

    void Set_Prismatic_Constraints(const TV_BOOL& constrain_input,const TV& min=TV(),const TV& max=TV())
    {constrain=constrain_input;prismatic_min=min;prismatic_max=max;}

private:
    TV_BOOL Equality_Constraint() const
    {return TV::Componentwise_And(constrain,TV::Componentwise_Greater_Equal(prismatic_min,prismatic_max));}
public:
//#####################################################################
    bool Has_Prismatic_Constraint() const PHYSBAM_OVERRIDE;
    void Constrain_Prismatically(TV& translation) const PHYSBAM_OVERRIDE; // This is in joint space
    void Constrain_Relative_Linear_Velocity(const FRAME<TV>& parent_frame,TV& relative_linear_velocity) const PHYSBAM_OVERRIDE; // This is in world space.
    TV Prismatic_Component_Translation() const PHYSBAM_OVERRIDE;
    void Prismatic_Constraint_Matrix(const FRAME<TV>& parent_frame,MATRIX_MXN<T>& constrained_matrix,MATRIX_MXN<T>* unconstrained_matrix=0) const PHYSBAM_OVERRIDE;
//#####################################################################
};

template<class T_input>
class PRISMATIC_TWIST_JOINT<VECTOR<T_input,1> >:public JOINT<VECTOR<T_input,1> >
{
    typedef T_input T;typedef JOINT<VECTOR<T,1> > BASE;
public:
    PRISMATIC_TWIST_JOINT()
    {PHYSBAM_FATAL_ERROR();}

    virtual ~PRISMATIC_TWIST_JOINT();
};
}
#endif

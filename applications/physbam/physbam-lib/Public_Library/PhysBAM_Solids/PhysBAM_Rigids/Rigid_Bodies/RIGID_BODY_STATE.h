//#####################################################################
// Copyright 2003-2006, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_STATE
//#####################################################################
#ifndef __RIGID_BODY_STATE__
#define __RIGID_BODY_STATE__

#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_STATE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_MASS.h>
namespace PhysBAM{

template<class TV>
class RIGID_BODY_STATE:public RIGID_GEOMETRY_STATE<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    typedef RIGID_GEOMETRY_STATE<TV> BASE;
public:
    using BASE::frame;using BASE::twist;
    T_SPIN angular_momentum;

    RIGID_BODY_STATE()
        :BASE(),angular_momentum()
    {}

    RIGID_BODY_STATE(const FRAME<TV>& frame_input)
       :BASE(frame_input),angular_momentum()
    {}

    RIGID_BODY_STATE(const FRAME<TV>& frame_input,const TWIST<TV>& twist_input)
       :BASE(frame_input,twist_input),angular_momentum()
    {}

    template<class TV2> explicit RIGID_BODY_STATE(const RIGID_BODY_STATE<TV2>& state_input)
        :BASE(state_input),angular_momentum(state_input.angular_momentum)
    {}

    MATRIX<T,0> World_Space_Inertia_Tensor(const MATRIX<T,0> moment_of_inertia) const // relative to the center of mass
    {STATIC_ASSERT(TV::m==1);return moment_of_inertia;}

    MATRIX<T,0> World_Space_Inertia_Tensor_Inverse(const MATRIX<T,0> moment_of_inertia) const // relative to the center of mass
    {STATIC_ASSERT(TV::m==1);return moment_of_inertia.Inverse();}

    MATRIX<T,1> World_Space_Inertia_Tensor(const MATRIX<T,1> moment_of_inertia) const // relative to the center of mass
    {STATIC_ASSERT(TV::m==2);return moment_of_inertia;}

    SYMMETRIC_MATRIX<T,3> World_Space_Inertia_Tensor(const DIAGONAL_MATRIX<T,3>& inertia_tensor) const // relative to the center of mass
    {STATIC_ASSERT(TV::m==3);MATRIX<T,3> object_to_world_transformation=frame.r.Rotation_Matrix();
    return SYMMETRIC_MATRIX<T,3>::Conjugate(object_to_world_transformation,inertia_tensor);}

    RIGID_BODY_MASS<TV,true> World_Space_Rigid_Mass(const RIGID_BODY_MASS<TV>& rigid_mass) const
    {return RIGID_BODY_MASS<TV,true>(rigid_mass.mass,World_Space_Inertia_Tensor(rigid_mass.inertia_tensor));}

    MATRIX<T,1> World_Space_Inertia_Tensor_Inverse(const MATRIX<T,1>& moment_of_inertia) const // relative to the center of mass
    {STATIC_ASSERT(TV::m==2);return moment_of_inertia.Inverse();}

    SYMMETRIC_MATRIX<T,3> World_Space_Inertia_Tensor_Inverse(const DIAGONAL_MATRIX<T,3>& inertia_tensor) const // relative to the center of mass
    {STATIC_ASSERT(TV::m==3);MATRIX<T,3> object_to_world_transformation=frame.r.Rotation_Matrix();
    return SYMMETRIC_MATRIX<T,3>::Conjugate(object_to_world_transformation,inertia_tensor.Inverse());}

    RIGID_BODY_MASS<TV,true> World_Space_Rigid_Mass_Inverse(const RIGID_BODY_MASS<TV>& rigid_mass) const
    {return RIGID_BODY_MASS<TV,true>(1/rigid_mass.mass,World_Space_Inertia_Tensor_Inverse(rigid_mass.inertia_tensor));}

    void Update_Angular_Velocity(const MATRIX<T,0> inertia_tensor) // needs to be called to keep the angular velocity valid
    {STATIC_ASSERT(TV::m==1);}

    void Update_Angular_Velocity(const MATRIX<T,1> moment_of_inertia) // needs to be called to keep the angular velocity valid
    {STATIC_ASSERT(TV::m==2);twist.angular=moment_of_inertia.Solve_Linear_System(angular_momentum);}

    void Update_Angular_Velocity(const DIAGONAL_MATRIX<T,3>& inertia_tensor) // needs to be called to keep the angular velocity valid
    {STATIC_ASSERT(TV::m==3);twist.angular=World_Space_Vector(inertia_tensor.Solve_Linear_System(Object_Space_Vector(angular_momentum)));}

    void Update_Angular_Momentum(const MATRIX<T,0> inertia_tensor) // needs to be called to keep the angular velocity valid
    {STATIC_ASSERT(TV::m==1);}

    void Update_Angular_Momentum(const MATRIX<T,1> moment_of_inertia) // assumes a valid angular_velocity
    {STATIC_ASSERT(TV::m==2);angular_momentum=moment_of_inertia*twist.angular;}

    void Update_Angular_Momentum(const DIAGONAL_MATRIX<T,3>& inertia_tensor) // assumes a valid angular_velocity
    {STATIC_ASSERT(TV::m==3);angular_momentum=World_Space_Vector(inertia_tensor*Object_Space_Vector(twist.angular));}

//#####################################################################
};
}
//#include <PhysBAM_Solids/PhysBAM_Rigids/Read_Write/Rigid_Bodies/READ_WRITE_RIGID_BODY_STATE.h>
#endif

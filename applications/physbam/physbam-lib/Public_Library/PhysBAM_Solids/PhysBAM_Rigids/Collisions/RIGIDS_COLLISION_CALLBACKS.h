//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_COLLISION_CALLBACKS
//##################################################################### 
#ifndef __RIGIDS_COLLISION_CALLBACKS__
#define __RIGIDS_COLLISION_CALLBACKS__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>

namespace PhysBAM{

template<class TV> class FRAME;
template<class TV> class TWIST;
template<class TV> class RIGID_BODY;
template<class T,int d> class VECTOR;

template<class TV>
class RIGIDS_COLLISION_CALLBACKS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
public:

    RIGIDS_COLLISION_CALLBACKS() {}
    virtual ~RIGIDS_COLLISION_CALLBACKS() {};

    void Swap_States(const int id_1,const int id_2)
    {Swap_State(id_1);Swap_State(id_2);}

//#####################################################################
    virtual void Restore_Size(const int size) {}
    virtual void Reevolve_Body_With_Saved_State(const int p,const T dt,const T time)=0;
    virtual void Restore_Position(const int p)=0;
    virtual void Restore_Positions()=0;
    virtual void Save_Position(const int p)=0;
    virtual void Restore_Velocity(const int p)=0;
    virtual void Save_Velocity(const int p)=0;
    virtual void Euler_Step_Position_With_New_Velocity(const int id,const T dt,const T time)=0;
    virtual void Euler_Step_Position(const int id,const T dt,const T time)=0;
    virtual void Swap_State(const int id)=0;
    virtual FRAME<TV> Saved_Particle_To_Levelset_Body_Transform(const int levelset_body,const int particle_body)=0;
    virtual void Exchange_Frame(const int id)=0;
    virtual TWIST<TV> Compute_Collision_Impulse(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const TV& location,const TV& normal,const TV& relative_velocity,const T coefficient_of_restitution,
        const T coefficient_of_friction,const bool clamp_friction_magnitude,const bool rolling_friction,const bool clamp_energy)=0;
    virtual void Subtract_Stored_Difference(TV& velocity,T_SPIN& momentum,const int particle_index)=0;
    virtual void Begin_Fracture(const int body_id)=0;
    virtual void End_Fracture(const int body_id,ARRAY<int>& added_bodies)=0;
    virtual void Begin_Asymmetric_Collisions(const int body_1,const int body_2)=0;
    virtual void End_Asymmetric_Collisions(const int body_1,const int body_2,VECTOR<ARRAY<int>,2>& added_bodies)=0;
//#####################################################################
};
}
#endif

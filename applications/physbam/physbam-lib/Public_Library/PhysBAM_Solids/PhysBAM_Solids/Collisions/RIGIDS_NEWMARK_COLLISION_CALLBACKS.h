//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_NEWMARK_COLLISION_CALLBACKS
//##################################################################### 
#ifndef __RIGIDS_NEWMARK_COLLISION_CALLBACKS__
#define __RIGIDS_NEWMARK_COLLISION_CALLBACKS__

#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGIDS_COLLISION_CALLBACKS.h>

namespace PhysBAM{

template<class TV> class NEWMARK_EVOLUTION;
template<class TV> class FRAME;
template<class TV> class TWIST;
template<class TV> class RIGID_BODY;

template<class TV>
class RIGIDS_NEWMARK_COLLISION_CALLBACKS:public RIGIDS_COLLISION_CALLBACKS<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
public:
    NEWMARK_EVOLUTION<TV>& evolution;
    TWIST<TV> old_stored_difference,asymmetric_collision_stored_difference;
    TV old_position;

    RIGIDS_NEWMARK_COLLISION_CALLBACKS(NEWMARK_EVOLUTION<TV>& evolution_input);
    virtual ~RIGIDS_NEWMARK_COLLISION_CALLBACKS();

//#####################################################################
    void Reevolve_Body_With_Saved_State(const int p,const T dt,const T time) PHYSBAM_OVERRIDE;
    void Restore_Positions() PHYSBAM_OVERRIDE;
    void Restore_Position(const int p) PHYSBAM_OVERRIDE;
    void Save_Position(const int p) PHYSBAM_OVERRIDE;
    void Restore_Velocity(const int p) PHYSBAM_OVERRIDE;
    void Save_Velocity(const int p) PHYSBAM_OVERRIDE;
    void Euler_Step_Position(const int id,const T dt,const T time);
    void Euler_Step_Position_With_New_Velocity(const int id,const T dt,const T time) PHYSBAM_OVERRIDE;
    void Swap_State(const int id) PHYSBAM_OVERRIDE;
    FRAME<TV> Saved_Particle_To_Levelset_Body_Transform(const int levelset_body,const int particle_body) PHYSBAM_OVERRIDE;
    void Exchange_Frame(const int id) PHYSBAM_OVERRIDE;
    TWIST<TV> Compute_Collision_Impulse(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const TV& location,const TV& normal,const TV& relative_velocity,const T coefficient_of_restitution,
        const T coefficient_of_friction,const bool clamp_friction_magnitude,const bool rolling_friction,const bool clamp_energy) PHYSBAM_OVERRIDE;
    void Subtract_Stored_Difference(TV& velocity,T_SPIN& momentum,const int particle_index) PHYSBAM_OVERRIDE;
    void Begin_Fracture(const int body_id) PHYSBAM_OVERRIDE;
    void End_Fracture(const int body_id,ARRAY<int>& added_bodies) PHYSBAM_OVERRIDE;
    void Begin_Asymmetric_Collisions(const int body_1,const int body_2) PHYSBAM_OVERRIDE;
    void End_Asymmetric_Collisions(const int body_1,const int body_2,VECTOR<ARRAY<int>,2>& added_bodies) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

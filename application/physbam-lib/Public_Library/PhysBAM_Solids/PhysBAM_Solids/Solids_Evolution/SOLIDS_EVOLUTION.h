//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_EVOLUTION
//#####################################################################
#ifndef __SOLIDS_EVOLUTION__
#define __SOLIDS_EVOLUTION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD_FORWARD.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_OBJECT_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigids_Evolution/RIGIDS_KINEMATIC_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION_CALLBACKS.h>
namespace PhysBAM{

template<class T_GRID> class RIGID_BODY_COLLISIONS;
template<class TV,bool world_space> class RIGID_BODY_MASS;
template<class TV> class RIGID_BODY_STATE;
template<class TV> class SOLID_BODY_COLLECTION;
template<class TV> class GRID;

template<class TV>
class SOLIDS_EVOLUTION:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    enum WORKAROUND {d=TV::m};
public:
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    SOLIDS_PARAMETERS<TV>& solids_parameters;
    RIGID_BODY_COLLISIONS<TV>* rigid_body_collisions;
    RIGID_DEFORMABLE_COLLISIONS<TV>* rigid_deformable_collisions;
    T time;
    SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callbacks;
protected:
    ARRAY<TV> F_full,R_full,B_full,S_full,AR_full;
    ARRAY<TWIST<TV> > rigid_F_full,rigid_R_full,rigid_S_full,rigid_B_full,rigid_AR_full;

    T energy_damped;
    T energy_lost_to_bodies;
    TV momentum_lost_to_bodies;
public:
    ARRAY<RIGID_BODY_MASS<TV,true> > world_space_rigid_mass;
    ARRAY<RIGID_BODY_MASS<TV,true> > world_space_rigid_mass_inverse;
private:
    static SOLIDS_EVOLUTION_CALLBACKS<TV> solids_evolution_callbacks_default;
public:
    RIGIDS_KINEMATIC_EVOLUTION<TV> kinematic_evolution;

    SOLIDS_EVOLUTION(SOLIDS_PARAMETERS<TV>& solids_parameters_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input);
    virtual ~SOLIDS_EVOLUTION();

    void Set_Solids_Evolution_Callbacks(SOLIDS_EVOLUTION_CALLBACKS<TV>& solids_evolution_callbacks_input)
    {solids_evolution_callbacks=&solids_evolution_callbacks_input;}

//#####################################################################
    virtual bool Use_CFL() const=0;
    virtual void Euler_Step_Position(const T dt,const T time);
    virtual void Advance_One_Time_Step_Position(const T dt,const T time,const bool solids)=0;
    virtual void Advance_One_Time_Step_Velocity(const T dt,const T time,const bool solids)=0;
    virtual void Initialize_Deformable_Objects(const T frame_rate,const bool restart);
    virtual void Initialize_Rigid_Bodies(const T frame_rate, const bool restart)=0;
    virtual bool Adjust_Velocity_For_Self_Repulsion_And_Self_Collisions(const T dt,const T time,int& repulsions_found,int& collisions_found,const bool exit_early=false);
    virtual void Postprocess_Frame(const int frame);
protected:
    void Save_Position(ARRAY<TV>& X,ARRAY<TV>& rigid_X,ARRAY<ROTATION<TV> >& rigid_rotation);
    void Restore_Position(ARRAY_VIEW<const TV> X,ARRAY_VIEW<const TV> rigid_X,ARRAY_VIEW<const ROTATION<TV> > rigid_Rotation);
    void Clamp_Velocities();
public:
    void Initialize_World_Space_Masses();
    void Restore_Position_After_Hypothetical_Position_Evolution(ARRAY<TV>& X_save,ARRAY<TV>& rigid_X_save,ARRAY<ROTATION<TV> >& rigid_rotation_save);
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time);
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time);
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time);
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time);
    void Euler_Step_Position(const T dt,const T time,const int p);
//#####################################################################
};
}
#endif

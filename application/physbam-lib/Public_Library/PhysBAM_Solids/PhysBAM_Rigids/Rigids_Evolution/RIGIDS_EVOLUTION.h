//#####################################################################
// Copyright 2006-2008, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_EVOLUTION
//#####################################################################
#ifndef __RIGIDS_EVOLUTION__
#define __RIGIDS_EVOLUTION__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Parallel_Computation/MPI_RIGIDS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigids_Evolution/RIGIDS_KINEMATIC_EVOLUTION.h>
namespace PhysBAM{

template<class TV> class RIGIDS_BACKWARD_EULER_SYSTEM;
template<class TV> class RIGIDS_ONLY_COLLISION_CALLBACKS;
template<class TV> class RIGIDS_EVOLUTION_CALLBACKS;
template<class TV> class RIGIDS_PARAMETERS;

template<class TV>
class RIGIDS_EVOLUTION
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    enum WORKAROUND {d=TV::m};
public:
    T time;
    RIGIDS_ONLY_COLLISION_CALLBACKS<TV>* rigids_collision_callbacks;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    RIGIDS_KINEMATIC_EVOLUTION<TV> kinematic_evolution;
    ARRAY<TWIST<TV> > rigid_velocity_save;
    ARRAY<T_SPIN> rigid_angular_momentum_save;
    ARRAY<TV> rigid_X_save;
    ARRAY<ROTATION<TV> > rigid_rotation_save;
    ARRAY<TV> rigid_velocity_difference;
    ARRAY<T_SPIN> rigid_angular_momentum_difference;
    bool use_existing_contact;

    RIGIDS_PARAMETERS<TV>& rigids_parameters;
    RIGID_BODY_COLLISIONS<TV>* rigid_body_collisions;
    ARRAY<TWIST<TV> > rigid_F_full,rigid_R_full,rigid_S_full,rigid_B_full,rigid_AR_full;

    ARRAY<RIGID_BODY_MASS<TV,true> > world_space_rigid_mass;
    ARRAY<RIGID_BODY_MASS<TV,true> > world_space_rigid_mass_inverse;
    RIGIDS_EVOLUTION_CALLBACKS<TV>* rigids_evolution_callbacks;
    MPI_RIGIDS<TV>* mpi_rigids;
private:
    ARRAY<TWIST<TV> > rigid_V_save,saved_pd;
    ARRAY<TV> rigid_X_save_for_constraints;
    ARRAY<ROTATION<TV> > rigid_rotation_save_for_constraints;
public:

    RIGIDS_EVOLUTION(RIGIDS_PARAMETERS<TV>& rigids_parameters_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input);
    ~RIGIDS_EVOLUTION();

    bool Use_CFL() const
    {return true;}

    void Set_Rigids_Evolution_Callbacks(RIGIDS_EVOLUTION_CALLBACKS<TV>& rigids_evolution_callbacks_input)
    {rigids_evolution_callbacks=&rigids_evolution_callbacks_input;}

//#####################################################################
    void Advance_One_Time_Step_Position(const T dt,const T time,const bool solids);
    void Advance_One_Time_Step_Velocity(const T dt,const T time,const bool solids);
    void Backward_Euler_Step_Velocity_Helper(const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update);
    void Apply_Constraints(const T dt,const T time);
    void Diagnostics(const T dt,const T time,const int velocity_time,const int position_time,int step,const char* description);
    void Print_Maximum_Velocities(const T time) const;
    void Update_Velocity_Using_Stored_Differences(const T dt,const T time,const int p);
    void Initialize_Rigid_Bodies(const T frame_rate, const bool restart);
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time);
    void Euler_Step_Position(const T dt,const T time,const int p);
    void Euler_Step_Position(const T dt,const T time);
    void Initialize_World_Space_Masses();
    void Clamp_Velocities();
    void CFL(const bool verbose);
    void Restore_Position(ARRAY_VIEW<const TV> rigid_X,ARRAY_VIEW<const ROTATION<TV> > rigid_rotation);
protected:
    void Average_And_Exchange_Position();
    void Trapezoidal_Step_Velocity(const T dt,const T time);
    void Prepare_Backward_Euler_System(RIGIDS_BACKWARD_EULER_SYSTEM<TV>& system,const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update);
    void Finish_Backward_Euler_Step(const T dt,const T current_position_time,const bool velocity_update);
    void Backward_Euler_Step_Velocity(const T dt,const T time);
    void Update_Positions_And_Apply_Contact_Forces(const T dt,const T time,const bool use_saved_pairs=false);
    void Update_Velocity_Using_Stored_Differences(const T dt,const T time);
    void Compute_Momentum_Differences();
    void Save_Velocity();
    void Save_Velocity(ARRAY<TWIST<TV> >& rigid_velocity_save,ARRAY<T_SPIN>& rigid_angular_momentum_save);
    void Restore_Velocity();
    void Restore_Velocity(ARRAY<TWIST<TV> >& rigid_velocity_save,ARRAY<T_SPIN>& rigid_angular_momentum_save);
    void Save_Position(ARRAY<TV>& rigid_X,ARRAY<ROTATION<TV> >& rigid_rotation);
    void Exchange_Velocity();
//#####################################################################
};
}
#endif

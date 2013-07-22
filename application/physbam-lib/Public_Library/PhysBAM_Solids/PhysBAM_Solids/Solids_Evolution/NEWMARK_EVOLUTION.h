//#####################################################################
// Copyright 2006-2008, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NEWMARK_EVOLUTION
//#####################################################################
#ifndef __NEWMARK_EVOLUTION__
#define __NEWMARK_EVOLUTION__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
namespace PhysBAM{

template<class T> class TRIANGLE_REPULSIONS;
template<class TV> class COLLISION_PARTICLE_STATE;
template<class TV> class BACKWARD_EULER_SYSTEM;
template<class TV> class ASYNCHRONOUS_EVOLUTION;
template<class TV> class RIGIDS_NEWMARK_COLLISION_CALLBACKS;
template<class T> class KRYLOV_SYSTEM_BASE;

template<class TV>
class NEWMARK_EVOLUTION:public SOLIDS_EVOLUTION<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    typedef SOLIDS_EVOLUTION<TV> BASE;
    enum WORKAROUND {d=TV::m};
protected:
    using BASE::F_full;using BASE::R_full;using BASE::S_full;using BASE::B_full;using BASE::rigid_F_full;using BASE::rigid_R_full;using BASE::rigid_S_full;
    using BASE::rigid_B_full;using BASE::rigid_AR_full;using BASE::AR_full;using BASE::Save_Position;using BASE::Restore_Position;
    using BASE::rigid_deformable_collisions;using BASE::world_space_rigid_mass;using BASE::world_space_rigid_mass_inverse;using BASE::time;
    using BASE::energy_damped;using BASE::energy_lost_to_bodies;using BASE::momentum_lost_to_bodies;
public:
    using BASE::solid_body_collection;using BASE::Euler_Step_Position;
    using BASE::Clamp_Velocities;using BASE::solids_evolution_callbacks;using BASE::rigid_body_collisions;using BASE::Initialize_World_Space_Masses;
    using BASE::solids_parameters;using BASE::kinematic_evolution;

    RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>& rigids_evolution_callbacks;
    TRIANGLE_REPULSIONS<TV>* repulsions;
    ARRAY<TV> X_save,V_save,V_original_save,X_original_save;
    ARRAY<TWIST<TV> > rigid_velocity_save;
    ARRAY<T_SPIN> rigid_angular_momentum_save;
    ARRAY<TV> rigid_X_save;
    ARRAY<ROTATION<TV> > rigid_rotation_save;
    ARRAY<TV> V_difference,rigid_velocity_difference;
    ARRAY<T_SPIN> rigid_angular_momentum_difference;
    bool use_existing_contact;
    ASYNCHRONOUS_EVOLUTION<TV>* asynchronous_evolution;
    bool print_matrix;

    T postelasticforces_PE;
    ARRAY<TV> V_postelasticforces;
    T postelasticforces_KE;

protected:
    ARRAY<TWIST<TV> > rigid_V_save,saved_pd;
    ARRAY<TV> X_save_for_constraints,rigid_X_save_for_constraints;
    ARRAY<ROTATION<TV> > rigid_rotation_save_for_constraints;
public:

    NEWMARK_EVOLUTION(SOLIDS_PARAMETERS<TV>& solids_parameters_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input);
    virtual ~NEWMARK_EVOLUTION();

//#####################################################################
    bool Use_CFL() const PHYSBAM_OVERRIDE;
    void Advance_One_Time_Step_Position(const T dt,const T time,const bool solids) PHYSBAM_OVERRIDE;
    void Advance_One_Time_Step_Velocity(const T dt,const T time,const bool solids) PHYSBAM_OVERRIDE;
    virtual void Backward_Euler_Step_Velocity_Helper(const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update);
    void Apply_Constraints(const T dt,const T time);
    void Diagnostics(const T dt,const T time,const int velocity_time,const int position_time,int step,const char* description);
    void Print_Maximum_Velocities(const T time) const;
    void Update_Velocity_Using_Stored_Differences(const T dt,const T time,const int p);
    void Initialize_Rigid_Bodies(const T frame_rate, const bool restart) PHYSBAM_OVERRIDE;
    void Write_Position_Update_Projection_Data(const STREAM_TYPE stream_type,const std::string& prefix);
    void Read_Position_Update_Projection_Data(const STREAM_TYPE stream_type,const std::string& prefix);
protected:
    void Average_And_Exchange_Position();
    virtual void Process_Collisions(const T dt,const T time,const bool advance_rigid_bodies);
    void Trapezoidal_Step_Velocity(const T dt,const T time);
    void Prepare_Backward_Euler_System(BACKWARD_EULER_SYSTEM<TV>& system,const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update);
    void Finish_Backward_Euler_Step(KRYLOV_SYSTEM_BASE<T>& system,const T dt,const T current_position_time,const bool velocity_update);
    void Backward_Euler_Step_Velocity(const T dt,const T time);
    void Make_Incompressible(const T dt,const bool correct_volume);
    void Update_Positions_And_Apply_Contact_Forces(const T dt,const T time,const bool use_saved_pairs=false);
    void Update_Velocity_Using_Stored_Differences(const T dt,const T time);
    void Compute_Momentum_Differences();
    void Save_Velocity();
    void Restore_Velocity() const;
    void Exchange_Velocity();
    void Apply_Projections_In_Position_Update(const T dt,const T time);
    void Set_Velocity_From_Positions_Position_Update(const T dt,const T time);
    void Set_Velocity_From_Positions_Velocity_Update(const T dt,const T time);
//#####################################################################
};
}
#endif

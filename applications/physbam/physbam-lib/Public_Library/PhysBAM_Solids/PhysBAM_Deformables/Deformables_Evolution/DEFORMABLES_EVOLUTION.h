//#####################################################################
// Copyright 2006-2008, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLES_EVOLUTION
//#####################################################################
#ifndef __DEFORMABLES_EVOLUTION__
#define __DEFORMABLES_EVOLUTION__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Geometry/Solids_Geometry_Evolution/KINEMATIC_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformables_Evolution/DEFORMABLES_EVOLUTION_CALLBACKS.h>
namespace PhysBAM{

template<class TV> class DEFORMABLES_BACKWARD_EULER_SYSTEM;
template<class TV> class DEFORMABLES_PARAMETERS;
template<class TV> class DEFORMABLE_BODY_COLLECTION;
template<class TV> class TRIANGLE_REPULSIONS;

template<class TV>
class DEFORMABLES_EVOLUTION
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    enum WORKAROUND {d=TV::m};
public:
    T time;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection;
    RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection;
    KINEMATIC_EVOLUTION<TV> kinematic_evolution;
    TRIANGLE_REPULSIONS<TV>* repulsions;
    ARRAY<TV> X_save,V_save,rigid_X_save,rigid_velocity_save;
    ARRAY<ROTATION<TV> > rigid_rotation_save;
    ARRAY<TV> V_difference,rigid_velocity_difference;

    DEFORMABLES_PARAMETERS<TV>& deformables_parameters;
    ARRAY<TV> F_full,R_full,B_full,S_full,AR_full;
    ARRAY<TWIST<TV> > rigid_F_full,rigid_R_full,rigid_S_full,rigid_B_full,rigid_AR_full;
    bool use_existing_contact;

    DEFORMABLES_EVOLUTION_CALLBACKS<TV>* deformables_evolution_callbacks;
private:
    ARRAY<TWIST<TV> > rigid_V_save;
    ARRAY<TV> X_save_for_constraints,rigid_X_save_for_constraints;    
    ARRAY<ROTATION<TV> > rigid_rotation_save_for_constraints;
public:

    DEFORMABLES_EVOLUTION(DEFORMABLES_PARAMETERS<TV>& deformables_parameters_input,DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection_input,RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection_input);
    ~DEFORMABLES_EVOLUTION();

    void Set_Deformables_Evolution_Callbacks(DEFORMABLES_EVOLUTION_CALLBACKS<TV>& deformables_evolution_callbacks_input)
    {deformables_evolution_callbacks=&deformables_evolution_callbacks_input;}

//#####################################################################
    void Advance_One_Time_Step_Position(const T dt,const T time);
    void Advance_One_Time_Step_Velocity(const T dt,const T time);
    bool Adjust_Velocity_For_Self_Repulsion_And_Self_Collisions(const T dt,const T time,int& repulsions_found,int& collisions_found,const bool exit_early=false);
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time);
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time);
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time);
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
protected:
    void Average_And_Exchange_Position();
    void Trapezoidal_Step_Velocity(const T dt,const T time);
    void Prepare_Backward_Euler_System(DEFORMABLES_BACKWARD_EULER_SYSTEM<TV>& system,const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update);
    void Finish_Backward_Euler_Step(const T dt,const T current_position_time,const bool velocity_update);
    void Backward_Euler_Step_Velocity(const T dt,const T time);
    void Update_Positions_And_Apply_Contact_Forces(const T dt,const T time,const bool use_saved_pairs=false);
    void Update_Velocity_Using_Stored_Differences(const T dt,const T time);
    void Compute_Momentum_Differences();
    void Save_Velocity();
    void Restore_Velocity() const;
    void Save_Position(ARRAY<TV>& X,ARRAY<TV>& rigid_X,ARRAY<ROTATION<TV> >& rigid_rotation);
    void Restore_Position(ARRAY_VIEW<const TV> X,ARRAY_VIEW<const TV> rigid_X,ARRAY_VIEW<const ROTATION<TV> > rigid_rotation);
    void Exchange_Velocity();
//#####################################################################
};
}
#endif

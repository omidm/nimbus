//#####################################################################
// Copyright 2008, Nipun Kwatra, Craig Schroeder, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ASYNCHRONOUS_EVOLUTION
//#####################################################################
#ifndef __ASYNCHRONOUS_EVOLUTION__
#define __ASYNCHRONOUS_EVOLUTION__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM{

template<class TV> class SOLID_BODY_COLLECTION;
template<class TV> class SOLIDS_PARAMETERS;
template<class TV> class SOLIDS_FORCES;
template<class TV> class SOLIDS_EVOLUTION;
template<class TV> class DEFORMABLES_FORCES;
template<class TV> class GENERALIZED_VELOCITY;

template<class TV>
struct BLOB
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;
    typedef typename MATRIX_POLICY<T_SPIN>::SYMMETRIC_MATRIX T_SYMMETRIC_MATRIX;

    ARRAY<int> blob_particles;
    T mass;
    T_SYMMETRIC_MATRIX inertia;
    TV center;
};

//For now this is a deformable concept
template<class TV>
class ASYNCHRONOUS_EVOLUTION:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;
    typedef typename MATRIX_POLICY<T_SPIN>::SYMMETRIC_MATRIX T_SYMMETRIC_MATRIX;
public:
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    SOLIDS_EVOLUTION<TV>* solids_evolution;
    
    enum ASYNCHRONOUS_MODE {FINE_SCALE,COARSE_SCALE};
    ASYNCHRONOUS_MODE asynchronous_mode;
    ARRAY<int> finescale_forces_indices,scaled_coarsescale_forces_indices;
    ARRAY<VECTOR<bool,3> > saved_force_settings;
    HASHTABLE<int> finescale_forces_particles_map,finescale_forces_rigid_body_particles_map,coarsescale_forces_particles_map,coarsescale_forces_rigid_body_particles_map;
    ARRAY<int> coarsescale_particle_indices,coarsescale_rigid_body_particle_indices,both_forces_particles_indices,both_forces_rigid_body_particles_indices,coarsescale_particles_pool;

    ARRAY<TV> X_n,V_n,X_np1,V_np1,V_last_finescale,V_save,X_last_finescale;
    ARRAY<TV> rigid_X_n,rigid_X_np1,rigid_X_last_finescale;
    ARRAY<ROTATION<TV> > rigid_rotation_n,rigid_rotation_np1,rigid_rotation_last_finescale;
    ARRAY<TV> rigid_V_n,rigid_V_np1,rigid_V_last_finescale,rigid_V_difference,rigid_V_save;
    ARRAY<T_SPIN> rigid_angular_velocity_n,rigid_angular_velocity_np1,rigid_angular_velocity_last_finescale,rigid_angular_velocity_difference,rigid_angular_velocity_save;
    ARRAY<T_SPIN> angular_momentum_n,angular_momentum_np1,angular_momentum_last_finescale,angular_momentum_difference;
    T time_n,time_np1,dt_epsilon;

    ARRAY<TV> finescale_force_on_particles,coarsescale_force_on_particles,accumulated_finescale_impulse_on_mixed_particles;
    ARRAY<TWIST<TV> > finescale_force_on_rigid_body_particles,coarsescale_force_on_rigid_body_particles;
    ARRAY<TV> F_full,B_full;
    ARRAY<TWIST<TV> > rigid_F_full,rigid_B_full;

    ARRAY<BLOB<TV> > blobs;

    bool initialized;
    bool use_velocity_averaging,use_adaptive,treat_mixed_particles_kinematic_in_finescale;
    bool use_squared_ratio_for_implicit_velocity_independent_scale; // True: use x_be=x+v*r*dts for BE update. False: use x_be=x+v*dts
    bool use_projection;
    T projection_rigidity;

    ASYNCHRONOUS_EVOLUTION(SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,SOLIDS_EVOLUTION<TV>* solids_evolution_input,const T cfl,const bool use_projection_input,
        const T projection_rigidity_input);
    ~ASYNCHRONOUS_EVOLUTION();

    bool In_Coarse_Mode()
    {return asynchronous_mode==COARSE_SCALE;}

    bool Take_Full_Backward_Euler_Step_For_Position_Update()
    {return (!use_projection && asynchronous_mode==COARSE_SCALE);}

//#####################################################################
    void Initialize();
    void Preprocess_Frame(const int frame,const T time_at_frame,const T time);
    void Preprocess_Substep(const T dt,const T time);
    void Postprocess_Substep(const T dt,const T time);
    void Setup_Finescale_Simulation();
    void Setup_Coarsescale_Simulation(const T dt,const T time);
    void Disable_Forces(const ARRAY<int>& disable_forces_indices,ARRAY<VECTOR<bool,3> >& saved_force_settings);
    void Enable_Forces(const ARRAY<int>& enable_forces_indices,const ARRAY<VECTOR<bool,3> >& saved_force_settings);
    void Rewind_Coarsescale_Particle_Positions();
    void Store_State(const T time);
    void Store_Time_np1_State();
    void Store_Time_Last_Finescale_State();
    void Add_Finescale_Force(DEFORMABLES_FORCES<TV>* force,const ARRAY<int>& affected_particle_indices_list,const ARRAY<int>& affected_rigid_body_particle_indices_list,const bool implicit_velocity_independent);
    void Add_Coarsescale_Force(DEFORMABLES_FORCES<TV>* force,const ARRAY<int>& affected_particle_indices_list,const ARRAY<int>& affected_rigid_body_particle_indices_list,const bool implicit_velocity_independent,const bool add_to_blobs);
    void Add_Force_Settings(DEFORMABLES_FORCES<TV>* force);
    void Setup_Particle_Indices_From_Maps();
    void Setup_Scaled_Forces_Rewind_Particles();
    void Set_Scaled_Coarsescale_Forces_Scale(const T dt_partial,const T dt_full,const bool position_update);
    void Set_Scaled_Coarsescale_Forces_Rewind_State(const bool rewind_positions,const bool rewind_velocities);
    void Update_Time_Varying_Material_Properties(const T time);
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > X,const T time); // To step asynchronous particles' positions only by remaining dt
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time); // To step asynchronous particles' positions only by remaining dt
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time);
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time);
    void Limit_Solids_Dt(T& dt,const T time);
    void Get_Active_Forces_Contributions(const T dt,const T time,ARRAY<TV>& force_on_particles,ARRAY<TWIST<TV> >& force_on_rigid_body_particles);
    void Average_Velocities();
    void Position_Velocity_Update(const T dt,const T time);
    void Accumulate_Finescale_Impulse_On_Mixed_Particles(const T dt, const T time);
    void Apply_Accumulated_Finescale_Impulse_To_Mixed_Particles();
    void Set_Asynchronous_Particles(const ARRAY<int>& async_list);
    void Set_Force_Active_Particles(DEFORMABLES_FORCES<TV>* force,const ARRAY<bool>& fine_list,bool is_fine);
    void Set_Asynchronous_Within_Distance(const T distance);
    void Add_Blob_From_Particles(const ARRAY<int>& blob_particles);
    void Compute_Blob_Properties();
    void Project(GENERALIZED_VELOCITY<TV>& V) const;
    void Boundary_Conditions(GENERALIZED_VELOCITY<TV>& V) const;
//#####################################################################
};
}
#endif

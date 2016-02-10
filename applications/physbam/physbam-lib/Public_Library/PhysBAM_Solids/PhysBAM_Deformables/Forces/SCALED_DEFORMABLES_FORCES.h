//#####################################################################
// Copyright 2008, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SCALED_DEFORMABLES_FORCES
//#####################################################################
#ifndef __SCALED_DEFORMABLES_FORCES__
#define __SCALED_DEFORMABLES_FORCES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV>
class SCALED_DEFORMABLES_FORCES:public DEFORMABLES_FORCES<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    typedef typename BASE::FREQUENCY_DATA DEFORMABLE_FREQUENCY_DATA;
    using BASE::particles;
    using BASE::use_velocity_independent_forces;using BASE::use_velocity_dependent_forces;using BASE::use_implicit_velocity_independent_forces;

    T velocity_dependent_force_scale,velocity_independent_force_scale,implicit_velocity_independent_force_scale;
    DEFORMABLES_FORCES<TV>* base_force;
private:
    bool rewind_positions; // positions should be rewound only during position update. Velocity update has positions at time n+1
    bool rewind_velocities; // velocities should be rewound in both position, velocity update. rewind_velocities=false only when simulation starts and then stays true
    const ARRAY<TV> &X_n;
    const ARRAY<TV> &V_n;
    mutable ARRAY<TV> X_rewind_particles_save,V_rewind_particles_save;
    mutable bool particle_positions_rewound,particle_velocities_rewound;

    ARRAY<int> rewind_to_time_n_particle_indices,rewind_to_time_n_rigid_body_particle_indices;
    mutable ARRAY<TV> F_temp;
    mutable ARRAY<TWIST<TV> > rigid_F_temp;

public:
    SCALED_DEFORMABLES_FORCES(DEFORMABLES_FORCES<TV>* base_force_input,PARTICLES<TV>& particles_input,const ARRAY<TV>& X_n_input,const ARRAY<TV>& V_n_input);

    virtual ~SCALED_DEFORMABLES_FORCES();

    void Set_Scale(const T velocity_dependent_force_scale_input,const T velocity_independent_force_scale_input,const T implicit_velocity_independent_force_scale_input)
    {velocity_dependent_force_scale=velocity_dependent_force_scale_input;
    velocity_independent_force_scale=velocity_independent_force_scale_input;
    implicit_velocity_independent_force_scale=implicit_velocity_independent_force_scale_input;}

    void Set_Rewind_To_Time_n_Particle_Indices(const ARRAY<int>& rewind_to_time_n_particle_indices_input)
    {rewind_to_time_n_particle_indices=rewind_to_time_n_particle_indices_input;
    X_rewind_particles_save.Resize(rewind_to_time_n_particle_indices.m);
    V_rewind_particles_save.Resize(rewind_to_time_n_particle_indices.m);}
    
    void Set_Rewind_State(const bool rewind_positions_input,const bool rewind_velocities_input)
    {rewind_positions=rewind_positions_input;rewind_velocities=rewind_velocities_input;}

//#####################################################################
private:
    void Rewind_Particle_Positions() const;
    void Restore_Rewound_Particle_Positions() const;
    void Rewind_Particle_Velocities() const;  // Need to get called in Add_Velocity_Independent_Forces only as it can depend on the time 'n' velocities. All other depend on the time 'n+1' velocity passed in by cg, so we dont need to rewind velocities there
    void Restore_Rewound_Particle_Velocities() const;
public:
    void Use_Rest_State_For_Strain_Rate(const bool use_rest_state_for_strain_rate_input=true) PHYSBAM_OVERRIDE;
    void Limit_Time_Step_By_Strain_Rate(const bool limit_time_step_by_strain_rate_input=true,const T max_strain_per_time_step_input=.1) PHYSBAM_OVERRIDE;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const PHYSBAM_OVERRIDE;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE;
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
    void Initialize_CFL(ARRAY_VIEW<DEFORMABLE_FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE;
    T Potential_Energy(const T time) const PHYSBAM_OVERRIDE;
    void Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name="") const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

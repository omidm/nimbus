//#####################################################################
// Copyright 2008, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SCALED_DEFORMABLES_FORCES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SCALED_DEFORMABLES_FORCES<TV>::
SCALED_DEFORMABLES_FORCES(DEFORMABLES_FORCES<TV>* base_force_input,PARTICLES<TV>& particles_input,const ARRAY<TV>& X_n_input,const ARRAY<TV>& V_n_input)
    :DEFORMABLES_FORCES<TV>(particles_input),base_force(base_force_input),
    rewind_positions(false),rewind_velocities(false),X_n(X_n_input),V_n(V_n_input),
    particle_positions_rewound(false),particle_velocities_rewound(false)
{
    use_velocity_independent_forces=base_force->use_velocity_independent_forces;
    use_velocity_dependent_forces=base_force->use_velocity_dependent_forces;
    use_implicit_velocity_independent_forces=base_force->use_implicit_velocity_independent_forces;
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SCALED_DEFORMABLES_FORCES<TV>::
~SCALED_DEFORMABLES_FORCES()
{
}
//#####################################################################
// Function Rewind_Particle_Positions
//#####################################################################
template<class TV> void SCALED_DEFORMABLES_FORCES<TV>::
Rewind_Particle_Positions() const
{
    if(rewind_positions) particle_positions_rewound=true;
    else{particle_positions_rewound=false;return;}

    for(int i=1;i<=rewind_to_time_n_particle_indices.m;i++){int p=rewind_to_time_n_particle_indices(i);
        X_rewind_particles_save(i)=particles.X(p);
        particles.X(p)=X_n(p);}
}
//#####################################################################
// Function Restore_Rewound_Particle_Positions
//#####################################################################
template<class TV> void SCALED_DEFORMABLES_FORCES<TV>::
Restore_Rewound_Particle_Positions() const
{
    if(!particle_positions_rewound) return;
    particle_positions_rewound=false;

    for(int i=1;i<=rewind_to_time_n_particle_indices.m;i++){int p=rewind_to_time_n_particle_indices(i);
        particles.X(p)=X_rewind_particles_save(i);}
}
//#####################################################################
// Function Rewind_Particle_Velocities
//#####################################################################
template<class TV> void SCALED_DEFORMABLES_FORCES<TV>::
Rewind_Particle_Velocities() const
{
    if(rewind_velocities) particle_velocities_rewound=true;
    else{particle_velocities_rewound=false;return;}

    for(int i=1;i<=rewind_to_time_n_particle_indices.m;i++){int p=rewind_to_time_n_particle_indices(i);
        V_rewind_particles_save(i)=particles.V(p);
        particles.V(p)=V_n(p);}
}
//#####################################################################
// Function Restore_Rewound_Particle_Velocities
//#####################################################################
template<class TV> void SCALED_DEFORMABLES_FORCES<TV>::
Restore_Rewound_Particle_Velocities() const
{
    if(!particle_velocities_rewound) return;
    particle_velocities_rewound=false;

    for(int i=1;i<=rewind_to_time_n_particle_indices.m;i++){int p=rewind_to_time_n_particle_indices(i);
        particles.V(p)=V_rewind_particles_save(i);}
}
//#####################################################################
// Function Use_Rest_State_For_Strain_Rate
//#####################################################################
template<class TV> void SCALED_DEFORMABLES_FORCES<TV>::
Use_Rest_State_For_Strain_Rate(const bool use_rest_state_for_strain_rate_input)
{
    Rewind_Particle_Positions();
    base_force->Use_Rest_State_For_Strain_Rate(use_rest_state_for_strain_rate_input);
    Restore_Rewound_Particle_Positions();
}
//#####################################################################
// Function Limit_Time_Step_By_Strain_Rate
//#####################################################################
template<class TV> void SCALED_DEFORMABLES_FORCES<TV>::
Limit_Time_Step_By_Strain_Rate(const bool limit_time_step_by_strain_rate_input,const T max_strain_per_time_step_input)
{
    Rewind_Particle_Positions();
    base_force->Limit_Time_Step_By_Strain_Rate(limit_time_step_by_strain_rate_input,max_strain_per_time_step_input);
    Restore_Rewound_Particle_Positions();
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void SCALED_DEFORMABLES_FORCES<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    base_force->Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void SCALED_DEFORMABLES_FORCES<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    base_force->Update_Mpi(particle_is_simulated,mpi_solids);
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void SCALED_DEFORMABLES_FORCES<TV>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    Rewind_Particle_Positions();
    base_force->Update_Position_Based_State(time,is_position_update);
    Restore_Rewound_Particle_Positions();
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void SCALED_DEFORMABLES_FORCES<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    Rewind_Particle_Positions();
    Rewind_Particle_Velocities(); // Need because Velocity_Independent_Forces can (e.g. wind) depend on the time 'n' velocities

    F_temp.Resize(F.Size());
    ARRAYS_COMPUTATIONS::Fill(F_temp,TV());
    base_force->Add_Velocity_Independent_Forces(F_temp,time);
    F+=velocity_independent_force_scale*F_temp;

    Restore_Rewound_Particle_Positions();
    Restore_Rewound_Particle_Velocities();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void SCALED_DEFORMABLES_FORCES<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    Rewind_Particle_Positions();

    ARRAYS_COMPUTATIONS::Fill(F_temp,TV());
    base_force->Add_Velocity_Dependent_Forces(V,F_temp,time);
    F+=velocity_dependent_force_scale*F_temp;

    Restore_Rewound_Particle_Positions();
}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class TV> void SCALED_DEFORMABLES_FORCES<TV>::
Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const
{
    Rewind_Particle_Positions();

    F_temp.Resize(dF.Size());ARRAYS_COMPUTATIONS::Fill(F_temp,TV());
    base_force->Add_Force_Differential(dX,F_temp,time);
    dF+=velocity_independent_force_scale*F_temp;

    Restore_Rewound_Particle_Positions();
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void SCALED_DEFORMABLES_FORCES<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    Rewind_Particle_Positions();

    ARRAYS_COMPUTATIONS::Fill(F_temp,TV());
    base_force->Add_Implicit_Velocity_Independent_Forces(V,F_temp,time);
    F+=implicit_velocity_independent_force_scale*F_temp;

    Restore_Rewound_Particle_Positions();
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void SCALED_DEFORMABLES_FORCES<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
    Rewind_Particle_Positions();
    base_force->Enforce_Definiteness(enforce_definiteness_input);
    Restore_Rewound_Particle_Positions();
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR SCALED_DEFORMABLES_FORCES<TV>::
CFL_Strain_Rate() const
{
    PHYSBAM_FATAL_ERROR("This should not be called");
    // TODO: Will we need to scale it?
    Rewind_Particle_Positions();
    T cfl_strain_rate=base_force->CFL_Strain_Rate();
    Restore_Rewound_Particle_Positions();
    return cfl_strain_rate;
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void SCALED_DEFORMABLES_FORCES<TV>::
Initialize_CFL(ARRAY_VIEW<DEFORMABLE_FREQUENCY_DATA> frequency)
{
    if(use_velocity_independent_forces || use_velocity_dependent_forces || use_implicit_velocity_independent_forces){
        Rewind_Particle_Positions();
        base_force->Initialize_CFL(frequency);
        Restore_Rewound_Particle_Positions();}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR SCALED_DEFORMABLES_FORCES<TV>::
Potential_Energy(const T time) const
{
    Rewind_Particle_Positions();
    T potential_energy=base_force->Potential_Energy(time);
    Restore_Rewound_Particle_Positions();
    return potential_energy;
}
//#####################################################################
// Function Add_Force_Data
//#####################################################################
template<class TV> void SCALED_DEFORMABLES_FORCES<TV>::
Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name) const
{
    Rewind_Particle_Positions();
    base_force->Add_Force_Data(force_data_list,"SCALED_DEFORMABLES_FORCES");
    Restore_Rewound_Particle_Positions();
}
template class SCALED_DEFORMABLES_FORCES<VECTOR<float,1> >;
template class SCALED_DEFORMABLES_FORCES<VECTOR<float,2> >;
template class SCALED_DEFORMABLES_FORCES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SCALED_DEFORMABLES_FORCES<VECTOR<double,1> >;
template class SCALED_DEFORMABLES_FORCES<VECTOR<double,2> >;
template class SCALED_DEFORMABLES_FORCES<VECTOR<double,3> >;
#endif

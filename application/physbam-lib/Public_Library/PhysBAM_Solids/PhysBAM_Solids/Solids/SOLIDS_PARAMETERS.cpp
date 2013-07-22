//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_PARAMETERS
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <climits>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLIDS_PARAMETERS<TV>::
SOLIDS_PARAMETERS()
    :triangle_collision_parameters(*new TRIANGLE_COLLISION_PARAMETERS<TV>),implicit_solve_parameters(*new IMPLICIT_SOLVE_PARAMETERS<TV>),
    rigid_body_collision_parameters(*new RIGID_BODY_COLLISION_PARAMETERS<TV>),rigid_body_evolution_parameters(*new RIGID_BODY_EVOLUTION_PARAMETERS<TV>),
    deformable_object_collision_parameters(*new DEFORMABLE_OBJECT_COLLISION_PARAMETERS<TV>),write_deformable_body(true),verbose(true),verbose_dt(false),cfl(10),min_dt(0),
    fracture(false),write_static_variables_every_frame(false),newton_tolerance((T)1e-3),newton_iterations(1),use_partially_converged_result(true),
    write_from_every_process(true),enforce_repulsions_in_cg(true),use_post_cg_constraints(true),use_rigid_deformable_contact(false),rigid_cluster_fracture_frequency(INT_MAX),
    use_trapezoidal_rule_for_velocities(true),enforce_poststabilization_in_cg(true),use_backward_euler_position_update(false),
    enforce_energy_conservation(false),energy_correction_iterations(20),no_contact_friction(false),use_projections_in_position_update(false),
    set_velocity_from_positions(false),set_velocity_from_positions_iterations(100),set_velocity_from_positions_tolerance((T)1e-6),set_velocity_from_positions_use_orig_force(true),
    set_velocity_from_positions_reset_alphas(false),set_velocity_from_positions_move_RE_to_KE_elastic(false),set_velocity_from_positions_move_RE_to_KE_damping(false),
    set_velocity_from_positions_damping(true),set_velocity_from_positions_physbam(false),set_velocity_from_positions_conserve_exactly(false),
    set_velocity_from_positions_percent_energy_recovered((T)1),allow_altitude_spring_change_between_updates(true)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLIDS_PARAMETERS<TV>::
~SOLIDS_PARAMETERS()
{
    delete &triangle_collision_parameters;
    delete &implicit_solve_parameters;
    delete &rigid_body_collision_parameters;
    delete &rigid_body_evolution_parameters;
    delete &deformable_object_collision_parameters;
}
//#####################################################################
template class SOLIDS_PARAMETERS<VECTOR<float,1> >;
template class SOLIDS_PARAMETERS<VECTOR<float,2> >;
template class SOLIDS_PARAMETERS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SOLIDS_PARAMETERS<VECTOR<double,1> >;
template class SOLIDS_PARAMETERS<VECTOR<double,2> >;
template class SOLIDS_PARAMETERS<VECTOR<double,3> >;
#endif
}

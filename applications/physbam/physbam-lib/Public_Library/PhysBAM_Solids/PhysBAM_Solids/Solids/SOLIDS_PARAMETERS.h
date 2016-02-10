//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_PARAMETERS
//#####################################################################
#ifndef __SOLIDS_PARAMETERS__
#define __SOLIDS_PARAMETERS__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

template<class TV> class TRIANGLE_COLLISION_PARAMETERS;
template<class TV> class IMPLICIT_SOLVE_PARAMETERS;
template<class TV> class RIGID_BODY_COLLISION_PARAMETERS;
template<class TV> class RIGID_BODY_EVOLUTION_PARAMETERS;
template<class TV> class DEFORMABLE_OBJECT_COLLISION_PARAMETERS;

template<class TV>
class SOLIDS_PARAMETERS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    TRIANGLE_COLLISION_PARAMETERS<TV>& triangle_collision_parameters;
    IMPLICIT_SOLVE_PARAMETERS<TV>& implicit_solve_parameters;
    RIGID_BODY_COLLISION_PARAMETERS<TV>& rigid_body_collision_parameters;
    RIGID_BODY_EVOLUTION_PARAMETERS<TV>& rigid_body_evolution_parameters;
    DEFORMABLE_OBJECT_COLLISION_PARAMETERS<TV>& deformable_object_collision_parameters;
    bool write_deformable_body;
    bool verbose,verbose_dt;
    T cfl;
    T min_dt;
    bool fracture;
    bool write_static_variables_every_frame;
    T newton_tolerance;
    int newton_iterations;
    bool use_partially_converged_result;
    bool write_from_every_process; // otherwise only root writes

    bool enforce_repulsions_in_cg;
    bool use_post_cg_constraints;

    // rigid/rigid and rigid/deformable collisions and contact
    bool use_rigid_deformable_contact;
    int rigid_cluster_fracture_frequency;

    bool use_trapezoidal_rule_for_velocities; // otherwise use backward euler from n->n+1
    bool enforce_poststabilization_in_cg;
    bool use_backward_euler_position_update;

    bool enforce_energy_conservation;
    int energy_correction_iterations;

    bool no_contact_friction;
    bool use_projections_in_position_update;
    bool set_velocity_from_positions;
    int set_velocity_from_positions_iterations;
    T set_velocity_from_positions_tolerance;
    bool set_velocity_from_positions_use_orig_force;
    bool set_velocity_from_positions_reset_alphas;
    T set_velocity_from_positions_move_RE_to_KE_elastic;
    T set_velocity_from_positions_move_RE_to_KE_damping;
    T set_velocity_from_positions_damping;
    bool set_velocity_from_positions_physbam;
    bool set_velocity_from_positions_conserve_exactly;
    T set_velocity_from_positions_percent_energy_recovered;

    bool allow_altitude_spring_change_between_updates;

    SOLIDS_PARAMETERS();
    virtual ~SOLIDS_PARAMETERS();
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_COLLISION_PARAMETERS
//#####################################################################
#ifndef __TRIANGLE_COLLISION_PARAMETERS__
#define __TRIANGLE_COLLISION_PARAMETERS__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_OBJECT_FORWARD.h>
namespace PhysBAM{

template<class TV>
class TRIANGLE_COLLISION_PARAMETERS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    bool perform_self_collision;
    bool initialize_collisions_without_objects;
    bool temporary_enable_collisions;
    int repulsion_pair_update_count,repulsion_pair_update_frequency;
    T self_collision_free_time;
    int steps_since_self_collision_free,total_collision_loops;
    int topological_hierarchy_build_count,topological_hierarchy_build_frequency;
    bool check_initial_mesh_for_self_intersection;
    bool check_mesh_for_self_intersection;
    bool turn_off_all_collisions;
    bool allow_intersections;
    T allow_intersections_tolerance;
    T collisions_small_number;
    T collisions_repulsion_thickness; // usually set equal to the cloth thickness
    bool clamp_repulsion_thickness;
    T collisions_repulsion_clamp_fraction;
    T collisions_collision_thickness;
    T repulsions_youngs_modulus;
    T collisions_final_repulsion_youngs_modulus;
    T repulsions_limiter_fraction; // fraction of the repulsion thickness that the elastic repulsion can take a point
    T collisions_final_repulsion_limiter_fraction; // fraction of the repulsion thickness that the final elastic repulsion can take a point
    T collisions_disable_repulsions_based_on_proximity_factor; // disables repulsions that would occur initially if distances were multiplied by this factor (0 disables nothing)
    bool collisions_output_repulsion_results,collisions_output_collision_results,collisions_output_number_checked; // diagnostics
    T self_collision_friction_coefficient;
    int collisions_nonrigid_collision_attempts;
    bool use_new_triangle_hierarchy,use_sort_for_new_triangle_hierarchy;
    bool perform_per_time_step_repulsions;
    bool perform_per_collision_step_repulsions;
    bool output_interaction_pairs;
    bool perform_repulsion_pair_attractions;
    T repulsion_pair_attractions_threshold;
    bool use_gauss_jacobi;

    TRIANGLE_COLLISION_PARAMETERS();
    virtual ~TRIANGLE_COLLISION_PARAMETERS();
//#####################################################################
};
}
#endif

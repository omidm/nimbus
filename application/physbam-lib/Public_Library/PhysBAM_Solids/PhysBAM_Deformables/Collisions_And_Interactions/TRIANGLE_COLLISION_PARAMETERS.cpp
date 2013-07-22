//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_COLLISION_PARAMETERS
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Collision_Detection/COLLISION_GEOMETRY_SPATIAL_PARTITION.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/STRUCTURE_INTERACTION_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> TRIANGLE_COLLISION_PARAMETERS<TV>::
TRIANGLE_COLLISION_PARAMETERS()
    :perform_self_collision(true),
    initialize_collisions_without_objects(false),temporary_enable_collisions(true),repulsion_pair_update_count(0),repulsion_pair_update_frequency(INT_MAX),
    self_collision_free_time(0),steps_since_self_collision_free(0),total_collision_loops(4),topological_hierarchy_build_count(0),topological_hierarchy_build_frequency(1),
    check_initial_mesh_for_self_intersection(false),check_mesh_for_self_intersection(true),turn_off_all_collisions(false),allow_intersections(false),
    allow_intersections_tolerance((T)1e-8),collisions_small_number((T)1e-8),collisions_repulsion_thickness((T)1e-2),clamp_repulsion_thickness(true),
    collisions_repulsion_clamp_fraction((T).9),collisions_collision_thickness((T)1e-6),repulsions_youngs_modulus((T)30),collisions_final_repulsion_youngs_modulus((T)30),
    repulsions_limiter_fraction((T).1),collisions_final_repulsion_limiter_fraction((T).1),
    collisions_disable_repulsions_based_on_proximity_factor((T)1.5),collisions_output_repulsion_results(false),collisions_output_collision_results(true),collisions_output_number_checked(true),
    self_collision_friction_coefficient((T).4),collisions_nonrigid_collision_attempts(3),use_new_triangle_hierarchy(false),
    use_sort_for_new_triangle_hierarchy(false),perform_per_time_step_repulsions(true),perform_per_collision_step_repulsions(false),output_interaction_pairs(false),
    perform_repulsion_pair_attractions(true),repulsion_pair_attractions_threshold((T)-2),use_gauss_jacobi(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> TRIANGLE_COLLISION_PARAMETERS<TV>::
~TRIANGLE_COLLISION_PARAMETERS()
{
}
//#####################################################################
template class TRIANGLE_COLLISION_PARAMETERS<VECTOR<float,1> >;
template class TRIANGLE_COLLISION_PARAMETERS<VECTOR<float,2> >;
template class TRIANGLE_COLLISION_PARAMETERS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TRIANGLE_COLLISION_PARAMETERS<VECTOR<double,1> >;
template class TRIANGLE_COLLISION_PARAMETERS<VECTOR<double,2> >;
template class TRIANGLE_COLLISION_PARAMETERS<VECTOR<double,3> >;
#endif
}

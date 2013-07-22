//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_COLLISION_PARAMETERS
//#####################################################################
#ifndef __RIGID_BODY_COLLISION_PARAMETERS__
#define __RIGID_BODY_COLLISION_PARAMETERS__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Collisions/COLLISIONS_GEOMETRY_FORWARD.h>
namespace PhysBAM{

template<class TV>
class RIGID_BODY_COLLISION_PARAMETERS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    bool rigid_collisions_use_particle_partition;
    bool rigid_collisions_use_particle_partition_center_phi_test;
    int rigid_collisions_particle_partition_size;
    bool rigid_collisions_use_triangle_hierarchy;
    bool rigid_collisions_use_triangle_hierarchy_center_phi_test;
    bool rigid_collisions_use_edge_intersection;
    bool rigid_collisions_print_interpenetration_statistics;
    SPATIAL_PARTITION_VOXEL_SIZE_HEURISTIC rigid_collisions_spatial_partition_voxel_size_heuristic;
    int rigid_collisions_spatial_partition_number_of_cells;
    T rigid_collisions_spatial_partition_voxel_size_scale_factor;
    T collision_body_thickness;
    T collision_bounding_box_thickness;
    bool use_push_out;
    bool use_legacy_push_out;
    bool use_projected_gauss_seidel_push_out;
    bool use_push_out_rotation;
    bool use_epsilon_scaling,use_epsilon_scaling_for_level;
    bool use_shock_propagation;
    bool use_persistant_contact;
    int collision_iterations;
    int contact_iterations;
    int contact_project_iterations;
    bool use_analytic_collisions;
    bool use_fracture_pattern;
    bool allow_refracturing;
    bool use_fracture_particle_optimization;
    bool enforce_rigid_rigid_contact_in_cg;
    bool perform_collisions;
    bool perform_contact;
    bool use_projected_gauss_seidel;
    bool use_ccd;
    bool thin_shells;
    T contact_proximity;
    T projected_gauss_seidel_tolerance;

    T collision_thickness;
    T contact_thickness;
    
    int threadid;

    RIGID_BODY_COLLISION_PARAMETERS();
    virtual ~RIGID_BODY_COLLISION_PARAMETERS();

    void Set_Collision_Body_Thickness(const T thickness=1e-6)
    {collision_body_thickness=thickness;}
//#####################################################################
};
}
#endif

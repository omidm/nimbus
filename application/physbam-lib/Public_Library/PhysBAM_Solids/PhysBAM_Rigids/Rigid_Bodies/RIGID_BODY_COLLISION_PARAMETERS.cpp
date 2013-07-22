//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_COLLISION_PARAMETERS
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_BODY_COLLISION_PARAMETERS<TV>::
RIGID_BODY_COLLISION_PARAMETERS()
    :rigid_collisions_use_particle_partition(true),rigid_collisions_use_particle_partition_center_phi_test(true),rigid_collisions_particle_partition_size(2),
    rigid_collisions_use_triangle_hierarchy(false),rigid_collisions_use_triangle_hierarchy_center_phi_test(false),rigid_collisions_use_edge_intersection(false),
    rigid_collisions_print_interpenetration_statistics(false),rigid_collisions_spatial_partition_voxel_size_heuristic(SPATIAL_PARTITION_MAX_BOX_SIZE),
    rigid_collisions_spatial_partition_number_of_cells(100),rigid_collisions_spatial_partition_voxel_size_scale_factor(8),collision_body_thickness(0),collision_bounding_box_thickness(0),
    use_push_out(true),use_legacy_push_out(false),use_projected_gauss_seidel_push_out(false),use_push_out_rotation(true),use_epsilon_scaling(true),use_epsilon_scaling_for_level(true),
    use_shock_propagation(true),use_persistant_contact(false),collision_iterations(5),contact_iterations(10),contact_project_iterations(5),use_analytic_collisions(false),
    use_fracture_pattern(false),allow_refracturing(false),use_fracture_particle_optimization(true),enforce_rigid_rigid_contact_in_cg(true),perform_collisions(true),perform_contact(true),
    use_projected_gauss_seidel(false),use_ccd(false),thin_shells(false),contact_proximity(0),projected_gauss_seidel_tolerance((T)1e-3),threadid(1)
{
    collision_thickness=(T)1e-2;
    contact_thickness=2*collision_thickness;

    Set_Collision_Body_Thickness();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_BODY_COLLISION_PARAMETERS<TV>::
~RIGID_BODY_COLLISION_PARAMETERS()
{
}
//#####################################################################
template class RIGID_BODY_COLLISION_PARAMETERS<VECTOR<float,1> >;
template class RIGID_BODY_COLLISION_PARAMETERS<VECTOR<float,2> >;
template class RIGID_BODY_COLLISION_PARAMETERS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_BODY_COLLISION_PARAMETERS<VECTOR<double,1> >;
template class RIGID_BODY_COLLISION_PARAMETERS<VECTOR<double,2> >;
template class RIGID_BODY_COLLISION_PARAMETERS<VECTOR<double,3> >;
#endif
}

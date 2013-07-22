//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace CONTACT_PAIRS
//##################################################################### 
#ifndef __SURFACE_CONTACT_PAIR__
#define __SURFACE_CONTACT_PAIR__

#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_2D.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_PARTICLE_INTERSECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_SKIP_COLLISION_CHECK.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_TRIANGLE_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_STRUCTURE_INTERACTION_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGIDS_COLLISION_CALLBACKS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/PARTICLES_IN_IMPLICIT_OBJECT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/SOLVE_CONTACT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
namespace PhysBAM{

namespace CONTACT_PAIRS
{
template<class TV>
bool Update_Surface_Contact_Pair(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,const int id_1,const int id_2,
    const bool correct_contact_energy,const int max_iterations,const typename TV::SCALAR epsilon_scale,const typename TV::SCALAR dt,const typename TV::SCALAR time,
    const bool use_triangle_hierarchy,const bool use_edge_intersection,const bool use_triangle_hierarchy_center_phi_test,const bool mpi_one_ghost)
{
    typedef typename TV::SCALAR T;

    bool need_another_iteration=true;int iteration=0;bool objects_moved=false;
    while(need_another_iteration && ++iteration<=max_iterations){need_another_iteration=false;
        T smallest_value;int smallest_index=0;TV collision_location,collision_normal,collision_relative_velocity;bool ignored_separating;
        bool found_intersection=rigid_body_collisions.Get_First_Intersection_Point(id_1,id_2,smallest_value,smallest_index,collision_location,collision_normal,
            collision_relative_velocity,true,0,ignored_separating,dt,true);
        if(!found_intersection){collision_callbacks.Swap_States(id_1,id_2);rigid_body_collisions.skip_collision_check.Set_Last_Checked(id_1,id_2);return objects_moved;}
        else{
            RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_body_cluster_bindings=rigid_body_collisions.rigid_body_cluster_bindings;
            int parent_id_1=rigid_body_cluster_bindings.Get_Parent_Index(id_1);if(!parent_id_1) parent_id_1=id_1;
            int parent_id_2=rigid_body_cluster_bindings.Get_Parent_Index(id_2);if(!parent_id_2) parent_id_2=id_2;
            collision_callbacks.Swap_States(id_1,id_2);collision_callbacks.Restore_Position(parent_id_1);collision_callbacks.Restore_Position(parent_id_2);
            collision_callbacks.Euler_Step_Position(parent_id_1,smallest_value,time);collision_callbacks.Euler_Step_Position(parent_id_2,smallest_value,time);
            SOLVE_CONTACT::Update_Contact_Pair_Helper<TV>(rigid_body_collisions,collision_callbacks,id_1,id_2,dt,time,epsilon_scale,
                collision_location,collision_normal,collision_relative_velocity,correct_contact_energy,rigid_body_collisions.rolling_friction,mpi_one_ghost);
            need_another_iteration=true;objects_moved=true;}}
    return objects_moved;
}
/*template<> 
bool Update_Surface_Contact_Pair<VECTOR<float,1> >(RIGID_BODY_COLLISIONS<VECTOR<float,1> >& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<VECTOR<float,1> >& collision_callbacks,const int id_1,const int id_2,
    const bool correct_contact_energy,const int max_iterations,const float epsilon_scale,const float dt,const float time,
    const bool use_triangle_hierarchy,const bool use_edge_intersection,const bool use_triangle_hierarchy_center_phi_test,const bool mpi_one_ghost)
{PHYSBAM_FATAL_ERROR();}
template<> 
bool Update_Surface_Contact_Pair<VECTOR<double,1> >(RIGID_BODY_COLLISIONS<VECTOR<double,1> >& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<VECTOR<double,1> >& collision_callbacks,const int id_1,const int id_2,
    const bool correct_contact_energy,const int max_iterations,const double epsilon_scale,const double dt,const double time,
    const bool use_triangle_hierarchy,const bool use_edge_intersection,const bool use_triangle_hierarchy_center_phi_test,const bool mpi_one_ghost)
{PHYSBAM_FATAL_ERROR();}*/

#define INSTANTIATION_HELPER(T,d) \
    template bool Update_Surface_Contact_Pair(RIGID_BODY_COLLISIONS<VECTOR<T,d> >& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks,\
        const int id_1,const int id_2,const bool correct_contact_energy,const int max_iterations,const T epsilon_scale,const T dt,const T time,const bool use_triangle_hierarchy, \
        const bool use_edge_intersection,const bool use_triangle_hierarchy_center_phi_test,const bool mpi_one_ghost);

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif
}
}
#endif

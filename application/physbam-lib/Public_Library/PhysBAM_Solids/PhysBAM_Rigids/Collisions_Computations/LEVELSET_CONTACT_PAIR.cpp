//#####################################################################
// Copyright 2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace CONTACT_PAIRS
//##################################################################### 
#ifndef __LEVELSET_CONTACT_PAIR__
#define __LEVELSET_CONTACT_PAIR__

#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_PARTICLE_INTERSECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_SKIP_COLLISION_CHECK.h>
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
bool Update_Levelset_Contact_Pair(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,const int id_1,const int id_2,
    const bool correct_contact_energy,const int max_iterations,const typename TV::SCALAR epsilon_scale,const typename TV::SCALAR dt,const typename TV::SCALAR time,
    const bool use_triangle_hierarchy,const bool use_edge_intersection,const bool use_triangle_hierarchy_center_phi_test,const bool mpi_one_ghost)
{
    typedef typename TV::SCALAR T;

    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=rigid_body_collisions.rigid_body_collection;
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_body_cluster_bindings=rigid_body_collisions.rigid_body_cluster_bindings;
    RIGID_BODY_INTERSECTIONS<TV> intersections(rigid_body_collection);

    int parent_id_1=rigid_body_cluster_bindings.Get_Parent(rigid_body_collection.Rigid_Body(id_1)).particle_index;
    int parent_id_2=rigid_body_cluster_bindings.Get_Parent(rigid_body_collection.Rigid_Body(id_2)).particle_index;
    T thickness=rigid_body_collisions.desired_separation_distance+max(rigid_body_collection.Rigid_Body(id_1).surface_roughness,rigid_body_collection.Rigid_Body(id_2).surface_roughness);
    ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> > particle_intersections;
    particle_intersections.Preallocate(100);
    bool need_another_iteration=true;int iteration=0;bool objects_moved=false;
    while(need_another_iteration && ++iteration<=max_iterations){need_another_iteration=false;
        particle_intersections.Remove_All();
        // bodies in x_{n+1} positions
        PARTICLES_IN_IMPLICIT_OBJECT::Append_All_Intersections(rigid_body_collection.Rigid_Body(id_1),rigid_body_collection.Rigid_Body(id_2),particle_intersections,thickness,
            use_triangle_hierarchy,use_edge_intersection,use_triangle_hierarchy_center_phi_test);
        if(!particle_intersections.m){rigid_body_collisions.skip_collision_check.Set_Last_Checked(id_1,id_2);continue;}
        // revert to the saved positions & save the proposed positions in rigid_frame_save - restore rigid_frame_save below
        collision_callbacks.Swap_States(id_1,id_2);
        T smallest_value=FLT_MAX;int smallest_index=0;TV collision_location,collision_normal,collision_relative_velocity;
        for(int i=1;i<=particle_intersections.m;i++){
            const RIGID_BODY_PARTICLE_INTERSECTION<TV>& intersection=particle_intersections(i);
            FRAME<TV> saved_transform=collision_callbacks.Saved_Particle_To_Levelset_Body_Transform(intersection.levelset_body,intersection.particle_body);
            T phi=(*rigid_body_collection.Rigid_Body(intersection.levelset_body).implicit_object->object_space_implicit_object)(saved_transform*intersection.particle_location);
            if(phi<smallest_value){ // inside in *proposed configuration*
                RIGID_BODY<TV> &body1=rigid_body_collection.Rigid_Body(intersection.particle_body),&body2=rigid_body_collection.Rigid_Body(intersection.levelset_body);
                RIGID_BODY<TV> &body_parent1=rigid_body_cluster_bindings.Get_Parent(body1),&body_parent2=rigid_body_cluster_bindings.Get_Parent(body2);
                TV location=body1.World_Space_Point(intersection.particle_location);
                // use extended normal so it's defined if we're outside the body in the old configuration
                TV normal=rigid_body_collisions.use_parent_normal?body_parent2.implicit_object->Extended_Normal(location):body2.implicit_object->Extended_Normal(location);
                TV relative_velocity=RIGID_GEOMETRY<TV>::Relative_Velocity_At_Geometry1_Particle(body_parent1,body_parent2,location,intersection.particle_index);

                if(TV::Dot_Product(relative_velocity,normal)<0){ // approaching in *old configuration*
                    //TODO: need to add pair to rigid_body_particle_intersections if separating?
                    smallest_value=phi;smallest_index=i;collision_location=location;collision_normal=normal;collision_relative_velocity=relative_velocity;}}}
        collision_callbacks.Swap_States(id_1,id_2);
        if(smallest_index){
            if(parent_id_1!=id_1) collision_callbacks.Exchange_Frame(id_1);
            if(parent_id_2!=id_2) collision_callbacks.Exchange_Frame(id_2);
            const RIGID_BODY_PARTICLE_INTERSECTION<TV>& intersection=particle_intersections(smallest_index);
            SOLVE_CONTACT::Update_Contact_Pair_Helper<TV>(rigid_body_collisions,collision_callbacks,intersection.particle_body,intersection.levelset_body,dt,time,epsilon_scale,
                collision_location,collision_normal,collision_relative_velocity,correct_contact_energy,rigid_body_collisions.rolling_friction,mpi_one_ghost);
            need_another_iteration=true;objects_moved=true;}}
    return objects_moved;
}

#define INSTANTIATION_HELPER(T,d) \
    template bool Update_Levelset_Contact_Pair(RIGID_BODY_COLLISIONS<VECTOR<T,d> >& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks,\
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

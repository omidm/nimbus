//#####################################################################
// Copyright 2007, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_DEFORMABLE_COLLISIONS
//#####################################################################
#ifndef __RIGID_DEFORMABLE_COLLISIONS__
#define __RIGID_DEFORMABLE_COLLISIONS__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_ID.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/SOLIDS_COLLISIONS_FORWARD.h>
namespace PhysBAM{

template<class TV> class ARTICULATED_RIGID_BODY;
template<class TV> class DEFORMABLE_OBJECT;
template<class TV> class RIGID_BODY_COLLISIONS;
template<class TV> class RIGID_BODY;
template<class TV> class RIGID_BODY_STATE;
template<class TV> struct RIGID_BODY_PARTICLE_INTERSECTION;
template<class TV> class SOLID_BODY_COLLECTION;
template<class TV> class RIGID_COLLISION_GEOMETRY;
template<class TV> class SOLIDS_PARAMETERS;
template<class TV>
class RIGID_DEFORMABLE_COLLISIONS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;
    typedef VECTOR<int,TV::dimension> ELEMENT;typedef typename MESH_POLICY<TV::dimension>::MESH T_MESH;typedef typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX T_SYMMETRIC_MATRIX;
    typedef typename RIGID_BODY_POLICY<TV>::WORLD_SPACE_INERTIA_TENSOR T_WORLD_SPACE_INERTIA_TENSOR;
public:
    struct PRECOMPUTE_CONTACT_PROJECTION
    {
        PRECOMPUTE_CONTACT_PROJECTION(const RIGID_BODY<TV>& rigid_body_input,const bool A_inverted_input)
            :rigid_body(rigid_body_input),A_inverted(A_inverted_input)
        {}

        const RIGID_BODY<TV>& rigid_body;
        ARRAY<int> particles;
        ARRAY<TV> V_rel_target;
        ARRAY<TV> N_over_NT_K_N,r,N;
        ARRAY<T_SPIN> rN;
        MATRIX<T,TV::dimension+T_SPIN::dimension,TV::dimension+T_SPIN::dimension> A;
        bool A_inverted; // if true, j=Ab, else Aj=b
    };

    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions;
    SOLIDS_PARAMETERS<TV>& solids_parameters;
    HASHTABLE<int,ARRAY<TETRAHEDRON_COLLISION_BODY<T>*> > particle_tetrahedron_candidates;
    HASHTABLE<int,ARRAY<RIGID_COLLISION_GEOMETRY<TV>*> > particle_rigid_body_candidates;
    ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID> tetrahedron_collision_bodies,rigid_collision_bodies;
    HASHTABLE<int,HASHTABLE<int> > particles_collided_with_rigid_body;
    HASHTABLE<int,HASHTABLE<int> > particles_contacting_rigid_body; // TODO: Make this MPI friendly
    bool fix_position_gap;
    ARRAY<PRECOMPUTE_CONTACT_PROJECTION*> precompute_contact_projections;
    ARRAY<JOINT_ID> contact_joints;
    const T normal_relative_velocity_threshold;

    RIGID_DEFORMABLE_COLLISIONS(SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions_input,SOLIDS_PARAMETERS<TV>& solids_parameters_input);
    virtual ~RIGID_DEFORMABLE_COLLISIONS();

//#####################################################################
    void Add_Elastic_Collisions(const T dt,const T time,ARRAY<ROTATION<TV> >& rigid_rotation_save,ARRAY<typename TV::SPIN>& rigid_angular_momentum_difference,
        ARRAY<TV>& rigid_velocity_difference,ARRAY<TV>& rigid_X_save,ARRAY<TWIST<TV> >& rigid_velocity_save,ARRAY<typename TV::SPIN> &rigid_angular_momentum_save,ARRAY<TV>& X_save,
        ARRAY<TV>& V_save);
    void Process_Contact(const T dt,const T time,ARTICULATED_RIGID_BODY<TV>* articulated_rigid_body,const bool use_saved_pairs,const bool use_existing_contact,
        ARRAY<TV>& rigid_X_save,ARRAY<ROTATION<TV> >& rigid_rotation_save,ARRAY<TV>& rigid_velocity_difference,ARRAY<typename TV::SPIN>& rigid_angular_momentum_difference,ARRAY<TV>& X_save,
        const T collision_body_thickness);
    void Process_Push_Out();
    void Initialize_All_Contact_Projections();
    void Set_Collision_Velocities(ARRAY_VIEW<TV> V,ARRAY_VIEW<TWIST<TV> > twist,ARRAY<TV>& X_save,ARRAY<TV>& rigid_X_save,ARRAY_VIEW<const ROTATION<TV> >& rigid_rotation_save,
        ARRAY<TWIST<TV> >& rigid_velocity_save,ARRAY<typename TV::SPIN>& rigid_angular_momentum_save,ARRAY<TV>& V_save);
    void Project_Contact_Pairs(ARRAY_VIEW<TV> V,ARRAY_VIEW<TWIST<TV> > twist);
    void Process_Precomputed_Contact_With_Kinematic_Rigid_Bodies();
    void Process_Precomputed_Contact_With_Rigid_Bodies();
    void Apply_Rigid_Deformable_Collision_Impulse(RIGID_BODY<TV>& rigid_body,const int particle,const TV& location,const TV& normal,const TV& relative_velocity,const T coefficient_of_restitution,
        const T coefficient_of_friction,const bool clamp_friction_magnitude,TV& impulse,bool allow_pull=false,bool apply_impulse=true);
private:
    void Get_Rigid_And_Tetrahedron_Collision_Bodies();
    bool Update_Rigid_Deformable_Collision_Pair(RIGID_BODY<TV>& rigid_body,const int particle_index,const T dt,const T time,ARRAY<TV>& X_save,ARRAY<TV>& V_save,ARRAY<TV>& rigid_X_save,
        ARRAY<ROTATION<TV> >& rigid_rotation_save,ARRAY<typename TV::SPIN>& rigid_angular_momentum_difference,ARRAY<TV>& rigid_velocity_difference);
    bool Update_Rigid_Deformable_Contact_Pair(RIGID_BODY<TV>& rigid_body,const int p,const T dt,const T time,const T epsilon_scale,ARRAY<TV>& X_save,ARRAY<TV>& rigid_X_save,
        ARRAY<ROTATION<TV> >& rigid_rotation_save,const T collision_body_thickness,const bool process_contact_unconditionally=false);
    bool Push_Out_From_Rigid_Body(RIGID_BODY<TV>& rigid_body,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const T move_fraction);
    void Get_Particles_Contacting_Rigid_Body(const RIGID_BODY<TV>& rigid_body,ARRAY<int>& particles,const bool include_soft_bound);
    bool Push_Out_From_Particle(const int particle);
    void Get_Particles_Intersecting_Rigid_Body(const RIGID_BODY<TV>& rigid_body,ARRAY<int>& particles,ARRAY<TV>& particle_distances) const;
    void Get_Rigid_Bodies_Intersecting_Rigid_Body(const int particle_index,ARRAY<int>& rigid_bodies,ARRAY<TV>& collision_location,ARRAY<TV>& body_distances,
        ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections) const;
    void Get_Objects_Intersecting_Particle(const int particle_index,ARRAY<ELEMENT>& triangles,ARRAY<TV>& weights,ARRAY<TV>& particle_distances,ARRAY<int>& bodies,
        ARRAY<TV>& body_distances);
    void Apply_Impulse(const int particle,RIGID_BODY<TV>& rigid_body,const TV& impulse);
    void Apply_Displacement_To_Particle(const int particle_index,const TV& particle_delta_X);
    bool Process_Deformable_Contact_With_Kinematic_Rigid_Body(RIGID_BODY<TV>& rigid_body,const T dt,const T time,ARRAY<TV>& rigid_X_save,ARRAY<ROTATION<TV> >& rigid_rotation_save,
        ARRAY<TV>& X_save);
    bool Apply_Rigid_Deformable_Contact_Projection(ARRAY_VIEW<const TV> X,ARRAY_VIEW<TV> V,TV& rigid_V,T_SPIN& rigid_angular_velocity,PRECOMPUTE_CONTACT_PROJECTION& precompute);
    // precompute.particles is the input list of particles to collide against
    void Initialize_Rigid_Deformable_Contact_Projection(PRECOMPUTE_CONTACT_PROJECTION& precompute,ARRAY_VIEW<const TV> X);
//#####################################################################
};
}
#endif

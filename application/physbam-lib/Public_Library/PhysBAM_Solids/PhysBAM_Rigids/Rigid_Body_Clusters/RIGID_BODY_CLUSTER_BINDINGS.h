//#####################################################################
// Copyright 2008-2009, Jon Gretarsson, Michael Lentine, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_CLUSTER_BINDINGS
//##################################################################### 
#ifndef __RIGID_BODY_CLUSTER_BINDINGS__
#define __RIGID_BODY_CLUSTER_BINDINGS__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_ID.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_CLUSTER_CONSTITUENT_ID.h>
namespace PhysBAM{

template<class TV> class JOINT;
template<class TV> class ARTICULATED_RIGID_BODY;
template<class TV> class GRID_BASED_COLLISION_GEOMETRY;
template<class TV> class GRID;

template<class TV>
class RIGID_BODY_CLUSTER_BINDINGS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename RIGID_BODY_POLICY<TV>::INERTIA_TENSOR T_INERTIA_TENSOR;
public:

    struct INTERNAL_JOINT
    {
        JOINT<TV>* joint;
        int parent;
        int child;
    };

    struct BOUNDARY_JOINT
    {
        FRAME<TV> original_frame;
        JOINT_ID joint_id;
        int boundary_body;
        bool is_parent;
    };

    struct CLUSTER
    {
        CLUSTER();
        bool active;
        bool stored_active;
        int infinite_body;
        int parent; // TODO: Make sure this is set up...
        T saved_mass;
        T_INERTIA_TENSOR saved_inertia_tensor;
        ARRAY<INTERNAL_JOINT> internal_joints; // joint, parent, child
        ARRAY<BOUNDARY_JOINT> boundary_joints; // joint, body, is_parent
        ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children;
        ARRAY<bool,RIGID_CLUSTER_CONSTITUENT_ID> kinematic_child;
        ARRAY<FRAME<TV>,RIGID_CLUSTER_CONSTITUENT_ID> child_to_parent;
        ARRAY<FRAME<TV>,RIGID_CLUSTER_CONSTITUENT_ID> saved_child_to_parent;
    };

    struct CALLBACKS
    {
        virtual ~CALLBACKS(){}

        virtual void Pre_Advance_Unclustered(const T dt,const T time)=0;
        virtual void Post_Advance_Unclustered(const T dt,const T time)=0;
        virtual void Compute_New_Clusters_Based_On_Unclustered_Strain()=0;
        virtual bool Create_New_Clusters()=0;
        virtual void Create_Cluster(const int parent)=0;
        virtual void Destroy_Cluster(const int parent)=0;
    };

    typedef typename HASHTABLE<int,CLUSTER* const>::ITERATOR T_REVERSE_BINDING_ITERATOR;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SIMPLICIAL_OBJECT;
    typedef typename TV::SPIN T_SPIN;

    RIGID_BODY_CLUSTER_BINDINGS(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body_input);

    CALLBACKS* callbacks;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body;
    bool collide_constituent_bodies;
    bool saved;
    bool clamp_kinematic_positions;
    
    ARRAY<ARRAY<PAIR<int,RIGID_CLUSTER_CONSTITUENT_ID> > > binding_index; // child_index -> list of (parent_particle_index,connection_id)
    HASHTABLE<int,CLUSTER*> reverse_bindings; // parent_particle_index -> { (child_particle_indices,child_to_parent_frame) }

    void Update_All_Joint_Structures();
    void Update_Joint_Structures(const int parent);
    void Merge_Bindings(const int parent1,const int parent2,const bool kinematic=false);
    void Append_To_Binding(const int parent,const int child_particle,const bool kinematic=false);
    void Append_To_Binding(const int parent,const ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID>& child_particles,const ARRAY<bool,RIGID_CLUSTER_CONSTITUENT_ID>& kinematic);
    void Append_Binding(const int particle,const PAIR<int,RIGID_CLUSTER_CONSTITUENT_ID> data);
    int Add_Binding(const ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID>& child_particles,const bool kinematic=false);
    void Delete_Binding(const int parent_particle);
    void Make_Active_Parent(const int parent_particle,ARRAY<PAIR<int,RIGID_CLUSTER_CONSTITUENT_ID> >& child_list);
    void Set_Binding_Active(const int parent_rigid_body_index,const bool active,GRID_BASED_COLLISION_GEOMETRY<GRID<TV> >* fluid_collision_body_list=0);
    void Deactivate_And_Return_Clusters(ARRAY<int>& active_bindings,GRID_BASED_COLLISION_GEOMETRY<GRID<TV> >* fluid_collision_body_list=0);
    void Reactivate_Bindings(const ARRAY<int>& active_bindings,GRID_BASED_COLLISION_GEOMETRY<GRID<TV> >* fluid_collision_body_list=0);
    void Save_Bindings_State();
    void Restore_Bindings_State();
    int Size() const{return reverse_bindings.Size();}
    bool Is_Parent(const int rigid_body_index) const{if(rigid_body_index<=binding_index.m && binding_index(rigid_body_index).m>0) return binding_index(rigid_body_index)(1).x==rigid_body_index; return false;}
    int Get_Parent_Index(const int rigid_body_index) const
    {if(rigid_body_index<=binding_index.m && binding_index(rigid_body_index).m>0 && reverse_bindings.Get(binding_index(rigid_body_index)(1).x)->active) return binding_index(rigid_body_index)(1).x; return 0;}
    RIGID_BODY<TV>& Get_Parent(RIGID_BODY<TV> &body) const{int parent=Get_Parent_Index(body.particle_index);return *(parent?&rigid_body_collection.Rigid_Body(parent):&body);}
    FRAME<TV> Get_Child_To_Parent_Frame(int particle_index) const;

    void Clamp_Particles_To_Embedded_Positions() const;
    void Clamp_Particles_To_Embedded_Velocities() const;
    void Clamp_Particles_To_Embedded_Velocities(ARRAY_VIEW<TWIST<TV> > twist) const;
    void Clamp_Particles_To_Embedded_Positions(const int parent) const;
    void Clamp_Particles_To_Embedded_Velocities(const int parent) const;
    void Clamp_Particles_To_Embedded_Velocities(const int parent,ARRAY_VIEW<TWIST<TV> >& twist) const;
    void Distribute_Force_To_Parents(ARRAY_VIEW<TWIST<TV> > wrench_full) const;
    void Distribute_Mass_To_All_Active_Parents();
    void Distribute_Mass_To_Parent(const int parent_rigid_body_index);
    void Build_Aggregate_Geometry(const int parent_rigid_body_index);
    void Update_Aggregate_Geometry(const int parent_rigid_body_index);
    void Clear_Hard_Bound_Particles(ARRAY<bool>& particle_is_simulated) const;

    int Binding(const int child_particle,FRAME<TV>& frame) const;
    void Clamp_Particle_To_Embedded_Position(const int child_particle) const;
    void Clamp_Particle_To_Embedded_Velocity(const int child_particle) const;
    bool Valid_Cluster_Collision_Helper(const int rigid_particle_index,const int parent);
    bool Valid_Cluster_Collision(const int rigid_body_1,const int rigid_body_2);

//##################################################################### 
};
}
#endif

//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_OBJECT_COLLISIONS
//#####################################################################
#ifndef __DEFORMABLE_OBJECT_COLLISIONS__
#define __DEFORMABLE_OBJECT_COLLISIONS__

#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_PARTICLE_STATE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLES_COLLISIONS_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
#endif
namespace PhysBAM{

template<class TV> class EMBEDDING;
template<class TV> class RIGID_DEFORMABLE_COLLISIONS;
template<class TV> class DEFORMABLE_BODY_COLLECTION;
template<class TV> class SOFT_BINDINGS;
template<class TV> class BINDING_LIST;

template<class TV>
class DEFORMABLE_OBJECT_COLLISIONS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    PARTICLES<TV>& particles;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection;
    ARRAY<STRUCTURE<TV>*>& deformable_object_structures;
    COLLISION_GEOMETRY_COLLECTION<TV>& collision_body_list;
    ARRAY<bool> check_collision;
    ARRAY<COLLISION_PARTICLE_STATE<TV> > particle_states; // TODO: remove this in exchange for sparser enforced_particles
    ARRAY<COLLISION_GEOMETRY_ID> particle_to_collision_body_id;
    ARRAY<PAIR<int,COLLISION_PARTICLE_STATE<TV> > > enforced_particles;
    ARRAY<ARRAY<COLLISION_GEOMETRY_ID> > protecting_bodies_of_nodes;
    T collision_tolerance;
    bool use_spatial_partition;
    bool disable_multiple_levelset_collisions;
    bool use_protectors;
    T maximum_levelset_collision_projection_velocity;
    T protection_thickness;
    bool ignore_priorities;
    ARRAY<STRUCTURE<TV>*> collision_structures;
    ARRAY<int> particle_to_structure;
    ARRAY<HASHTABLE<COLLISION_GEOMETRY_ID> > structure_collide_collision_body;  // whether each structure collides with each collision body
    bool collisions_on;
    ARRAY<ARRAY<int>,COLLISION_GEOMETRY_ID> collision_body_candidate_nodes;
    ARRAY<int> ignored_nodes; // must be set before Initialize_Object_Collisions is called
    ARRAY<T>* collision_tolerances; // must be set before Initialize_Object_Collisions is called
    HASHTABLE<int,T> *thickness_table,*friction_table;

    OPERATION_HASH<COLLISION_GEOMETRY_ID> get_potential_collisions_already_added;
private:
    bool use_structure_collide_collision_body;
public:

    DEFORMABLE_OBJECT_COLLISIONS(PARTICLES<TV>& particles,DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection,ARRAY<STRUCTURE<TV>*>& deformable_object_structures,
        COLLISION_GEOMETRY_COLLECTION<TV>& collision_body_list);

    virtual ~DEFORMABLE_OBJECT_COLLISIONS();

public:
    void Activate_Collisions(const bool collisions_on_input=true)
    {collisions_on=collisions_on_input;}

//#####################################################################
    void Initialize_Object_Collisions(const bool collide_with_interior=false,const T collision_tolerance=(T)1e-6,const bool use_spatial_partition_for_levelset_collisions=true,
        const bool disable_multiple_levelset_collisions=true,const T maximum_levelset_collision_projection_velocity=FLT_MAX);
    void Use_Structure_Collide_Collision_Body(const bool value=true);
    void Compute_Candidate_Nodes_For_Collision_Body_Collisions(const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>& bodies);
    int Adjust_Nodes_For_Collision_Body_Collisions(BINDING_LIST<TV>& binding_list,SOFT_BINDINGS<TV>& soft_bindings,ARRAY_VIEW<const TV> X_old,const T dt,
        const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>* bodies);
    int Adjust_Existing_Nodes_For_Collision_Body_Collisions(BINDING_LIST<TV>& binding_list,SOFT_BINDINGS<TV>& soft_bindings,ARRAY_VIEW<const TV> X_old,const T dt,
        const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>* bodies);
    void Set_Collision_Velocities(ARRAY_VIEW<TV> V); // for external forces and velocities
    void Zero_Out_Collision_Velocities(ARRAY_VIEW<TV> V); // for external forces and velocities
    void Reset_Object_Collisions();
    void Update_Simulated_Particles();
private:
    template<class T_MESH> void Add_Collision_Mesh(T_MESH& mesh,const bool collide_with_interior);
//#####################################################################
};
}
#endif

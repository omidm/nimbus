//#####################################################################
// Copyright 2011, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __THREADED_RIGIDS__
#define __THREADED_RIGIDS__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Parallel_Computation/PARTITION_ID.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_PACKAGE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Collision_Detection/COLLISION_GEOMETRY_SPATIAL_PARTITION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#ifdef USE_PTHREADS
    #include <pthread.h>
    #include <PhysBAM_Tools/Parallel_Computation/PTHREAD.h>
#endif

namespace PhysBAM{

template<class TV>
class THREADED_RIGIDS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename TV::SPIN T_SPIN;

public:
    ARRAY<THREAD_PACKAGE>& buffers; // global memory for exchanges
    int tid;
    int number_of_processors;

    ARRAY<ARRAY<int>,PARTITION_ID> particles_of_partition;
    ARRAY<PARTITION_ID> partition_id_from_particle_index;
    ARRAY<ARRAY<T> > boundaries;
    ARRAY<int,VECTOR<int,TV::dimension> > partition_coordinates_to_rank;
    TV_INT processes_per_dimension;
    RANGE<TV> global_domain;
    RANGE<TV> local_domain;

#ifdef USE_PTHREADS
    mutable pthread_mutex_t* lock;
    mutable pthread_barrier_t* barr;
#endif

    THREADED_RIGIDS(ARRAY<THREAD_PACKAGE>& buffers_input,const int tid_input,const int number_of_threads);
    ~THREADED_RIGIDS();

    bool Is_Real_Body(int p) const
    {return partition_id_from_particle_index(p)==PARTITION_ID(tid);}

    bool Is_Dynamic_Ghost_Body(RIGID_BODY<TV>& body) const
    {return (partition_id_from_particle_index(body.particle_index)!=PARTITION_ID(tid)) &&
            (partition_id_from_particle_index(body.particle_index)!=PARTITION_ID(particles_of_partition.m)) && body.Is_Simulated();}

//#####################################################################
    void Initialize(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,TV_INT processes_per_dimension_input,
        COLLISION_GEOMETRY_SPATIAL_PARTITION<COLLISION_GEOMETRY<TV>,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>& spatial_partition);
    void Broadcast_Positions(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input);
    void Update_Partitions(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,
        COLLISION_GEOMETRY_SPATIAL_PARTITION<COLLISION_GEOMETRY<TV>,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>& spatial_partition);
    RANGE<TV> Split_Range(TV_INT& processes_per_dimension,const RANGE<TV>& global_range);
    VECTOR<int,1> Split_Range(const VECTOR<T,1>& global_range,const VECTOR<int,1>& processes_per_dimension);
    VECTOR<int,2> Split_Range(const VECTOR<T,2>& global_range,const VECTOR<int,2>& processes_per_dimension);
    VECTOR<int,3> Split_Range(const VECTOR<T,3>& global_range,const VECTOR<int,3>& processes_per_dimension);
    void Split_Dimension(const T x,const int processes,ARRAY<T>& boundaries);
    void Simple_Partition(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,
        COLLISION_GEOMETRY_SPATIAL_PARTITION<COLLISION_GEOMETRY<TV>,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>& spatial_partition);
    void Clear_Impulse_Accumulators(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input);
    void Exchange_All_Impulses(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<TWIST<TV> >& mpi_rigid_velocity_save,ARRAY<T_SPIN>& mpi_rigid_angular_momentum_save,
        RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,const bool euler_step_with_new_velocity,const T dt,const T time);
    void Exchange_All_Pushes(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<TV>& mpi_rigid_X_save,ARRAY<ROTATION<TV> >& mpi_rigid_rotation_save,
        RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions);
    void Prune_And_Exchange_Impulses(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<bool>& need_to_reevolve);
    void Exchange_Bounding_Box_Collision_Pairs(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<VECTOR<int,2> >& pairs,const bool contact_pairs);
    void Synchronize_Dt(T& dt) const;
    template<class T_DATA> T_DATA Reduce_Min(const T_DATA local_value) const;
    template<class T_DATA> T_DATA Reduce_Max(const T_DATA local_value) const;
    template<class T_DATA> T_DATA Reduce_Add(const T_DATA local_value) const;
//#####################################################################
};
}
#endif

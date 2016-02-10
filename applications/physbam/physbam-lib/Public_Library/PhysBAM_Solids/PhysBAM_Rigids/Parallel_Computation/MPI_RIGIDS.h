//#####################################################################
// Copyright 2010, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPI_RIGIDS__
#define __MPI_RIGIDS__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Parallel_Computation/PARTITION_ID.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_PACKAGE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Collision_Detection/COLLISION_GEOMETRY_SPATIAL_PARTITION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGIDS_ONLY_COLLISION_CALLBACKS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Parallel_Computation/THREADED_RIGIDS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>

namespace MPI{class Group;class Intracomm;class Request;class Status;class Op;}
namespace PhysBAM{

template<class TV>
class MPI_RIGIDS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename TV::SPIN T_SPIN;
    
public:
    ARRAY<ARRAY<int>,PARTITION_ID> particles_of_partition;
    ARRAY<PARTITION_ID> partition_id_from_particle_index;
    ARRAY<ARRAY<T> > boundaries;
    ARRAY<int,VECTOR<int,TV::dimension> > partition_coordinates_to_rank;
    TV_INT processes_per_dimension;
    RANGE<TV> global_domain;
    RANGE<TV> local_domain;

    MPI::Intracomm* comm;
    MPI::Group* group;
    mutable int current_tag;
public:
    int rank; // my 0-indexed processor
    int number_of_processors;
    bool recompute;

    THREADED_RIGIDS<TV>* threaded_rigids;

    MPI_RIGIDS();
    MPI_RIGIDS(ARRAY<THREAD_PACKAGE>& buffers_input,const int tid_input,const int nthreads);
    ~MPI_RIGIDS();

    int Get_Unique_Tag() const
    {if(threaded_rigids) PHYSBAM_NOT_IMPLEMENTED();
    else return current_tag=max(32,(current_tag+1)&((1<<15)-1));}

    bool Is_Real_Body(int p) const
    {if(threaded_rigids) return threaded_rigids->Is_Real_Body(p);
    else return partition_id_from_particle_index(p)==PARTITION_ID(rank+1);}

    bool Is_Dynamic_Ghost_Body(RIGID_BODY<TV>& body) const
    {if(threaded_rigids) return threaded_rigids->Is_Dynamic_Ghost_Body(body);
    else return (partition_id_from_particle_index(body.particle_index)!=PARTITION_ID(rank+1)) &&
            (partition_id_from_particle_index(body.particle_index)!=PARTITION_ID(particles_of_partition.m)) && body.Is_Simulated();}

//#####################################################################

    void Initialize(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,TV_INT processes_per_dimension_input,
        COLLISION_GEOMETRY_SPATIAL_PARTITION<COLLISION_GEOMETRY<TV>,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>& spatial_partition)
    {
        if(threaded_rigids) threaded_rigids->Initialize(rigid_body_collection_input,processes_per_dimension_input,spatial_partition);
        else Initialize_Helper(rigid_body_collection_input,processes_per_dimension_input,spatial_partition);
    }

    void Broadcast_Positions(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
    {
        if(threaded_rigids) threaded_rigids->Broadcast_Positions(rigid_body_collection_input);
        else Broadcast_Positions_Helper(rigid_body_collection_input);
    }

    void Update_Partitions(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,
        COLLISION_GEOMETRY_SPATIAL_PARTITION<COLLISION_GEOMETRY<TV>,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>& spatial_partition)
    {
        if(threaded_rigids) threaded_rigids->Update_Partitions(rigid_body_collection_input,spatial_partition);
        else Update_Partitions_Helper(rigid_body_collection_input,spatial_partition);
    }

    RANGE<TV> Split_Range(TV_INT& processes_per_dimension,const RANGE<TV>& global_range)
    {
        if(threaded_rigids) return threaded_rigids->Split_Range(processes_per_dimension,global_range);
        else return Split_Range_Helper(processes_per_dimension,global_range);
    }

    VECTOR<int,1> Split_Range(const VECTOR<T,1>& global_range,const VECTOR<int,1>& processes_per_dimension)
    {
        if(threaded_rigids) return threaded_rigids->Split_Range(global_range,processes_per_dimension);
        else return Split_Range_Helper(global_range,processes_per_dimension);
    }

    VECTOR<int,2> Split_Range(const VECTOR<T,2>& global_range,const VECTOR<int,2>& processes_per_dimension)
    {
        if(threaded_rigids) return threaded_rigids->Split_Range(global_range,processes_per_dimension);
        else return Split_Range_Helper(global_range,processes_per_dimension);
    }

    VECTOR<int,3> Split_Range(const VECTOR<T,3>& global_range,const VECTOR<int,3>& processes_per_dimension)
    {
        if(threaded_rigids) return threaded_rigids->Split_Range(global_range,processes_per_dimension);
        else return Split_Range_Helper(global_range,processes_per_dimension);
    }

    void Split_Dimension(const T x,const int processes,ARRAY<T>& boundaries)
    {
        if(threaded_rigids) threaded_rigids->Split_Dimension(x,processes,boundaries);
        else Split_Dimension_Helper(x,processes,boundaries);
    }

    void Simple_Partition(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,
        COLLISION_GEOMETRY_SPATIAL_PARTITION<COLLISION_GEOMETRY<TV>,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>& spatial_partition)
    {
        if(threaded_rigids) threaded_rigids->Simple_Partition(rigid_body_collection_input,spatial_partition);
        else Simple_Partition_Helper(rigid_body_collection_input,spatial_partition);
    }

    void Clear_Impulse_Accumulators(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
    {
        if(threaded_rigids) threaded_rigids->Clear_Impulse_Accumulators(rigid_body_collection_input);
        else Clear_Impulse_Accumulators_Helper(rigid_body_collection_input);
    }

    void Exchange_All_Impulses(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<TWIST<TV> >& mpi_rigid_velocity_save,ARRAY<T_SPIN>& mpi_rigid_angular_momentum_save,
        RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,const bool euler_step_with_new_velocity,const T dt,const T time)
    {
        if(threaded_rigids) threaded_rigids->Exchange_All_Impulses(rigid_body_collection,mpi_rigid_velocity_save,mpi_rigid_angular_momentum_save,rigid_body_collisions,euler_step_with_new_velocity,dt,time);
        else Exchange_All_Impulses_Helper(rigid_body_collection,mpi_rigid_velocity_save,mpi_rigid_angular_momentum_save,rigid_body_collisions,euler_step_with_new_velocity,dt,time);
    }

    void Exchange_All_Pushes(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<TV>& mpi_rigid_X_save,ARRAY<ROTATION<TV> >& mpi_rigid_rotation_save,
        RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions)
    {
        if(threaded_rigids) threaded_rigids->Exchange_All_Pushes(rigid_body_collection,mpi_rigid_X_save,mpi_rigid_rotation_save,rigid_body_collisions);
        else Exchange_All_Pushes_Helper(rigid_body_collection,mpi_rigid_X_save,mpi_rigid_rotation_save,rigid_body_collisions);
    }

    void Prune_And_Exchange_Impulses(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<bool>& need_to_reevolve)
    {
        if(threaded_rigids) threaded_rigids->Prune_And_Exchange_Impulses(rigid_body_collection,need_to_reevolve);
        else Prune_And_Exchange_Impulses_Helper(rigid_body_collection,need_to_reevolve);
    }

    void Exchange_Bounding_Box_Collision_Pairs(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<VECTOR<int,2> >& pairs,const bool contact_pairs)
    {
        if(threaded_rigids) threaded_rigids->Exchange_Bounding_Box_Collision_Pairs(rigid_body_collection,pairs,contact_pairs);
        else Exchange_Bounding_Box_Collision_Pairs_Helper(rigid_body_collection,pairs,contact_pairs);
    }

    void Synchronize_Dt(T& dt) const
    {
        if(threaded_rigids) threaded_rigids->Synchronize_Dt(dt);
        else Synchronize_Dt_Helper(dt);
    }

    template<class T_DATA> T_DATA Reduce_Min(const T_DATA local_value) const
    {
        if(threaded_rigids) return threaded_rigids->Reduce_Min(local_value);
        else return Reduce_Min_Helper(local_value);
    }

    template<class T_DATA> T_DATA Reduce_Max(const T_DATA local_value) const
    {
        if(threaded_rigids) return threaded_rigids->Reduce_Max(local_value);
        else return Reduce_Max_Helper(local_value);
    }

    template<class T_DATA> T_DATA Reduce_Add(const T_DATA local_value) const
    {
        if(threaded_rigids) return threaded_rigids->Reduce_Add(local_value);
        else return Reduce_Add_Helper(local_value);
    }

//#####################################################################
    void Initialize_Helper(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,TV_INT processes_per_dimension_input,
        COLLISION_GEOMETRY_SPATIAL_PARTITION<COLLISION_GEOMETRY<TV>,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>& spatial_partition);
    void Broadcast_Positions_Helper(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input);
    void Update_Partitions_Helper(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,
        COLLISION_GEOMETRY_SPATIAL_PARTITION<COLLISION_GEOMETRY<TV>,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>& spatial_partition);
    RANGE<TV> Split_Range_Helper(TV_INT& processes_per_dimension,const RANGE<TV>& global_range);
    VECTOR<int,1> Split_Range_Helper(const VECTOR<T,1>& global_range,const VECTOR<int,1>& processes_per_dimension);
    VECTOR<int,2> Split_Range_Helper(const VECTOR<T,2>& global_range,const VECTOR<int,2>& processes_per_dimension);
    VECTOR<int,3> Split_Range_Helper(const VECTOR<T,3>& global_range,const VECTOR<int,3>& processes_per_dimension);
    void Split_Dimension_Helper(const T x,const int processes,ARRAY<T>& boundaries);
    void Simple_Partition_Helper(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,
        COLLISION_GEOMETRY_SPATIAL_PARTITION<COLLISION_GEOMETRY<TV>,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>& spatial_partition);
    void Clear_Impulse_Accumulators_Helper(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input);
    void Exchange_All_Impulses_Helper(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<TWIST<TV> >& mpi_rigid_velocity_save,ARRAY<T_SPIN>& mpi_rigid_angular_momentum_save,
        RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,const bool euler_step_with_new_velocity,const T dt,const T time);
    void Exchange_All_Pushes_Helper(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<TV>& mpi_rigid_X_save,ARRAY<ROTATION<TV> >& mpi_rigid_rotation_save,
        RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions);
    void Prune_And_Exchange_Impulses_Helper(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<bool>& need_to_reevolve);
    void Exchange_Bounding_Box_Collision_Pairs_Helper(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<VECTOR<int,2> >& pairs,const bool contact_pairs);
    void Synchronize_Dt_Helper(T& dt) const;
    template<class T_DATA> T_DATA Reduce_Min_Helper(const T_DATA local_value) const;
    template<class T_DATA> T_DATA Reduce_Max_Helper(const T_DATA local_value) const;
    template<class T_DATA> T_DATA Reduce_Add_Helper(const T_DATA local_value) const;
//#####################################################################
};
}
#endif

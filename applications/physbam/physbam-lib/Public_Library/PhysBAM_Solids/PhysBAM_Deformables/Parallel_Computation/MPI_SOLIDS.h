//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPI_SOLIDS__
#define __MPI_SOLIDS__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Parallel_Computation/PARTITION_ID.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#ifdef Status // Splendid choice for a macro.
#undef Status
#endif
namespace MPI{class Group;class Intracomm;class Request;class Status;class Op;}
namespace PhysBAM{

class SEGMENT_MESH;
class MPI_PACKAGE;
template<class T> class TRIANGLE_REPULSIONS;
template<class TV> class DEFORMABLE_BODY_COLLECTION;
template<class TV> class RIGID_GEOMETRY_COLLECTION;

template<class TV>
class MPI_SOLIDS:public NONCOPYABLE
{
public:
    DEFORMABLE_BODY_COLLECTION<TV>* deformable_body_collection;
    ARRAY<ARRAY<int>,PARTITION_ID> particles_of_partition;
    ARRAY<PARTITION_ID> partition_id_from_particle_index;

    struct MPI_PARTITION{
        ARRAY<ARRAY<int>,PARTITION_ID> ghost_dynamic_particles;
        ARRAY<ARRAY<int>,PARTITION_ID> boundary_dynamic_particles;
        ARRAY<PARTITION_ID> neighbor_partitions;
    };
    MPI_PARTITION mpi_partition_binding,mpi_partition_force;

    // boundary_dynamic_particles, neighbor_partitions, ghost_dynamic_particles, neighbor_ranks, comm

    MPI::Intracomm* comm;
    MPI::Group* group;
    mutable int current_tag;
public:
    int rank; // my 0-indexed processor
    int number_of_ranks;

    static PARTITION_ID Rank_To_Partition(int rank_input)
    {return PARTITION_ID(rank_input+1);}

    static int Partition_To_Rank(PARTITION_ID partition)
    {return Value(partition)-1;}

    PARTITION_ID Partition() const
    {return Rank_To_Partition(rank);}

    PARTITION_ID Number_Of_Partitions() const
    {return PARTITION_ID(number_of_ranks);}

    int Number_Of_Processors() const
    {return number_of_ranks;}

    MPI_SOLIDS();
    ~MPI_SOLIDS();

    int Get_Unique_Tag() const
    {return current_tag=max(32,(current_tag+1)&((1<<15)-1));}

    template<class T1,class T2> void
    Gather_Data(ARRAY<T1>& array1,ARRAY<T2>& array2)
    {Gather_Data(ARRAY_VIEW<T1>(array1),ARRAY_VIEW<T2>(array2));}
    
    template<class T1,class T2> void
    Scatter_Data(ARRAY<T1>& array1,ARRAY<T2>& array2)
    {Scatter_Data(ARRAY_VIEW<T1>(array1),ARRAY_VIEW<T2>(array2));}

    template<class T1> void
    Scatter_Data(ARRAY<T1>& array1)
    {Scatter_Data(ARRAY_VIEW<T1>(array1));}

    template<class T1,class T2> void
    Broadcast_Data(ARRAY<T1>& array1,ARRAY<T2>& array2)
    {Broadcast_Data(ARRAY_VIEW<T1>(array1),ARRAY_VIEW<T2>(array2));}

    template<class T1> void
    Broadcast_Data(ARRAY<T1>& array1)
    {Broadcast_Data(ARRAY_VIEW<T1>(array1));}

    template<class T1,class T2> void
    Exchange_Data(ARRAY<T1>& array1,ARRAY<T2>& array2)
    {Exchange_Data(ARRAY_VIEW<T1>(array1),ARRAY_VIEW<T2>(array2));}

    template<class T1,class T2> void
    Gather_Data(ARRAY_VIEW<T1> array1,ARRAY_VIEW<T2> array2)
    {ARRAY<MPI_PACKAGE*> packages;packages.Append(Package(array1));packages.Append(Package(array2));Gather_Helper(packages);}

    // Scatters all data from the root to the proper processing partitions on other processors
    template<class T1,class T2> void
    Scatter_Data(ARRAY_VIEW<T1> array1,ARRAY_VIEW<T2> array2)
    {ARRAY<MPI_PACKAGE*> packages;packages.Append(Package(array1));packages.Append(Package(array2));Scatter_Helper(packages);}

    template<class T1> void
    Scatter_Data(ARRAY_VIEW<T1> array1)
    {ARRAY<MPI_PACKAGE*> packages;packages.Append(Package(array1));Scatter_Helper(packages);}

    // Broadcasts all data to every processor
    template<class T1,class T2> void
    Broadcast_Data(ARRAY_VIEW<T1> array1,ARRAY_VIEW<T2> array2)
    {ARRAY<MPI_PACKAGE*> packages;packages.Append(Package(array1));packages.Append(Package(array2));Broadcast_Helper(packages);}

    template<class T1> void
    Broadcast_Data(ARRAY_VIEW<T1> array1)
    {ARRAY<MPI_PACKAGE*> packages;packages.Append(Package(array1));Broadcast_Helper(packages);}

    // Exchanges nodes required for force computation
    template<class T1,class T2> void
    Exchange_Data(ARRAY_VIEW<T1> array1,ARRAY_VIEW<T2> array2)
    {ARRAY<MPI_PACKAGE*> packages;packages.Append(Package(array1));packages.Append(Package(array2));Exchange_Helper(packages,mpi_partition_force);}

    template<class T>
    void Exchange_Force_Boundary_Data(ARRAY_VIEW<T> array) const
    {Exchange_Boundary_Data(array,mpi_partition_force);}

    template<class T>
    void Exchange_Binding_Boundary_Data(ARRAY_VIEW<T> array) const
    {Exchange_Boundary_Data(array,mpi_partition_binding);}

    template<class T> void Exchange_Binding_Boundary_Data_Global(ARRAY<T>& array) const
    {Exchange_Binding_Boundary_Data_Global(ARRAY_VIEW<T>(array));}

    template<class T> void Exchange_Binding_Boundary_Data(ARRAY<T>& array) const
    {Exchange_Binding_Boundary_Data(ARRAY_VIEW<T>(array));}

    template<class T> void Exchange_Force_Boundary_Data_Global(ARRAY<T>& array) const
    {Exchange_Force_Boundary_Data_Global(ARRAY_VIEW<T>(array));}

    template<class T> void Exchange_Force_Boundary_Data(ARRAY<T>& array) const
    {Exchange_Force_Boundary_Data(ARRAY_VIEW<T>(array));}

    void KD_Tree_Partition(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection_input,RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection_input,const ARRAY<TV>& X)
    {KD_Tree_Partition(deformable_body_collection_input,rigid_geometry_collection_input,ARRAY_VIEW<const TV>(X));}

//#####################################################################
    void Barrier() const;
    void Update(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_connection_input,const SEGMENT_MESH& force_connectivity,const SEGMENT_MESH& binding_connectivity);
    void Compute_Mpi_Partition(MPI_PARTITION& mpi_partition,const SEGMENT_MESH& connectivity);
    template<class T> void Exchange_Boundary_Data(ARRAY_VIEW<T> array,const MPI_PARTITION& mpi_partition) const;
    template<class T> void Exchange_Boundary_Data_Global(ARRAY_VIEW<T> array,const MPI_PARTITION& mpi_partition) const;
    template<class T> void Exchange_Force_Boundary_Data_Global(ARRAY_VIEW<T> array) const;
    template<class T> void Exchange_Binding_Boundary_Data_Global(ARRAY_VIEW<T> array) const;
    template<class T> void Gather_Data(ARRAY_VIEW<T> array) const;
    template<class T> void Broadcast(T* data,const int length) const;
    void All_Gather_Particles(ARRAY_VIEW<TV> X,ARRAY_VIEW<TV> V) const;
    template<class T_PAIR1,class T_PAIR2> void Gather_Interaction_Pairs(ARRAY<T_PAIR1>& point_triangle_pairs,ARRAY<T_PAIR2>& edge_edge_pairs) const;
    void Broadcast_Collision_Modified_Data(ARRAY_VIEW<bool> modified,ARRAY_VIEW<bool> recently_modified,ARRAY_VIEW<TV> X,ARRAY_VIEW<TV> V) const;
    void Gather_Collision_Modified_Data(ARRAY_VIEW<bool> modified,ARRAY_VIEW<bool> recently_modified,ARRAY_VIEW<TV> X,ARRAY_VIEW<TV> V) const;
    template<int d1,int d2> void All_Gather_Intersecting_Pairs(HASHTABLE<VECTOR<int,d1> >& intersecting_point_face_pairs,HASHTABLE<VECTOR<int,d2> >& intersecting_edge_edge_pairs);
    void All_Scatter_Adhesion_Pairs(ARRAY<ARRAY<PAIR<VECTOR<int,2>,VECTOR<typename TV::SCALAR,2> > >,PARTITION_ID>& pairs_to_scatter,
        ARRAY<ARRAY<PAIR<VECTOR<int,2>,VECTOR<typename TV::SCALAR,2> > >,PARTITION_ID>& pairs_received);
private:
    template<class T_DATA> MPI_PACKAGE* Package(ARRAY_VIEW<T_DATA> array);
    void Gather_Helper(const ARRAY<MPI_PACKAGE*>& packages) const;
    void Scatter_Helper(const ARRAY<MPI_PACKAGE*>& packages) const;
    void Broadcast_Helper(const ARRAY<MPI_PACKAGE*>& packages) const;
    void Exchange_Helper(const ARRAY<MPI_PACKAGE*>& packages,MPI_PARTITION& mpi_partition) const;
public:
    template<class T_DATA> T_DATA Reduce_Add_Global(const T_DATA local_value) const;
    template<class T_DATA> T_DATA Reduce_Max_Global(const T_DATA local_value) const;
    template<class T_DATA> T_DATA Reduce_Min_Global(const T_DATA local_value) const;
    template<class T_DATA> T_DATA Reduce_Add(const T_DATA local_value) const;
    template<class T_DATA> T_DATA Reduce_Max(const T_DATA local_value) const;
    template<class T_DATA> T_DATA Reduce_Min(const T_DATA local_value) const;
    template<class T_DATA> T_DATA Debug_Reduce_Array_Sum(ARRAY_VIEW<const T_DATA> values) const; // sends all values to one processor and adds them up in order (very slow)
    template<class T_ARRAY_PAIR> void Distribute_Repulsion_Pairs(const T_ARRAY_PAIR& pairs,ARRAY<ARRAY<int>,PARTITION_ID>& send_particles,
        ARRAY<ARRAY<int>,PARTITION_ID>& receive_particles,T_ARRAY_PAIR& my_boundary_pairs,T_ARRAY_PAIR& my_internal_pairs) const;
    void Gather_Repulsion_Inputs(ARRAY_VIEW<TV> X_self_collision_free,ARRAY_VIEW<TV> X,ARRAY_VIEW<TV> V,ARRAY<ARRAY<int>,PARTITION_ID>& send_particles,
        ARRAY<ARRAY<int>,PARTITION_ID>& receive_particles) const;
    void Scatter_Repulsion_Outputs(ARRAY_VIEW<TV> X_self_collision_free,ARRAY_VIEW<TV> X,ARRAY_VIEW<TV> V,ARRAY<ARRAY<int>,PARTITION_ID>& send_particles,
        ARRAY<ARRAY<int>,PARTITION_ID>& receive_particles) const;
    void Simple_Partition(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection_input,RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection_input,ARRAY_VIEW<const TV> X,const VECTOR<int,TV::m>& counts);
    void KD_Tree_Partition(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection_input,RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection_input,ARRAY_VIEW<const TV> X);
    template<class ID> void KD_Tree_Partition_Subset(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection_input,RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection_input,const ARRAY<int,ID>& nodes,ARRAY_VIEW<const TV> X);
    void Print_Diagnostics() const;
private:
    const MPI::Intracomm& Comm() const;
//#####################################################################
};
}
#endif

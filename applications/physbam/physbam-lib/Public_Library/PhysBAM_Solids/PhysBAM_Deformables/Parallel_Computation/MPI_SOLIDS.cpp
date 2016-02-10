//#####################################################################
// Copyright 2006-2009, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Arrays_Computations/SUMMATIONS.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Data_Structures/KD_TREE.h>
#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/HAIR_ID.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <stdexcept>
#ifdef USE_MPI
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_REPULSION_UTILITIES.h>
#endif
using namespace PhysBAM;

#ifdef USE_MPI

//#####################################################################
// Function MPI_SOLIDS
//#####################################################################
template<class TV> MPI_SOLIDS<TV>::
MPI_SOLIDS()
    :deformable_body_collection(0),current_tag(0)
{
    LOG::SCOPE scope("MPI INITIALIZE","Initializing MPI");
    number_of_ranks=MPI::COMM_WORLD.Get_size();
    {std::stringstream ss;ss<<"number of ranks = "<<number_of_ranks<<std::endl;LOG::filecout(ss.str());}

    comm=new MPI::Intracomm;
    group=new MPI::Group;
    *comm=MPI::COMM_WORLD.Dup();
    *group=comm->Get_group();
    rank=comm->Get_rank();
}
//#####################################################################
// Function ~MPI_SOLIDS
//#####################################################################
template<class TV> MPI_SOLIDS<TV>::
~MPI_SOLIDS()
{
    if(comm){comm->Free();delete comm;}
    if(group){group->Free();delete group;}
}
//#####################################################################
// Function Barrier
//#####################################################################
template<class TV> void MPI_SOLIDS<TV>::
Barrier() const
{comm->Barrier();}
//#####################################################################
// Function Update
//#####################################################################
template<class TV> void MPI_SOLIDS<TV>::
Update(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection_input,const SEGMENT_MESH& force_connectivity,const SEGMENT_MESH& binding_connectivity)
{
    deformable_body_collection=&deformable_body_collection_input;
    Compute_Mpi_Partition(mpi_partition_force,force_connectivity);
    Compute_Mpi_Partition(mpi_partition_binding,binding_connectivity);
}
//#####################################################################
// Function Compute_Mpi_Partition
//#####################################################################
template<class TV> void MPI_SOLIDS<TV>::
Compute_Mpi_Partition(MPI_PARTITION& mpi_partition,const SEGMENT_MESH& connectivity)
{
    PHYSBAM_ASSERT(connectivity.neighbor_nodes);
    mpi_partition.boundary_dynamic_particles.Remove_All();
    mpi_partition.boundary_dynamic_particles.Resize(Number_Of_Partitions());
    mpi_partition.ghost_dynamic_particles.Remove_All();
    mpi_partition.ghost_dynamic_particles.Resize(Number_Of_Partitions());
    mpi_partition.neighbor_partitions.Remove_All();

    ARRAY<bool> done(partition_id_from_particle_index.Size());
    HASHTABLE<PAIR<int,PARTITION_ID> > pair_done;
    for(int j=1;j<=particles_of_partition(Partition()).m;j++){int p=particles_of_partition(Partition())(j);
        if(p<=deformable_body_collection->particles.array_collection->Size()){
            const ARRAY<int>& neighbor_particles=(*connectivity.neighbor_nodes)(p);
            for(int i=1;i<=neighbor_particles.m;i++){int n=neighbor_particles(i);
                PARTITION_ID np=partition_id_from_particle_index(n);
                if(np==Partition()) continue;
                if(pair_done.Set(PAIR<int,PARTITION_ID>(p,np)))
                    mpi_partition.boundary_dynamic_particles(np).Append(p);
                if(done(n)) continue;
                done(n)=true;
                mpi_partition.ghost_dynamic_particles(np).Append(n);}}}
    for(PARTITION_ID p(1);p<=mpi_partition.ghost_dynamic_particles.m;p++){
        if(mpi_partition.ghost_dynamic_particles(p).m)
            mpi_partition.neighbor_partitions.Append(p);
        Sort(mpi_partition.ghost_dynamic_particles(p));
        Sort(mpi_partition.boundary_dynamic_particles(p));}
}
//#####################################################################
// Function Exchange_Boundary_Data
//#####################################################################
template<class TV> template<class T_DATA> void MPI_SOLIDS<TV>::
Exchange_Boundary_Data(ARRAY_VIEW<T_DATA> array,const MPI_PARTITION& mpi_partition) const
{
    int tag=1;
    ARRAY<MPI_PACKAGE> packages;ARRAY<MPI::Request> requests;
    // send
    for(int i=1;i<=mpi_partition.neighbor_partitions.m;i++){int neighbor_rank=Partition_To_Rank(mpi_partition.neighbor_partitions(i));
        MPI_PACKAGE package(array.Subset(mpi_partition.boundary_dynamic_particles(mpi_partition.neighbor_partitions(i))));
        packages.Append(package);requests.Append(package.Isend(*comm,neighbor_rank,tag));}
    // receive
    for(int i=1;i<=mpi_partition.neighbor_partitions.m;i++){int neighbor_rank=Partition_To_Rank(mpi_partition.neighbor_partitions(i));
        MPI_PACKAGE package(array.Subset(mpi_partition.ghost_dynamic_particles(mpi_partition.neighbor_partitions(i))));
        packages.Append(package);requests.Append(package.Irecv(*comm,neighbor_rank,tag));}
    // finish
    MPI_UTILITIES::Wait_All(requests);MPI_PACKAGE::Free_All(packages);
}
//#####################################################################
// Function Exchange_Boundary_Data
//#####################################################################
template<class TV> template<class T_DATA> void MPI_SOLIDS<TV>::
Exchange_Binding_Boundary_Data_Global(ARRAY_VIEW<T_DATA> array) const // TODO: I am broken.
{
    Exchange_Boundary_Data(array,mpi_partition_binding);
}
//#####################################################################
// Function Exchange_Boundary_Data
//#####################################################################
template<class TV> template<class T_DATA> void MPI_SOLIDS<TV>::
Exchange_Force_Boundary_Data_Global(ARRAY_VIEW<T_DATA> array) const // TODO: I am broken.
{
    Exchange_Boundary_Data(array,mpi_partition_force);
}
//#####################################################################
// Function Gather_Data
//#####################################################################
template<class TV> template<class T_DATA> void MPI_SOLIDS<TV>::
Gather_Data(ARRAY_VIEW<T_DATA> array) const
{
    ARRAY<T_DATA> data_of_partition(array.Subset(particles_of_partition(Partition())));
    if(rank==0){
        // Count others
        ARRAY<int,PARTITION_ID> counts(particles_of_partition.m),offsets(particles_of_partition.m);
        for(PARTITION_ID p(1);p<=particles_of_partition.m;p++){
            counts(p)=particles_of_partition(p).m;
            if(p>PARTITION_ID(1)) offsets(p)=offsets(p-1)+counts(p-1);}

        ARRAY<T_DATA> all_data(offsets(offsets.Size())+counts(counts.Size()),false);
        comm->Gatherv(data_of_partition.Get_Array_Pointer(),data_of_partition.m,MPI_UTILITIES::Datatype<T_DATA>(),all_data.Get_Array_Pointer(),counts.Get_Array_Pointer(),
            offsets.Get_Array_Pointer(),MPI_UTILITIES::Datatype<T_DATA>(),0);

        for(PARTITION_ID p(1);p<=particles_of_partition.m;p++){
            array.Subset(particles_of_partition(p))=all_data.Subset(IDENTITY_ARRAY<>(counts(p))+offsets(p));}}
    else comm->Gatherv(data_of_partition.Get_Array_Pointer(),data_of_partition.m,MPI_UTILITIES::Datatype<T_DATA>(),0,0,0,MPI_UTILITIES::Datatype<T_DATA>(),0);
}
//#####################################################################
// Function Broadcast
//#####################################################################
template<class TV> template<class T_DATA> void MPI_SOLIDS<TV>::
Broadcast(T_DATA* data,const int length) const
{
    comm->Bcast(data,length,MPI_UTILITIES::Datatype<T_DATA>(),0);
}
//#####################################################################
// Function All_Gather_Particles
//#####################################################################
template<class TV> void MPI_SOLIDS<TV>::
All_Gather_Particles(ARRAY_VIEW<TV> X,ARRAY_VIEW<TV> V) const
{
    ARRAY<TV> data_of_partition(X.Subset(particles_of_partition(Partition())));
    data_of_partition.Append_Elements(V.Subset(particles_of_partition(Partition())));

    ARRAY<int,PARTITION_ID> counts(particles_of_partition.m),offsets(particles_of_partition.m);
    for(PARTITION_ID p(1);p<=particles_of_partition.m;p++){
        counts(p)=particles_of_partition(p).m*2;
        if(p>PARTITION_ID(1)) offsets(p)=offsets(p-1)+counts(p-1);}

    ARRAY<TV> all_data(offsets(offsets.Size())+counts(counts.Size()),false);
    comm->Allgatherv(data_of_partition.Get_Array_Pointer(),data_of_partition.m,MPI_UTILITIES::Datatype<TV>(),all_data.Get_Array_Pointer(),counts.Get_Array_Pointer(),
        offsets.Get_Array_Pointer(),MPI_UTILITIES::Datatype<TV>());

    for(PARTITION_ID p(1);p<=particles_of_partition.m;p++){
        int half_count=counts(p)/2;
        X.Subset(particles_of_partition(p))=all_data.Subset(IDENTITY_ARRAY<>(half_count)+offsets(p));
        V.Subset(particles_of_partition(p))=all_data.Subset(IDENTITY_ARRAY<>(half_count)+(offsets(p)+half_count));}
}
//#####################################################################
// Function Gather_Interaction_Pairs
//#####################################################################
namespace{
template<class T,class T_PAIR1,class T_PAIR2> void Gather_Interaction_Pairs_Helper(const MPI_SOLIDS<VECTOR<T,1> >& mpi_solids,ARRAY<T_PAIR1>& point_triangle_pairs,ARRAY<T_PAIR2>& edge_edge_pairs)
{PHYSBAM_NOT_IMPLEMENTED();}

template<class TV,class T_PAIR1,class T_PAIR2> void Gather_Interaction_Pairs_Helper(const MPI_SOLIDS<TV>& mpi_solids,ARRAY<T_PAIR1>& point_triangle_pairs,ARRAY<T_PAIR2>& edge_edge_pairs)
{
    int tag=mpi_solids.Get_Unique_Tag();
    if(mpi_solids.rank!=0){
        ARRAY<char> buffer(1+MPI_UTILITIES::Pack_Size(point_triangle_pairs,edge_edge_pairs,*mpi_solids.comm),false);int position=0; // TODO: remove +1 when LAM fixes assert
        MPI_UTILITIES::Pack(point_triangle_pairs,edge_edge_pairs,buffer,position,*mpi_solids.comm);
        mpi_solids.comm->Send(buffer.Get_Array_Pointer(),position,MPI::PACKED,0,tag);}
    else{
        for(int received=0;received<mpi_solids.number_of_ranks-1;received++){
            MPI::Status probe_status;mpi_solids.comm->Probe(MPI::ANY_SOURCE,tag,probe_status);
            ARRAY<char> buffer(probe_status.Get_count(MPI::PACKED));int position=0;
            mpi_solids.comm->Recv(buffer.Get_Array_Pointer(),buffer.m,MPI::PACKED,probe_status.Get_source(),probe_status.Get_tag());
            ARRAY<T_PAIR1> temp_point_triangle_pairs;ARRAY<T_PAIR2> temp_edge_edge_pairs;
            MPI_UTILITIES::Unpack(temp_point_triangle_pairs,temp_edge_edge_pairs,buffer,position,*mpi_solids.comm); // TODO: Make append unpack
            point_triangle_pairs.Append_Elements(temp_point_triangle_pairs);edge_edge_pairs.Append_Elements(temp_edge_edge_pairs);}}
}
};
template<class TV> template<class T_PAIR1,class T_PAIR2> void MPI_SOLIDS<TV>::
Gather_Interaction_Pairs(ARRAY<T_PAIR1>& point_triangle_pairs,ARRAY<T_PAIR2>& edge_edge_pairs) const
{
    Gather_Interaction_Pairs_Helper(*this,point_triangle_pairs,edge_edge_pairs);
}
//#####################################################################
// Function Broadcast_Collision_Modified_Data
//#####################################################################
template<class TV> void MPI_SOLIDS<TV>::
Broadcast_Collision_Modified_Data(ARRAY_VIEW<bool> modified,ARRAY_VIEW<bool> recently_modified,ARRAY_VIEW<TV> X,ARRAY_VIEW<TV> V) const
{
    ARRAY<int> recently_modified_indices;
    if(rank!=0){
        int count;comm->Bcast(&count,1,MPI_UTILITIES::Datatype<int>(),0);
        ARRAY<char> buffer(count,false);int position=0;
        comm->Bcast(buffer.Get_Array_Pointer(),buffer.m,MPI::PACKED,0);
        MPI_UTILITIES::Unpack(recently_modified_indices,buffer,position,*comm);
        INDIRECT_ARRAY<ARRAY_VIEW<bool> > modified_subset(modified,recently_modified_indices);INDIRECT_ARRAY<ARRAY_VIEW<bool> > recently_modified_subset(recently_modified,recently_modified_indices);
        ARRAYS_COMPUTATIONS::Fill(modified_subset,true);ARRAYS_COMPUTATIONS::Fill(recently_modified_subset,true);
        INDIRECT_ARRAY<ARRAY_VIEW<TV> > X_modified(X,recently_modified_indices),V_modified(V,recently_modified_indices);
        MPI_UTILITIES::Unpack(X_modified,V_modified,buffer,position,*comm);}
    else{
        for(int i=1;i<=recently_modified.Size();i++) if(recently_modified(i)) recently_modified_indices.Append(i);
        // TODO: remove +1 when LAM fixes assert
        ARRAY<char> buffer(1+MPI_UTILITIES::Pack_Size(recently_modified_indices,X.Subset(recently_modified_indices),V.Subset(recently_modified_indices),*comm),false);int position=0;
        MPI_UTILITIES::Pack(recently_modified_indices,X.Subset(recently_modified_indices),V.Subset(recently_modified_indices),buffer,position,*comm);
        comm->Bcast(&position,1,MPI_UTILITIES::Datatype<int>(),0);
        comm->Bcast(buffer.Get_Array_Pointer(),position,MPI::PACKED,0);}
}
//#####################################################################
// Function Gather_Collision_Modified_Data
//#####################################################################
template<class TV> void MPI_SOLIDS<TV>::
Gather_Collision_Modified_Data(ARRAY_VIEW<bool> modified,ARRAY_VIEW<bool> recently_modified,ARRAY_VIEW<TV> X,ARRAY_VIEW<TV> V) const
{
    ARRAY<int> recently_modified_indices;int tag=Get_Unique_Tag();
    if(rank==0){
        for(int received=0;received<number_of_ranks-1;received++){
            MPI::Status probe_status;comm->Probe(MPI::ANY_SOURCE,tag,probe_status);
            ARRAY<char> buffer(probe_status.Get_count(MPI::PACKED));int position=0;
            comm->Recv(buffer.Get_Array_Pointer(),buffer.m,MPI::PACKED,probe_status.Get_source(),probe_status.Get_tag());
            MPI_UTILITIES::Unpack(recently_modified_indices,buffer,position,*comm);
            INDIRECT_ARRAY<ARRAY_VIEW<bool> > modified_subset(modified,recently_modified_indices);INDIRECT_ARRAY<ARRAY_VIEW<bool> > recently_modified_subset(recently_modified,recently_modified_indices);
            ARRAYS_COMPUTATIONS::Fill(modified_subset,true);ARRAYS_COMPUTATIONS::Fill(recently_modified_subset,true);
            INDIRECT_ARRAY<ARRAY_VIEW<TV> > X_modified(X,recently_modified_indices),V_modified(V,recently_modified_indices);
            MPI_UTILITIES::Unpack(X_modified,V_modified,buffer,position,*comm);}}
    else{
        for(int i=1;i<=recently_modified.Size();i++) if(recently_modified(i) && partition_id_from_particle_index(i)==Partition()) recently_modified_indices.Append(i);
        // TODO: remove +1 when LAM fixes assert
        ARRAY<char> buffer(1+MPI_UTILITIES::Pack_Size(recently_modified_indices,X.Subset(recently_modified_indices),V.Subset(recently_modified_indices),*comm),false);int position=0;
        MPI_UTILITIES::Pack(recently_modified_indices,X.Subset(recently_modified_indices),V.Subset(recently_modified_indices),buffer,position,*comm);
        comm->Send(buffer.Get_Array_Pointer(),buffer.m,MPI::PACKED,0,tag);}
}
//#####################################################################
// Function All_Gather_Intersecting_Pairs
//#####################################################################
template<class TV> template<int d1,int d2> void MPI_SOLIDS<TV>::
All_Gather_Intersecting_Pairs(HASHTABLE<VECTOR<int,d1> >& intersecting_point_face_pairs,HASHTABLE<VECTOR<int,d2> >& intersecting_edge_edge_pairs)
{
    // put local pairs into arrays
    ARRAY<VECTOR<int,d1> > local_point_face_pairs(intersecting_point_face_pairs.Size(),false);
    {HASHTABLE_ITERATOR<VECTOR<int,d1> > iterator(intersecting_point_face_pairs);
    for(int i=1;iterator.Valid();iterator.Next(),i++) local_point_face_pairs(i)=iterator.Key();}
    ARRAY<VECTOR<int,d2> > local_edge_edge_pairs(intersecting_edge_edge_pairs.Size(),false);
    {HASHTABLE_ITERATOR<VECTOR<int,d2> > iterator(intersecting_edge_edge_pairs);
    for(int i=1;iterator.Valid();iterator.Next(),i++) local_edge_edge_pairs(i)=iterator.Key();}

    // gather counts
    VECTOR<int,2> local_counts(intersecting_point_face_pairs.Size(),intersecting_edge_edge_pairs.Size());
    ARRAY<VECTOR<int,2> > global_counts(number_of_ranks,false);
    comm->Allgather(&local_counts[1],2,MPI::INT,global_counts.Get_Array_Pointer(),2,MPI::INT);

    // gather point face pairs
    {int total_point_face_count=0;
    ARRAY<int> point_face_counts(number_of_ranks,false),point_face_displacements(number_of_ranks,false);
    for(int k=1;k<=number_of_ranks;k++){
        point_face_displacements(k)=total_point_face_count;
        point_face_counts(k)=global_counts(k)[1];
        total_point_face_count+=global_counts(k)[1];}
    if(total_point_face_count){
        MPI::Datatype point_face_type=MPI_UTILITIES::Datatype<VECTOR<int,d1> >();
        ARRAY<VECTOR<int,d1> > global_point_face_pairs(total_point_face_count,false);
        comm->Allgatherv(local_point_face_pairs.Get_Array_Pointer(),local_point_face_pairs.Size(),point_face_type,
            global_point_face_pairs.Get_Array_Pointer(),point_face_counts.Get_Array_Pointer(),point_face_displacements.Get_Array_Pointer(),point_face_type);
        // fill hash table
        for(int i=1;i<=global_point_face_pairs.Size();i++) intersecting_point_face_pairs.Set(global_point_face_pairs(i));}}

    // gather edge edge pairs
    {int total_edge_edge_count=0;
    ARRAY<int> edge_edge_counts(number_of_ranks,false),edge_edge_displacements(number_of_ranks,false);
    for(int k=1;k<=number_of_ranks;k++){
        edge_edge_displacements(k)=total_edge_edge_count;
        edge_edge_counts(k)=global_counts(k)[2];
        total_edge_edge_count+=global_counts(k)[2];}
    if(total_edge_edge_count){
        MPI::Datatype edge_edge_type=MPI_UTILITIES::Datatype<VECTOR<int,d2> >();
        ARRAY<VECTOR<int,d2> > global_edge_edge_pairs(total_edge_edge_count,false);
        comm->Allgatherv(local_edge_edge_pairs.Get_Array_Pointer(),local_edge_edge_pairs.Size(),edge_edge_type,
            global_edge_edge_pairs.Get_Array_Pointer(),edge_edge_counts.Get_Array_Pointer(),edge_edge_displacements.Get_Array_Pointer(),edge_edge_type);
        // fill hash table
        for(int i=1;i<=global_edge_edge_pairs.Size();i++) intersecting_edge_edge_pairs.Set(global_edge_edge_pairs(i));}}
}
//#####################################################################
// Function All_Scatter_Adhesion_Pairs
//#####################################################################
template<class TV> void MPI_SOLIDS<TV>::
All_Scatter_Adhesion_Pairs(ARRAY<ARRAY<PAIR<VECTOR<int,2>,VECTOR<typename TV::SCALAR,2> > >,PARTITION_ID>& pairs_to_scatter,
        ARRAY<ARRAY<PAIR<VECTOR<int,2>,VECTOR<typename TV::SCALAR,2> > >,PARTITION_ID>& pairs_received)
{
    int tag=Get_Unique_Tag();
    ARRAY<MPI::Request> requests;ARRAY<ARRAY<char>,PARTITION_ID> send_buffers(Number_Of_Partitions());
    for(PARTITION_ID partition(1);partition<=pairs_to_scatter.Size();partition++){
        if(partition!=Partition()){
            send_buffers(partition).Resize(1+MPI_UTILITIES::Pack_Size(pairs_to_scatter(partition),*comm));int position=0;
            MPI_UTILITIES::Pack(pairs_to_scatter(partition),send_buffers(partition),position,*comm);
            requests.Append(comm->Isend(&send_buffers(partition)(1),position,MPI::PACKED,Partition_To_Rank(partition),tag));}}
    pairs_received.Resize(Number_Of_Partitions());
    for(PARTITION_ID dummy_partition(2);dummy_partition<=pairs_received.Size();dummy_partition++){ // n-1 dummy indices
        MPI::Status probe_status;comm->Probe(MPI::ANY_SOURCE,tag,probe_status);
        ARRAY<char> buffer(probe_status.Get_count(MPI::PACKED));int position=0;
        comm->Recv(buffer.Get_Array_Pointer(),buffer.m,MPI::PACKED,probe_status.Get_source(),probe_status.Get_tag());
        MPI_UTILITIES::Unpack(pairs_received(Rank_To_Partition(probe_status.Get_source())),buffer,position,*comm);}
    MPI_UTILITIES::Wait_All(requests);
}
//#####################################################################
// Function Package
//#####################################################################
template<class TV> template<class T_DATA> MPI_PACKAGE* MPI_SOLIDS<TV>::
Package(ARRAY_VIEW<T_DATA> array)
{
    return new MPI_PACKAGE(array);
}
//#####################################################################
// Function Gather_Helper
//#####################################################################
template<class TV> void MPI_SOLIDS<TV>::
Gather_Helper(const ARRAY<MPI_PACKAGE*>& packages) const
{
    int tag=Get_Unique_Tag();
    ARRAY<MPI_PACKAGE> union_packages;ARRAY<MPI::Request> requests;
    if(rank!=0){
        ARRAY<MPI_PACKAGE> packages_to_union;
        for(int i=1;i<=packages.m;i++) packages_to_union.Append(packages(i)->Subset(particles_of_partition(Partition())));
        union_packages.Append(MPI_PACKAGE::Union(packages_to_union));
        requests.Append(union_packages.Last().Isend(*comm,0,tag));
        MPI_PACKAGE::Free_All(packages_to_union);}
    else for(int other_rank=1;other_rank<number_of_ranks;other_rank++){
        ARRAY<MPI_PACKAGE> packages_to_union;
        for(int i=1;i<=packages.m;i++) packages_to_union.Append(packages(i)->Subset(particles_of_partition(Rank_To_Partition(other_rank))));
        union_packages.Append(MPI_PACKAGE::Union(packages_to_union));
        requests.Append(union_packages.Last().Irecv(*comm,other_rank,tag));
        MPI_PACKAGE::Free_All(packages_to_union);}
    MPI_UTILITIES::Wait_All(requests);MPI_PACKAGE::Free_All(union_packages);
    for(int i=1;i<=packages.m;i++){packages(i)->Free();delete packages(i);}
}
//#####################################################################
// Function Scatter_Helper
//#####################################################################
template<class TV> void MPI_SOLIDS<TV>::
Scatter_Helper(const ARRAY<MPI_PACKAGE*>& packages) const
{
    int tag=Get_Unique_Tag();
    ARRAY<MPI_PACKAGE> union_packages;ARRAY<MPI::Request> requests;
    if(rank!=0){
        ARRAY<MPI_PACKAGE> packages_to_union;
        for(int i=1;i<=packages.m;i++) packages_to_union.Append(packages(i)->Subset(particles_of_partition(Partition())));
        union_packages.Append(MPI_PACKAGE::Union(packages_to_union));
        requests.Append(union_packages.Last().Irecv(*comm,0,tag));
        MPI_PACKAGE::Free_All(packages_to_union);}
    else for(int other_rank=1;other_rank<number_of_ranks;other_rank++){
        ARRAY<MPI_PACKAGE> packages_to_union;
        for(int i=1;i<=packages.m;i++) packages_to_union.Append(packages(i)->Subset(particles_of_partition(Rank_To_Partition(other_rank))));
        union_packages.Append(MPI_PACKAGE::Union(packages_to_union));
        requests.Append(union_packages.Last().Isend(*comm,other_rank,tag));
        MPI_PACKAGE::Free_All(packages_to_union);}
    MPI_UTILITIES::Wait_All(requests);MPI_PACKAGE::Free_All(union_packages);
    for(int i=1;i<=packages.m;i++){packages(i)->Free();delete packages(i);}
}
//#####################################################################
// Function Broadcast_Helper
//#####################################################################
template<class TV> void MPI_SOLIDS<TV>::
Broadcast_Helper(const ARRAY<MPI_PACKAGE*>& packages) const
{
    ARRAY<MPI_PACKAGE> union_packages;ARRAY<MPI_PACKAGE> packages_to_union;
    for(int i=1;i<=packages.m;i++) packages_to_union.Append(*packages(i));
    union_packages.Append(MPI_PACKAGE::Union(packages_to_union));
    union_packages(1).Broadcast(*comm);
    MPI_PACKAGE::Free_All(union_packages);
    for(int i=1;i<=packages.m;i++){packages(i)->Free();delete packages(i);}
}
//#####################################################################
// Function Exchange_Helper
//#####################################################################
template<class TV> void MPI_SOLIDS<TV>::
Exchange_Helper(const ARRAY<MPI_PACKAGE*>& packages,MPI_PARTITION& mpi_partition) const
{
    int tag=1;
    ARRAY<MPI_PACKAGE> union_packages;ARRAY<MPI::Request> requests;
    // send
    for(int i=1;i<=mpi_partition.neighbor_partitions.m;i++){PARTITION_ID neighbor_partition=mpi_partition.neighbor_partitions(i);int neighbor_rank=Partition_To_Rank(neighbor_partition);
        ARRAY<MPI_PACKAGE> packages_to_union;
        for(int j=1;j<=packages.m;j++) packages_to_union.Append(packages(j)->Subset(mpi_partition.boundary_dynamic_particles(neighbor_partition)));
        union_packages.Append(MPI_PACKAGE::Union(packages_to_union));
        requests.Append(union_packages.Last().Isend(*comm,neighbor_rank,tag));
        MPI_PACKAGE::Free_All(packages_to_union);}
    // receive
    for(int i=1;i<=mpi_partition.neighbor_partitions.m;i++){PARTITION_ID neighbor_partition=mpi_partition.neighbor_partitions(i);int neighbor_rank=Partition_To_Rank(neighbor_partition);
        ARRAY<MPI_PACKAGE> packages_to_union;
        for(int j=1;j<=packages.m;j++) packages_to_union.Append(packages(j)->Subset(mpi_partition.ghost_dynamic_particles(neighbor_partition)));
        union_packages.Append(MPI_PACKAGE::Union(packages_to_union));
        requests.Append(union_packages.Last().Irecv(*comm,neighbor_rank,tag));
        MPI_PACKAGE::Free_All(packages_to_union);}
    // finish
    MPI_UTILITIES::Wait_All(requests);MPI_PACKAGE::Free_All(union_packages);
    for(int i=1;i<=packages.m;i++){packages(i)->Free();delete packages(i);}
}
//#####################################################################
// Function Reduce_Global
//#####################################################################
template<class TV> template<class T_DATA> T_DATA MPI_SOLIDS<TV>::
Reduce_Add_Global(const T_DATA local_value) const
{
    T_DATA global_value;
    MPI_UTILITIES::Reduce(local_value,global_value,MPI::SUM,*comm);
    return global_value;
}
//#####################################################################
// Function Reduce_Global
//#####################################################################
template<class TV> template<class T_DATA> T_DATA MPI_SOLIDS<TV>::
Reduce_Max_Global(const T_DATA local_value) const
{
    T_DATA global_value;
    MPI_UTILITIES::Reduce(local_value,global_value,MPI::MAX,*comm);
    return global_value;
}
//#####################################################################
// Function Reduce_Global
//#####################################################################
template<class TV> template<class T_DATA> T_DATA MPI_SOLIDS<TV>::
Reduce_Min_Global(const T_DATA local_value) const
{
    T_DATA global_value;
    MPI_UTILITIES::Reduce(local_value,global_value,MPI::MIN,*comm);
    return global_value;
}
//#####################################################################
// Function Reduce_Add
//#####################################################################
template<class TV> template<class T_DATA> T_DATA MPI_SOLIDS<TV>::
Reduce_Add(const T_DATA local_value) const
{
    if(number_of_ranks==1) return local_value;
    T_DATA global_value;
    MPI_UTILITIES::Reduce(local_value,global_value,MPI::SUM,*comm);
    return global_value;
}
//#####################################################################
// Function Reduce_Max
//#####################################################################
template<class TV> template<class T_DATA> T_DATA MPI_SOLIDS<TV>::
Reduce_Max(const T_DATA local_value) const
{
    if(number_of_ranks==1) return local_value;
    T_DATA global_value;
    MPI_UTILITIES::Reduce(local_value,global_value,MPI::MAX,*comm);
    return global_value;
}
//#####################################################################
// Function Reduce_Min
//#####################################################################
template<class TV> template<class T_DATA> T_DATA MPI_SOLIDS<TV>::
Reduce_Min(const T_DATA local_value) const
{
    if(number_of_ranks==1) return local_value;
    T_DATA global_value;
    MPI_UTILITIES::Reduce(local_value,global_value,MPI::MIN,*comm);
    return global_value;
}
//#####################################################################
// Function Debug_Reduce_Array_Sum
//#####################################################################
template<class TV> template<class T_DATA> T_DATA MPI_SOLIDS<TV>::
Debug_Reduce_Array_Sum(ARRAY_VIEW<const T_DATA> local_values) const // sends all values to one processor and adds them up in order (very slow)
{
    const int root=0;
    if(number_of_ranks==1) return ARRAYS_COMPUTATIONS::Sum(local_values);
    MPI::Datatype type=MPI_UTILITIES::Datatype<T_DATA>();
    int local_count=local_values.Size();
    if(comm->Get_rank()!=root){ // slave
        comm->Gather(&local_count,1,MPI::INT,0,1,MPI::INT,root);
        comm->Gatherv(local_values.Get_Array_Pointer(),local_count,type,0,0,0,type,root);
        T_DATA sum;
        comm->Bcast(&sum,1,type,root);
        return sum;}
    else{ // master
        ARRAY<int> counts(comm->Get_size(),false);
        comm->Gather(&local_count,1,MPI::INT,counts.Get_Array_Pointer(),1,MPI::INT,root);
        int total_count=0;
        ARRAY<int> displacements(counts.m,false);
        for(int k=1;k<=counts.m;k++){
            displacements(k)=total_count;
            total_count+=counts(k);}
        ARRAY<T_DATA> global_values(total_count,false);
        comm->Gatherv(local_values.Get_Array_Pointer(),local_count,type,
            global_values.Get_Array_Pointer(),counts.Get_Array_Pointer(),displacements.Get_Array_Pointer(),type,root);
        {std::stringstream ss;ss<<"global values hash = "<<Hash(global_values)<<std::endl;LOG::filecout(ss.str());}
        T_DATA sum=ARRAYS_COMPUTATIONS::Sum(global_values);
        comm->Bcast(&sum,1,type,root);
        return sum;}
}
//#####################################################################
// Function Simple_Partition
//#####################################################################
template<class TV> void MPI_SOLIDS<TV>::
Simple_Partition(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection_input,RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection_input,ARRAY_VIEW<const TV> X,const VECTOR<int,TV::m>& counts)
{
    PHYSBAM_ASSERT(counts.Product()==number_of_ranks);
    particles_of_partition.Resize(Number_Of_Partitions());
    RANGE<TV> box=RANGE<TV>::Bounding_Box(X);box.Change_Size((typename TV::SCALAR)1e-6);GRID<TV> grid(counts+1,box);
    partition_id_from_particle_index.Remove_All();
    partition_id_from_particle_index.Resize(deformable_body_collection_input.particles.array_collection->Size()+rigid_geometry_collection_input.particles.array_collection->Size());
    for(int p=1;p<=X.Size();p++){
        VECTOR<int,TV::m> cell=grid.Clamp_To_Cell(X(p));
        int cell_number=cell[1]-1;for(int i=2;i<=TV::m;i++) cell_number=cell_number*counts[i]+cell[i]-1;
        partition_id_from_particle_index(p)=Rank_To_Partition(cell_number);
        particles_of_partition(Rank_To_Partition(cell_number)).Append(p);}
}
//#####################################################################
// Function Simple_Partition
//#####################################################################
template<class T_ARRAY,class T> void KD_Tree_Partition_Helper(const T_ARRAY& local_p_to_global_p,KD_TREE_NODE<T>& node,const int depth,const int target_depth,PARTITION_ID& partition,ARRAY<ARRAY<int>,PARTITION_ID>& particle_of_partition)
{
    if(node.split_axis==0) particle_of_partition(partition).Append(local_p_to_global_p(typename T_ARRAY::INDEX(node.node_index)));
    else{
        if(depth==target_depth) partition++;
        KD_Tree_Partition_Helper(local_p_to_global_p,*node.left,depth+1,target_depth,partition,particle_of_partition);
        KD_Tree_Partition_Helper(local_p_to_global_p,*node.right,depth+1,target_depth,partition,particle_of_partition);}
}
template<class TV> void MPI_SOLIDS<TV>::
KD_Tree_Partition(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection_input,RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection_input,ARRAY_VIEW<const TV> X)
{
    int height=integer_log(number_of_ranks);
    if(1<<height!=number_of_ranks)
        throw std::runtime_error("KD_Tree_Partition requires number of processors to be a power of two");
    particles_of_partition.Resize(Number_Of_Partitions());
    KD_TREE<TV> kd_tree(false);kd_tree.Create_KD_Tree(X);PARTITION_ID partition;
    KD_Tree_Partition_Helper(IDENTITY_ARRAY<int>(X.Size()),*kd_tree.root_node,0,height,partition,particles_of_partition);

    std::stringstream ss;
    ss<<"Partition Statistics: ";
    for(PARTITION_ID i(1);i<=particles_of_partition.Size();i++) ss<<particles_of_partition(i).Size()<<", ";
    ss<<std::endl;
    LOG::filecout(ss.str());
    partition_id_from_particle_index.Remove_All();
    partition_id_from_particle_index.Resize(deformable_body_collection_input.particles.array_collection->Size()+rigid_geometry_collection_input.particles.array_collection->Size());
    for(PARTITION_ID i(1);i<=particles_of_partition.Size();i++){
        INDIRECT_ARRAY<ARRAY<PARTITION_ID> > partition_subset(partition_id_from_particle_index,particles_of_partition(i));
        if(partition_id_from_particle_index.Size()) ARRAYS_COMPUTATIONS::Fill(partition_subset,i);}
}
//#####################################################################
// Function KD_Tree_Partition_Subset
//#####################################################################
template<class TV> template<class ID> void MPI_SOLIDS<TV>::
KD_Tree_Partition_Subset(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection_input,RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection_input,const ARRAY<int,ID>& nodes,ARRAY_VIEW<const TV> X)
{
    int height=integer_log(number_of_ranks);
    if(1<<height!=number_of_ranks)
        throw std::runtime_error("KD_Tree_Partition requires number of processors to be a power of two");
    particles_of_partition.Resize(Number_Of_Partitions());
    ARRAY<TV> X_nodes(Value(nodes.Size()),false);
    for(typename ARRAY<int,ID>::INDEX i(1);i<=nodes.Size();i++) X_nodes(Value(i))=X(nodes(i));
    KD_TREE<TV> kd_tree(false);kd_tree.Create_KD_Tree(X_nodes);PARTITION_ID partition;
    KD_Tree_Partition_Helper(nodes,*kd_tree.root_node,0,height,partition,particles_of_partition);
    std::stringstream ss;
    ss<<"Partition Statistics: ";
    for(PARTITION_ID i(1);i<=particles_of_partition.Size();i++) ss<<particles_of_partition(i).Size()<<", ";
    ss<<std::endl;LOG::filecout(ss.str());
    partition_id_from_particle_index.Remove_All();
    partition_id_from_particle_index.Resize(deformable_body_collection_input.particles.array_collection->Size()+rigid_geometry_collection_input.particles.array_collection->Size());
    for(PARTITION_ID i(1);i<=particles_of_partition.Size();i++){
        INDIRECT_ARRAY<ARRAY<PARTITION_ID> > partition_subset(partition_id_from_particle_index,particles_of_partition(i));
        ARRAYS_COMPUTATIONS::Fill(partition_subset,i);}
}
//#####################################################################
// Function Distribute_Repulsion_Pairs
//#####################################################################
namespace{
template<class T,class T_ARRAY_PAIR> void Distribute_Repulsion_Pairs_Helper(const MPI_SOLIDS<VECTOR<T,1> >& mpi_solids,const T_ARRAY_PAIR& pairs,ARRAY<ARRAY<int>,PARTITION_ID>& send_particles,ARRAY<ARRAY<int>,PARTITION_ID>& receive_particles,
                                                                     T_ARRAY_PAIR& my_boundary_pairs,T_ARRAY_PAIR& my_internal_pairs)
{PHYSBAM_NOT_IMPLEMENTED();}
template<class TV,class T_ARRAY_PAIR> void Distribute_Repulsion_Pairs_Helper(const MPI_SOLIDS<TV>& mpi_solids,const T_ARRAY_PAIR& pairs,ARRAY<ARRAY<int>,PARTITION_ID>& send_particles,ARRAY<ARRAY<int>,PARTITION_ID>& receive_particles,
                                                                     T_ARRAY_PAIR& my_boundary_pairs,T_ARRAY_PAIR& my_internal_pairs)
{
    typedef typename T_ARRAY_PAIR::ELEMENT T_PAIR;
    typedef typename T_PAIR::SCALAR T_DATA;
    typedef typename IF<IS_SAME<T_PAIR,POINT_FACE_REPULSION_PAIR<VECTOR<T_DATA,T_PAIR::d> > >::value,VECTOR<PARTITION_ID,T_PAIR::d+1>,VECTOR<PARTITION_ID,T_PAIR::d*2-2> >::TYPE PROCESSOR_VECTOR;
    STATIC_ASSERT(sizeof(int)==4);PHYSBAM_ASSERT(mpi_solids.particles_of_partition.Size()<=PARTITION_ID(32));

    LOG::SCOPE scope("Repulsion Pair Distribution","Repulsion Pair Distribution");
    int tag=mpi_solids.Get_Unique_Tag();
    if(mpi_solids.rank==0){
        // Color nodes
        UNION_FIND<> union_find(mpi_solids.deformable_body_collection->particles.array_collection->Size());
        ARRAY<unsigned int> processor_masks(pairs.m,false);
        ARRAY<ARRAY<int>,PARTITION_ID> internal_pairs(mpi_solids.particles_of_partition.Size(),false);
        int boundary_count=0;
        for(int i=1;i<=pairs.m;i++){const T_PAIR& pair=pairs(i);
            PROCESSOR_VECTOR processors(mpi_solids.partition_id_from_particle_index.Subset(pair.nodes));
            processor_masks(i)=0;for(int j=1;j<=processors.Size();j++) processor_masks(i)|=1<<(mpi_solids.Partition_To_Rank(processors[j]));
            if(!power_of_two(processor_masks(i))){union_find.Union(pair.nodes);boundary_count++;}else internal_pairs(processors[1]).Append(i);}
        // Build connected components
        ARRAY<ARRAY<int> > components;ARRAY<unsigned int> component_processors;
        ARRAY<int> particle_to_component(mpi_solids.partition_id_from_particle_index.m);
        for(int i=1;i<=pairs.m;i++){const T_PAIR& pair=pairs(i);unsigned int processor_mask=processor_masks(i);
            if(!power_of_two(processor_mask)){
                int root=union_find.Find(pair.nodes[1]);
                if(!particle_to_component(root)){particle_to_component(root)=components.Append(ARRAY<int>());component_processors.Append(0);}
                int component_id=particle_to_component(root);
                INDIRECT_ARRAY<ARRAY<int>,VECTOR<int,T_PAIR::count>&> particle_subset(particle_to_component,pair.nodes);
                ARRAYS_COMPUTATIONS::Fill(particle_subset,component_id);
                components(component_id).Append(i);
                component_processors(component_id)|=processor_mask;}}
        // choose processors to dump on and make list of boundary pairs processed on each processor
        ARRAY<PARTITION_ID> component_processor(components.m,false);ARRAY<ARRAY<int>,PARTITION_ID> processor_pair_indices(mpi_solids.particles_of_partition.Size());
        for(int i=1;i<=components.m;i++){
            PARTITION_ID processor=mpi_solids.Rank_To_Partition(integer_log_exact(rightmost_bit(component_processors(i))));component_processor(i)=processor;
            processor_pair_indices(processor).Append_Elements(components(i));}
        // construct particle list to send
        ARRAY<ARRAY<ARRAY<int>,PARTITION_ID>,PARTITION_ID> send_matrix(mpi_solids.particles_of_partition.Size()),receive_matrix(mpi_solids.particles_of_partition.Size());
        for(PARTITION_ID i(1);i<=mpi_solids.particles_of_partition.Size();i++){
            send_matrix(i).Resize(mpi_solids.particles_of_partition.Size());receive_matrix(i).Resize(mpi_solids.particles_of_partition.Size());}
        for(int p=1;p<=particle_to_component.m;p++){int component=particle_to_component(p);
            if(component){
                PARTITION_ID source=mpi_solids.partition_id_from_particle_index(p),destination=component_processor(component);
                if(source!=destination){send_matrix(source)(destination).Append(p);receive_matrix(destination)(source).Append(p);}}}
        // send pairs
        ARRAY<MPI::Request> requests;ARRAY<ARRAY<char>,PARTITION_ID> buffers(mpi_solids.particles_of_partition.Size());
        for(PARTITION_ID processor(2);processor<=mpi_solids.particles_of_partition.Size();processor++){
            INDIRECT_ARRAY<const T_ARRAY_PAIR> boundary(pairs,processor_pair_indices(processor)),internal(pairs,internal_pairs(processor));
            buffers(processor).Resize(1+MPI_UTILITIES::Pack_Size(send_matrix(processor),receive_matrix(processor),boundary,internal,*mpi_solids.comm)); // TODO: remove +1 when LAM fixes assert
            int position=0;
            MPI_UTILITIES::Pack(send_matrix(processor),receive_matrix(processor),boundary,internal,buffers(processor),position,*mpi_solids.comm);
            requests.Append(mpi_solids.comm->Isend(&buffers(processor)(1),position,MPI::PACKED,mpi_solids.Partition_To_Rank(processor),tag));}
        send_particles=send_matrix(PARTITION_ID(1));receive_particles=receive_matrix(PARTITION_ID(1));
        my_boundary_pairs=pairs.Subset(processor_pair_indices(PARTITION_ID(1)));my_internal_pairs=pairs.Subset(internal_pairs(PARTITION_ID(1)));
        MPI_UTILITIES::Wait_All(requests);}
    else{
        MPI::Status probe_status;
        mpi_solids.comm->Probe(0,tag,probe_status);
        ARRAY<char> buffer(probe_status.Get_count(MPI::PACKED));int position=0;
        mpi_solids.comm->Recv(buffer.Get_Array_Pointer(),buffer.m,MPI::PACKED,probe_status.Get_source(),probe_status.Get_tag());
        MPI_UTILITIES::Unpack(send_particles,receive_particles,my_boundary_pairs,my_internal_pairs,buffer,position,*mpi_solids.comm);}
}
};

template<class TV> template<class T_ARRAY_PAIR> void MPI_SOLIDS<TV>::
Distribute_Repulsion_Pairs(const T_ARRAY_PAIR& pairs,ARRAY<ARRAY<int>,PARTITION_ID>& send_particles,ARRAY<ARRAY<int>,PARTITION_ID>& receive_particles,T_ARRAY_PAIR& my_boundary_pairs,
    T_ARRAY_PAIR& my_internal_pairs) const
{
    Distribute_Repulsion_Pairs_Helper(*this,pairs,send_particles,receive_particles,my_boundary_pairs,my_internal_pairs);
}
//#####################################################################
// Function Gather_Repulsion_Inputs
//#####################################################################
template<class TV> void MPI_SOLIDS<TV>::
Gather_Repulsion_Inputs(ARRAY_VIEW<TV> X_self_collision_free,ARRAY_VIEW<TV> X,ARRAY_VIEW<TV> V,ARRAY<ARRAY<int>,PARTITION_ID>& send_particles,
    ARRAY<ARRAY<int>,PARTITION_ID>& receive_particles) const
{
    int tag=Get_Unique_Tag();
    ARRAY<MPI_PACKAGE> packages;ARRAY<MPI::Request> requests;
    // send
    for(PARTITION_ID other_proc(1);other_proc<=send_particles.Size();other_proc++) if(other_proc!=Partition() && send_particles(other_proc).Size()>0){
        ARRAY<MPI_PACKAGE> packages_to_union;
        //packages_to_union.Append(X_self_collision_free.Subset(send_particles(other_proc)));
        packages_to_union.Append(MPI_PACKAGE(X.Subset(send_particles(other_proc))));
        packages_to_union.Append(MPI_PACKAGE(V.Subset(send_particles(other_proc))));
        MPI_PACKAGE package=MPI_PACKAGE::Union(packages_to_union);
        packages.Append(package);
        requests.Append(package.Isend(*comm,Partition_To_Rank(other_proc),tag));
        MPI_PACKAGE::Free_All(packages_to_union);}
    // receive
    for(PARTITION_ID other_proc(1);other_proc<=receive_particles.Size();other_proc++) if(other_proc!=Partition() && receive_particles(other_proc).Size()>0){
        ARRAY<MPI_PACKAGE> packages_to_union;
        //packages_to_union.Append(X_self_collision_free.Subset(receive_particles(other_proc)));
        packages_to_union.Append(MPI_PACKAGE(X.Subset(receive_particles(other_proc))));
        packages_to_union.Append(MPI_PACKAGE(V.Subset(receive_particles(other_proc))));
        MPI_PACKAGE package=MPI_PACKAGE::Union(packages_to_union);
        packages.Append(package);
        requests.Append(package.Irecv(*comm,Partition_To_Rank(other_proc),tag));
        MPI_PACKAGE::Free_All(packages_to_union);}
    MPI_UTILITIES::Wait_All(requests);MPI_PACKAGE::Free_All(packages);
}
//#####################################################################
// Function Scatter_Repulsion_Results
//#####################################################################
template<class TV> void MPI_SOLIDS<TV>::
Scatter_Repulsion_Outputs(ARRAY_VIEW<TV> X_self_collision_free,ARRAY_VIEW<TV> X,ARRAY_VIEW<TV> V,ARRAY<ARRAY<int>,PARTITION_ID>& send_particles,
    ARRAY<ARRAY<int>,PARTITION_ID>& receive_particles) const
{
    Gather_Repulsion_Inputs(X_self_collision_free,X,V,receive_particles,send_particles);
}
//#####################################################################
// Function Print_Diagnostics
//#####################################################################
template<class TV> void MPI_SOLIDS<TV>::
Print_Diagnostics() const
{
    std::stringstream ss;
    ss<<"======================================================================"<<std::endl;
    ss<<"MPI_SOLIDS Diagnostics - Rank="<<rank<<", number_of_processes="<<Number_Of_Processors()<<std::endl;
    ss<<"======================================================================"<<std::endl;
    ss<<"particles_of_partition.m="<<particles_of_partition.Size()<<std::endl;
    ss<<"mpi_partition_binding.neighbor_partitions : "<<mpi_partition_binding.neighbor_partitions<<std::endl;
    ss<<"mpi_partition_binding.ghost_dynamic_particles : "<<mpi_partition_binding.ghost_dynamic_particles<<std::endl;
    ss<<"mpi_partition_binding.boundary_dynamic_particles : "<<mpi_partition_binding.boundary_dynamic_particles<<std::endl;
    ss<<"mpi_partition_force.neighbor_partitions : "<<mpi_partition_force.neighbor_partitions<<std::endl;
    ss<<"mpi_partition_force.ghost_dynamic_particles : "<<mpi_partition_force.ghost_dynamic_particles<<std::endl;
    ss<<"mpi_partition_force.boundary_dynamic_particles : "<<mpi_partition_force.boundary_dynamic_particles<<std::endl;
    LOG::filecout(ss.str());
}
//#####################################################################
// Function Comm
//#####################################################################
template<class TV> const MPI::Intracomm& MPI_SOLIDS<TV>::
Comm() const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################

#else

//#####################################################################
template<class TV> MPI_SOLIDS<TV>::MPI_SOLIDS(){PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> MPI_SOLIDS<TV>::~MPI_SOLIDS(){}
template<class TV> void MPI_SOLIDS<TV>::Update(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection_input,const SEGMENT_MESH& force_connectivity,const SEGMENT_MESH& binding_connectivity) {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> void MPI_SOLIDS<TV>::Compute_Mpi_Partition(MPI_PARTITION& mpi_partition,const SEGMENT_MESH& connectivity) {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> template<class T_DATA> void MPI_SOLIDS<TV>::Exchange_Boundary_Data(ARRAY_VIEW<T_DATA> array,const MPI_PARTITION& mpi_partitions) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> template<class T_DATA> void MPI_SOLIDS<TV>::Exchange_Force_Boundary_Data_Global(ARRAY_VIEW<T_DATA> array) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> template<class T_DATA> void MPI_SOLIDS<TV>::Exchange_Binding_Boundary_Data_Global(ARRAY_VIEW<T_DATA> array) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> template<class T_DATA> void MPI_SOLIDS<TV>::Gather_Data(ARRAY_VIEW<T_DATA> array) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> template<class T_DATA> void MPI_SOLIDS<TV>::Broadcast(T_DATA* data,const int length) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> template<class T_DATA> T_DATA MPI_SOLIDS<TV>::Reduce_Add_Global(const T_DATA local_value) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> template<class T_DATA> T_DATA MPI_SOLIDS<TV>::Reduce_Max_Global(const T_DATA local_value) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> template<class T_DATA> T_DATA MPI_SOLIDS<TV>::Reduce_Min_Global(const T_DATA local_value) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> template<class T_DATA> T_DATA MPI_SOLIDS<TV>::Reduce_Add(const T_DATA local_value) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> template<class T_DATA> T_DATA MPI_SOLIDS<TV>::Reduce_Max(const T_DATA local_value) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> template<class T_DATA> T_DATA MPI_SOLIDS<TV>::Reduce_Min(const T_DATA local_value) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> template<class T_DATA> T_DATA MPI_SOLIDS<TV>::Debug_Reduce_Array_Sum(ARRAY_VIEW<const T_DATA> local_values) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> void MPI_SOLIDS<TV>::Gather_Repulsion_Inputs(ARRAY_VIEW<TV> X_self_collision_free,ARRAY_VIEW<TV> X,ARRAY_VIEW<TV> V,ARRAY<ARRAY<int>,PARTITION_ID>& send_particles,
    ARRAY<ARRAY<int>,PARTITION_ID>& receive_particles) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> void MPI_SOLIDS<TV>::Scatter_Repulsion_Outputs(ARRAY_VIEW<TV> X_self_collision_free,ARRAY_VIEW<TV> X,ARRAY_VIEW<TV> V,ARRAY<ARRAY<int>,PARTITION_ID>& send_particles,
    ARRAY<ARRAY<int>,PARTITION_ID>& receive_particles) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> void MPI_SOLIDS<TV>::Simple_Partition(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection_input,RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection_input,ARRAY_VIEW<const TV> X,const VECTOR<int,TV::m>& counts){PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> void MPI_SOLIDS<TV>::KD_Tree_Partition(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection_input,RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection_input,ARRAY_VIEW<const TV> X){PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> template<class ID> void MPI_SOLIDS<TV>::KD_Tree_Partition_Subset(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection_input,RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection_input,const ARRAY<int,ID>& nodes,ARRAY_VIEW<const TV> X){PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> template<class T_ARRAY_PAIR> void MPI_SOLIDS<TV>::Distribute_Repulsion_Pairs(const T_ARRAY_PAIR& pairs,ARRAY<ARRAY<int>,PARTITION_ID>& send_particles,
    ARRAY<ARRAY<int>,PARTITION_ID>& receive_particles,T_ARRAY_PAIR& my_boundary_pairs,T_ARRAY_PAIR& my_internal_pairs) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> void MPI_SOLIDS<TV>::Print_Diagnostics() const{PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> template<class T_DATA> MPI_PACKAGE* MPI_SOLIDS<TV>::Package(ARRAY_VIEW<T_DATA> array){PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> void MPI_SOLIDS<TV>::Gather_Helper(const ARRAY<MPI_PACKAGE*>& packages) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> void MPI_SOLIDS<TV>::Scatter_Helper(const ARRAY<MPI_PACKAGE*>& packages) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> void MPI_SOLIDS<TV>::Exchange_Helper(const ARRAY<MPI_PACKAGE*>& packages,MPI_PARTITION& mpi_partition) const 
    {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> const MPI::Intracomm& MPI_SOLIDS<TV>::Comm() const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> void MPI_SOLIDS<TV>::All_Gather_Particles(ARRAY_VIEW<TV> X,ARRAY_VIEW<TV> V) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> template<class T_PAIR1,class T_PAIR2> void MPI_SOLIDS<TV>::Gather_Interaction_Pairs(ARRAY<T_PAIR1>& point_triangle_pairs,ARRAY<T_PAIR2>& edge_edge_pairs) const
    {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> void MPI_SOLIDS<TV>::Broadcast_Collision_Modified_Data(ARRAY_VIEW<bool> modified,ARRAY_VIEW<bool> recently_modified,ARRAY_VIEW<TV> X,ARRAY_VIEW<TV> V) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> void MPI_SOLIDS<TV>::Gather_Collision_Modified_Data(ARRAY_VIEW<bool> modified,ARRAY_VIEW<bool> recently_modified,ARRAY_VIEW<TV> X,ARRAY_VIEW<TV> V) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> void MPI_SOLIDS<TV>::Broadcast_Helper(const ARRAY<MPI_PACKAGE*>& packages) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> template<int d1,int d2> void MPI_SOLIDS<TV>::All_Gather_Intersecting_Pairs(HASHTABLE<VECTOR<int,d1> >& intersecting_point_face_pairs,HASHTABLE<VECTOR<int,d2> >& intersecting_edge_edge_pairs)
    {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> void MPI_SOLIDS<TV>::All_Scatter_Adhesion_Pairs(ARRAY<ARRAY<PAIR<VECTOR<int,2>,VECTOR<typename TV::SCALAR,2> > >,PARTITION_ID>& pairs_to_scatter,
    ARRAY<ARRAY<PAIR<VECTOR<int,2>,VECTOR<typename TV::SCALAR,2> > >,PARTITION_ID>& pairs_received)
    {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> void MPI_SOLIDS<TV>::Barrier() const {PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################

#endif

//#####################################################################
// TODO: Clean up this mess.
#define P(...) __VA_ARGS__
#define INSTANTIATION_HELPER2(V,T,d)

#define INSTANTIATION_HELPER(V,T) \
    INSTANTIATION_HELPER2(P(V),T,1); \
    INSTANTIATION_HELPER2(P(V),T,2); \
    INSTANTIATION_HELPER2(P(V),T,3); \
    template void MPI_SOLIDS<P(V)>::Gather_Data(ARRAY_VIEW<T>) const; \
    template void MPI_SOLIDS<P(V)>::Gather_Data(ARRAY_VIEW<VECTOR<T,2> >) const; \
    template void MPI_SOLIDS<P(V)>::Gather_Data(ARRAY_VIEW<VECTOR<T,3> >) const; \
    template void MPI_SOLIDS<P(V)>::Broadcast(T* data,const int length) const; \
    template MPI_PACKAGE* MPI_SOLIDS<P(V)>::Package(ARRAY_VIEW<T>); \
    template void MPI_SOLIDS<P(V)>::Exchange_Boundary_Data(ARRAY_VIEW<T>,const MPI_PARTITION&) const; \
    template void MPI_SOLIDS<P(V)>::Exchange_Force_Boundary_Data_Global(ARRAY_VIEW<T>) const; \
    template void MPI_SOLIDS<P(V)>::Exchange_Binding_Boundary_Data_Global(ARRAY_VIEW<T>) const;

#define INSTANTIATION_HELPER_ALL2(V,T,d,e) \
    template void MPI_SOLIDS<P(V)>::Gather_Interaction_Pairs(ARRAY<VECTOR<T,d> >& point_triangle_pairs,ARRAY<VECTOR<T,e> >& edge_edge_pairs) const;

#define INSTANTIATION_HELPER_ALL3(V,T,d,e) \
    template void MPI_SOLIDS<P(V)>::Gather_Interaction_Pairs(ARRAY<POINT_FACE_REPULSION_PAIR<VECTOR<T,d> > >& point_triangle_pairs, \
        ARRAY<EDGE_EDGE_REPULSION_PAIR<VECTOR<T,e> > >& edge_edge_pairs) const;

#define INSTANTIATION_HELPER_ALL_NEED1(V,T) \
    template int MPI_SOLIDS<P(V)>::Reduce_Add(const int) const; \
    template float MPI_SOLIDS<P(V)>::Reduce_Add(const float) const; \
    template double MPI_SOLIDS<P(V)>::Reduce_Add(const double) const; \
    template int MPI_SOLIDS<P(V)>::Reduce_Add_Global(const int) const; \
    template float MPI_SOLIDS<P(V)>::Reduce_Add_Global(const float) const; \
    template double MPI_SOLIDS<P(V)>::Reduce_Add_Global(const double) const; \
    template int MPI_SOLIDS<P(V)>::Reduce_Max_Global(const int) const; \
    template float MPI_SOLIDS<P(V)>::Reduce_Max_Global(const float) const; \
    template double MPI_SOLIDS<P(V)>::Reduce_Max_Global(const double) const; \
    template int MPI_SOLIDS<P(V)>::Reduce_Min_Global(const int) const; \
    template float MPI_SOLIDS<P(V)>::Reduce_Min_Global(const float) const; \
    template double MPI_SOLIDS<P(V)>::Reduce_Min_Global(const double) const; \
    template int MPI_SOLIDS<P(V)>::Reduce_Max(const int) const; \
    template float MPI_SOLIDS<P(V)>::Reduce_Max(const float) const; \
    template double MPI_SOLIDS<P(V)>::Reduce_Max(const double) const; \
    template int MPI_SOLIDS<P(V)>::Reduce_Min(const int) const; \
    template float MPI_SOLIDS<P(V)>::Reduce_Min(const float) const; \
    template double MPI_SOLIDS<P(V)>::Reduce_Min(const double) const; \
    template void MPI_SOLIDS<P(V)>::Exchange_Boundary_Data(ARRAY_VIEW<VECTOR<T,1> >,const MPI_PARTITION&) const; \
    template void MPI_SOLIDS<P(V)>::Exchange_Boundary_Data(ARRAY_VIEW<VECTOR<T,2> >,const MPI_PARTITION&) const; \
    template void MPI_SOLIDS<P(V)>::Exchange_Boundary_Data(ARRAY_VIEW<VECTOR<T,3> >,const MPI_PARTITION&) const; \
    template void MPI_SOLIDS<P(V)>::Exchange_Force_Boundary_Data_Global(ARRAY_VIEW<VECTOR<T,1> >) const; \
    template void MPI_SOLIDS<P(V)>::Exchange_Force_Boundary_Data_Global(ARRAY_VIEW<VECTOR<T,2> >) const; \
    template void MPI_SOLIDS<P(V)>::Exchange_Force_Boundary_Data_Global(ARRAY_VIEW<VECTOR<T,3> >) const; \
    template void MPI_SOLIDS<P(V)>::Exchange_Binding_Boundary_Data_Global(ARRAY_VIEW<VECTOR<T,1> >) const; \
    template void MPI_SOLIDS<P(V)>::Exchange_Binding_Boundary_Data_Global(ARRAY_VIEW<VECTOR<T,2> >) const; \
    template void MPI_SOLIDS<P(V)>::Exchange_Binding_Boundary_Data_Global(ARRAY_VIEW<VECTOR<T,3> >) const; \
    template MPI_PACKAGE* MPI_SOLIDS<P(V)>::Package(ARRAY_VIEW<VECTOR<T,1> >); \
    template MPI_PACKAGE* MPI_SOLIDS<P(V)>::Package(ARRAY_VIEW<VECTOR<T,2> >); \
    template MPI_PACKAGE* MPI_SOLIDS<P(V)>::Package(ARRAY_VIEW<VECTOR<T,3> >);

#define INSTANTIATION_HELPER_1D_ONLY(T) \
    template void MPI_SOLIDS<VECTOR<T,1> >::Gather_Interaction_Pairs(ARRAY<POINT_FACE_REPULSION_PAIR<VECTOR<T,1> > >& point_triangle_pairs, \
        ARRAY<EDGE_EDGE_REPULSION_PAIR<VECTOR<T,1> > >& edge_edge_pairs) const; \
    template void MPI_SOLIDS<VECTOR<T,1> >::Distribute_Repulsion_Pairs(const ARRAY<POINT_FACE_REPULSION_PAIR<VECTOR<T,1> > >&,ARRAY<ARRAY<int>,PARTITION_ID>&, \
        ARRAY<ARRAY<int>,PARTITION_ID>&,ARRAY<POINT_FACE_REPULSION_PAIR<VECTOR<T,1> > >&,ARRAY<POINT_FACE_REPULSION_PAIR<VECTOR<T,1> > >&) const; \
    template void MPI_SOLIDS<VECTOR<T,1> >::Distribute_Repulsion_Pairs(const ARRAY<EDGE_EDGE_REPULSION_PAIR<VECTOR<T,1> > >&,ARRAY<ARRAY<int>,PARTITION_ID>&, \
        ARRAY<ARRAY<int>,PARTITION_ID>&,ARRAY<EDGE_EDGE_REPULSION_PAIR<VECTOR<T,1> > >&,ARRAY<EDGE_EDGE_REPULSION_PAIR<VECTOR<T,1> > >&) const;

#define INSTANTIATION_HELPER_ALL(V,T) \
    template double MPI_SOLIDS<P(V)>::Debug_Reduce_Array_Sum(ARRAY_VIEW<const double>) const; \
    template void MPI_SOLIDS<P(V)>::Distribute_Repulsion_Pairs(const ARRAY<POINT_FACE_REPULSION_PAIR<P(V)> >&,ARRAY<ARRAY<int>,PARTITION_ID>&, \
        ARRAY<ARRAY<int>,PARTITION_ID>&,ARRAY<POINT_FACE_REPULSION_PAIR<P(V)> >&,ARRAY<POINT_FACE_REPULSION_PAIR<P(V)> >&) const; \
    template void MPI_SOLIDS<P(V)>::Distribute_Repulsion_Pairs(const ARRAY<EDGE_EDGE_REPULSION_PAIR<P(V)> >&,ARRAY<ARRAY<int>,PARTITION_ID>&, \
        ARRAY<ARRAY<int>,PARTITION_ID>&,ARRAY<EDGE_EDGE_REPULSION_PAIR<P(V)> >&,ARRAY<EDGE_EDGE_REPULSION_PAIR<P(V)> >&) const; \
    template void MPI_SOLIDS<P(V)>::All_Gather_Intersecting_Pairs<P(V::dimension+1),P(V::dimension*2-2)>(HASHTABLE<VECTOR<int,P(V::dimension+1)>,void>&, \
        HASHTABLE<VECTOR<int,P(V::dimension*2-2)>,void>&); \
    template void MPI_SOLIDS<P(V)>::KD_Tree_Partition_Subset(DEFORMABLE_BODY_COLLECTION<P(V)>& deformable_body_collection_input,RIGID_GEOMETRY_COLLECTION<P(V)>& rigid_geometry_collection_input,const ARRAY<int,HAIR_ID>& nodes,ARRAY_VIEW<const P(V)> X); \
    INSTANTIATION_HELPER(P(V),T); \
    INSTANTIATION_HELPER(P(V),int); \
    INSTANTIATION_HELPER_ALL2(P(V),int,2,2); \
    INSTANTIATION_HELPER_ALL2(P(V),int,3,2); \
    INSTANTIATION_HELPER_ALL2(P(V),int,3,3); \
    INSTANTIATION_HELPER_ALL2(P(V),int,4,3); \
    INSTANTIATION_HELPER_ALL2(P(V),int,4,4); \
    INSTANTIATION_HELPER_ALL2(P(V),T,2,2); \
    INSTANTIATION_HELPER_ALL2(P(V),T,3,2); \
    INSTANTIATION_HELPER_ALL2(P(V),T,3,3); \
    INSTANTIATION_HELPER_ALL3(P(V),T,2,2); \
    INSTANTIATION_HELPER_ALL3(P(V),T,3,2); \
    INSTANTIATION_HELPER_ALL3(P(V),T,3,3); \
    INSTANTIATION_HELPER_ALL_NEED1(P(V),T);

template class MPI_SOLIDS<VECTOR<float,1> >;
template class MPI_SOLIDS<VECTOR<float,2> >;
template class MPI_SOLIDS<VECTOR<float,3> >;
INSTANTIATION_HELPER_1D_ONLY(float);
INSTANTIATION_HELPER_ALL(P(VECTOR<float,2>),float);
INSTANTIATION_HELPER_ALL(P(VECTOR<float,3>),float);
INSTANTIATION_HELPER_ALL_NEED1(P(VECTOR<float,1>),float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MPI_SOLIDS<VECTOR<double,1> >;
template class MPI_SOLIDS<VECTOR<double,2> >;
template class MPI_SOLIDS<VECTOR<double,3> >;
INSTANTIATION_HELPER_1D_ONLY(double);
INSTANTIATION_HELPER_ALL(P(VECTOR<double,2>),double);
INSTANTIATION_HELPER_ALL(P(VECTOR<double,3>),double);
INSTANTIATION_HELPER_ALL_NEED1(P(VECTOR<double,1>),double);
#endif

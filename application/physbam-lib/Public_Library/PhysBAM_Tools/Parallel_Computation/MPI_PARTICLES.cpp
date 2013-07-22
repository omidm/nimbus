//#####################################################################
// Copyright 2011, Wen Zheng
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_PARTICLES
//#####################################################################
#ifdef USE_MPI
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Math_Tools/READ_WRITE_RANGE.h>
#endif
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PARTICLES.h>
using namespace PhysBAM;

#ifdef USE_MPI
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> MPI_PARTICLES<T_GRID>::
MPI_PARTICLES(T_MPI_GRID& mpi_grid_input)
    :mpi_grid(mpi_grid_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> MPI_PARTICLES<T_GRID>::
~MPI_PARTICLES()
{
}
//#####################################################################
// Function Update_Boundary_Particle_Mapping
//#####################################################################
template<class T_GRID> template<class T_ARRAY> void MPI_PARTICLES<T_GRID>::
Update_Boundary_Particle_Mapping(T_ARRAY& array_collection,ARRAY_VIEW<TV> positions,const int bandwidth,const ARRAY<int>* indices,const bool include_ghost_regions,const bool include_corners)
{
    const T_GRID& local_grid=mpi_grid.local_grid;
    RANGE<TV_INT> sentinels=RANGE<TV_INT>::Zero_Box();
    neighbor_ranks=include_corners?mpi_grid.all_neighbor_ranks:mpi_grid.side_neighbor_ranks;
    neighbor_directions=include_corners?mpi_grid.all_neighbor_directions:mpi_grid.side_neighbor_directions;
    ARRAY<MPI::Request> requests;
    mpi_grid.Find_Boundary_Regions(send_regions,sentinels,false,RANGE<VECTOR<int,1> >(include_ghost_regions?1-bandwidth:0,bandwidth-1),include_corners,include_ghost_regions,local_grid);
    // send
    send_particles.Clean_Memory();send_particles.Resize(send_regions.m);
    send_particle_counts.Resize(send_regions.m);ARRAYS_COMPUTATIONS::Fill(send_particle_counts,0);
    for(int n=1;n<=send_regions.m;n++)if(neighbor_ranks(n)!=MPI::PROC_NULL){
        RANGE<TV> send_region_domain(local_grid.Node(send_regions(n).min_corner),local_grid.Node(send_regions(n).max_corner+1));
        if(!indices){for(int i=1;i<=positions.m;i++) if(Inside_Local_Domain(positions(i),send_region_domain)) send_particles(n).Append(i);}
        else{for(int i=1;i<=indices->m;i++) if(Inside_Local_Domain(positions((*indices)(i)),send_region_domain)) send_particles(n).Append((*indices)(i));}
        send_particle_counts(n)=send_particles(n).m;
        requests.Append(mpi_grid.comm->Isend(&send_particle_counts(n),1,MPI::INT,neighbor_ranks(n),mpi_grid.Get_Send_Tag(neighbor_directions(n))));}
    // receive
    recv_particle_counts.Resize(send_regions.m);ARRAYS_COMPUTATIONS::Fill(recv_particle_counts,0);
    for(int n=1;n<=send_regions.m;n++)if(neighbor_ranks(n)!=MPI::PROC_NULL){
        requests.Append(mpi_grid.comm->Irecv(&recv_particle_counts(n),1,MPI::INT,neighbor_ranks(n),mpi_grid.Get_Recv_Tag(neighbor_directions(n))));}
    // finish
    MPI_UTILITIES::Wait_All(requests);
    // create space for receiving particles
    recv_particle_range.Clean_Memory();recv_particle_range.Resize(send_regions.m);
    for(int n=1;n<=send_regions.m;n++)if(neighbor_ranks(n)!=MPI::PROC_NULL && recv_particle_counts(n)>0){
        recv_particle_range(n).min_corner=VECTOR<int,1>(array_collection.Size()+1);
        array_collection.Resize(array_collection.Size()+recv_particle_counts(n));
        recv_particle_range(n).max_corner=VECTOR<int,1>(array_collection.Size());}
}
//#####################################################################
// Function Output_Send_And_Recv_Indices
//#####################################################################
template<class T_GRID> void MPI_PARTICLES<T_GRID>::
Output_Send_And_Recv_Indices(bool verbose)const
{
    for(int n=1;n<=send_regions.m;n++)if(neighbor_ranks(n)!=MPI::PROC_NULL){
        LOG::cout<<"Send "<<mpi_grid.rank<<"->"<<neighbor_ranks(n)<<" num="<<send_particles(n).m;
        if(verbose) LOG::cout<<": "<<send_particles(n)<<std::endl;
        else LOG::cout<<std::endl;
        LOG::cout<<"Recv "<<mpi_grid.rank<<"<-"<<neighbor_ranks(n)<<" num="<<(recv_particle_range(n).Empty()?0:(recv_particle_range(n).max_corner(1)-recv_particle_range(n).min_corner(1)+1));
        if(verbose) LOG::cout<<": "<<recv_particle_range(n).min_corner(1)<<"-"<<recv_particle_range(n).max_corner(1)<<std::endl;
        else LOG::cout<<std::endl;}
}
//#####################################################################
// Function Exchange_Boundary_Particles
//#####################################################################
template<class T_GRID> template<class T2> void MPI_PARTICLES<T_GRID>::
Exchange_Boundary_Particle_Values(ARRAY_VIEW<T2> values) const
{
    ARRAY<MPI_PACKAGE> packages;ARRAY<MPI::Request> requests;
    // send
    ARRAY<INDIRECT_ARRAY<ARRAY_VIEW<T2> >* > send_values(send_regions.m);
    for(int n=1;n<=send_regions.m;n++)if(neighbor_ranks(n)!=MPI::PROC_NULL && send_particles(n).m>0){
        send_values(n)=new INDIRECT_ARRAY<ARRAY_VIEW<T2> >(values,send_particles(n));
        MPI_PACKAGE package(*send_values(n));
        packages.Append(package);requests.Append(package.Isend(*mpi_grid.comm,neighbor_ranks(n),mpi_grid.Get_Send_Tag(neighbor_directions(n))));}
    // receive
    for(int n=1;n<=send_regions.m;n++)if(neighbor_ranks(n)!=MPI::PROC_NULL && !recv_particle_range(n).Empty()){
        MPI_PACKAGE package(values,recv_particle_range(n));
        packages.Append(package);requests.Append(package.Irecv(*mpi_grid.comm,neighbor_ranks(n),mpi_grid.Get_Recv_Tag(neighbor_directions(n))));}
    // finish
    MPI_UTILITIES::Wait_All(requests);MPI_PACKAGE::Free_All(packages);
    send_values.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Backwards_Exchange_Boundary_Particle_Values
//#####################################################################
template<class T_GRID> template<class T2> void MPI_PARTICLES<T_GRID>::
Backwards_Exchange_Boundary_Particle_Values(ARRAY_VIEW<T2> send_values,ARRAY<ARRAY_VIEW<T2>*> recv_values) const
{
    ARRAY<MPI_PACKAGE> packages;ARRAY<MPI::Request> requests;
    // send
    for(int n=1;n<=send_regions.m;n++)if(neighbor_ranks(n)!=MPI::PROC_NULL && !recv_particle_range(n).Empty()){
        MPI_PACKAGE package(send_values,recv_particle_range(n));
        packages.Append(package);requests.Append(package.Isend(*mpi_grid.comm,neighbor_ranks(n),mpi_grid.Get_Send_Tag(neighbor_directions(n))));}
    // receive
    ARRAY<INDIRECT_ARRAY<ARRAY_VIEW<T2> >* > recv_values_indexed(send_regions.m);
    for(int n=1;n<=send_regions.m;n++)if(neighbor_ranks(n)!=MPI::PROC_NULL && send_particles(n).m>0){
        recv_values_indexed(n)=new INDIRECT_ARRAY<ARRAY_VIEW<T2> >(*recv_values(n),send_particles(n));
        MPI_PACKAGE package(*recv_values_indexed(n));
        packages.Append(package);requests.Append(package.Irecv(*mpi_grid.comm,neighbor_ranks(n),mpi_grid.Get_Recv_Tag(neighbor_directions(n))));}
    // finish
    MPI_UTILITIES::Wait_All(requests);MPI_PACKAGE::Free_All(packages);
    recv_values_indexed.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Remove_Particles_Outside_Local_Domain
//#####################################################################
template<class T_GRID> void MPI_PARTICLES<T_GRID>::
Remove_Particles_Outside_Local_Domain(ARRAY_COLLECTION& array_collection,ARRAY_VIEW<TV> positions,const ARRAY<int>* indices)
{
    if(!indices){for(int i=1;i<=positions.m;i++) if(!Inside_Local_Domain(positions(i))) array_collection.Add_To_Deletion_List(i);}
    else{for(int i=1;i<=indices->m;i++) if(!Inside_Local_Domain(positions((*indices)(i)))) array_collection.Add_To_Deletion_List((*indices)(i));}
    if(array_collection.deletion_list.m>0) array_collection.Delete_Elements_On_Deletion_List();
}
//#####################################################################
// Function Synchronize_Flag
//#####################################################################
template<class T_GRID> void MPI_PARTICLES<T_GRID>::
Synchronize_Flag(bool& flag) const
{
    int flag_local=flag,flag_global=flag;
    MPI_UTILITIES::Reduce(flag_local,flag_global,MPI::MAX,*mpi_grid.comm);
    flag=(flag_global>0);
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d)                                        \
    template class MPI_PARTICLES<GRID<VECTOR<T,d> > >;                  \
    template void MPI_PARTICLES<GRID<VECTOR<T,d> > >::Update_Boundary_Particle_Mapping(ARRAY_COLLECTION&,ARRAY_VIEW<VECTOR<T,d> >,const int,const ARRAY<int>*,const bool,const bool); \
    template void MPI_PARTICLES<GRID<VECTOR<T,d> > >::Update_Boundary_Particle_Mapping(ARRAY<VECTOR<int,d+1> >&,ARRAY_VIEW<VECTOR<T,d> >,const int,const ARRAY<int>*,const bool,const bool); \
    template void MPI_PARTICLES<GRID<VECTOR<T,d> > >::Update_Boundary_Particle_Mapping(ARRAY<VECTOR<int,d> >&,ARRAY_VIEW<VECTOR<T,d> >,const int,const ARRAY<int>*,const bool,const bool); \
    template void MPI_PARTICLES<GRID<VECTOR<T,d> > >::Exchange_Boundary_Particle_Values(ARRAY_VIEW<VECTOR<int,d> >) const; \
    template void MPI_PARTICLES<GRID<VECTOR<T,d> > >::Exchange_Boundary_Particle_Values(ARRAY_VIEW<VECTOR<T,d> >) const; \
    template void MPI_PARTICLES<GRID<VECTOR<T,d> > >::Exchange_Boundary_Particle_Values(ARRAY_VIEW<VECTOR<T,d+1> >) const; \
    template void MPI_PARTICLES<GRID<VECTOR<T,d> > >::Exchange_Boundary_Particle_Values(ARRAY_VIEW<T>) const; \
    template void MPI_PARTICLES<GRID<VECTOR<T,d> > >::Exchange_Boundary_Particle_Values(ARRAY_VIEW<int>) const; \
    template void MPI_PARTICLES<GRID<VECTOR<T,d> > >::Exchange_Boundary_Particle_Values(ARRAY_VIEW<bool>) const; \
    template void MPI_PARTICLES<GRID<VECTOR<T,d> > >::Exchange_Boundary_Particle_Values(ARRAY_VIEW<VECTOR<int,d+1> >) const; \
    template void MPI_PARTICLES<GRID<VECTOR<T,d> > >::Backwards_Exchange_Boundary_Particle_Values(ARRAY_VIEW<VECTOR<int,d> >,ARRAY<ARRAY_VIEW<VECTOR<int,d> >*>) const; \
    template void MPI_PARTICLES<GRID<VECTOR<T,d> > >::Backwards_Exchange_Boundary_Particle_Values(ARRAY_VIEW<VECTOR<T,d> >,ARRAY<ARRAY_VIEW<VECTOR<T,d> >*>) const; \
    template void MPI_PARTICLES<GRID<VECTOR<T,d> > >::Backwards_Exchange_Boundary_Particle_Values(ARRAY_VIEW<VECTOR<T,d+1> >,ARRAY<ARRAY_VIEW<VECTOR<T,d+1> >*>) const; \
    template void MPI_PARTICLES<GRID<VECTOR<T,d> > >::Backwards_Exchange_Boundary_Particle_Values(ARRAY_VIEW<T>,ARRAY<ARRAY_VIEW<T>*>) const; \
    template void MPI_PARTICLES<GRID<VECTOR<T,d> > >::Backwards_Exchange_Boundary_Particle_Values(ARRAY_VIEW<int>,ARRAY<ARRAY_VIEW<int>*>) const;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif
#endif

#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Particles/VORTICITY_PARTICLES.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_RLE_PARTICLES.h>
#ifdef USE_MPI
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_BLOCK_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_BLOCK_3D.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_RLE.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#endif
namespace PhysBAM{

#ifdef USE_MPI

//#####################################################################
// Function ISend_Particles
//#####################################################################
template<class T_GRID,class T_PARTICLES,class T> void
ISend_Particles(const MPI_RLE_GRID<T_GRID>& mpi_grid,ARRAY<T_PARTICLES*>* particles_in_short_cells,T_PARTICLES* particles_in_long_cells,const int bandwidth,const T ghost_distance,const int tag,
    ARRAY<ARRAY<char> >& buffers,ARRAY<MPI::Request>& requests)
{
    typedef typename T_GRID::HORIZONTAL_GRID T_HORIZONTAL_GRID;
    typedef typename T_HORIZONTAL_GRID::VECTOR_T TV_HORIZONTAL;
    typedef typename T_HORIZONTAL_GRID::VECTOR_INT TV_HORIZONTAL_INT;
    typedef typename T_GRID::BOX_HORIZONTAL T_BOX_HORIZONTAL;
    typedef typename T_GRID::BOX_HORIZONTAL_INT T_BOX_HORIZONTAL_INT;
    typedef typename T_GRID::BLOCK_ITERATOR BLOCK_ITERATOR;
    ARRAY<T_BOX_HORIZONTAL_INT> send_regions;
    mpi_grid.Find_Boundary_Regions(send_regions,BLOCK_ITERATOR::Sentinels(),false,RANGE<VECTOR<int,1> >(0,bandwidth-1),true,false);
    // figure out domains of neighboring processors (TODO: clean up)
    T_BOX_HORIZONTAL domain=mpi_grid.local_grid.horizontal_grid.Domain();
    ARRAY<T_BOX_HORIZONTAL> neighbor_domains(T_HORIZONTAL_GRID::number_of_one_ring_neighbors_per_cell);
    for(int n=1;n<=neighbor_domains.m;n++){
        TV_HORIZONTAL neighbor_direction=TV_HORIZONTAL(T_HORIZONTAL_GRID::One_Ring_Neighbor(TV_HORIZONTAL_INT(),n));
        neighbor_domains(n)=(domain+neighbor_direction*domain.Edge_Lengths()).Thickened((T).01*mpi_grid.local_grid.Minimum_Edge_Length());}
    if(ghost_distance){
        domain.Change_Size(-ghost_distance);
        for(int n=1;n<=neighbor_domains.m;n++)neighbor_domains(n).Change_Size(ghost_distance);}
    // send particles that have exited the domain
    buffers.Resize(send_regions.m);
    ARRAY<ARRAY<PAIR<T_PARTICLES*,int> > > exchange_particles(send_regions.m);
    // TODO: this is inefficient because it does an entire box check even if only some sides are needed, and doesn't locally delete sent particles
    for(int n=1;n<=send_regions.m;n++)if(mpi_grid.all_neighbor_ranks(n)!=MPI::PROC_NULL){
        exchange_particles(n).Preallocate(100);
        if(particles_in_short_cells)
            for(BLOCK_ITERATOR block(mpi_grid.local_grid,send_regions(n));block;block++){int b=block.Block();if((*particles_in_short_cells)(b)){
                T_PARTICLES& cell_particles=*(*particles_in_short_cells)(b);
                for(int i=1;i<=cell_particles.array_collection->Size();i++){
                    TV_HORIZONTAL X_horizontal=cell_particles.X(i).Horizontal_Vector();
                    if(!domain.Lazy_Inside(X_horizontal) && neighbor_domains(n).Lazy_Inside(X_horizontal))
                        exchange_particles(n).Append(PAIR<T_PARTICLES*,int>(&cell_particles,i));}}}
        if(particles_in_long_cells)
            for(int i=1;i<=particles_in_long_cells->array_collection->Size();i++){
                TV_HORIZONTAL X_horizontal=particles_in_long_cells->X(i).Horizontal_Vector();
                if(!domain.Lazy_Inside(X_horizontal) && neighbor_domains(n).Lazy_Inside(X_horizontal))
                    exchange_particles(n).Append(PAIR<T_PARTICLES*,int>(particles_in_long_cells,i));}
        requests.Append(ISend_Particles(mpi_grid,exchange_particles(n),mpi_grid.all_neighbor_ranks(n),mpi_grid.all_neighbor_directions(n),tag,buffers(n)));}
}
//#####################################################################
// Function ISend_Particles
//#####################################################################
template<class T_GRID,class T_PARTICLES,class TV_HORIZONTAL_INT> MPI::Request
ISend_Particles(const MPI_RLE_GRID<T_GRID>& mpi_grid,const ARRAY<PAIR<T_PARTICLES*,int> >& particles,const int destination_rank,const TV_HORIZONTAL_INT& destination_direction,
    const int tag,ARRAY<char>& buffer)
{
    MPI::Comm& comm=*mpi_grid.comm;
    int position=0;
    buffer.Resize(MPI_UTILITIES::Pack_Size(destination_direction,comm)+MPI_UTILITIES::Pack_Size(particles.m,comm)
        +(particles.m?particles.m*MPI_UTILITIES::Pack_Size(*particles(1).x,comm):0));
    MPI_UTILITIES::Pack(destination_direction,buffer,position,comm);
    MPI_UTILITIES::Pack(particles.m,buffer,position,comm);
    for(int i=1;i<=particles.m;i++) MPI_UTILITIES::Pack(*particles(i).x,particles(i).y,buffer,position,comm);
    return comm.Isend(buffer.Get_Array_Pointer(),position,MPI::PACKED,destination_rank,tag);
}
//#####################################################################
// Function Recv_Particles
//#####################################################################
template<class T_GRID,class T_PARTICLES> void
Recv_Particles(const MPI_RLE_GRID<T_GRID>& mpi_grid,T_PARTICLES& particles,const int tag)
{
    typedef typename T_GRID::VECTOR_T TV;
    MPI::Comm& comm=*mpi_grid.comm;
    MPI::Status probe_status;
    comm.Probe(MPI::ANY_SOURCE,tag,probe_status);
    ARRAY<char> buffer(probe_status.Get_count(MPI::PACKED));int position=0;
    MPI::Status status;
    comm.Recv(buffer.Get_Array_Pointer(),buffer.m,MPI::PACKED,probe_status.Get_source(),tag,status);
    typename T_GRID::HORIZONTAL_GRID::VECTOR_INT direction;MPI_UTILITIES::Unpack(direction,buffer,position,comm);
    TV wrap_offset=-mpi_grid.Wrap_Offset(-direction).Insert(0,2);
    int m;MPI_UTILITIES::Unpack(m,buffer,position,comm);
    int number=particles.array_collection->Size();particles.array_collection->Add_Elements(m);
    for(int i=1;i<=m;i++){
        MPI_UTILITIES::Unpack(particles,++number,buffer,position,comm);
        particles.X(number)+=wrap_offset;}
}
//#####################################################################
// Function Recv_Particles
//#####################################################################
template<class T_GRID,class T_PARTICLES> void
Recv_Particles(const MPI_RLE_GRID<T_GRID>& mpi_grid,PARTICLE_LEVELSET_RLE<T_GRID>& particle_levelset,ARRAY<T_PARTICLES*>& particles,T_PARTICLES* particles_in_long_cells,const int tag)
{
    typedef typename T_GRID::VECTOR_T TV;
    MPI::Comm& comm=*mpi_grid.comm;
    MPI::Status probe_status;
    comm.Probe(MPI::ANY_SOURCE,tag,probe_status);
    ARRAY<char> buffer(probe_status.Get_count(MPI::PACKED));int position=0;
    MPI::Status status;
    comm.Recv(buffer.Get_Array_Pointer(),buffer.m,MPI::PACKED,probe_status.Get_source(),tag,status);
    typename T_GRID::HORIZONTAL_GRID::VECTOR_INT direction;MPI_UTILITIES::Unpack(direction,buffer,position,comm);
    TV wrap_offset=-mpi_grid.Wrap_Offset(-direction).Insert(0,2);
    int m;MPI_UTILITIES::Unpack(m,buffer,position,comm);
    T_PARTICLES recv_particles; // message will unpack into this particle object
    for(int i=1;i<=m;i++){
        if(!recv_particles.array_collection->Size()) recv_particles.array_collection->Add_Element();
        MPI_UTILITIES::Unpack(recv_particles,1,buffer,position,comm);
        recv_particles.X(1)+=wrap_offset;
        int final_b=mpi_grid.local_grid.Clamped_Block_Index(recv_particles.X(1),1);
        if(!final_b){
            if(particles_in_long_cells) particle_levelset.Move_Particle(recv_particles,*particles_in_long_cells,1);}
        else{
            if(!particles(final_b)) particles(final_b)=particle_levelset.template Allocate_Particle<T_PARTICLES>();
            particle_levelset.Move_Particle(recv_particles,*particles(final_b),1);}}
}
//#####################################################################
// Function Exchange_Boundary_Particles
//#####################################################################
template<class T_GRID,class T_PARTICLES> void
Exchange_Boundary_Particles(MPI_RLE_GRID<T_GRID>& mpi_grid,T_PARTICLES& particles,const typename T_GRID::SCALAR ghost_distance)
{
    int tag=mpi_grid.Get_Unique_Tag();
    ARRAY<MPI::Request> requests;
    ARRAY<ARRAY<char> > buffers;
    ISend_Particles(mpi_grid,(ARRAY<T_PARTICLES*>*)0,&particles,0,ghost_distance,tag,buffers,requests);
    for(int message=1;message<=requests.m;message++) Recv_Particles(mpi_grid,particles,tag);
    MPI_UTILITIES::Wait_All(requests);
}
//#####################################################################
// Function Exchange_Boundary_Particles
//#####################################################################
template<class T_GRID,class T_PARTICLES> void
Exchange_Boundary_Particles(MPI_RLE_GRID<T_GRID>& mpi_grid,PARTICLE_LEVELSET_RLE<T_GRID>& particle_levelset,ARRAY<T_PARTICLES*>& particles,T_PARTICLES* particles_in_long_cells,
    const int bandwidth)
{
    int tag=mpi_grid.Get_Unique_Tag();
    ARRAY<MPI::Request> requests;
    ARRAY<ARRAY<char> > buffers;
    ISend_Particles(mpi_grid,&particles,particles_in_long_cells,bandwidth,0,tag,buffers,requests);
    for(int message=1;message<=requests.m;message++) Recv_Particles(mpi_grid,particle_levelset,particles,particles_in_long_cells,tag);
    MPI_UTILITIES::Wait_All(requests);
}
//#####################################################################

#else

//#####################################################################
template<class T_GRID,class T_PARTICLES> void Exchange_Boundary_Particles(MPI_RLE_GRID<T_GRID>& mpi_grid,T_PARTICLES&,const typename T_GRID::SCALAR){PHYSBAM_NOT_IMPLEMENTED();}
template<class T_GRID,class T_PARTICLES> void Exchange_Boundary_Particles(MPI_RLE_GRID<T_GRID>& mpi_grid,PARTICLE_LEVELSET_RLE<T_GRID>&,ARRAY<T_PARTICLES*>&,
    T_PARTICLES*,const int){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################

#endif

//#####################################################################
#define INSTANTIATION_HELPER_T(T_GRID,T,d) \
    template void Exchange_Boundary_Particles(MPI_RLE_GRID<T_GRID >&,VORTICITY_PARTICLES<VECTOR<T,d> >&,const T); \
    template void Exchange_Boundary_Particles(MPI_RLE_GRID<T_GRID >&,PARTICLE_LEVELSET_RLE<T_GRID >&,ARRAY<PARTICLE_LEVELSET_PARTICLES<VECTOR<T,d> >*>&, \
        PARTICLE_LEVELSET_PARTICLES<VECTOR<T,d> >*,const int); \
    template void Exchange_Boundary_Particles(MPI_RLE_GRID<T_GRID >&,PARTICLE_LEVELSET_RLE<T_GRID >&,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<T,d> >*>&, \
        PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<T,d> >*,const int);
#define INSTANTIATION_HELPER(T) INSTANTIATION_HELPER_T(RLE_GRID_2D< T >,T,2);INSTANTIATION_HELPER_T(RLE_GRID_3D< T >,T,3);
INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif
}
#endif

//#####################################################################
// Copyright 2005-2008, Geoffrey Irving, Eran Guendelman, Frank Losasso, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Particles/VORTICITY_PARTICLES.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_UNIFORM_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#ifdef USE_MPI
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#endif
namespace PhysBAM{

#ifdef USE_MPI

//#####################################################################
// Function ISend_Particles
//#####################################################################
template<class T_GRID,class T_PARTICLES> MPI::Request
ISend_Particles(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,const ARRAY<PAIR<T_PARTICLES*,int> >& particles,const int destination_rank,
    const typename T_GRID::VECTOR_INT& destination_direction,const int tag,ARRAY<char>& buffer)
{
    int position=0;
    MPI::Comm& comm=*mpi_grid.comm;
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
template<class T_GRID,class T_PARTICLES,class T_ARRAYS_PARTICLES> void
Recv_Particles(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int tag,
    const MPI::Status& probe_status,PARTICLE_LEVELSET<T_GRID>& particle_levelset)
{
    typedef typename T_GRID::VECTOR_T TV;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    MPI::Comm& comm=*mpi_grid.comm;
    T_PARTICLES* recv_particles(template_particles.Clone()); // message will unpack into this particle object
    recv_particles->array_collection->Add_Element();
    ARRAY<char> buffer(probe_status.Get_count(MPI::PACKED));int position=0;
    MPI::Status status;
    comm.Recv(buffer.Get_Array_Pointer(),buffer.m,MPI::PACKED,probe_status.Get_source(),tag,status);
    TV_INT direction;MPI_UTILITIES::Unpack(direction,buffer,position,comm);
    TV wrap_offset=-mpi_grid.Wrap_Offset(-direction);
    int m;MPI_UTILITIES::Unpack(m,buffer,position,comm);
    RANGE<TV> domain=mpi_grid.local_grid.Domain();
    for(int i=1;i<=m;i++){
        MPI_UTILITIES::Unpack(*recv_particles,1,buffer,position,comm);
        TV& X=recv_particles->X(1);X+=wrap_offset;
        if(!domain.Lazy_Inside(X)) continue;
        TV_INT final_b=mpi_grid.local_grid.Block_Index(X,1); // TODO: check whether this is a good block
        if(!particles(final_b)) particles(final_b)=particle_levelset.Allocate_Particles(template_particles);
        particle_levelset.Copy_Particle(*recv_particles,*particles(final_b),1);}
    delete recv_particles;
}
//#####################################################################
// Function Recv_Block_Particles
//#####################################################################
template<class T_GRID,class T_PARTICLES,class T_ARRAYS_PARTICLES> void
Recv_Block_Particles(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int tag,
    const MPI::Status& probe_status,PARTICLE_LEVELSET<T_GRID>& particle_levelset)
{
    typedef typename T_GRID::VECTOR_T TV;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    MPI::Comm& comm=*mpi_grid.comm;
    T_PARTICLES* recv_particles(template_particles.Clone()); // message will unpack into this particle object
    recv_particles->array_collection->Add_Element();
    ARRAY<char> buffer(probe_status.Get_count(MPI::PACKED));int position=0;
    MPI::Status status;
    comm.Recv(buffer.Get_Array_Pointer(),buffer.m,MPI::PACKED,probe_status.Get_source(),tag,status);
    TV_INT direction;MPI_UTILITIES::Unpack(direction,buffer,position,comm);
    TV wrap_offset=-mpi_grid.Wrap_Offset(-direction);
    int m;MPI_UTILITIES::Unpack(m,buffer,position,comm);
    RANGE<TV> domain=mpi_grid.local_grid.Domain();
    for(int i=1;i<=m;i++){
        MPI_UTILITIES::Unpack(*recv_particles,1,buffer,position,comm);
        TV& X=recv_particles->X(1);X+=wrap_offset;
        if(domain.Lazy_Inside(X)) continue;
        TV_INT final_b=mpi_grid.local_grid.Block_Index(X,1); // TODO: check whether this is a good block
        if(!particles(final_b)) particles(final_b)=particle_levelset.Allocate_Particles(template_particles);
        particle_levelset.Copy_Particle(*recv_particles,*particles(final_b),1);}
    delete recv_particles;
}
//#####################################################################
// Function Recv_Ghost_Particles
//#####################################################################
template<class T_GRID,class T_PARTICLES,class T_ARRAYS_PARTICLES> void
Recv_Ghost_Particles(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int tag,
    const MPI::Status& probe_status,const int bandwidth,PARTICLE_LEVELSET<T_GRID>& particle_levelset)
{
    typedef typename T_GRID::VECTOR_T TV;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    MPI::Comm& comm=*mpi_grid.comm;
    T_PARTICLES* recv_particles(template_particles.Clone()); // message will unpack into this particle object
    recv_particles->array_collection->Add_Element();
    ARRAY<char> buffer(probe_status.Get_count(MPI::PACKED));int position=0;
    MPI::Status status;
    comm.Recv(buffer.Get_Array_Pointer(),buffer.m,MPI::PACKED,probe_status.Get_source(),tag,status);
    TV_INT direction;MPI_UTILITIES::Unpack(direction,buffer,position,comm);
    TV wrap_offset=-mpi_grid.Wrap_Offset(-direction);
    int m;MPI_UTILITIES::Unpack(m,buffer,position,comm);
    RANGE<TV> domain=mpi_grid.local_grid.Domain();
    for(int i=1;i<=m;i++){
        MPI_UTILITIES::Unpack(*recv_particles,1,buffer,position,comm);
        TV& X=recv_particles->X(1);X+=wrap_offset;
        if(domain.Lazy_Inside(X)) continue;
        TV_INT final_b=mpi_grid.local_grid.Block_Index(X,bandwidth-1); // TODO: check whether this is a good block
        if(!particles(final_b)) particles(final_b)=particle_levelset.Allocate_Particles(template_particles);
        particle_levelset.Copy_Particle(*recv_particles,*particles(final_b),1);}
    delete recv_particles;
}
//#####################################################################
// Function ISend_Particles
//#####################################################################
template<class T_GRID,class T_PARTICLES,class T> void
ISend_Particles(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,T_PARTICLES& particles,const T ghost_distance,const int tag,ARRAY<ARRAY<char> >& buffers,ARRAY<MPI::Request>& requests)
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;
    RANGE<TV> domain=mpi_grid.local_grid.Domain();
    ARRAY<RANGE<TV> > neighbor_domains(T_GRID::number_of_one_ring_neighbors_per_cell);
    for(int n=1;n<=neighbor_domains.m;n++){
        TV neighbor_direction=TV(T_GRID::One_Ring_Neighbor(TV_INT(),n));
        neighbor_domains(n)=(domain+neighbor_direction*domain.Edge_Lengths()).Thickened((T).01*mpi_grid.local_grid.Minimum_Edge_Length());}
    if(ghost_distance){
        domain.Change_Size(-ghost_distance);
        for(int n=1;n<=neighbor_domains.m;n++)neighbor_domains(n).Change_Size(ghost_distance);}
    // send particles that have exited the domain
    buffers.Resize(T_GRID::number_of_one_ring_neighbors_per_cell);
    ARRAY<ARRAY<PAIR<T_PARTICLES*,int> > > exchange_particles(T_GRID::number_of_one_ring_neighbors_per_cell);
    // TODO: this is inefficient because it does an entire box check even if only some sides are needed, and doesn't locally delete sent particles
    for(int n=1;n<=T_GRID::number_of_one_ring_neighbors_per_cell;n++)if(mpi_grid.all_neighbor_ranks(n)!=MPI::PROC_NULL){
        exchange_particles(n).Preallocate(100);
        for(int i=1;i<=particles.array_collection->Size();i++){
            if(!domain.Lazy_Inside(particles.X(i)) && neighbor_domains(n).Lazy_Inside(particles.X(i)))
                exchange_particles(n).Append(PAIR<T_PARTICLES*,int>(&particles,i));}
        requests.Append(ISend_Particles(mpi_grid,exchange_particles(n),mpi_grid.all_neighbor_ranks(n),mpi_grid.all_neighbor_directions(n),tag,buffers(n)));}
}
//#####################################################################
// Function Recv_Particles
//#####################################################################
template<class T_GRID,class T_PARTICLES> void
Recv_Particles(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,T_PARTICLES& particles,const int tag)
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;
    MPI::Comm& comm=*mpi_grid.comm;
    MPI::Status probe_status;
    comm.Probe(MPI::ANY_SOURCE,tag,probe_status);
    ARRAY<char> buffer(probe_status.Get_count(MPI::PACKED));int position=0;
    MPI::Status status;
    comm.Recv(buffer.Get_Array_Pointer(),buffer.m,MPI::PACKED,probe_status.Get_source(),tag,status);
    TV_INT direction;MPI_UTILITIES::Unpack(direction,buffer,position,comm);
    TV wrap_offset=-mpi_grid.Wrap_Offset(-direction);
    int m;MPI_UTILITIES::Unpack(m,buffer,position,comm);
    int number=particles.array_collection->Size();particles.array_collection->Add_Elements(m);
    for(int i=1;i<=m;i++){
        MPI_UTILITIES::Unpack(particles,++number,buffer,position,comm);
        particles.X(number)+=wrap_offset;}
}
//#####################################################################
// Function Exchange_Boundary_Particles
//#####################################################################
template<class T_GRID,class T_PARTICLES,class T_ARRAYS_PARTICLES> void
Exchange_Boundary_Particles_Threaded(const THREADED_UNIFORM_GRID<T_GRID>& threaded_grid,const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int bandwidth,PARTICLE_LEVELSET<T_GRID>& particle_levelset)
{
#ifdef USE_PTHREADS
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    STATIC_ASSERT((IS_SAME<T_PARTICLES,typename REMOVE_POINTER<typename T_ARRAYS_PARTICLES::ELEMENT>::TYPE>::value));
    ARRAY<RANGE<TV_INT> > send_regions;
    // this way Find_Boundary_Regions will return block indices which for uniform grids are the node index not minimum corner cell
    RANGE<TV_INT> sentinels=RANGE<TV_INT>(TV_INT(),TV_INT::All_Ones_Vector());
    threaded_grid.Find_Boundary_Regions(send_regions,sentinels,false,RANGE<VECTOR<int,1> >(0,bandwidth-1),true,false);
    // send particles that have exited the domain
    ARRAY<ARRAY<PAIR<T_PARTICLES*,int> > > exchange_particles(T_GRID::number_of_one_ring_neighbors_per_cell);
    // TODO: this is inefficient because it does an entire box check even if only some sides are needed, and send a lot more corner particles than it should
    RANGE<TV> domain=threaded_grid.local_grid.Domain();
    for(int n=1;n<=send_regions.m;n++) if(threaded_grid.all_neighbor_ranks(n)!=-1){
        exchange_particles(n).Preallocate(100);
        for(NODE_ITERATOR iterator(threaded_grid.local_grid,send_regions(n));iterator.Valid();iterator.Next())if(particles(iterator.Node_Index())){
            T_PARTICLES* cell_particles=particles(iterator.Node_Index());
            while(cell_particles){for(int i=1;i<=cell_particles->array_collection->Size();i++)if(!domain.Lazy_Inside(cell_particles->X(i))) // could be optimized further
                exchange_particles(n).Append(PAIR<T_PARTICLES*,int>(cell_particles,i));cell_particles=cell_particles->next;}} // TODO: delete the particle locally?
        if(!exchange_particles(n).m) continue;
        THREAD_PACKAGE pack((sizeof(int)+sizeof(T_PARTICLES*))*exchange_particles(n).m+sizeof(int));pack.send_tid=threaded_grid.tid;pack.recv_tid=threaded_grid.all_neighbor_ranks(n);
        int position=0;*(int*)(&pack.buffer(position+1))=exchange_particles(n).m;position+=sizeof(int);
        for(int i=1;i<=exchange_particles(n).m;i++){
            *(T_PARTICLES**)(&pack.buffer(position+1))=exchange_particles(n)(i).x;position+=sizeof(T_PARTICLES*);
            *(int*)(&pack.buffer(position+1))=exchange_particles(n)(i).y;position+=sizeof(int);}
        pthread_mutex_lock(threaded_grid.lock);
        threaded_grid.buffers.Append(pack);
        pthread_mutex_unlock(threaded_grid.lock);}
    pthread_barrier_wait(threaded_grid.barr);
    // probe and receive
    for(int buf=1;buf<=threaded_grid.buffers.m;buf++){if(threaded_grid.tid!=threaded_grid.buffers(buf).recv_tid) continue;
        THREAD_PACKAGE& pack=threaded_grid.buffers(buf);int position=0;int size=*(int*)(&pack.buffer(position+1));position+=sizeof(int);
        for(int i=1;i<=size;i++){
            T_PARTICLES* recv_particles=*(T_PARTICLES**)(&pack.buffer(position+1));position+=sizeof(T_PARTICLES*);
            int index=*(int*)(&pack.buffer(position+1));position+=sizeof(int);
            TV& X=recv_particles->X(index);
            TV_INT final_b=threaded_grid.local_grid.Block_Index(X,bandwidth-1); // TODO: check whether this is a good block
            if(!particles(final_b)) particles(final_b)=particle_levelset.Allocate_Particles(template_particles);
            particle_levelset.Copy_Particle(*recv_particles,*particles(final_b),index);}}
    // wait for sends to complete
    pthread_barrier_wait(threaded_grid.barr);
    if(threaded_grid.tid==1) threaded_grid.buffers.m=0;
    pthread_barrier_wait(threaded_grid.barr);
#endif
}
template<class T_GRID,class T_PARTICLES,class T_ARRAYS_PARTICLES> void
Exchange_Boundary_Particles(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int bandwidth,PARTICLE_LEVELSET<T_GRID>& particle_levelset)
{
    if(mpi_grid.threaded_grid){Exchange_Boundary_Particles_Threaded(*mpi_grid.threaded_grid,template_particles,particles,bandwidth,particle_levelset);return;}
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    STATIC_ASSERT((IS_SAME<T_PARTICLES,typename REMOVE_POINTER<typename T_ARRAYS_PARTICLES::ELEMENT>::TYPE>::value));
    int tag=mpi_grid.Get_Unique_Tag();
    ARRAY<RANGE<TV_INT> > send_regions;
    // this way Find_Boundary_Regions will return block indices which for uniform grids are the node index not minimum corner cell
    RANGE<TV_INT> sentinels=RANGE<TV_INT>(TV_INT(),TV_INT::All_Ones_Vector());
    mpi_grid.Find_Boundary_Regions(send_regions,sentinels,false,RANGE<VECTOR<int,1> >(0,bandwidth-1),true,false);
    // send particles that have exited the domain
    ARRAY<MPI::Request> requests;
    ARRAY<ARRAY<char> > buffers(T_GRID::number_of_one_ring_neighbors_per_cell);
    ARRAY<ARRAY<PAIR<T_PARTICLES*,int> > > exchange_particles(T_GRID::number_of_one_ring_neighbors_per_cell);
    // TODO: this is inefficient because it does an entire box check even if only some sides are needed, and send a lot more corner particles than it should
    RANGE<TV> domain=mpi_grid.local_grid.Domain();
    for(int n=1;n<=send_regions.m;n++) if(mpi_grid.all_neighbor_ranks(n)!=MPI::PROC_NULL){
        exchange_particles(n).Preallocate(100);
        for(NODE_ITERATOR iterator(mpi_grid.local_grid,send_regions(n));iterator.Valid();iterator.Next())if(particles(iterator.Node_Index())){
            T_PARTICLES* cell_particles=particles(iterator.Node_Index());
            while(cell_particles){for(int i=1;i<=cell_particles->array_collection->Size();i++)if(!domain.Lazy_Inside(cell_particles->X(i))) // could be optimized further
                exchange_particles(n).Append(PAIR<T_PARTICLES*,int>(cell_particles,i));cell_particles=cell_particles->next;}} // TODO: delete the particle locally?
        requests.Append(ISend_Particles(mpi_grid,exchange_particles(n),mpi_grid.all_neighbor_ranks(n),mpi_grid.all_neighbor_directions(n),tag,buffers(n)));}
    // probe and receive
    for(int message=1;message<=requests.m;message++){
        MPI::Status probe_status;
        mpi_grid.comm->Probe(MPI::ANY_SOURCE,tag,probe_status);
        Recv_Particles(mpi_grid,template_particles,particles,tag,probe_status,particle_levelset);}
    // wait for sends to complete
    MPI_UTILITIES::Wait_All(requests);
}
//#####################################################################
// Function Exchange_Overlapping_Block_Particles
//#####################################################################
template<class T_GRID,class T_PARTICLES,class T_ARRAYS_PARTICLES> void
Exchange_Overlapping_Block_Particles_Threaded(const THREADED_UNIFORM_GRID<T_GRID>& threaded_grid,const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int bandwidth,PARTICLE_LEVELSET<T_GRID>& particle_levelset)
{
#ifdef USE_PTHREADS
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    STATIC_ASSERT((IS_SAME<T_PARTICLES,typename REMOVE_POINTER<typename T_ARRAYS_PARTICLES::ELEMENT>::TYPE>::value));
    ARRAY<RANGE<TV_INT> > send_regions;
    // this way Find_Boundary_Regions will return block indices which for uniform grids are the node index not minimum corner cell
    RANGE<TV_INT> sentinels=RANGE<TV_INT>(TV_INT(),TV_INT::All_Ones_Vector());
    threaded_grid.Find_Boundary_Regions(send_regions,sentinels,false,RANGE<VECTOR<int,1> >(0,bandwidth-1),true,false);
    // send particles that have exited the domain
    ARRAY<ARRAY<PAIR<T_PARTICLES*,int> > > exchange_particles(T_GRID::number_of_one_ring_neighbors_per_cell);
    // TODO: this is inefficient because it does an entire box check even if only some sides are needed, and send a lot more corner particles than it should
    RANGE<TV> block_domain=threaded_grid.local_grid.Domain();block_domain.Change_Size(-(T).5*threaded_grid.local_grid.DX());
    for(int n=1;n<=send_regions.m;n++) if(threaded_grid.all_neighbor_ranks(n)!=-1){
        exchange_particles(n).Preallocate(100);
        for(NODE_ITERATOR iterator(threaded_grid.local_grid,send_regions(n));iterator.Valid();iterator.Next())if(particles(iterator.Node_Index())){
            T_PARTICLES* cell_particles=particles(iterator.Node_Index());
            while(cell_particles){for(int i=1;i<=cell_particles->array_collection->Size();i++) if(!block_domain.Thickened(3*threaded_grid.local_grid.Minimum_Edge_Length()).Lazy_Inside(cell_particles->X(i))) 
                exchange_particles(n).Append(PAIR<T_PARTICLES*,int>(cell_particles,i));cell_particles=cell_particles->next;}} // TODO: delete the particle locally?
        if(!exchange_particles(n).m) continue;
        THREAD_PACKAGE pack((sizeof(int)+sizeof(T_PARTICLES*))*exchange_particles(n).m+sizeof(int));pack.send_tid=threaded_grid.tid;pack.recv_tid=threaded_grid.all_neighbor_ranks(n);
        int position=0;*(int*)(&pack.buffer(position+1))=exchange_particles(n).m;position+=sizeof(int);
        for(int i=1;i<=exchange_particles(n).m;i++){
            *(T_PARTICLES**)(&pack.buffer(position+1))=exchange_particles(n)(i).x;position+=sizeof(T_PARTICLES*);
            *(int*)(&pack.buffer(position+1))=exchange_particles(n)(i).y;position+=sizeof(int);}
        pthread_mutex_lock(threaded_grid.lock);
        threaded_grid.buffers.Append(pack);
        pthread_mutex_unlock(threaded_grid.lock);}
    // probe and receive
    pthread_barrier_wait(threaded_grid.barr);
    for(int buf=1;buf<=threaded_grid.buffers.m;buf++){if(threaded_grid.tid!=threaded_grid.buffers(buf).recv_tid) continue;
        THREAD_PACKAGE& pack=threaded_grid.buffers(buf);int position=0;int size=*(int*)(&pack.buffer(position+1));position+=sizeof(int);
        for(int i=1;i<=size;i++){
            T_PARTICLES* recv_particles=*(T_PARTICLES**)(&pack.buffer(position+1));position+=sizeof(T_PARTICLES*);
            int index=*(int*)(&pack.buffer(position+1));position+=sizeof(int);
            TV& X=recv_particles->X(index);
            TV_INT final_b=threaded_grid.local_grid.Block_Index(X,bandwidth-1); // TODO: check whether this is a good block
            if(!particles(final_b)) particles(final_b)=particle_levelset.Allocate_Particles(template_particles);
            particle_levelset.Copy_Particle(*recv_particles,*particles(final_b),index);}}
    // wait for sends to complete
    pthread_barrier_wait(threaded_grid.barr);
    if(threaded_grid.tid==1) threaded_grid.buffers.m=0;
    pthread_barrier_wait(threaded_grid.barr);
#endif
}
template<class T_GRID,class T_PARTICLES,class T_ARRAYS_PARTICLES> void
Exchange_Overlapping_Block_Particles(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int bandwidth,PARTICLE_LEVELSET<T_GRID>& particle_levelset)
{
    if(mpi_grid.threaded_grid){Exchange_Overlapping_Block_Particles_Threaded(*mpi_grid.threaded_grid,template_particles,particles,bandwidth,particle_levelset);return;}
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    STATIC_ASSERT((IS_SAME<T_PARTICLES,typename REMOVE_POINTER<typename T_ARRAYS_PARTICLES::ELEMENT>::TYPE>::value));
    int tag=mpi_grid.Get_Unique_Tag();
    ARRAY<RANGE<TV_INT> > send_regions;
    // this way Find_Boundary_Regions will return block indices which for uniform grids are the node index not minimum corner cell
    RANGE<TV_INT> sentinels=RANGE<TV_INT>(TV_INT(),TV_INT::All_Ones_Vector());
    mpi_grid.Find_Boundary_Regions(send_regions,sentinels,false,RANGE<VECTOR<int,1> >(0,bandwidth-1),true,false);
    // send particles that have exited the domain
    ARRAY<MPI::Request> requests;
    ARRAY<ARRAY<char> > buffers(T_GRID::number_of_one_ring_neighbors_per_cell);
    ARRAY<ARRAY<PAIR<T_PARTICLES*,int> > > exchange_particles(T_GRID::number_of_one_ring_neighbors_per_cell);
    // TODO: this is inefficient because it does an entire box check even if only some sides are needed, and send a lot more corner particles than it should
    RANGE<TV> block_domain=mpi_grid.local_grid.Domain();block_domain.Change_Size(-(T).5*mpi_grid.local_grid.DX());
    for(int n=1;n<=send_regions.m;n++) if(mpi_grid.all_neighbor_ranks(n)!=MPI::PROC_NULL){
        exchange_particles(n).Preallocate(100);
        for(NODE_ITERATOR iterator(mpi_grid.local_grid,send_regions(n));iterator.Valid();iterator.Next())if(particles(iterator.Node_Index())){
            T_PARTICLES* cell_particles=particles(iterator.Node_Index());
            while(cell_particles){for(int i=1;i<=cell_particles->array_collection->Size();i++) if(!block_domain.Thickened(3*mpi_grid.local_grid.Minimum_Edge_Length()).Lazy_Inside(cell_particles->X(i))) 
                exchange_particles(n).Append(PAIR<T_PARTICLES*,int>(cell_particles,i));cell_particles=cell_particles->next;}} // TODO: delete the particle locally?
        requests.Append(ISend_Particles(mpi_grid,exchange_particles(n),mpi_grid.all_neighbor_ranks(n),mpi_grid.all_neighbor_directions(n),tag,buffers(n)));}
    // probe and receive
    for(int message=1;message<=requests.m;message++){
        MPI::Status probe_status;
        mpi_grid.comm->Probe(MPI::ANY_SOURCE,tag,probe_status);
        Recv_Block_Particles(mpi_grid,template_particles,particles,tag,probe_status,particle_levelset);}
    // wait for sends to complete
    MPI_UTILITIES::Wait_All(requests);
}
//#####################################################################
// Function Exchange_Ghost_Particles
//#####################################################################
template<class T_GRID,class T_PARTICLES,class T_ARRAYS_PARTICLES> void
Exchange_Ghost_Particles(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int bandwidth,PARTICLE_LEVELSET<T_GRID>& particle_levelset)
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    STATIC_ASSERT((IS_SAME<T_PARTICLES,typename REMOVE_POINTER<typename T_ARRAYS_PARTICLES::ELEMENT>::TYPE>::value));
    int tag=mpi_grid.Get_Unique_Tag();
    ARRAY<RANGE<TV_INT> > send_regions;
    // this way Find_Boundary_Regions will return block indices which for uniform grids are the node index not minimum corner cell
    RANGE<TV_INT> sentinels=RANGE<TV_INT>(TV_INT(),TV_INT::All_Ones_Vector());
    mpi_grid.Find_Boundary_Regions(send_regions,sentinels,false,RANGE<VECTOR<int,1> >(0,bandwidth-1),true,false);
    // send particles that are in the ghost cells of an adjacent processor
    ARRAY<MPI::Request> requests;
    ARRAY<ARRAY<char> > buffers(T_GRID::number_of_one_ring_neighbors_per_cell);
    ARRAY<ARRAY<PAIR<T_PARTICLES*,int> > > exchange_particles(T_GRID::number_of_one_ring_neighbors_per_cell);
    // TODO: this is inefficient because it does an entire box check even if only some sides are needed, and send a lot more corner particles than it should
    RANGE<TV> domain=mpi_grid.local_grid.Domain();
    for(int n=1;n<=send_regions.m;n++) if(mpi_grid.all_neighbor_ranks(n)!=MPI::PROC_NULL){
        exchange_particles(n).Preallocate(100);
        for(NODE_ITERATOR iterator(mpi_grid.local_grid,send_regions(n));iterator.Valid();iterator.Next())if(particles(iterator.Node_Index())){
            T_PARTICLES* cell_particles=particles(iterator.Node_Index());
            while(cell_particles){for(int i=1;i<=cell_particles->array_collection->Size();i++) if(domain.Lazy_Inside(cell_particles->X(i))) // could be optimized further
                exchange_particles(n).Append(PAIR<T_PARTICLES*,int>(cell_particles,i));cell_particles=cell_particles->next;}} // TODO: delete the particle locally?
        requests.Append(ISend_Particles(mpi_grid,exchange_particles(n),mpi_grid.all_neighbor_ranks(n),mpi_grid.all_neighbor_directions(n),tag,buffers(n)));}
    // probe and receive
    for(int message=1;message<=requests.m;message++){
        MPI::Status probe_status;
        mpi_grid.comm->Probe(MPI::ANY_SOURCE,tag,probe_status);
        Recv_Ghost_Particles(mpi_grid,template_particles,particles,tag,probe_status,bandwidth,particle_levelset);}
    // wait for sends to complete
    MPI_UTILITIES::Wait_All(requests);
}
//#####################################################################
// Function Exchange_Boundary_Particles
//#####################################################################
template<class T_GRID,class T_PARTICLES> void
Exchange_Boundary_Particles_Flat(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,T_PARTICLES& particles,const typename T_GRID::SCALAR ghost_distance)
{
    int tag=mpi_grid.Get_Unique_Tag();
    ARRAY<MPI::Request> requests;
    ARRAY<ARRAY<char> > buffers;
    ISend_Particles(mpi_grid,particles,ghost_distance,tag,buffers,requests);
    for(int message=1;message<=requests.m;message++)Recv_Particles(mpi_grid,particles,tag);
    MPI_UTILITIES::Wait_All(requests);
}
//#####################################################################

#else

//#####################################################################
template<class T_GRID,class T_PARTICLES> void Exchange_Boundary_Particles_Flat(const MPI_UNIFORM_GRID<T_GRID>&,T_PARTICLES&,const typename T_GRID::SCALAR)
{PHYSBAM_NOT_IMPLEMENTED();}
template<class T_GRID,class T_PARTICLES,class T_ARRAYS_PARTICLES> void Exchange_Boundary_Particles(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,
    const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int bandwidth,PARTICLE_LEVELSET<T_GRID>& particle_levelset)
{PHYSBAM_NOT_IMPLEMENTED();}
template<class T_GRID,class T_PARTICLES,class T_ARRAYS_PARTICLES> void Exchange_Overlapping_Block_Particles(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,
    const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int bandwidth,PARTICLE_LEVELSET<T_GRID>& particle_levelset)
{PHYSBAM_NOT_IMPLEMENTED();}
template<class T_GRID,class T_PARTICLES,class T_ARRAYS_PARTICLES> void Exchange_Ghost_Particles(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,const T_PARTICLES& template_particles,
    T_ARRAYS_PARTICLES& particles,const int bandwidth,PARTICLE_LEVELSET<T_GRID>& particle_levelset)
{PHYSBAM_NOT_IMPLEMENTED();}

//#####################################################################

#endif

//#####################################################################
#define INSTANTIATION_HELPER(T,T_GRID,d) \
    template void Exchange_Boundary_Particles(const MPI_UNIFORM_GRID<T_GRID >&,const PARTICLE_LEVELSET_PARTICLES<T_GRID::VECTOR_T>&, \
        GRID_ARRAYS_POLICY<T_GRID >::ARRAYS_SCALAR::REBIND<PARTICLE_LEVELSET_PARTICLES<T_GRID::VECTOR_T>*>::TYPE&,const int,PARTICLE_LEVELSET<T_GRID>&); \
    template void Exchange_Overlapping_Block_Particles(const MPI_UNIFORM_GRID<T_GRID >&,const PARTICLE_LEVELSET_PARTICLES<T_GRID::VECTOR_T>&, \
        GRID_ARRAYS_POLICY<T_GRID >::ARRAYS_SCALAR::REBIND<PARTICLE_LEVELSET_PARTICLES<T_GRID::VECTOR_T>*>::TYPE&,const int,PARTICLE_LEVELSET<T_GRID>&); \
    template void Exchange_Boundary_Particles(const MPI_UNIFORM_GRID<T_GRID >&,const PARTICLE_LEVELSET_REMOVED_PARTICLES<T_GRID::VECTOR_T>&, \
        GRID_ARRAYS_POLICY<T_GRID >::ARRAYS_SCALAR::REBIND<PARTICLE_LEVELSET_REMOVED_PARTICLES<T_GRID::VECTOR_T>*>::TYPE&,const int,PARTICLE_LEVELSET<T_GRID>&); \
    template void Exchange_Boundary_Particles_Flat(const MPI_UNIFORM_GRID<T_GRID >&,VORTICITY_PARTICLES<T_GRID::VECTOR_T>& particles,const T_GRID::SCALAR ghost_distance); \
    template void Exchange_Ghost_Particles(const MPI_UNIFORM_GRID<T_GRID >&,const PARTICLE_LEVELSET_PARTICLES<T_GRID::VECTOR_T>&, \
        GRID_ARRAYS_POLICY<T_GRID >::ARRAYS_SCALAR::REBIND<PARTICLE_LEVELSET_PARTICLES<T_GRID::VECTOR_T>*>::TYPE&,const int,PARTICLE_LEVELSET<T_GRID>&); \
    template void Exchange_Ghost_Particles(const MPI_UNIFORM_GRID<T_GRID >&,const PARTICLE_LEVELSET_REMOVED_PARTICLES<T_GRID::VECTOR_T>&, \
        GRID_ARRAYS_POLICY<T_GRID >::ARRAYS_SCALAR::REBIND<PARTICLE_LEVELSET_REMOVED_PARTICLES<T_GRID::VECTOR_T>*>::TYPE&,const int,PARTICLE_LEVELSET<T_GRID>&);
#define P(...) __VA_ARGS__
INSTANTIATION_HELPER(float,P(GRID<VECTOR<float,1> >),1);
INSTANTIATION_HELPER(float,P(GRID<VECTOR<float,2> >),2);
INSTANTIATION_HELPER(float,P(GRID<VECTOR<float,3> >),3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,P(GRID<VECTOR<double,1> >),1);
INSTANTIATION_HELPER(double,P(GRID<VECTOR<double,2> >),2);
INSTANTIATION_HELPER(double,P(GRID<VECTOR<double,3> >),3);
#endif
}

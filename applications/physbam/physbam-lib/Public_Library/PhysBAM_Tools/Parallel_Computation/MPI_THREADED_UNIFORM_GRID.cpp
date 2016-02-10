//#####################################################################
// Copyright 2011, Linhai Qiu, Valeriya Nikolaenko
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Parallel_Computation/MPI_THREADED_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> MPI_THREADED_UNIFORM_GRID<T_GRID>::
MPI_THREADED_UNIFORM_GRID(T_GRID& local_grid_input,const int number_of_ghost_cells_input,const int number_of_threads_input,const int tid_input,ARRAY<THREAD_PACKAGE>* buffers_input,const bool skip_initialization,const TV_INT& mpi_processes_per_dimension,const TV_INT& threaded_processes_per_dimension,const TV_BOOL& periodic_input,MPI::Group* group_input)
    :MPI_GRID<T_GRID>(local_grid_input,number_of_ghost_cells_input,true,mpi_processes_per_dimension*threaded_processes_per_dimension,periodic_input,group_input),number_of_threads(number_of_threads_input),tid(tid_input),buffers(buffers_input),mpi_packages_buffers(0),mpi_requests_buffers(0),status(0)
{
#ifdef USE_MPI
#ifdef USE_PTHREADS
    // set up lock, barrier, and mpi_buffers that are common for all threads in one mpi process
    if(number_of_threads>=1){
        if(tid==1){
            lock=new pthread_mutex_t;
            pthread_mutex_init(lock,NULL);
            THREAD_PACKAGE pack(sizeof(pthread_mutex_t*));
            *(pthread_mutex_t**)(&pack.buffer(1))=lock;
            buffers->Append(pack);}
        else{
            while(buffers->m==0) PROCESS_UTILITIES::PB_Sleep(1);
            lock=*(pthread_mutex_t**)(&(*buffers)(1).buffer(1));}
        if(tid==1){
            barr=new pthread_barrier_t;
            pthread_barrier_init(barr,NULL,number_of_threads);
            THREAD_PACKAGE pack(sizeof(pthread_barrier_t*));
            *(pthread_barrier_t**)(&pack.buffer(1))=barr;
            buffers->Append(pack);}
        else{
            while(buffers->m<=1) PROCESS_UTILITIES::PB_Sleep(1);
            barr=*(pthread_barrier_t**)(&(*buffers)(2).buffer(1));}
        if(tid==1){
            mpi_packages_buffers=new ARRAY<MPI_PACKAGE>;
            THREAD_PACKAGE pack(sizeof(ARRAY<MPI_PACKAGE>*));
            *(ARRAY<MPI_PACKAGE>**)(&pack.buffer(1))=mpi_packages_buffers;
            buffers->Append(pack);}
        else{
            while(buffers->m<=2) PROCESS_UTILITIES::PB_Sleep(1);
            mpi_packages_buffers=*(ARRAY<MPI_PACKAGE>**)(&(*buffers)(3).buffer(1));}
        if(tid==1){
            mpi_requests_buffers=new ARRAY<MPI::Request>;
            THREAD_PACKAGE pack(sizeof(ARRAY<MPI::Request>*));
            *(ARRAY<MPI::Request>**)(&pack.buffer(1))=mpi_requests_buffers;
            buffers->Append(pack);}
        else{
            while(buffers->m<=3) PROCESS_UTILITIES::PB_Sleep(1);
            mpi_requests_buffers=*(ARRAY<MPI::Request>**)(&(*buffers)(4).buffer(1));}
        pthread_barrier_wait(barr);
        if(tid==1) buffers->m=0;
        pthread_barrier_wait(barr);}
    // Initialize
    global_grid=local_grid_input.Get_Regular_Grid();
    pthread_mutex_lock(lock);
    status=new MPI::Status;
    if(!group) number_of_mpi_processes=MPI::COMM_WORLD.Get_size();
    else number_of_mpi_processes=group->Get_size();
    pthread_mutex_unlock(lock);
    number_of_processes=number_of_mpi_processes*number_of_threads;
    // setup mpi_grid and threaded_grid: skip initialization in both
    mpi_grid=new MPI_UNIFORM_GRID<T_GRID>(local_grid,number_of_ghost_cells_input,true,mpi_processes_per_dimension,periodic_input,group_input);
    if(number_of_threads>=1) threaded_grid=new THREADED_UNIFORM_GRID<T_GRID>(*buffers_input,tid_input,number_of_threads,local_grid,number_of_ghost_cells_input,true,threaded_processes_per_dimension,periodic_input);
    // setup topology of mpi_grid
    Initialize_MPI_Grid(mpi_processes_per_dimension);
    // setup topology of threaded_grid
    Initialize_Threaded_Grid(threaded_processes_per_dimension);
    // setup topology of the whole mpi_threaded_grid
    assert(&mpi_grid->local_grid==&threaded_grid->local_grid && &local_grid==&threaded_grid->local_grid);
    Initialize_MPI_Threaded_Grid();
#endif
#endif
}
//#####################################################################
// Function Initialize_MPI_Grid
//#####################################################################
template<class T_GRID> void MPI_THREADED_UNIFORM_GRID<T_GRID>::
Initialize_MPI_Grid(const TV_INT& mpi_processes_per_dimension)
{
#ifdef USE_MPI
#ifdef USE_PTHREADS
    // extract global grid and divide among processes
    mpi_grid->global_grid=global_grid;
    mpi_grid->number_of_processes=number_of_mpi_processes;
    mpi_grid->Split_Grid(mpi_processes_per_dimension);
    // setup communicator and topology
    Initialize_MPI_Communicator(true,group);
    if(tid==1){
        assert(!buffers->m);
        THREAD_PACKAGE pack(sizeof(MPI::Intracomm*));
        *(MPI::Intracomm**)(&pack.buffer(1))=mpi_grid->comm;
        buffers->Append(pack);}
    else{
        while(buffers->m==0) PROCESS_UTILITIES::PB_Sleep(1);
        mpi_grid->comm=*(MPI::Intracomm**)(&(*buffers)(1).buffer(1));
        pthread_mutex_lock(lock);
        mpi_grid->rank=mpi_grid->comm->Get_rank();
        pthread_mutex_unlock(lock);
        mpi_grid->coordinates=mpi_grid->all_coordinates(mpi_grid->rank+1);}
    pthread_barrier_wait(barr);
    if(tid==1) buffers->m=0;
    pthread_barrier_wait(barr);

    mpi_grid->side_neighbor_ranks.Resize(T_GRID::number_of_neighbors_per_node);
    mpi_grid->side_neighbor_directions.Resize(T_GRID::number_of_neighbors_per_node);
    mpi_grid->all_neighbor_ranks.Resize(T_GRID::number_of_one_ring_neighbors_per_cell);
    mpi_grid->all_neighbor_directions.Resize(T_GRID::number_of_one_ring_neighbors_per_cell);
    for(int n=1;n<=T_GRID::number_of_neighbors_per_node;n++){
        mpi_grid->side_neighbor_ranks(n)=mpi_grid->process_ranks(T_GRID::Node_Neighbor(mpi_grid->coordinates,n));
        mpi_grid->side_neighbor_directions(n)=T_GRID::Node_Neighbor(TV_INT(),n);}
    for(int n=1;n<=T_GRID::number_of_one_ring_neighbors_per_cell;n++){
        mpi_grid->all_neighbor_ranks(n)=mpi_grid->process_ranks(T_GRID::One_Ring_Neighbor(mpi_grid->coordinates,n));
        mpi_grid->all_neighbor_directions(n)=T_GRID::One_Ring_Neighbor(TV_INT(),n);}
    if(!group){
        if(tid==1){
            group=new MPI::Group(mpi_grid->comm->Get_group());
            mpi_grid->group=group;
            THREAD_PACKAGE pack(sizeof(MPI::Group*));
            *(MPI::Group**)(&pack.buffer(1))=group;
            buffers->Append(pack);}
        else{
            while(buffers->m==0) PROCESS_UTILITIES::PB_Sleep(1);
            group=*(MPI::Group**)(&(*buffers)(1).buffer(1));
            mpi_grid->group=group;}
        pthread_barrier_wait(barr);
        if(tid==1) buffers->m=0;
        pthread_barrier_wait(barr);}
    
    local_grid=mpi_grid->Restrict_Grid(mpi_grid->coordinates);
    // initialize offset
    TV_INT start_index;TV_INT end_index;
    for(int axis=1;axis<=T_GRID::dimension;axis++)
        start_index[axis]=mpi_grid->boundaries(axis)(mpi_grid->coordinates[axis]);
    mpi_grid->local_to_global_offset=start_index-TV_INT::All_Ones_Vector();
#endif
#endif
}
//#####################################################################
// Function Initialize_MPI_Communicator
//#####################################################################
template<class T_GRID> void MPI_THREADED_UNIFORM_GRID<T_GRID>::
Initialize_MPI_Communicator(const bool manual,MPI::Group* group)
{
#ifdef USE_MPI
    mpi_grid->process_ranks.Resize(mpi_grid->process_grid.Domain_Indices(1));
    ARRAYS_COMPUTATIONS::Fill(mpi_grid->process_ranks.array,MPI::PROC_NULL);
    TV_INT extents=mpi_grid->process_grid.Domain_Indices().Maximum_Corner();
    if(tid==1) mpi_grid->comm=new MPI::Intracomm;
    // setup communicator manually to guarantee reasonable topology for SMP clusters
    // sort axes in decreasing order of how much we have to communicate along them
    ARRAY<int> axes(T_GRID::dimension);ARRAY<T> axis_lengths(T_GRID::dimension);
    for(int axis=1;axis<=T_GRID::dimension;axis++){axes(axis)=axis;axis_lengths(axis)=(T)mpi_grid->global_grid.Domain_Indices().Maximum_Corner()[axis]/extents[axis];}
    Sort(axes,Indirect_Comparison(axis_lengths));
    // lay out process ranks on grid
    Fill_Process_Ranks(mpi_grid->process_grid,mpi_grid->process_ranks,axes);
    // fill in ghost process_ranks for periodic domains
    if(mpi_grid->periodic!=TV_BOOL()) for(NODE_ITERATOR iterator(mpi_grid->process_grid,1,T_GRID::GHOST_REGION);iterator.Valid();iterator.Next()){
        TV_INT node=iterator.Node_Index(),wrapped_node=node;
        for(int axis=1;axis<=T_GRID::dimension;axis++) if(mpi_grid->periodic[axis]) wrapped_node[axis]=(node[axis]+mpi_grid->process_grid.Counts()[axis]-1)%mpi_grid->process_grid.Counts()[axis]+1;
        mpi_grid->process_ranks(node)=mpi_grid->process_ranks(wrapped_node);}
    // allocate communicator
    if(tid==1) *mpi_grid->comm=group?MPI::COMM_WORLD.Create(*group):MPI::COMM_WORLD.Dup();
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    LOG::filecout("After Create comm\n");
#endif
    mpi_grid->all_coordinates.Resize(mpi_grid->number_of_processes);
    for(NODE_ITERATOR iterator(mpi_grid->process_grid);iterator.Valid();iterator.Next())
        mpi_grid->all_coordinates(mpi_grid->process_ranks(iterator.Node_Index())+1)=iterator.Node_Index();
    if(tid==1) mpi_grid->rank=mpi_grid->comm->Get_rank();
    if(tid==1) mpi_grid->coordinates=mpi_grid->all_coordinates(mpi_grid->rank+1);
#endif
}
#ifdef USE_MPI
template<class T> static void Fill_Process_Ranks(GRID<VECTOR<T,1> >& process_grid,ARRAY<int,VECTOR<int,1> >& process_ranks,ARRAY<int>& axes,bool fill_threaded_ranks=false)
{
    VECTOR<int,1> extents=process_grid.Domain_Indices().Maximum_Corner();
    for(int i=1;i<=extents.x;i++)process_ranks(i)=i-1;
}
template<class T> static void Fill_Process_Ranks(GRID<VECTOR<T,2> >& process_grid,ARRAY<int,VECTOR<int,2> >& process_ranks,ARRAY<int>& axes,bool fill_threaded_ranks=false)
{
    int null_process=fill_threaded_ranks?-1:MPI::PROC_NULL;
    int next_rank=0;
    VECTOR<int,2> extents=process_grid.Domain_Indices().Maximum_Corner(),half_extents=extents/2;
    for(int i=1;i<=half_extents[axes(1)];i++)for(int j=1;j<=half_extents[axes(2)];j++)
        for(int ii=0;ii<2;ii++)for(int jj=0;jj<2;jj++){
            VECTOR<int,2> permuted_index(2*i+ii-1,2*j+jj-1),index;
            for(int a=1;a<=2;a++)index[axes(a)]=permuted_index[a];
            process_ranks(index)=next_rank++;}
    for(int i=1;i<=extents.x;i++)for(int j=1;j<=extents.y;j++)if(process_ranks(i,j)==null_process) process_ranks(i,j)=next_rank++;
}
template<class T> static void Fill_Process_Ranks(GRID<VECTOR<T,3> >& process_grid,ARRAY<int,VECTOR<int,3> >& process_ranks,ARRAY<int>& axes,bool fill_threaded_ranks=false)
{
    int null_process=fill_threaded_ranks?-1:MPI::PROC_NULL;
    int next_rank=0;
    VECTOR<int,3> extents=process_grid.Domain_Indices().Maximum_Corner(),half_extents=extents/2;
    for(int i=1;i<=half_extents[axes(1)];i++)for(int j=1;j<=half_extents[axes(2)];j++)for(int ij=1;ij<=half_extents[axes(3)];ij++)
        for(int ii=0;ii<2;ii++)for(int jj=0;jj<2;jj++)for(int ijij=0;ijij<2;ijij++){
            VECTOR<int,3> permuted_index(2*i+ii-1,2*j+jj-1,2*ij+ijij-1),index;
            for(int a=1;a<=3;a++)index[axes(a)]=permuted_index[a];
            process_ranks(index)=next_rank++;}
    for(int i=1;i<=extents.x;i++)for(int j=1;j<=extents.y;j++)for(int ij=1;ij<=extents.z;ij++)if(process_ranks(i,j,ij)==null_process) process_ranks(i,j,ij)=next_rank++;
}
#endif
//#####################################################################
// Function Initialize_Threaded_Grid
//#####################################################################
template<class T_GRID> void MPI_THREADED_UNIFORM_GRID<T_GRID>::
Initialize_Threaded_Grid(TV_INT threaded_processes_per_dimension)
{
#ifdef USE_PTHREADS
#ifdef USE_MPI
    // set up topology of threaded_grid
    threaded_grid->lock=lock;threaded_grid->barr=barr;
    threaded_grid->number_of_processes=number_of_threads;
    TV_INT first_mpi_coordinate;first_mpi_coordinate.Fill(1);
    T_GRID global_grid_first_mpi=(mpi_grid->Restrict_Grid(first_mpi_coordinate)).Get_Regular_Grid();
    if(global_grid_first_mpi.counts!=(mpi_grid->local_grid).Get_Regular_Grid().counts){
        threaded_grid->global_grid=global_grid_first_mpi;
        threaded_grid->Split_Grid(threaded_processes_per_dimension);
        for(int axis=1;axis<=TV::dimension;axis++) threaded_processes_per_dimension(axis)=threaded_grid->boundaries(axis).Size()-1;
        threaded_grid->global_grid=(mpi_grid->local_grid).Get_Regular_Grid();
        threaded_grid->Split_Grid(threaded_processes_per_dimension);}
    else{
        threaded_grid->global_grid=(mpi_grid->local_grid).Get_Regular_Grid();
        threaded_grid->Split_Grid(threaded_processes_per_dimension);} 
    threaded_grid->process_ranks.Resize(threaded_grid->process_grid.Domain_Indices(1));
    ARRAYS_COMPUTATIONS::Fill(threaded_grid->process_ranks.array,-1);
    TV_INT extents=threaded_grid->process_grid.Domain_Indices().Maximum_Corner();
    // sort axes global grid of mpi_threaded_grid
    ARRAY<int> axes(T_GRID::dimension);ARRAY<T> axis_lengths(T_GRID::dimension);
    for(int axis=1;axis<=T_GRID::dimension;axis++){axes(axis)=axis;axis_lengths(axis)=(T)global_grid_first_mpi.Domain_Indices().Maximum_Corner()[axis]/extents[axis];}
    Sort(axes,Indirect_Comparison(axis_lengths));
    // lay out process ranks on grid
    Fill_Process_Ranks(threaded_grid->process_grid,threaded_grid->process_ranks,axes,true);
    threaded_grid->all_coordinates.Resize(number_of_threads);
    for(NODE_ITERATOR iterator(threaded_grid->process_grid);iterator.Valid();iterator.Next())
        threaded_grid->all_coordinates(threaded_grid->process_ranks(iterator.Node_Index())+1)=iterator.Node_Index();
    threaded_grid->coordinates=threaded_grid->all_coordinates(tid);
    local_grid=threaded_grid->Restrict_Grid(threaded_grid->coordinates);
    // determine neighbor threads
    threaded_grid->side_neighbor_ranks.Resize(T_GRID::number_of_neighbors_per_node);
    threaded_grid->side_neighbor_directions.Resize(T_GRID::number_of_neighbors_per_node);
    threaded_grid->all_neighbor_ranks.Resize(T_GRID::number_of_one_ring_neighbors_per_cell);
    threaded_grid->all_neighbor_directions.Resize(T_GRID::number_of_one_ring_neighbors_per_cell);
    for(int n=1;n<=T_GRID::number_of_neighbors_per_node;n++){
        threaded_grid->side_neighbor_ranks(n)=threaded_grid->process_ranks(T_GRID::Node_Neighbor(threaded_grid->coordinates,n));
        threaded_grid->side_neighbor_directions(n)=T_GRID::Node_Neighbor(TV_INT(),n);}
    for(int n=1;n<=T_GRID::number_of_one_ring_neighbors_per_cell;n++){
        threaded_grid->all_neighbor_ranks(n)=threaded_grid->process_ranks(T_GRID::One_Ring_Neighbor(threaded_grid->coordinates,n));
        threaded_grid->all_neighbor_directions(n)=T_GRID::One_Ring_Neighbor(TV_INT(),n);}
#endif
#endif
}
//#####################################################################
// Initialize_MPI_Threaded_Grid
//#####################################################################
template<class T_GRID> void MPI_THREADED_UNIFORM_GRID<T_GRID>::
Initialize_MPI_Threaded_Grid()
{
#ifdef USE_MPI
#ifdef USE_PTHREADS
    all_neighbor_directions=mpi_grid->all_neighbor_directions;
    side_neighbor_directions=mpi_grid->side_neighbor_directions;
    // modify neighbor ranks of mpi grid according to the actual mpi grid in each direction due to the existence of multiple threaded grids
    for(int n=1;n<=T_GRID::number_of_neighbors_per_node;n++)
        if(threaded_grid->side_neighbor_ranks(n)>=0)
            mpi_grid->side_neighbor_ranks(n)=MPI::PROC_NULL;
    for(int n=1;n<=T_GRID::number_of_one_ring_neighbors_per_cell;n++){
        if(threaded_grid->all_neighbor_ranks(n)>=0)
            mpi_grid->all_neighbor_ranks(n)=MPI::PROC_NULL;
        if(threaded_grid->all_neighbor_ranks(n)<0 && all_neighbor_directions(n).L1_Norm()>1){
            TV_INT current_direction=all_neighbor_directions(n);
            for(int axis=1;axis<=TV::dimension;axis++)
                if(current_direction(axis)!=0){
                    TV_INT axis_vector;axis_vector(axis)=current_direction(axis);
                    if(threaded_grid->process_ranks(threaded_grid->coordinates+axis_vector)>=0) current_direction(axis)=0;}
            mpi_grid->all_neighbor_ranks(n)=mpi_grid->process_ranks(mpi_grid->coordinates+current_direction);}}
    //initialize process ranks of mpi/threaded grid with rank=number_of_threads*mpi_rank+tid-1
    process_grid=T_GRID(mpi_grid->process_grid.counts*threaded_grid->process_grid.counts,RANGE<TV>::Centered_Box());
    process_ranks.Resize(process_grid.Domain_Indices(1));ARRAYS_COMPUTATIONS::Fill(process_ranks.array,-1);
    for(int mpi_rank=0;mpi_rank<=mpi_grid->number_of_processes-1;mpi_rank++)
        for(int threaded_rank=0;threaded_rank<=threaded_grid->number_of_processes-1;threaded_rank++){
            TV_INT global_coordinates;
            for(int axis=1;axis<=TV::dimension;axis++)
                global_coordinates(axis)=threaded_grid->process_grid.counts(axis)*(mpi_grid->all_coordinates(mpi_rank+1)(axis)-1)+threaded_grid->all_coordinates(threaded_rank+1)(axis);
            process_ranks(global_coordinates)=mpi_rank*number_of_threads+threaded_rank;}
    all_coordinates.Resize(number_of_processes);
    for(NODE_ITERATOR iterator(process_grid);iterator.Valid();iterator.Next())
        all_coordinates(process_ranks(iterator.Node_Index())+1)=iterator.Node_Index();
    rank=mpi_grid->rank*number_of_threads+tid-1;
    coordinates=all_coordinates(rank+1);
    //initialize neighbor ranks according to the new mpi/threaded grid
    side_neighbor_ranks.Resize(T_GRID::number_of_neighbors_per_node);
    all_neighbor_ranks.Resize(T_GRID::number_of_one_ring_neighbors_per_cell);
    for(int n=1;n<=T_GRID::number_of_neighbors_per_node;n++)
        side_neighbor_ranks(n)=process_ranks(T_GRID::Node_Neighbor(coordinates,n));
    for(int n=1;n<=T_GRID::number_of_one_ring_neighbors_per_cell;n++)
        all_neighbor_ranks(n)=process_ranks(T_GRID::One_Ring_Neighbor(coordinates,n));
#endif
#endif
}
//#####################################################################
// Get_Send_Tag
//#####################################################################
template<class T_GRID> int MPI_THREADED_UNIFORM_GRID<T_GRID>::
Get_Send_Tag(TV_INT direction) const
{
    STATIC_ASSERT(TV_INT::m<=3);int tag=0;
    int recv_rank=process_ranks(coordinates+direction);
    tag=rank*number_of_processes+recv_rank;
    return tag;
}
//#####################################################################
// Get_Recv_Tag
//#####################################################################
template<class T_GRID> int MPI_THREADED_UNIFORM_GRID<T_GRID>::
Get_Recv_Tag(TV_INT direction) const
{
    STATIC_ASSERT(TV_INT::m<=3);int tag=0;
    int send_rank=process_ranks(coordinates+direction);
    tag=send_rank*number_of_processes+rank;
    return tag;
}
//#####################################################################
// Exchange_Boundary_Cell_Data
//#####################################################################
template<class T_GRID> template<class T2> void MPI_THREADED_UNIFORM_GRID<T_GRID>::
Exchange_Boundary_Cell_Data(ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> >& data,const int bandwidth,const bool include_corners) const
{
#ifdef USE_MPI
#ifdef USE_PTHREADS
    if(threaded_grid) threaded_grid->Exchange_Boundary_Cell_Data(data,bandwidth,include_corners);
    if(threaded_grid && mpi_grid){
        RANGE<TV_INT> sentinels=RANGE<TV_INT>::Zero_Box();
        const ARRAY<int>& neighbor_ranks=include_corners?mpi_grid->all_neighbor_ranks:mpi_grid->side_neighbor_ranks;
        const ARRAY<TV_INT>& neighbor_directions=include_corners?all_neighbor_directions:side_neighbor_directions;
        // send
        ARRAY<RANGE<TV_INT> > send_regions;
        Find_Boundary_Regions(send_regions,sentinels,false,RANGE<VECTOR<int,1> >(0,bandwidth-1),include_corners,true,mpi_grid->local_grid);
        for(int n=1;n<=send_regions.m;n++)if(neighbor_ranks(n)!=MPI::PROC_NULL){
            pthread_mutex_lock(lock);
            MPI_PACKAGE package=mpi_grid->Package_Cell_Data(data,send_regions(n));
            mpi_packages_buffers->Append(package);
            mpi_requests_buffers->Append(package.Isend(*mpi_grid->comm,neighbor_ranks(n),Get_Send_Tag(neighbor_directions(n))));
            pthread_mutex_unlock(lock);}
        // receive
        ARRAY<RANGE<TV_INT> > recv_regions;
        Find_Boundary_Regions(recv_regions,sentinels,false,RANGE<VECTOR<int,1> >(-bandwidth,-1),include_corners,true,mpi_grid->local_grid);
        for(int n=1;n<=recv_regions.m;n++)if(neighbor_ranks(n)!=MPI::PROC_NULL){
            pthread_mutex_lock(lock);
            MPI_PACKAGE package=mpi_grid->Package_Cell_Data(data,recv_regions(n));
            mpi_packages_buffers->Append(package);
            mpi_requests_buffers->Append(package.Irecv(*mpi_grid->comm,neighbor_ranks(n),Get_Recv_Tag(neighbor_directions(n))));
            pthread_mutex_unlock(lock);}
        pthread_barrier_wait(barr);
        if(tid==1){MPI_UTILITIES::Wait_All(*mpi_requests_buffers);MPI_PACKAGE::Free_All(*mpi_packages_buffers);}
        pthread_barrier_wait(barr);}
#endif
#endif
}
//#####################################################################
// Exchange_Boundary_Face_Data
//#####################################################################
template<class T_GRID> template<class T2> void MPI_THREADED_UNIFORM_GRID<T_GRID>::
Exchange_Boundary_Face_Data(ARRAY<T2,FACE_INDEX<TV::dimension> >& data,const int bandwidth) const
{
#ifdef USE_MPI
#ifdef USE_PTHREADS
    if(threaded_grid) threaded_grid->Exchange_Boundary_Face_Data(data,bandwidth);
    if(mpi_grid && threaded_grid){
        RANGE<VECTOR<int,1> > boundary_band(0,bandwidth-1),ghost_band(-bandwidth,-1);
        // send
        ARRAY<ARRAY<RANGE<TV_INT> > > send_regions(GRID_T::dimension);
        for(int axis=1;axis<=GRID_T::dimension;axis++)Find_Boundary_Regions(send_regions(axis),mpi_grid->Face_Sentinels(axis),true,boundary_band,true,true,mpi_grid->local_grid);
        for(int n=1;n<=send_regions(1).m;n++)if(mpi_grid->all_neighbor_ranks(n)!=MPI::PROC_NULL){
            ARRAY<RANGE<TV_INT> > send_regions_n(GRID_T::dimension);for(int axis=1;axis<=GRID_T::dimension;axis++)send_regions_n(axis)=send_regions(axis)(n);
            pthread_mutex_lock(lock);
            MPI_PACKAGE package=mpi_grid->Package_Face_Data(data,send_regions_n);
            mpi_packages_buffers->Append(package);
            mpi_requests_buffers->Append(package.Isend(*mpi_grid->comm,mpi_grid->all_neighbor_ranks(n),Get_Send_Tag(all_neighbor_directions(n))));
            pthread_mutex_unlock(lock);}
        // receive
        ARRAY<ARRAY<RANGE<TV_INT> > > recv_regions(GRID_T::dimension);
        for(int axis=1;axis<=GRID_T::dimension;axis++)Find_Boundary_Regions(recv_regions(axis),mpi_grid->Face_Sentinels(axis),true,ghost_band,true,true,mpi_grid->local_grid);
        for(int n=1;n<=recv_regions(1).m;n++)if(mpi_grid->all_neighbor_ranks(n)!=MPI::PROC_NULL){
            ARRAY<RANGE<TV_INT> > recv_regions_n(GRID_T::dimension);for(int axis=1;axis<=GRID_T::dimension;axis++)recv_regions_n(axis)=recv_regions(axis)(n);
            pthread_mutex_lock(lock);
            MPI_PACKAGE package=mpi_grid->Package_Face_Data(data,recv_regions_n);
            mpi_packages_buffers->Append(package);
            mpi_requests_buffers->Append(package.Irecv(*mpi_grid->comm,mpi_grid->all_neighbor_ranks(n),Get_Recv_Tag(all_neighbor_directions(n))));
            pthread_mutex_unlock(lock);}
        // finish
        pthread_barrier_wait(barr);
        if(tid==1){MPI_UTILITIES::Wait_All(*mpi_requests_buffers);MPI_PACKAGE::Free_All(*mpi_packages_buffers);}
        pthread_barrier_wait(barr);}
#endif
#endif
}
//#####################################################################
// Union_Common_Face_Data
//#####################################################################
template<class T_GRID> template<class T_ARRAYS> void MPI_THREADED_UNIFORM_GRID<T_GRID>::
Union_Common_Face_Data(T_ARRAYS& data) const
{
    if(mpi_grid) mpi_grid->Union_Common_Face_Data(data);
}
//#####################################################################
// Average_Common_Face_Data
//#####################################################################
template<class T_GRID> template<class T_FACE_ARRAYS_2> void MPI_THREADED_UNIFORM_GRID<T_GRID>::
Average_Common_Face_Data(T_FACE_ARRAYS_2& data) const
{
#ifdef USE_MPI
#ifdef USE_PTHREADS
    if(threaded_grid) threaded_grid->Average_Common_Face_Data(data);
    if(mpi_grid && threaded_grid){
        ARRAY<ARRAY<RANGE<TV_INT> > > regions(T_GRID::dimension);
        for(int axis=1;axis<=T_GRID::dimension;axis++)Find_Boundary_Regions(regions(axis),mpi_grid->Parallel_Face_Sentinels(axis),false,RANGE<VECTOR<int,1> >(0,0),false);
        ARRAY<ARRAY<T> > buffers(regions(1).m);
        pthread_mutex_lock(lock);
        ARRAY<MPI_PACKAGE> packages(regions(1).m);
        MPI::Datatype T_type=MPI_UTILITIES::Datatype<T>();
        pthread_mutex_unlock(lock);
        // send and receive into temporary buffers
        for(int n=1;n<=regions(1).m;n++)if(mpi_grid->side_neighbor_ranks(n)!=MPI::PROC_NULL){
            int axis=(n-1)/2+1;
            pthread_mutex_lock(lock);
            packages(n)=mpi_grid->Package_Common_Face_Data(data,axis,regions(axis)(n));
            mpi_requests_buffers->Append(packages(n).Isend(*mpi_grid->comm,mpi_grid->side_neighbor_ranks(n),Get_Send_Tag(side_neighbor_directions(n))));
            buffers(n).Resize(packages(n).Size()/sizeof(T));
            mpi_requests_buffers->Append(mpi_grid->comm->Irecv(buffers(n).Get_Array_Pointer(),buffers(n).m,T_type,mpi_grid->side_neighbor_ranks(n),Get_Recv_Tag(side_neighbor_directions(n))));
            pthread_mutex_unlock(lock);}
        // wait
        pthread_barrier_wait(barr);
        if(tid==1) MPI_UTILITIES::Wait_All(*mpi_requests_buffers);
        pthread_barrier_wait(barr);
        // average received data with local data (TODO: find a cleaner general way to do this)
        for(int n=1;n<=regions(1).m;n++)if(mpi_grid->side_neighbor_ranks(n)!=MPI::PROC_NULL){
            pthread_mutex_lock(lock);
            ARRAY<char> pack_buffer(packages(n).Pack_Size(*mpi_grid->comm));packages(n).Pack(pack_buffer,*mpi_grid->comm);
            ARRAY<T> local_buffer(buffers(n).m);
            pthread_mutex_unlock(lock);
            int position=0;
            pthread_mutex_lock(lock);
            T_type.Unpack(pack_buffer.Get_Array_Pointer(),pack_buffer.m,local_buffer.Get_Array_Pointer(),local_buffer.m,position,*mpi_grid->comm);
            pthread_mutex_unlock(lock);
            ARRAY<T>::Copy((T).5,buffers(n),(T).5,local_buffer,local_buffer); // average
            position=0;
            pthread_mutex_lock(lock);
            T_type.Pack(local_buffer.Get_Array_Pointer(),local_buffer.m,pack_buffer.Get_Array_Pointer(),pack_buffer.m,position,*mpi_grid->comm);
            packages(n).Unpack(pack_buffer,*mpi_grid->comm);
            pthread_mutex_unlock(lock);}
        pthread_mutex_lock(lock);
        MPI_PACKAGE::Free_All(packages);
        pthread_mutex_unlock(lock);}
#endif
#endif
}
//#####################################################################
// Sum_Common_Face_Data
//#####################################################################
template<class T_GRID> template<class T_FACE_ARRAYS_2> void MPI_THREADED_UNIFORM_GRID<T_GRID>::
Sum_Common_Face_Data(T_FACE_ARRAYS_2& data) const
{
    if(threaded_grid && !mpi_grid) PHYSBAM_FATAL_ERROR("Sum_Common_Face_Data not implemented for threading!");
    if(mpi_grid) mpi_grid->Sum_Common_Face_Data(data);
}
//#####################################################################
// Copy_Common_Face_Data
//#####################################################################
template<class T_GRID> template<class T_FACE_ARRAYS_2> void MPI_THREADED_UNIFORM_GRID<T_GRID>::
Copy_Common_Face_Data(T_FACE_ARRAYS_2& data) const
{
    if(mpi_grid) mpi_grid->Copy_Common_Face_Data(data);
}
//#####################################################################
// Assert_Common_Face_Data
//#####################################################################
template<class T_GRID> template<class T_FACE_ARRAYS_2> void MPI_THREADED_UNIFORM_GRID<T_GRID>::
Assert_Common_Face_Data(T_FACE_ARRAYS_2& data,const T tolerance) const
{
    if(threaded_grid) threaded_grid->Assert_Common_Face_Data(data,tolerance);
    if(mpi_grid) mpi_grid->Assert_Common_Face_Data(data,tolerance);
}
//#####################################################################
// Synchronize_Dt
//#####################################################################
template<class T_GRID> void MPI_THREADED_UNIFORM_GRID<T_GRID>::
Synchronize_Dt(T& dt) const
{
#ifdef USE_MPI
#ifdef USE_PTHREADS
    if(threaded_grid && mpi_grid){
        threaded_grid->Synchronize_Dt(dt);
        if(tid==1){
            mpi_grid->Synchronize_Dt(dt);
            THREAD_PACKAGE pack(sizeof(T));*(T*)(&pack.buffer(1))=dt;
            buffers->Append(pack);}
        pthread_barrier_wait(barr);
        dt=*(T*)(&(*buffers)(1).buffer(1));
        pthread_barrier_wait(barr);
        if(tid==1) buffers->m=0;
        pthread_barrier_wait(barr);}
    else{
        if(threaded_grid) threaded_grid->Synchronize_Dt(dt);
        if(mpi_grid) mpi_grid->Synchronize_Dt(dt);}
#endif
#endif
}
//#####################################################################
// Synchronize_Ghost_Cells
//#####################################################################
template<class T_GRID> void MPI_THREADED_UNIFORM_GRID<T_GRID>::
Synchronize_Ghost_Cells(int& ghost_cells) const
{
#ifdef USE_MPI
#ifdef USE_PTHREADS
    if(threaded_grid && mpi_grid){
        threaded_grid->Synchronize_Ghost_Cells(ghost_cells);
        if(tid==1){
            mpi_grid->Synchronize_Ghost_Cells(ghost_cells);
            THREAD_PACKAGE pack(sizeof(int));*(int*)(&pack.buffer(1))=ghost_cells;
            buffers->Append(pack);}
        pthread_barrier_wait(barr);
        ghost_cells=*(int*)(&(*buffers)(1).buffer(1));
        pthread_barrier_wait(barr);
        if(tid==1) buffers->m=0;
        pthread_barrier_wait(barr);}
    else{
        if(threaded_grid) threaded_grid->Synchronize_Ghost_Cells(ghost_cells);
        if(mpi_grid) mpi_grid->Synchronize_Ghost_Cells(ghost_cells);}
#endif
#endif
}
//#####################################################################
// Initialize
//#####################################################################
template<class T_GRID> void MPI_THREADED_UNIFORM_GRID<T_GRID>::
Initialize(VECTOR<VECTOR<bool,2>,T_GRID::dimension>& domain_walls)
{
    if(mpi_grid) mpi_grid->Initialize(domain_walls);
    if(threaded_grid) threaded_grid->Initialize(domain_walls);
}
//#####################################################################
// Neighbor
//#####################################################################
template<class T_GRID> bool MPI_THREADED_UNIFORM_GRID<T_GRID>::
Neighbor(const int axis,const int axis_side) const
{
    int side=2*(axis-1)+axis_side;
    return side_neighbor_ranks(side)!=-1;
}
//#####################################################################
// Allgather
//#####################################################################
template<class T_GRID> void MPI_THREADED_UNIFORM_GRID<T_GRID>::
Allgather(int sent_data,ARRAY<int>& data)
{
#ifdef USE_MPI
#ifdef USE_PTHREADS
    ARRAY<int> gathered_data_threaded(threaded_grid->number_of_processes);
    gathered_data_threaded(tid)=sent_data;
    threaded_grid->Allgather(gathered_data_threaded);
    if(tid==1){
        mpi_grid->comm->Allgather(gathered_data_threaded.base_pointer,threaded_grid->number_of_processes,MPI_UTILITIES::Datatype<int>(),data.base_pointer,threaded_grid->number_of_processes,MPI_UTILITIES::Datatype<int>());
        THREAD_PACKAGE pack(sizeof(ARRAY<int>*));*(ARRAY<int>**)(&pack.buffer(1))=&data;
        buffers->Append(pack);}
    pthread_barrier_wait(barr);
    if(tid!=1){
        ARRAY<int>* data_tmp=*(ARRAY<int>**)(&(*buffers)(1).buffer(1));
        data=*data_tmp;}
    pthread_barrier_wait(barr);
    if(tid==1) buffers->m=0;
    pthread_barrier_wait(barr);
#endif
#endif
}
//#####################################################################
#define P(...) __VA_ARGS__
#define INSTANTIATION_HELPER_LENGTH(T_GRID,length)                      \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Exchange_Boundary_Cell_Data(T_ARRAYS_BASE::REBIND<VECTOR<T_GRID::SCALAR,length> >::TYPE&,const int,const bool) const;
#define INSTANTIATION_HELPER(T,T_GRID,d)                                 \
    template class MPI_THREADED_UNIFORM_GRID<T_GRID >; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Exchange_Boundary_Cell_Data(T_ARRAYS_BASE&,const int,const bool) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Exchange_Boundary_Cell_Data(T_ARRAYS_BASE::REBIND<int>::TYPE&,const int,const bool) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Exchange_Boundary_Cell_Data(T_ARRAYS_BASE::REBIND<bool>::TYPE&,const int,const bool) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Exchange_Boundary_Cell_Data(T_ARRAYS_BASE::REBIND<MATRIX<T_GRID::SCALAR,1> >::TYPE&,const int,const bool) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Exchange_Boundary_Cell_Data(T_ARRAYS_BASE::REBIND<SYMMETRIC_MATRIX<T,2> >::TYPE&,const int,const bool) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Exchange_Boundary_Cell_Data(T_ARRAYS_BASE::REBIND<SYMMETRIC_MATRIX<T,3> >::TYPE&,const int,const bool) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Exchange_Boundary_Face_Data(T_FACE_ARRAYS&,const int) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Exchange_Boundary_Face_Data(T_FACE_ARRAYS::REBIND<bool>::TYPE&,const int) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Exchange_Boundary_Face_Data(T_FACE_ARRAYS::REBIND<VECTOR<T,d> >::TYPE&,const int) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Exchange_Boundary_Face_Data(T_FACE_ARRAYS::REBIND<VECTOR<T,d+2> >::TYPE&,const int) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Exchange_Boundary_Face_Data(T_FACE_ARRAYS::REBIND<MATRIX<T,1> >::TYPE&,const int) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Exchange_Boundary_Face_Data(T_FACE_ARRAYS::REBIND<SYMMETRIC_MATRIX<T,2> >::TYPE&,const int) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Exchange_Boundary_Face_Data(T_FACE_ARRAYS::REBIND<SYMMETRIC_MATRIX<T,3> >::TYPE&,const int) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Union_Common_Face_Data(T_FACE_ARRAYS::REBIND<bool>::TYPE&) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Average_Common_Face_Data(T_FACE_ARRAYS&) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Average_Common_Face_Data(T_FACE_ARRAYS::REBIND<bool>::TYPE&) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Average_Common_Face_Data(T_FACE_ARRAYS::REBIND<VECTOR<T,d> >::TYPE&) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Average_Common_Face_Data(T_FACE_ARRAYS::REBIND<VECTOR<T,d+2> >::TYPE&) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Average_Common_Face_Data(T_FACE_ARRAYS::REBIND<MATRIX<T,1> >::TYPE&) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Average_Common_Face_Data(T_FACE_ARRAYS::REBIND<SYMMETRIC_MATRIX<T,2> >::TYPE&) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Average_Common_Face_Data(T_FACE_ARRAYS::REBIND<SYMMETRIC_MATRIX<T,3> >::TYPE& ) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Sum_Common_Face_Data(T_FACE_ARRAYS&) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Copy_Common_Face_Data(T_FACE_ARRAYS&) const; \
    template void MPI_THREADED_UNIFORM_GRID<T_GRID >::Assert_Common_Face_Data(T_FACE_ARRAYS&,const T tolerance) const; \
    INSTANTIATION_HELPER_LENGTH(P(T_GRID),1);INSTANTIATION_HELPER_LENGTH(P(T_GRID),2); \
    INSTANTIATION_HELPER_LENGTH(P(T_GRID),3);INSTANTIATION_HELPER_LENGTH(P(T_GRID),4);INSTANTIATION_HELPER_LENGTH(P(T_GRID),5);

#define INSTANTIATION_HELPER_ALL(T) \
    INSTANTIATION_HELPER(T,P(GRID<VECTOR<T,1> >),1);         \
    INSTANTIATION_HELPER(T,P(GRID<VECTOR<T,2> >),2);         \
    INSTANTIATION_HELPER(T,P(GRID<VECTOR<T,3> >),3);

INSTANTIATION_HELPER_ALL(float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER_ALL(double);
#endif

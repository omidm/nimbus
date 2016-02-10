//#####################################################################
// Copyright 2011. Linhai Qiu, Valeriya Nikolaenko
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_THREADED_UNIFORM_GRID
//#####################################################################
#ifndef __MPI_THREADED_UNIFORM_GRID__
#define __MPI_THREADED_UNIFORM_GRID__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/THREADED_UNIFORM_GRID.h>

namespace PhysBAM{

template<class T_GRID> class MPI_UNIFORM_GRID;

template<class T_GRID>
class MPI_THREADED_UNIFORM_GRID:public MPI_GRID<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS;
    typedef typename TV::template REBIND<bool>::TYPE TV_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<RANGE<TV_INT> >::TYPE T_ARRAYS_BOX_INT;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;

public:
    typedef T_GRID GRID_T;

    typedef MPI_GRID<T_GRID> BASE;
    using BASE::number_of_processes;using BASE::rank;using BASE::process_grid;using BASE::process_ranks;using BASE::coordinates;using BASE::all_coordinates;using BASE::boundaries;
    using BASE::global_grid;using BASE::local_grid;using BASE::side_neighbor_ranks;using BASE::side_neighbor_directions;using BASE::all_neighbor_ranks;
    using BASE::all_neighbor_directions;using BASE::group;using BASE::Find_Boundary_Regions;using BASE::Restrict_Grid;using BASE::Split_Grid;

    MPI_UNIFORM_GRID<T_GRID>* mpi_grid;
    THREADED_UNIFORM_GRID<T_GRID>* threaded_grid;
    int number_of_mpi_processes;
    int number_of_threads;
    int tid;
    ARRAY<THREAD_PACKAGE>* buffers; //global memory
    ARRAY<MPI_PACKAGE>* mpi_packages_buffers;
    ARRAY<MPI::Request>* mpi_requests_buffers;
    MPI::Status* status;
#ifdef USE_PTHREADS
    pthread_mutex_t* lock;
    pthread_barrier_t* barr;
#endif
    MPI_THREADED_UNIFORM_GRID(T_GRID& local_grid_input,const int number_of_ghost_cells_input,const int number_of_threads,const int tid_input=1,ARRAY<THREAD_PACKAGE>* buffers_input=0,const bool skip_initialization=false,const TV_INT& mpi_processes_per_dimension=TV_INT(),const TV_INT& threaded_processes_per_dimension=TV_INT(),const TV_BOOL& periodic_input=TV_BOOL(),MPI::Group* group_input=0);

//##############################################################
    void Initialize_MPI_Grid(const TV_INT& mpi_processes_per_dimension);
    void Initialize_MPI_Communicator(const bool manual,MPI::Group* group);
    void Initialize_Threaded_Grid(TV_INT threaded_processes_per_dimension);
    void Initialize_MPI_Threaded_Grid();
    int Get_Send_Tag(TV_INT direction) const;
    int Get_Recv_Tag(TV_INT direction) const;
    template<class T2> void Exchange_Boundary_Cell_Data(ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> >& data,const int bandwidth,const bool include_corners=true) const;
    template<class T2> void Exchange_Boundary_Face_Data(ARRAY<T2,FACE_INDEX<TV::dimension> >& data,const int bandwidth) const;
    template<class T_ARRAYS> void Union_Common_Face_Data(T_ARRAYS& data) const;
    template<class T_FACE_ARRAYS_2> void Average_Common_Face_Data(T_FACE_ARRAYS_2& data) const;
    template<class T_FACE_ARRAYS_2> void Sum_Common_Face_Data(T_FACE_ARRAYS_2& data) const;
    template<class T_FACE_ARRAYS_2> void Copy_Common_Face_Data(T_FACE_ARRAYS_2& data) const;
    template<class T_FACE_ARRAYS_2> void Assert_Common_Face_Data(T_FACE_ARRAYS_2& data,const T tolerance=0) const;
    void Synchronize_Dt(T& dt) const;
    void Synchronize_Ghost_Cells(int& ghost_cells) const;
    void Initialize(VECTOR<VECTOR<bool,2>,T_GRID::dimension>& domain_walls);
    bool Neighbor(const int axis,const int axis_side) const;
    void Allgather(int sent_data,ARRAY<int>& data);
//#############################################################
};
}
#endif


//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SMOKE_EXAMPLE__
#define __SMOKE_EXAMPLE__
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include "application/smoke/cache_data_include.h"
#include "application/smoke/cache_face_array.h"
#include "application/smoke/cache_options.h"
#include "application/smoke/nimbus_thread_queue.h"
#include "application/smoke/options.h"
#include "application/smoke/projection/laplace_solver_wrapper.h"
#include "data/physbam/translator_physbam_old.h"
#include "shared/nimbus.h"
namespace PhysBAM{

template<class TV>
class SMOKE_EXAMPLE
{
    typedef application::DataConfig DataConfig;
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    enum workaround1{d=TV::m};

    // data types
    typedef ARRAY<T,FACE_INDEX<TV::dimension> > T_FACE_ARRAY;
    typedef ARRAY<bool,FACE_INDEX<TV::dimension> > BOOL_FACE_ARRAY;
    typedef ARRAY<T,TV_INT> T_SCALAR_ARRAY;
    typedef ARRAY<bool,TV_INT> BOOL_SCALAR_ARRAY;
    typedef ARRAY<int,TV_INT> INT_SCALAR_ARRAY;

public:
    T dt_buffer;
    nimbus::NimbusThreadQueue* nimbus_thread_queue;
    nimbus::int_dimension_t kScale;
    GeometricRegion local_region;
    GeometricRegion relative_region;
    DataConfig data_config;
    STREAM_TYPE stream_type;
    T initial_time;
    int first_frame,last_frame;
    T frame_rate;
    std::string frame_title;
    int write_substeps_level;
    // bool write_output_files;
    std::string output_directory;
    int restart;
    bool write_debug_data;
    int number_of_ghost_cells;

    T cfl;

    GRID<TV> mac_grid;
    MPI_UNIFORM_GRID<GRID<TV> > *mpi_grid;
    PROJECTION_UNIFORM<GRID<TV> > projection;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities_ghost;
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<GRID<TV>,T, AVERAGING_UNIFORM<GRID<TV>, FACE_LOOKUP_UNIFORM<GRID<TV> > >,LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T,FACE_LOOKUP_UNIFORM<GRID<TV> > > > advection_scalar;
    BOUNDARY_UNIFORM<GRID<TV>,T> boundary_scalar;
    BOUNDARY_UNIFORM<GRID<TV>,T> *boundary;
    ARRAY<T, TV_INT> density;
    ARRAY<T, TV_INT> density_ghost;
    VECTOR<VECTOR<bool,2>,TV::dimension> domain_boundary;
    RANGE<TV> source;
    pthread_mutex_t lock;

    LaplaceSolverWrapper laplace_solver_wrapper;

    nimbus::TranslatorPhysBAMOld<TV> translator;

    // cache objects
    bool use_cache;
    typedef typename application::CacheFaceArray<T> TCacheFaceArray;
    typedef typename application::CacheFaceArray<bool> BoolCacheFaceArray;
    typedef typename application::CacheScalarArray<T> TCacheScalarArray;
    typedef typename application::CacheScalarArray<int> IntCacheScalarArray;
    typedef typename application::CacheScalarArray<bool> BoolCacheScalarArray;
    // typedef typename application::CacheParticleLevelsetEvolution<float> TCachePLE;
    TCacheFaceArray *cache_fv;
    TCacheFaceArray *cache_fvg;
    BoolCacheFaceArray *cache_psi_n;
    TCacheScalarArray *cache_dens, *cache_dens_ghost;
    // TCacheScalarArray *cache_phi3, *cache_phi7, *cache_phi8;
    TCacheScalarArray *cache_pressure, *cache_divergence;
    IntCacheScalarArray *cache_colors;
    BoolCacheScalarArray *cache_psi_d;
    // TCachePLE *cache_ple;
    // bool create_destroy_ple;

    SMOKE_EXAMPLE(const STREAM_TYPE stream_type_input,
                  bool use_threading,
                  int core_quota);
    SMOKE_EXAMPLE(const STREAM_TYPE stream_type_input,
                  application::AppCacheObjects *cache,
                  bool use_threading,
                  int core_quota);
    virtual ~SMOKE_EXAMPLE();

    T CFL(ARRAY<T, FACE_INDEX<TV::dimension> >& face_velocities);
    void CFL_Threaded(RANGE<TV_INT>& domain, ARRAY<T, FACE_INDEX<TV::dimension> >& face_velocities, T& dt);
    
    T Time_At_Frame(const int frame) const
    {return initial_time+(frame-first_frame)/frame_rate;}

    void Initialize_Grid(TV_INT counts,RANGE<TV> domain)
    {mac_grid.Initialize(counts, domain, true);}

    void Initialize_Fields()
    {for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()) face_velocities(iterator.Full_Index())=0;
    for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()) density(iterator.Cell_Index())=0;}

    void Get_Scalar_Field_Sources(const T time)
    {for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next())
	if(source.Lazy_Inside(iterator.Location())) density(iterator.Cell_Index())=1;}

    void Write_Output_Files(const int frame, int rank = -1);
    void Read_Output_Files(const int frame);
    void Set_Boundary_Conditions(const T time);

    void Save_To_Nimbus_No_Cache(const nimbus::Job *job, const nimbus::DataArray &da, const int frame);
    void Save_To_Nimbus(const nimbus::Job *job, const nimbus::DataArray &da, const int frame);
    void Load_From_Nimbus_No_Cache(const nimbus::Job *job, const nimbus::DataArray &da, const int frame);
    void Load_From_Nimbus(const nimbus::Job *job, const nimbus::DataArray &da, const int frame);

//#####################################################################
};
}
#endif

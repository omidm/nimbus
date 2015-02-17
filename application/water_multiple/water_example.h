//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __WATER_EXAMPLE__
#define __WATER_EXAMPLE__
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Boundaries/GEOMETRY_BOUNDARY_POLICY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Boundaries/BOUNDARY_PHI_WATER.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include "application/water_multiple/app_data_include.h"
#include "application/water_multiple/app_data_face_array.h"
#include "application/water_multiple/app_data_options.h"
#include "application/water_multiple/nimbus_thread_queue.h"
#include "application/water_multiple/options.h"
#include "application/water_multiple/projection/laplace_solver_wrapper.h"
#include "data/physbam/translator_physbam_old.h"
#include "shared/nimbus.h"
namespace PhysBAM{

template<class T_GRID> class LEVELSET_MULTIPLE_UNIFORM;

//NOTE: Should eventually derive off of an incompressible project
template<class TV>
class WATER_EXAMPLE:public LEVELSET_CALLBACKS<GRID<TV> >
{
    typedef application::DataConfig DataConfig;
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;
    typedef typename LEVELSET_POLICY<GRID<TV> >::PARTICLE_LEVELSET T_PARTICLE_LEVELSET;
    typedef typename GEOMETRY_BOUNDARY_POLICY<GRID<TV> >::BOUNDARY_PHI_WATER T_BOUNDARY_PHI_WATER;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<TV> >::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    enum workaround1{d=TV::m};

    // data types
    typedef ARRAY<T,FACE_INDEX<TV::dimension> > T_FACE_ARRAY;
    typedef ARRAY<bool,FACE_INDEX<TV::dimension> > T_FACE_ARRAY_BOOL;
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
    bool write_output_files;
    std::string output_directory;
    int number_of_ghost_cells;

    T cfl;

    GRID<TV> mac_grid;
    MPI_UNIFORM_GRID<GRID<TV> > *mpi_grid;
    PROJECTION_DYNAMICS_UNIFORM<GRID<TV> >& projection;
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<TV> >  &particle_levelset_evolution;
    INCOMPRESSIBLE_UNIFORM<GRID<TV> > incompressible;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities_ghost;
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T> advection_scalar;
    BOUNDARY_UNIFORM<GRID<TV>,T> boundary_scalar;
    BOUNDARY_UNIFORM<GRID<TV>,T> *boundary,*phi_boundary;
    T_BOUNDARY_PHI_WATER phi_boundary_water;
    VECTOR<VECTOR<bool,2>,TV::dimension> domain_boundary;
    T_GRID_BASED_COLLISION_GEOMETRY* collision_bodies_affecting_fluid;
    ARRAY<IMPLICIT_OBJECT<TV>*> sources;
    LaplaceSolverWrapper laplace_solver_wrapper;

    T_FACE_ARRAY u_interface_dummy;
    T_FACE_ARRAY force_dummy;
    T_FACE_ARRAY_BOOL valid_mask_dummy;
    T_FACE_ARRAY t_face_dummy;
    T_SCALAR_ARRAY t_scalar_dummy;
    BOOL_FACE_ARRAY b_face_dummy;
    BOOL_SCALAR_ARRAY b_scalar_dummy;
    INT_SCALAR_ARRAY i_scalar_dummy;

    nimbus::TranslatorPhysBAMOld<TV> translator;

    ARRAY<T, TV_INT> phi_ghost_bandwidth_seven;
    ARRAY<T, TV_INT> phi_ghost_bandwidth_eight;

    typedef application::StaticConfigValidMask StaticConfigValidMask;
    typedef application::StaticConfigUInterface StaticConfigUInterface;
    typedef application::StaticConfigForce StaticConfigForce;
    typedef application::StaticConfigCollisionBody StaticConfigCollisionBody;
    StaticConfigValidMask* static_config_valid_mask;
    StaticConfigUInterface* static_config_u_interface;
    StaticConfigForce* static_config_force;
    StaticConfigCollisionBody* static_config_collision_body;
    // app_data objects
    typedef typename application::AppDataFaceArray<T> TAppDataFaceArray;
    typedef typename application::AppDataFaceArray<bool> BoolAppDataFaceArray;
    typedef typename application::AppDataScalarArray<T> TAppDataScalarArray;
    typedef typename application::AppDataScalarArray<int> IntAppDataScalarArray;
    typedef typename application::AppDataScalarArray<bool> BoolAppDataScalarArray;
    typedef typename application::AppDataParticleLevelsetEvolution<float> TAppDataPLE;
    typedef application::AppDataSparseMatrix TAppDataSparseMatrix;
    typedef application::AppDataArrayM2C TAppDataArrayM2C;
    TAppDataFaceArray *app_data_fv;
    TAppDataFaceArray *app_data_fvg;
    BoolAppDataFaceArray *app_data_psi_n;
    TAppDataScalarArray *app_data_phi3, *app_data_phi7, *app_data_phi8;
    TAppDataScalarArray *app_data_pressure, *app_data_divergence;
    IntAppDataScalarArray *app_data_colors;
    application::AppDataRawGridArray *app_data_index_c2m;
    BoolAppDataScalarArray *app_data_psi_d;
    TAppDataPLE *app_data_ple;
    TAppDataSparseMatrix *app_data_matrix_a;
    TAppDataArrayM2C *app_data_index_m2c;
    application::AppDataVector *app_data_vector_b;
    bool create_destroy_ple;

    WATER_EXAMPLE(StaticConfigCollisionBody* app_data_collision_body,
                  const STREAM_TYPE stream_type_input,
                  nimbus::TaskThreadPool::TaskThreadList* allocated_threads);
    WATER_EXAMPLE(StaticConfigCollisionBody* app_data_collision_body,
                  const STREAM_TYPE stream_type_input,
                  application::AppAppObjects *app_data,
                  nimbus::TaskThreadPool::TaskThreadList* allocated_threads);
    WATER_EXAMPLE(StaticConfigCollisionBody* app_data_collision_body,
                  const STREAM_TYPE stream_type_input,
                  application::AppAppObjects *app_data,
                  PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<TV> > *ple,
                  nimbus::TaskThreadPool::TaskThreadList* allocated_threads);
    virtual ~WATER_EXAMPLE();
    
    T Time_At_Frame(const int frame) const
    {return initial_time+(frame-first_frame)/frame_rate;}

    void Get_Levelset_Velocity(const GRID<TV>& grid,T_LEVELSET& levelset,ARRAY<T,FACE_INDEX<TV::dimension> >& V_levelset,const T time) const PHYSBAM_OVERRIDE
    {V_levelset=face_velocities;}

    void Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time);
    void Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET_MULTIPLE_UNIFORM<GRID<TV> >& levelset_multiple,ARRAY<T,FACE_INDEX<TV::dimension> >& V_levelset,const T time) const PHYSBAM_OVERRIDE {}
    void Initialize_Grid(TV_INT counts,RANGE<TV> range);
    void Set_Boundary_Conditions(const T time);
    void Adjust_Phi_With_Sources(const T time);

    void Sum_Phi_3();

    void Initialize_Phi();

    void Write_Output_Files(const int frame, int rank = -1);
    void Read_Output_Files(const int frame);

    void Save_To_Nimbus_No_AppData(const nimbus::Job *job, const nimbus::DataArray &da, const int frame);
    void Save_To_Nimbus(const nimbus::Job *job, const nimbus::DataArray &da, const int frame);
    void Load_From_Nimbus(const nimbus::Job *job, const nimbus::DataArray &da, const int frame);

//#####################################################################
};
}
#endif

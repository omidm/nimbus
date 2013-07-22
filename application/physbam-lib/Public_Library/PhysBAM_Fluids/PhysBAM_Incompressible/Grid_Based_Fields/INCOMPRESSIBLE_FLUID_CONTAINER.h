//#####################################################################
// Copyright 2011, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_FLUID_CONTAINER
//#####################################################################
#ifndef __INCOMPRESSIBLE_FLUID_CONTAINER__
#define __INCOMPRESSIBLE_FLUID_CONTAINER__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/THREADED_UNIFORM_GRID.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GRID.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Boundaries/GEOMETRY_BOUNDARY_POLICY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Grid_Based_Fields/DENSITY_CONTAINER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Grid_Based_Fields/TEMPERATURE_CONTAINER.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID>
class INCOMPRESSIBLE_FLUID_CONTAINER:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_ARRAYS_SLIP T_FACE_ARRAYS_SLIP_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
public:
    RIGID_GRID<T_GRID> rigid_grid;
    T_FACE_ARRAYS_SCALAR face_velocities;
    T_FACE_ARRAYS_BOOL psi_N;
    T_ARRAYS_BOOL psi_D;
    T_ARRAYS_SCALAR viscosity;
    T_ARRAYS_SCALAR pressure;
    T_ARRAYS_SCALAR phi;
    DENSITY_CONTAINER<T_GRID> density_container;
    TEMPERATURE_CONTAINER<T_GRID> temperature_container;
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID> particle_levelset_evolution;
    T_GRID_BASED_COLLISION_GEOMETRY collision_bodies_affecting_fluid;
    T_FACE_ARRAYS_BOOL face_velocities_valid_mask;
    int number_of_ghost_cells;

    INCOMPRESSIBLE_FLUID_CONTAINER(RIGID_GRID_COLLECTION<T_GRID>& rigid_grid_collection,int number_of_ghost_cells_input=3);
    INCOMPRESSIBLE_FLUID_CONTAINER(RIGID_GRID<T_GRID> rigid_grid_input=RIGID_GRID<T_GRID>());
    ~INCOMPRESSIBLE_FLUID_CONTAINER();

//#####################################################################
    void Write_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const;
    void Read_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame);
    void Initialize_Grids(const bool initialize_phi=false,const bool initialize_viscosity=false,const bool initialize_density=false,const bool initialize_temperature=false,const bool use_ghost_cells=false);
    void Sync_Data(INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& fluid_container,THREADED_UNIFORM_GRID<T_GRID>& threaded_grid);
    void Distribute_Data(INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& fluid_container,THREADED_UNIFORM_GRID<T_GRID>& threaded_grid);
//#####################################################################
};
}
#endif

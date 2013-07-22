//#####################################################################
// Copyright 2011, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMPRESSIBLE_EXAMPLE
//#####################################################################
#ifndef __COMPRESSIBLE_EXAMPLE__
#define __COMPRESSIBLE_EXAMPLE__

#include <PhysBAM_Dynamics/Coupled_Evolution/COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_EIGENSYSTEM.h>
// #include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Solvers/COMPRESSIBLE_ADVECTION.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_CALLBACKS.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/EXAMPLE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM {

template<class TV>
class COMPRESSIBLE_EXAMPLE:public EXAMPLE<TV>, public EXAMPLE_FORCES_AND_VELOCITIES<TV>, public CONSERVATION_CALLBACKS<GRID<TV>,VECTOR<typename TV::SCALAR,TV::dimension+2> >
{
    typedef EXAMPLE<TV> BASE;typedef typename TV::SCALAR T;typedef GRID<TV> T_GRID;
    typedef VECTOR<T,TV::dimension+2> TV_DIMENSION;typedef VECTOR<int,TV::dimension> TV_INT;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY_SCALAR;
    typedef typename REBIND_LENGTH<T_BOUNDARY_SCALAR,TV::dimension+2>::TYPE T_BOUNDARY_DIMENSION_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_FACE_ARRAYS_DIMENSION_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    enum {d=TV::dimension+2};

  public:
    using BASE::output_directory;using BASE::first_frame;using BASE::restart;using BASE::Write_Frame_Title;using BASE::stream_type;using BASE::parse_args;

    T_GRID grid;
    MPI_UNIFORM_GRID<T_GRID> *mpi_grid;
    T_ARRAYS_DIMENSION_SCALAR U,U_save,U_ghost;
    T_FACE_ARRAYS_SCALAR face_velocities;
    T_ARRAYS_BOOL psi;
    T_FACE_ARRAYS_BOOL psi_N;
    VECTOR<VECTOR<bool,2>,TV::dimension> domain_walls;

    T cfl,global_min_dt;
    int resolution,eno_order,rk_order;
    bool use_rf,solid_affects_fluid,adaptive_time_step,clamp_fluxes,write_debug_data;
    EOS<T>* eos;
    T_BOUNDARY_DIMENSION_SCALAR* boundary;
    T_BOUNDARY_DIMENSION_SCALAR* boundary_scalar;
    VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension> cfl_eigensystem;
    CONSERVATION<T_GRID,d>* conservation_law_solver;
    COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES<TV>* kinematic_coupling_utilities;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>* collision_bodies_affecting_fluid;
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;

    COMPRESSIBLE_EXAMPLE(const STREAM_TYPE& stream_type,const int array_collection_type=0);
    virtual ~COMPRESSIBLE_EXAMPLE();

//#####################################################################
    virtual void Preprocess_Frame(const int frame) {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Postprocess_Frame(const int frame) {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Preprocess_Substep(const int frame,const int substep) {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Postprocess_Substep(const T dt,const T time) {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Initialize_Grids();
    virtual void Initialize_Bodies() {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Initialize_Fluid_State() = 0;
    virtual void Initialize_Boundaries() = 0;
    virtual void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE;
    virtual void Clamp_Dt_Adaptively(T_GRID& grid,const T_ARRAYS_DIMENSION_SCALAR& rhs,const T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_BOOL& psi,T_ARRAYS_SCALAR& rho_dt,T_ARRAYS_SCALAR& e_dt,const T& dt,T& clamp_rho,T& clamp_e);
    virtual void Clamp_Fluxes(T_GRID& grid,T_FACE_ARRAYS_DIMENSION_SCALAR& fluxes,T_ARRAYS_DIMENSION_SCALAR& rhs,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const T_ARRAYS_BOOL& psi,const T& dt,T& clamp_rho,T& clamp_e);
    void Compute_Dt_For_Internal_Energy(T_GRID& grid,T_FACE_ARRAYS_DIMENSION_SCALAR& fluxes,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,TV_INT& face_index,TV_INT& cell_index,const bool first,T& clamp_e,T& rho_dt,T& e_dt,const int axis);
    void Save_State();
    void Restore_State();
    void Set_Eigensystems();
    void Log_Parameters() const PHYSBAM_OVERRIDE;
    void Register_Options() PHYSBAM_OVERRIDE;
    void Parse_Options() PHYSBAM_OVERRIDE;
    void Parse_Late_Options() PHYSBAM_OVERRIDE;
//#####################################################################
    void Read_Output_Files(const int frame);
    void Limit_Dt(T& dt,const T time);
//#####################################################################
};
}
#endif // __COMPRESSIBLE_EXAMPLE__

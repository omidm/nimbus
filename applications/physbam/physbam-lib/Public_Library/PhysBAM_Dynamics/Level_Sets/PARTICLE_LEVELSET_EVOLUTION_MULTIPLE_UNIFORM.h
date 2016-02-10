//#####################################################################
// Copyright 2005-2006, Ron Fedkiw, Eran Guendelman, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM
//##################################################################### 
#ifndef __PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM__
#define __PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM__

#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_MULTIPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_POLICY.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_MULTIPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_ARRAY> class RUNGEKUTTA;

template<class T_GRID>
class PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM:public PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<PARTICLE_LEVELSET_PARTICLES<TV>*>::TYPE T_ARRAYS_PARTICLE_LEVELSET_PARTICLES;
    typedef typename T_ARRAYS_SCALAR::template REBIND<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>::TYPE T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES;
    typedef typename LEVELSET_POLICY<T_GRID>::FAST_LEVELSET_T T_FAST_LEVELSET;
    typedef typename LEVELSET_ADVECTION_POLICY<T_GRID>::LEVELSET_ADVECTION_MULTIPLE T_LEVELSET_ADVECTION_MULTIPLE;
    typedef typename LEVELSET_ADVECTION_POLICY<T_GRID>::FAST_LEVELSET_ADVECTION_T T_FAST_LEVELSET_ADVECTION;
public:
    typedef PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID> BASE;
    using BASE::track_mass;using BASE::runge_kutta_order_levelset;using BASE::runge_kutta_order_particles;using BASE::use_particle_levelset;using BASE::use_frozen_velocity;
    using BASE::cfl_number;using BASE::reseeding_frequency;using BASE::time;using BASE::use_fmm;using BASE::use_reinitialization;using BASE::V;using BASE::grid;

    ARRAY<T_ARRAYS_SCALAR> phis;
    PARTICLE_LEVELSET_MULTIPLE_UNIFORM<T_GRID> particle_levelset_multiple;
    ARRAY<T> initial_mass;
    ARRAY<RUNGEKUTTA<T_ARRAYS_SCALAR>*> rungekutta_phis;

    T_LEVELSET_ADVECTION_MULTIPLE levelset_advection_multiple;

    PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM(const T_GRID& grid_input,const int number_of_ghost_cells_input);
    virtual ~PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM();

//#####################################################################
    virtual PARTICLE_LEVELSET_MULTIPLE_UNIFORM<T_GRID>& Particle_Levelset_Multiple();
    virtual LEVELSET_MULTIPLE_UNIFORM<T_GRID>& Levelset_Multiple();
    virtual PARTICLE_LEVELSET_UNIFORM<T_GRID>& Particle_Levelset(const int i) PHYSBAM_OVERRIDE;
    virtual T_FAST_LEVELSET& Levelset(const int i) PHYSBAM_OVERRIDE;
    virtual T_FAST_LEVELSET_ADVECTION& Levelset_Advection(const int i) PHYSBAM_OVERRIDE;
    void Use_Semi_Lagrangian_Advection() PHYSBAM_OVERRIDE;
    void Use_Hamilton_Jacobi_Weno_Advection() PHYSBAM_OVERRIDE;
    void Use_Hamilton_Jacobi_Eno_Advection(const int order) PHYSBAM_OVERRIDE;
    void Track_Mass(const bool track_mass_input=true) PHYSBAM_OVERRIDE;
    void Initialize_Domain(const T_GRID& grid_input,const int number_of_regions);  // don't call up to the base class here because we don't need those variables initialized OVERRIDE PROBLEM
    void Initialize_Domain(const T_GRID& grid_input) PHYSBAM_OVERRIDE;
    void Make_Signed_Distance() PHYSBAM_OVERRIDE;
    void Make_Signed_Distance(const T stopping_distance) PHYSBAM_OVERRIDE;
    void Set_Number_Particles_Per_Cell(const int number_particles_per_cell) PHYSBAM_OVERRIDE;
    void Set_Levelset_Callbacks(LEVELSET_CALLBACKS<T_GRID>& levelset_callbacks) PHYSBAM_OVERRIDE;
    void Initialize_FMM_Initialization_Iterative_Solver(const bool refine_fmm_initialization_with_iterative_solver_input=true,const int fmm_initialization_iterations_input=10,
        const T fmm_initialization_iterative_tolerance_input=1e-2,const T fmm_initialization_iterative_drift_fraction_input=.1);
    void Bias_Towards_Negative_Particles(const bool bias_towards_negative_particles) PHYSBAM_OVERRIDE;
    void Set_Seed(const int seed) PHYSBAM_OVERRIDE;
    void Seed_Particles(const T time) PHYSBAM_OVERRIDE;
    void Delete_Particles_Outside_Grid() PHYSBAM_OVERRIDE;
    void Set_CFL_Number(const T cfl_number_input) PHYSBAM_OVERRIDE;
    void Initialize_Runge_Kutta() PHYSBAM_OVERRIDE;
    void Advance_To_Time(T_FACE_ARRAYS_SCALAR* face_velocities,const T stopping_time,const bool verbose=true) PHYSBAM_OVERRIDE;
    T Time_Step(const T stopping_time,bool& limited_by_stopping_time) PHYSBAM_OVERRIDE;
    T CFL(const bool need_to_get_velocity=true,const bool analytic_test=false) PHYSBAM_OVERRIDE;
    void Advance_One_Time_Step(T_FACE_ARRAYS_SCALAR* face_velocities,const T dt) PHYSBAM_OVERRIDE;
    void Advance_Levelset(const T dt) PHYSBAM_OVERRIDE;
    void Advance_Particles(const T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const bool analytic_test=false) PHYSBAM_OVERRIDE;
    T Advance_Particles(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T input_time) PHYSBAM_OVERRIDE;
    T Advance_Particles(T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T input_time) PHYSBAM_OVERRIDE;
    void Modify_Levelset_And_Particles(T_FACE_ARRAYS_SCALAR* face_velocities) PHYSBAM_OVERRIDE;
    void Modify_Levelset_And_Particles(T_FACE_ARRAYS_SCALAR* face_velocities,const T stopping_distance) PHYSBAM_OVERRIDE;
    void Reseed_Particles(const T time,const int time_step=0,T_ARRAYS_BOOL* cell_centered_mask=0,const bool verbose=true) PHYSBAM_OVERRIDE;
    void Fill_Levelset_Ghost_Cells(const T time) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

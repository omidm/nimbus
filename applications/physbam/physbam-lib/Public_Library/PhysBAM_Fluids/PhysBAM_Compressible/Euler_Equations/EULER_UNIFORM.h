//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_UNIFORM
//#####################################################################
//
// Input an EOS class, that is a class that inherits the virtual base class EOS.
// Input a GRID_1D class.
// Input U as 3 by (1,m) for mass, momentum, and energy.
//
// Use Set_Up_Cut_Out_Grid(psi) to define a cut out grid with psi as (1,m).
// When psi=true, solve the equaitions.
// When psi=false, do NOT solve the equations.
//
//#####################################################################
#ifndef __EULER_UNIFORM__
#define __EULER_UNIFORM__

#include <PhysBAM_Tools/Parallel_Computation/MPI_GRID_POLICY.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_CAVITATION_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_EIGENSYSTEM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_PROJECTION_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID>
class EULER_UNIFORM:public EULER<T_GRID>
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef VECTOR<T,T_GRID::dimension+2> TV_DIMENSION;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename REBIND<T_ARRAYS_SCALAR,bool>::TYPE T_ARRAYS_BOOL;typedef typename REBIND<T_ARRAYS_SCALAR,TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_FACE_ARRAYS_DIMENSION_SCALAR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename MPI_GRID_POLICY<T_GRID>::MPI_GRID T_MPI_GRID;
    typedef EULER<T_GRID> BASE;
    typedef VECTOR<bool,T_GRID::dimension> TV_BOOL;
    typedef typename T_ARRAYS_DIMENSION_SCALAR::ELEMENT T_ARRAYS_ELEMENT;
    typedef BOUNDARY_UNIFORM<T_GRID,TV_DIMENSION> T_BOUNDARY;
protected:
    using BASE::max_time_step;using BASE::cut_out_grid;using BASE::gravity;using BASE::downward_direction;
public:
    using BASE::boundary;using BASE::conservation;using BASE::eos;using BASE::Get_Velocity;using BASE::e;
    using BASE::Set_Max_Time_Step;using BASE::Set_Custom_Conservation;using BASE::Set_CFL_Number;using BASE::open_boundaries;
    using BASE::use_force;

    T_GRID grid;
    T_MPI_GRID* mpi_grid;
    T_ARRAYS_DIMENSION_SCALAR U,U_save; // mass, momentum, and energy
    const T_ARRAYS_DIMENSION_SCALAR& U_ghost;
    T_ARRAYS_BOOL* psi_pointer; // defines cut out grid
    T_ARRAYS_BOOL psi;
    VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension> eigensystems,eigensystems_default,eigensystems_pressureonly;
    bool timesplit,use_sound_speed_for_cfl,perform_rungekutta_for_implicit_part,compute_pressure_fluxes,thinshell;
    bool use_sound_speed_based_dt_multiple_for_cfl; // if set, dt will be set to multiplication_factor_for_sound_speed_based_dt*dt_based_on_c whenever this number is less than dt_based_on_u
    T multiplication_factor_for_sound_speed_based_dt;
    EULER_PROJECTION_UNIFORM<T_GRID> euler_projection;
    bool apply_cavitation_correction;
    EULER_CAVITATION_UNIFORM<TV> euler_cavitation_density;
    EULER_CAVITATION_UNIFORM<TV> euler_cavitation_internal_energy;
    T e_min,last_dt;
    TV_INT pressure_flux_dimension_indices;
    T_FACE_ARRAYS_SCALAR force;
    TV_DIMENSION initial_total_conserved_quantity,accumulated_boundary_flux;
private:
    mutable T_ARRAYS_DIMENSION_SCALAR U_ghost_private; // Forces us to only touch U_ghost through the const reference and let Fill_Ghost_Cells be the only thing modifying the variable;
    T_FACE_ARRAYS_DIMENSION_SCALAR fluxes_pressure;
    T_ARRAYS_SCALAR added_internal_energy;
    mutable bool ghost_cells_valid;
    mutable int ghost_cells_valid_ring;
public:
    bool need_to_remove_added_internal_energy,need_to_remove_added_internal_energy_save;

    EULER_UNIFORM(const T_GRID& grid_input);
    ~EULER_UNIFORM();

    T Get_Temperature(const TV_INT& cell)
    {
        T density=Get_Density(U,cell);
        T internal_energy=e(U,cell);
        return eos->T(density,internal_energy);
    }

//#####################################################################
    void Set_Up_Cut_Out_Grid(T_ARRAYS_BOOL& psi_input);
    void Set_Custom_Equation_Of_State(EOS<T>& eos_input);
    void Set_Custom_Boundary(T_BOUNDARY& boundary_input);
    void Set_Body_Force(const bool use_force_input=true);
    void Initialize_Domain(const T_GRID& grid_input);
    void Save_State(T_ARRAYS_DIMENSION_SCALAR& U_save,T_FACE_ARRAYS_SCALAR& face_velocities_save,bool& need_to_remove_added_internal_energy_save);
    void Restore_State(T_ARRAYS_DIMENSION_SCALAR& U_save,T_FACE_ARRAYS_SCALAR& face_velocities_save,bool& need_to_remove_added_internal_energy_save);
    void Get_Cell_Velocities(const T dt,const T time,const int ghost_cells,T_ARRAYS_VECTOR& centered_velocities);
    void Compute_Total_Conserved_Quantity(const bool update_boundary_flux,const T dt,TV_DIMENSION& total_conserved_quantity);
    void Invalidate_Ghost_Cells();
    void Warn_For_Low_Internal_Energy() const;
    bool Equal_Real_Data(const T_ARRAYS_DIMENSION_SCALAR& U1,const T_ARRAYS_DIMENSION_SCALAR& U2) const;
    void Fill_Ghost_Cells(const T dt,const T time,const int ghost_cells) const;
    void Get_Dirichlet_Boundary_Conditions(const T dt,const T time);
    void Advance_One_Time_Step_Forces(const T dt,const T time);
    void Compute_Cavitation_Velocity(T_ARRAYS_SCALAR& rho_n, T_FACE_ARRAYS_SCALAR& face_velocities_n, T_ARRAYS_DIMENSION_SCALAR& momentum_n);
    void Advance_One_Time_Step_Explicit_Part(const T dt,const T time,const int rk_substep,const int rk_order);
    void Advance_One_Time_Step_Implicit_Part(const T dt,const T time);
    void Clamp_Internal_Energy(const T dt,const T time); // TODO(kwatra): Do we really need dt, time here?
    void Clamp_Internal_Energy_Ghost(T_ARRAYS_DIMENSION_SCALAR& U_ghost,const int number_of_ghost_cells) const;
    void Remove_Added_Internal_Energy(const T dt,const T time); // TODO(kwatra): Do we really need dt, time here?
    void Euler_Step(const T dt,const T time,const bool is_time_n);
    T CFL_Using_Sound_Speed() const;
    T CFL(const T time) const;
    void Set_Eigensystems(const bool advection_only);
    void Log_Parameters() const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

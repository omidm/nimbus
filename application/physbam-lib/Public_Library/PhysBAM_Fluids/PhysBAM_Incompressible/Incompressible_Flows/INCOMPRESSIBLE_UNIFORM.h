//#####################################################################
// Copyright 2002-2010, Mridul Aanjaneya, Doug Enright, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Frank Losasso, Duc Nguyen, Nick Rasmussen, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_UNIFORM  
//#####################################################################
#ifndef __INCOMPRESSIBLE_UNIFORM__
#define __INCOMPRESSIBLE_UNIFORM__

#include <PhysBAM_Tools/Parallel_Computation/MPI_GRID_POLICY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_COLLIDABLE_POLICY_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/FLUID_STRAIN_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
namespace PhysBAM{

template<class TV> class BOUNDARY_CONDITIONS_CALLBACKS;

template<class T_GRID>
class INCOMPRESSIBLE_UNIFORM:public INCOMPRESSIBLE<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename ADVECTION_COLLIDABLE_POLICY<T_GRID>::ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;typedef typename T_ARRAYS_BASE::template REBIND<bool>::TYPE T_ARRAYS_BOOL_BASE;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_FACE_ARRAYS_BOOL::template REBIND<VECTOR<bool,T_GRID::dimension> >::TYPE T_FACE_ARRAYS_BOOL_DIMENSION;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    typedef EXTRAPOLATION_UNIFORM<T_GRID,T> T_EXTRAPOLATION_SCALAR;typedef typename MPI_GRID_POLICY<T_GRID>::MPI_GRID T_MPI_GRID;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef typename TV::SPIN T_SPIN;typedef typename T_ARRAYS_SCALAR::template REBIND<T_SPIN>::TYPE T_ARRAYS_SPIN;
public:
    typedef INCOMPRESSIBLE<T_GRID> BASE;
    using BASE::use_force;using BASE::surface_tension;using BASE::use_variable_surface_tension;using BASE::viscosity;using BASE::use_variable_viscosity;
    using BASE::use_variable_vorticity_confinement;using BASE::dt_old;using BASE::gravity;using BASE::nonzero_surface_tension;using BASE::nonzero_viscosity;
    using BASE::downward_direction;using BASE::Set_Custom_Advection;using BASE::valid_mask;
    using BASE::use_explicit_part_of_implicit_viscosity;using BASE::vorticity_confinement;using BASE::max_time_step;using BASE::advection;
    using BASE::maximum_implicit_viscosity_iterations;using BASE::conserve_kinetic_energy;
    
    T_GRID grid;
    T_MPI_GRID* mpi_grid;
    BOUNDARY_UNIFORM<T_GRID,T>* boundary;
    PROJECTION_DYNAMICS_UNIFORM<T_GRID>& projection;
    T_ARRAYS_SCALAR variable_surface_tension;
    T_ARRAYS_SCALAR variable_viscosity;
    T_FACE_ARRAYS_SCALAR force;
    T_FACE_ARRAYS_SCALAR kinetic_energy,potential_energy;
    T_ARRAYS_SCALAR variable_vorticity_confinement;
    FLUID_STRAIN_UNIFORM<T_GRID>* strain;
    T_GRID_BASED_COLLISION_GEOMETRY* collision_body_list;
    bool momentum_conserving_vorticity;    
    T allowed_energy_gained,analytic_energy,analytic_potential;
    bool use_analytic_energy;
    bool use_vorticity_weights;
    T_FACE_ARRAYS_SCALAR vorticity_weights;
    T energy_clamp;
    int vc_projection_direction;
    T buoyancy_constant;
    THREAD_QUEUE* thread_queue;
protected:               
    BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<T_GRID>& boundary_default;
    ADVECTION_MACCORMACK_UNIFORM<T_GRID,T,ADVECTION<T_GRID,T> >* advection_maccormack;
public:

    INCOMPRESSIBLE_UNIFORM(const T_GRID& grid_input,PROJECTION_DYNAMICS_UNIFORM<T_GRID>& projection_input,THREAD_QUEUE* thread_queue_input=0);
    virtual ~INCOMPRESSIBLE_UNIFORM();

    void Set_Custom_Boundary(BOUNDARY_UNIFORM<T_GRID,T>& boundary_input)
    {boundary=&boundary_input;}

    void Set_Collision_Body_List(T_GRID_BASED_COLLISION_GEOMETRY& collision_body_list_input)
    {collision_body_list=&collision_body_list_input;}

    void Conserve_Kinetic_Energy(bool conserve=true)
    {conserve_kinetic_energy=conserve;if(conserve){kinetic_energy.Resize(grid);potential_energy.Resize(grid);}else{kinetic_energy.Clean_Memory();potential_energy.Clean_Memory();}}

//#####################################################################
    void Initialize_Grids(const T_GRID& grid_input) PHYSBAM_OVERRIDE;
    void Set_Body_Force(const bool use_force_input=true);
    void Set_Variable_Surface_Tension(const bool use_variable_surface_tension_input=true);
    void Set_Variable_Viscosity(const bool use_variable_viscosity_input=true);
    void Use_Variable_Vorticity_Confinement(const bool use_variable_vorticity_confinement_input=true);
    void Use_Variable_Vorticity_Confinement(T_GRID& grid,const bool use_variable_vorticity_confinement_input=true);
    void Use_Strain();
    void Apply_Pressure_Kinetic_Energy(T_FACE_ARRAYS_SCALAR& face_velocities_new,T_FACE_ARRAYS_SCALAR& face_velocities_old,const T dt,const T time);
    void Calculate_Kinetic_Energy(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Update_Kinetic_Energy(T_FACE_ARRAYS_SCALAR& face_velocities_new,T_FACE_ARRAYS_SCALAR& face_velocities_old,const T dt,const T time);
    void Update_Potential_Energy(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Update_Potential_Energy(T_FACE_ARRAYS_SCALAR& face_velocities_new,T_FACE_ARRAYS_SCALAR& face_velocities_old,const T dt,const T time);
    void Correct_Kinetic_Energy(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Advect_With_Vorticity(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,const int iterations,const int number_of_ghost_cells);
    void Add_Energy_With_Vorticity(T_FACE_ARRAYS_SCALAR& face_velocities,const VECTOR<VECTOR<bool,2>,TV::dimension>& domain_boundary,const T dt,const T time,const int number_of_ghost_cells,T_LEVELSET* lsv=0,T_ARRAYS_SCALAR* density=0);
    void Advance_One_Time_Step_Convection(const T dt,const T time,const T_FACE_ARRAYS_SCALAR& advecting_face_velocities,T_FACE_ARRAYS_SCALAR& face_velocities_to_advect,const int number_of_ghost_cells);
    void Advance_One_Time_Step_Forces(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,const bool implicit_viscosity,const T_ARRAYS_SCALAR* phi_ghost,const int number_of_ghost_cells);
    void Add_Gravity_Threaded(RANGE<TV_INT>&domain,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,int axis);
    void Add_Body_Force_Threaded(RANGE<TV_INT>&domain,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,int axis);
    virtual void Advance_One_Time_Step_Implicit_Part(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,const bool implicit_viscosity=false,BOUNDARY_UNIFORM<T_GRID,T>* projection_boundary=0,
        bool use_levelset_viscosity=false,BOUNDARY_CONDITIONS_CALLBACKS<TV>* bc_callbacks=0,bool print_viscosity_matrix=false);
    int Real_CFL(T_FACE_ARRAYS_SCALAR& face_velocities,const bool inviscid,const bool viscous_only,T input_dt) const;
    T CFL(T_FACE_ARRAYS_SCALAR& face_velocities,const bool inviscid=false,const bool viscous_only=false) const;
    void Apply_Vorticity_Confinement_Force(T_FACE_ARRAYS_SCALAR& face_velocities,T_ARRAYS_VECTOR& F);
    virtual void Compute_Vorticity_Confinement_Force(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_ARRAYS_VECTOR& F);
    void Extrapolate_Velocity_Across_Interface(T_FACE_ARRAYS_SCALAR& face_velocities_input,T_ARRAYS_SCALAR& phi_ghost,const bool enforce_divergence_free=false,const T band_width=3,
        const T damping=0,const TV& air_speed=TV(),const T_FACE_ARRAYS_BOOL_DIMENSION* face_neighbors_visible=0,const T_FACE_ARRAYS_BOOL* fixed_faces_input=0);
    void Set_Dirichlet_Boundary_Conditions(const T_ARRAYS_SCALAR* phi,const T pressure);
    void Set_Dirichlet_Boundary_Conditions(const T_ARRAYS_SCALAR* phi,const T_ARRAYS_SCALAR& pressure);
    void Set_Boundary_Conditions_With_Extrapolation(const T_ARRAYS_SCALAR* phi,const T pressure,const T_ARRAYS_BOOL &extrapolated_psi_D);
    void Add_Surface_Tension(T_LEVELSET& levelset,const T time);
    void Use_Maccormack_Advection(const T_ARRAYS_BOOL* node_mask, const T_ARRAYS_BOOL* cell_mask, const T_FACE_ARRAYS_BOOL* face_mask);
    bool Consistent_Boundary_Conditions(T_FACE_ARRAYS_SCALAR& face_velocities) const;
    void Implicit_Viscous_Update(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
//#####################################################################
};
}
#endif


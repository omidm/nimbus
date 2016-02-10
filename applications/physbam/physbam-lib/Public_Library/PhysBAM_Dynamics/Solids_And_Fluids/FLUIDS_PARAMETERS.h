//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUIDS_PARAMETERS
//#####################################################################
#ifndef __FLUIDS_PARAMETERS__
#define __FLUIDS_PARAMETERS__

#include <PhysBAM_Tools/Fourier_Transforms_Calculations/TURBULENCE_POLICY.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Boundaries/GEOMETRY_BOUNDARY_POLICY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Grid_Based_Fields/INCOMPRESSIBLE_FLUID_CONTAINER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_POLICY.h>
namespace PhysBAM{

template<class T> class FLUIDS_PARAMETERS_CALLBACKS;
template<class T> class EOS;
template<class T_GRID,int d> class CONSERVATION;
template<class TV> class GRID;

template <class T_GRID>
class FLUIDS_PARAMETERS:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename TV::template REBIND<bool>::TYPE TV_BOOL;typedef VECTOR<int,TV::dimension> TV_INT;
    typedef typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX T_SYMMETRIC_MATRIX;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef typename ADVECTION_POLICY<T_GRID>::ADVECTION_SEMI_LAGRANGIAN_SCALAR T_ADVECTION_SEMI_LAGRANGIAN_SCALAR;
    typedef typename REBIND<T_ADVECTION_SEMI_LAGRANGIAN_SCALAR,TV>::TYPE T_ADVECTION_SEMI_LAGRANGIAN_VECTOR;
    typedef typename ADVECTION_POLICY<T_GRID>::ADVECTION_HAMILTON_JACOBI_WENO_SCALAR T_ADVECTION_HAMILTON_JACOBI_WENO_SCALAR;
    typedef typename TURBULENCE_POLICY<TV>::TURBULENCE T_TURBULENCE;typedef typename LEVELSET_POLICY<T_GRID>::PARTICLE_LEVELSET_EVOLUTION T_PARTICLE_LEVELSET_EVOLUTION;
    typedef typename LEVELSET_POLICY<T_GRID>::PARTICLE_LEVELSET T_PARTICLE_LEVELSET;typedef typename INCOMPRESSIBLE_POLICY<T_GRID>::INCOMPRESSIBLE T_INCOMPRESSIBLE;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY_SCALAR;typedef typename GEOMETRY_BOUNDARY_POLICY<T_GRID>::BOUNDARY_REFLECTION T_BOUNDARY_REFLECTION;
    typedef typename GEOMETRY_BOUNDARY_POLICY<T_GRID>::BOUNDARY_PHI_WATER T_BOUNDARY_PHI_WATER;
    typedef typename GEOMETRY_BOUNDARY_POLICY<T_GRID>::BOUNDARY_MAC_GRID_SOLID_WALL_SLIP T_BOUNDARY_MAC_GRID_SOLID_WALL_SLIP;
    typedef typename REBIND<T_BOUNDARY_SCALAR,T_SYMMETRIC_MATRIX>::TYPE T_BOUNDARY_SYMMETRIC_MATRIX;
    typedef typename REBIND_LENGTH<T_BOUNDARY_SCALAR,T_GRID::dimension+2>::TYPE T_BOUNDARY_DIMENSION_SCALAR;
public:
    const bool smoke,fire,water,sph,compressible;
    bool quadtree,octree;
    int number_of_ghost_cells;
    T cfl;
    T gravity;
    TV gravity_direction;
    T_GRID* grid;
protected:
    bool need_destroy_grid;
public:
    int maximum_tree_depth;
    VECTOR<VECTOR<bool,2>,TV::dimension> domain_walls;
    T levelset_refinement_bandwidth;
    bool minimal_air_bandwidth;
    T_BOUNDARY_SCALAR* phi_boundary;
    ARRAY<T_BOUNDARY_SCALAR*> phi_boundary_multiphase;
    T_BOUNDARY_REFLECTION& phi_boundary_reflection;
    T_BOUNDARY_PHI_WATER& phi_boundary_water;
    T particle_half_bandwidth;
    int reseeding_frame_rate;
    int reinitialize_geometry_frame_rate;
    bool bias_towards_negative_particles;
    bool use_particle_levelset;
    int number_particles_per_cell,positive_particles_buffer_size,negative_particles_buffer_size;
    bool use_removed_positive_particles,use_removed_negative_particles;
    int removed_positive_particles_buffer_size,removed_negative_particles_buffer_size;
    bool store_particle_ids;
    bool reincorporate_removed_particle_velocity;
    T removed_particle_mass_scaling;
    bool use_sph_for_removed_negative_particles;
    T normal_flame_speed,curvature_flame_speed;
    T_BOUNDARY_SCALAR* fluid_boundary;
    T_BOUNDARY_SCALAR& fluid_boundary_water;
    T_BOUNDARY_MAC_GRID_SOLID_WALL_SLIP& boundary_mac_slip;
    T incompressible_tolerance;
    int incompressible_iterations;
    bool show_residual;
    int cg_restart_iterations;
    bool use_body_force;
    T density,outside_density,density_fuel;
    T surface_tension;
    bool variable_surface_tension;
    T viscosity,viscosity_fuel;
    bool variable_viscosity;
    bool implicit_viscosity,use_explicit_part_of_implicit_viscosity;
    int implicit_viscosity_iterations;
    T implicit_viscosity_tolerance;
    bool use_coupled_implicit_viscosity;
    bool use_vorticity_confinement,use_vorticity_confinement_fuel,use_variable_vorticity_confinement;
    T confinement_parameter,confinement_parameter_fuel;
    bool use_non_zero_divergence;
    bool second_order_cut_cell_method,enforce_divergence_free_extrapolation;
    bool solve_neumann_regions;
    bool solve_single_cell_neumann_regions;
    KRYLOV_SOLVER_TYPE evolution_solver_type;
    T kolmogorov;
    T turbulence_lowest_angular_frequency;
    int turbulence_update_frame_rate;
    GRID<TV> turbulence_grid;
    T_TURBULENCE& turbulence;
    bool use_external_velocity;
    bool use_soot,use_density,use_temperature;
    bool use_fixed_soot_boundary,use_fixed_density_boundary,use_fixed_temperature_boundary;
    int soot_advection_order;
    T ambient_soot,ambient_density,ambient_temperature;
    DENSITY_CONTAINER<T_GRID> soot_container;
    DENSITY_CONTAINER<T_GRID> soot_fuel_container;
    DENSITY_CONTAINER<T_GRID>& density_container;
    TEMPERATURE_CONTAINER<T_GRID>& temperature_container;
    bool use_soot_fuel_combustion;
    T burn_temperature_threshold,burn_rate,soot_fuel_calorific_value;
    T_BOUNDARY_SCALAR *soot_boundary,*density_boundary,*temperature_boundary;
    T temperature_fuel,temperature_products;
    T density_buoyancy_threshold,density_buoyancy_constant,temperature_buoyancy_constant;
    T removed_positive_particle_buoyancy_constant;
    T rho_bottom,rho_top;
    bool use_strain;
    T elastic_modulus;
    T plasticity_alpha,plasticity_gamma;
    T adhesion_coefficient,adhesion_normal_strain,adhesion_half_bandwidth;
    T_BOUNDARY_SYMMETRIC_MATRIX* strain_boundary;
    T_BOUNDARY_SYMMETRIC_MATRIX& strain_boundary_default;
    T object_friction; // for solid object (0 = no friction, 1 = full friction)
    bool delete_fluid_inside_objects; // both phi and the particles
    bool move_grid,move_grid_explicitly;
    int moving_grid_number_of_cells;
    bool monitor_mass;
    T mass;
    bool write_levelset,write_particles,write_removed_positive_particles,write_removed_negative_particles,write_flattened_particles;
    bool write_velocity,write_strain;
    bool write_debug_data,write_restart_data,write_ghost_values;
    int restart_data_write_rate; // defaults to writing restart data every frame, but increase to avoid writing large data each frame
    T_ADVECTION_SEMI_LAGRANGIAN_SCALAR& semi_lagrangian;
    T_ADVECTION_HAMILTON_JACOBI_WENO_SCALAR& hamilton_jacobi_weno;
    bool simulate;
    T min_collision_distance_factor,max_collision_distance_factor;
    bool solid_affects_fluid,fluid_affects_solid;
    bool thin_shells_refine_near_objects;
    TV_BOOL periodic;
    bool use_separation_inside_water;
    T separation_velocity_tolerance;
    FLUIDS_PARAMETERS_CALLBACKS<T_GRID>* callbacks;
    enum TYPE {NONE,SMOKE,FIRE,WATER,SPH,COMPRESSIBLE};
    bool refine_fmm_initialization_with_iterative_solver;
    bool modify_wall_tangential_velocities;
    T_GRID_BASED_COLLISION_GEOMETRY* collision_bodies_affecting_fluid;
protected:
    bool need_destroy_collision_bodies_affecting_fluid;
public:
    T collidable_contour_value;
    T collidable_phi_replacement_value;
    bool flood_fill_for_bubbles;
    bool use_maccormack_semi_lagrangian_advection;
    bool use_maccormack_compute_mask;
    bool use_maccormack_for_level_set;
    bool use_maccormack_for_incompressible;
    T bandwidth_without_maccormack_near_interface;
    bool mass_conservation;
    int mass_conservation_minimum_refinement_depth;
    int mass_conservation_maximum_refinement_depth;
    bool analytic_test;
    T_BOUNDARY_DIMENSION_SCALAR* compressible_boundary;
    T_BOUNDARY_SCALAR* compressible_pressure_boundary;
    EOS<T>* compressible_eos;
    CONSERVATION<T_GRID,T_GRID::dimension+2>* compressible_conservation_method;
    bool compressible_set_max_time_step;
    T compressible_max_time_step;
    int compressible_spatial_order;
    int compressible_rungekutta_order;
    T compressible_tolerance;
    int compressible_iterations;
    bool compressible_timesplit;
    bool compressible_apply_isobaric_fix;
    bool compressible_monitor_conservation_error;
    bool compressible_perform_rungekutta_for_implicit_part;
    bool compressible_use_sound_speed_for_cfl;
    bool compressible_use_sound_speed_based_dt_multiple_for_cfl;
    T compressible_multiplication_factor_for_sound_speed_based_dt;
    bool compressible_apply_cavitation_correction;
    bool compressible_adaptive_time_step;
    bool compressible_log_extremas;
    bool use_poisson;
    bool solids_override_source_velocities;
    bool use_slip;
    bool use_slip_constraints_across_non_occluded_faces;
    bool use_preconditioner_for_slip_system;
    bool stokes_flow;
    THREAD_QUEUE* thread_queue;
    int number_of_threads;
    bool use_trapezoid_rule;

    FLUIDS_PARAMETERS(const TYPE type,INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& incompressible_fluid_container);
    virtual ~FLUIDS_PARAMETERS();

    void Set_Fluids_Parameters_Callbacks(FLUIDS_PARAMETERS_CALLBACKS<T_GRID>& callbacks_input)
    {callbacks=&callbacks_input;}

    void Set_Default_Number_Particles_Per_Cell(const VECTOR<T,1>&)
    {number_particles_per_cell=4;}

    void Set_Default_Number_Particles_Per_Cell(const VECTOR<T,2>&)
    {number_particles_per_cell=16;}

    void Set_Default_Number_Particles_Per_Cell(const VECTOR<T,3>&)
    {number_particles_per_cell=64;}

    static GRID<TV> Default_Turbulence_Grid()
    {return GRID<TV>(TV_INT()+64,RANGE<TV>::Unit_Box());}

//#####################################################################
    virtual void Use_Fluid_Coupling_Defaults()=0;
    virtual void Use_No_Fluid_Coupling_Defaults()=0;
    void Initialize_Turbulence(const T time,const T frame_rate);
    void Initialize_Domain_Boundary_Conditions();
    void Initialize_Soot(const T time);
    void Evolve_Soot(const T dt,const T time);
    void Initialize_Density_And_Temperature(const T time);
    void Evolve_Density_And_Temperature(const T dt,const T time);
    virtual void Log_Parameters() const;
//#####################################################################
};
}
#endif

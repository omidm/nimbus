//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Avi Robinson-Mosher, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUIDS_PARAMETERS_UNIFORM
//#####################################################################
#ifndef __FLUIDS_PARAMETERS_UNIFORM__    
#define __FLUIDS_PARAMETERS_UNIFORM__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/THREADED_UNIFORM_GRID.h>
#include <PhysBAM_Geometry/Grids_PDE_Linear/LAPLACE_COLLIDABLE_POLICY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_FORWARD.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM_FORWARD.h>
#include <PhysBAM_Dynamics/Particles/PARTICLES_FORWARD.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS.h>
namespace PhysBAM{

template<class T_GRID> class MPI_UNIFORM_GRID;
template<class T_GRID> class PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM;
template<class T_GRID> class EULER_UNIFORM;
template<class TV> class SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES;
template<class TV> class COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES;
template<class T_GRID> class LAPLACE_UNIFORM;

template<class T_GRID>
class FLUIDS_PARAMETERS_UNIFORM:public FLUIDS_PARAMETERS<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef VECTOR<T,T_GRID::dimension+2> TV_DIMENSION;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef INCOMPRESSIBLE_UNIFORM<T_GRID> T_INCOMPRESSIBLE;typedef typename LEVELSET_POLICY<T_GRID>::EXTRAPOLATION_SCALAR T_EXTRAPOLATION_SCALAR;
    typedef typename LEVELSET_POLICY<T_GRID>::EXTRAPOLATION_VECTOR T_EXTRAPOLATION_VECTOR;typedef typename TURBULENCE_POLICY<TV>::TURBULENCE T_TURBULENCE;
    typedef typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX T_SYMMETRIC_MATRIX;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;typedef typename REBIND<T_ARRAYS_SCALAR,bool>::TYPE T_ARRAYS_BOOL;
    typedef typename REBIND<T_ARRAYS_SCALAR,TV>::TYPE T_ARRAYS_VECTOR;typedef typename LEVELSET_POLICY<T_GRID>::PARTICLE_LEVELSET_EVOLUTION T_PARTICLE_LEVELSET_EVOLUTION;
    typedef typename REBIND<T_ARRAYS_SCALAR,PARTICLE_LEVELSET_PARTICLES<TV>*>::TYPE T_ARRAYS_PARTICLE_LEVELSET_PARTICLES;
    typedef typename REBIND<T_ARRAYS_SCALAR,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>::TYPE T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename REBIND<T_FACE_ARRAYS_SCALAR,bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;typedef typename LEVELSET_POLICY<T_GRID>::FAST_LEVELSET_T T_FAST_LEVELSET;
    typedef typename ADVECTION_POLICY<T_GRID>::ADVECTION_SEMI_LAGRANGIAN_SCALAR T_ADVECTION_SEMI_LAGRANGIAN_SCALAR;
    typedef typename REBIND<T_ADVECTION_SEMI_LAGRANGIAN_SCALAR,T_SYMMETRIC_MATRIX>::TYPE T_ADVECTION_SEMI_LAGRANGIAN_SYMMETRIC_MATRIX;
    typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_LOOKUP_COLLIDABLE T_FACE_LOOKUP_COLLIDABLE;typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename LEVELSET_POLICY<T_GRID>::PARTICLE_LEVELSET T_PARTICLE_LEVELSET;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET_MULTIPLE T_LEVELSET_MULTIPLE;
    typedef typename REBIND<T_ARRAYS_SCALAR,T_SYMMETRIC_MATRIX>::TYPE T_ARRAYS_SYMMETRIC_MATRIX;
public:
    typedef FLUIDS_PARAMETERS<T_GRID> BASE;
    using BASE::smoke;using BASE::fire;using BASE::water;using BASE::use_body_force;using BASE::grid;
    using BASE::soot_container;using BASE::soot_fuel_container;using BASE::density_container;using BASE::temperature_container;
    using BASE::use_soot_fuel_combustion;using BASE::burn_temperature_threshold;using BASE::burn_rate;
    using BASE::soot_fuel_calorific_value;
    using BASE::domain_walls;
    using BASE::callbacks;using BASE::gravity;using BASE::gravity_direction;using BASE::phi_boundary;using BASE::fluid_boundary;using BASE::fluid_boundary_water;
    using BASE::phi_boundary_water;using BASE::boundary_mac_slip;using BASE::normal_flame_speed;using BASE::curvature_flame_speed;using BASE::surface_tension;
    using BASE::variable_surface_tension;using BASE::viscosity;using BASE::viscosity_fuel;using BASE::variable_viscosity;using BASE::implicit_viscosity;
    using BASE::implicit_viscosity_iterations;using BASE::use_vorticity_confinement;using BASE::confinement_parameter;
    using BASE::use_vorticity_confinement_fuel;using BASE::use_variable_vorticity_confinement;using BASE::use_strain;using BASE::elastic_modulus;using BASE::plasticity_alpha;
    using BASE::plasticity_gamma;using BASE::adhesion_coefficient;using BASE::confinement_parameter_fuel;using BASE::use_explicit_part_of_implicit_viscosity;
    using BASE::levelset_refinement_bandwidth;using BASE::kolmogorov;using BASE::use_external_velocity;
    using BASE::use_soot;using BASE::use_density;using BASE::use_temperature;using BASE::density;using BASE::soot_advection_order;
    using BASE::density_fuel;using BASE::temperature_products;using BASE::temperature_fuel;using BASE::number_particles_per_cell;using BASE::turbulence_lowest_angular_frequency;
    using BASE::turbulence_update_frame_rate;using BASE::use_non_zero_divergence;using BASE::solve_neumann_regions;using BASE::solve_single_cell_neumann_regions;
    using BASE::temperature_buoyancy_constant;using BASE::object_friction;
    using BASE::move_grid;using BASE::move_grid_explicitly;using BASE::moving_grid_number_of_cells;using BASE::write_velocity;using BASE::write_levelset;using BASE::write_particles;
    using BASE::write_debug_data;using BASE::restart_data_write_rate;using BASE::write_removed_positive_particles;using BASE::write_removed_negative_particles;using BASE::write_strain;
    using BASE::write_ghost_values;using BASE::solid_affects_fluid;using BASE::fluid_affects_solid;using BASE::Initialize_Density_And_Temperature;using BASE::separation_velocity_tolerance;using BASE::write_restart_data;
    using BASE::use_separation_inside_water;using BASE::adhesion_normal_strain;using BASE::modify_wall_tangential_velocities;using BASE::store_particle_ids;using BASE::turbulence;
    using BASE::collision_bodies_affecting_fluid;using BASE::collidable_contour_value;using BASE::collidable_phi_replacement_value;
    using BASE::semi_lagrangian;
    using BASE::phi_boundary_reflection;using BASE::phi_boundary_multiphase;using BASE::flood_fill_for_bubbles;using BASE::adhesion_half_bandwidth;using BASE::sph;
    using BASE::use_maccormack_semi_lagrangian_advection;using BASE::use_maccormack_for_incompressible;using BASE::use_maccormack_for_level_set;using BASE::number_of_ghost_cells;
    using BASE::cfl;using BASE::density_buoyancy_constant;using BASE::rho_top;using BASE::rho_bottom;using BASE::density_buoyancy_threshold;using BASE::use_sph_for_removed_negative_particles;
    using BASE::mass_conservation;using BASE::analytic_test;using BASE::hamilton_jacobi_weno;using BASE::compressible;using BASE::compressible_boundary;using BASE::compressible_pressure_boundary;
    using BASE::compressible_eos;using BASE::compressible_conservation_method;using BASE::compressible_set_max_time_step;using BASE::compressible_max_time_step;
    using BASE::compressible_spatial_order;using BASE::compressible_rungekutta_order;using BASE::compressible_timesplit;using BASE::compressible_apply_isobaric_fix;using BASE::compressible_apply_cavitation_correction;using BASE::compressible_adaptive_time_step;
    using BASE::write_flattened_particles;using BASE::use_poisson;using BASE::simulate;using BASE::use_slip;using BASE::thread_queue;using BASE::number_of_threads;using BASE::removed_positive_particle_buoyancy_constant;
    using BASE::bandwidth_without_maccormack_near_interface;

    MPI_UNIFORM_GRID<T_GRID>* mpi_grid;
    T_GRID p_grid;

    T_PARTICLE_LEVELSET_EVOLUTION* particle_levelset_evolution;
    INCOMPRESSIBLE_UNIFORM<T_GRID>* incompressible;
    PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>* particle_levelset_evolution_multiple;
    INCOMPRESSIBLE_MULTIPHASE_UNIFORM<T_GRID>* incompressible_multiphase;
    SPH_EVOLUTION_UNIFORM<T_GRID>* sph_evolution;
    T_ARRAYS_BOOL& maccormack_node_mask;
    T_ARRAYS_BOOL& maccormack_cell_mask;
    T_FACE_ARRAYS_BOOL& maccormack_face_mask;
    ADVECTION_MACCORMACK_UNIFORM<T_GRID,T,T_ADVECTION_SEMI_LAGRANGIAN_SCALAR>& maccormack_semi_lagrangian;
    EULER_UNIFORM<T_GRID>* euler;
    SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>* euler_solid_fluid_coupling_utilities;
    COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<TV>* compressible_incompressible_coupling_utilities;
    ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>* advection_conservative;
    PROJECTION_DYNAMICS_UNIFORM<T_GRID>* projection;

    // multiphase parameters
    ARRAY<T> masses; // used to keep track of mass loss in the driver
    ARRAY<T> densities;
    ARRAY<T> viscosities;
    ARRAY<T,VECTOR<int,2> > surface_tensions;
    ARRAY<T> confinement_parameters;
    ARRAY<bool> dirichlet_regions;
    ARRAY<bool> pseudo_dirichlet_regions;
    ARRAY<bool> fuel_region; // is this region a fuel region
    bool use_reacting_flow;
    ARRAY<T,VECTOR<int,2> > normal_flame_speeds; // specify 0 for non reacting flow.  must be symmetric
    ARRAY<T,VECTOR<int,2> > curvature_flame_speeds; // specify 0 for non reacting flow.  must be skew symmetric
    int number_of_regions;
    bool use_flame_speed_multiplier;
    bool use_dsd;
    ARRAY<bool> use_multiphase_strain;
    ARRAY<T> elastic_moduli;
    ARRAY<T> plasticity_alphas;
    ARRAY<T> plasticity_gammas;
    bool use_psi_R;
    bool use_levelset_viscosity;
    bool print_viscosity_matrix;
    bool use_second_order_pressure;
    bool use_conservative_advection;
    bool use_modified_projection;
    bool use_surface_solve;
    int projection_scale;
    VECTOR<bool,T_GRID::dimension> periodic_boundary;

    FLUIDS_PARAMETERS_UNIFORM(const int number_of_regions,const typename FLUIDS_PARAMETERS<T_GRID>::TYPE type,INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& incompressible_fluid_container);
    virtual ~FLUIDS_PARAMETERS_UNIFORM();

//#####################################################################
    virtual void Initialize_Grids();
    void Initialize_Fluid_Evolution(T_FACE_ARRAYS_SCALAR& face_velocities);
    void Use_Fluid_Coupling_Defaults() PHYSBAM_OVERRIDE;
    void Use_No_Fluid_Coupling_Defaults() PHYSBAM_OVERRIDE;
    void Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time);
    void Delete_Particles_Inside_Objects(const T time);
    void Initialize_Number_Of_Regions(const int number_of_regions_input);
private:
    template<class T_PARTICLES> void Delete_Particles_Inside_Objects(typename T_ARRAYS_SCALAR::template REBIND<T_PARTICLES*>::TYPE& particles,
        const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time);
public:
    void Set_Projection(PROJECTION_DYNAMICS_UNIFORM<T_GRID>* projection_input);
    void Update_Fluid_Parameters(const T dt,const T time);
    void Get_Body_Force(T_FACE_ARRAYS_SCALAR& force,const T dt,const T time);
    void Apply_Isobaric_Fix(const T dt,const T time);
    void Get_Neumann_And_Dirichlet_Boundary_Conditions(LAPLACE_UNIFORM<T_GRID>* elliptic_solver,
            T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Set_Domain_Boundary_Conditions(LAPLACE_UNIFORM<T_GRID>& elliptic_solver,T_FACE_ARRAYS_SCALAR& face_velocities,const T time);
    void Blend_In_External_Velocity(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Move_Grid(T_FACE_ARRAYS_SCALAR& face_velocities,const T time);
    void Move_Grid(T_FACE_ARRAYS_SCALAR& face_velocities,const TV_INT& shift_domain,const T time);
    void Adjust_Strain_For_Object(T_LEVELSET& levelset_object,T_ARRAYS_SYMMETRIC_MATRIX& e_ghost,const T time);
    void Combustion(const T dt,const T time);
    void Evolve_Soot(const T dt,const T time);
    void Sync_Parameters(FLUIDS_PARAMETERS_UNIFORM<T_GRID>& single_parameters,THREADED_UNIFORM_GRID<T_GRID>& threaded_grid);
    void Distribute_Parameters(FLUIDS_PARAMETERS_UNIFORM<T_GRID>& single_parameters,THREADED_UNIFORM_GRID<T_GRID>& threaded_grid);
    template<class T_ARRAYS_PARTICLES> int Total_Number_Of_Particles(const T_ARRAYS_PARTICLES& particles) const;
    template<class T_ARRAYS_PARTICLES> void Write_Particles(const STREAM_TYPE stream_type,const POINT_CLOUD<TV>& template_particles,const T_ARRAYS_PARTICLES& particles,
        const std::string& output_directory,const std::string& prefix,const int frame) const;
    template<class T_PARTICLES,class T_ARRAYS_PARTICLES> void Read_Particles(const STREAM_TYPE stream_type,const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,
        const std::string& output_directory,const std::string& prefix,const int frame);
    void Read_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame);
    void Write_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int first_frame,const int frame) const;
    void Log_Parameters() const PHYSBAM_OVERRIDE;
//#####################################################################
};      
}
#endif

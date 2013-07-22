//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUIDS_PARAMETERS_DYADIC
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __FLUIDS_PARAMETERS_DYADIC__    
#define __FLUIDS_PARAMETERS_DYADIC__

#include <PhysBAM_Tools/Fourier_Transforms_Calculations/TURBULENCE_POLICY.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_OCTREE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_QUADTREE.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_DYADIC.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS.h>
namespace PhysBAM{

template<class T_GRID>
class FLUIDS_PARAMETERS_DYADIC:public FLUIDS_PARAMETERS<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::MAP_MESH MAP_MESH;typedef typename T_GRID::CELL T_CELL;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FLOOD_FILL FLOOD_FILL;
    typedef typename T_GRID::UNIFORM_GRID UNIFORM_GRID;typedef typename UNIFORM_GRID::CELL_ITERATOR UNIFORM_CELL_ITERATOR;typedef typename UNIFORM_GRID::NODE_ITERATOR UNIFORM_NODE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_UNIFORM_ARRAYS;typedef typename LEVELSET_POLICY<T_GRID>::PARTICLE_LEVELSET T_PARTICLE_LEVELSET;
    typedef typename TURBULENCE_POLICY<TV>::TURBULENCE T_TURBULENCE;
    typedef typename INCOMPRESSIBLE_POLICY<T_GRID>::INCOMPRESSIBLE T_INCOMPRESSIBLE;typedef typename LEVELSET_POLICY<T_GRID>::EXTRAPOLATION_SCALAR EXTRAPOLATION_SCALAR;
    typedef typename T_UNIFORM_ARRAYS::template REBIND<TV>::TYPE T_UNIFORM_ARRAYS_VECTOR;typedef typename MATRIX_POLICY<TV>::TRANSFORMATION_MATRIX TRANSFORMATION_MATRIX;
    typedef typename LEVELSET_POLICY<UNIFORM_GRID>::LEVELSET UNIFORM_LEVELSET;typedef typename LEVELSET_POLICY<T_GRID>::EXTRAPOLATION_VECTOR EXTRAPOLATION_VECTOR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_HELPER T_LINEAR_INTERPOLATION_HELPER;typedef typename FIRE_INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP_FIRE T_FACE_LOOKUP_FIRE;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_HELPER T_LINEAR_INTERPOLATION_DYADIC_HELPER;typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename LEVELSET_POLICY<T_GRID>::PARTICLE_LEVELSET_EVOLUTION T_PARTICLE_LEVELSET_EVOLUTION;
    typedef AVERAGING_DYADIC<T_GRID,T_FACE_LOOKUP_FIRE> T_AVERAGING_FIRE;typedef LINEAR_INTERPOLATION_DYADIC<T_GRID,T,T_FACE_LOOKUP_FIRE> T_LINEAR_INTERPOLATION_FIRE;
public: 
    typedef FLUIDS_PARAMETERS<T_GRID> BASE;
    using BASE::smoke;using BASE::fire;using BASE::water;using BASE::use_body_force;using BASE::grid;using BASE::density_container;using BASE::temperature_container;using BASE::callbacks;
    using BASE::minimal_air_bandwidth;using BASE::domain_walls;
    using BASE::gravity;using BASE::gravity_direction;using BASE::phi_boundary;using BASE::fluid_boundary;using BASE::phi_boundary_water;using BASE::fluid_boundary_water;
    using BASE::normal_flame_speed;using BASE::curvature_flame_speed;using BASE::viscosity;using BASE::viscosity_fuel;using BASE::variable_viscosity;using BASE::implicit_viscosity;
    using BASE::implicit_viscosity_iterations;using BASE::use_vorticity_confinement;using BASE::confinement_parameter;using BASE::use_vorticity_confinement_fuel;using BASE::semi_lagrangian;
    using BASE::use_variable_vorticity_confinement;using BASE::use_strain;using BASE::elastic_modulus;using BASE::plasticity_alpha;using BASE::plasticity_gamma;using BASE::adhesion_coefficient;
    using BASE::confinement_parameter_fuel;using BASE::use_explicit_part_of_implicit_viscosity;using BASE::levelset_refinement_bandwidth;using BASE::kolmogorov;using BASE::use_external_velocity;
    using BASE::use_density;using BASE::use_temperature;using BASE::density;using BASE::density_fuel;using BASE::temperature_products;using BASE::temperature_fuel;
    using BASE::number_particles_per_cell;using BASE::turbulence_lowest_angular_frequency;using BASE::turbulence_update_frame_rate;using BASE::use_non_zero_divergence;
    using BASE::solve_neumann_regions;using BASE::temperature_buoyancy_constant;using BASE::object_friction;using BASE::move_grid;using BASE::move_grid_explicitly;using BASE::moving_grid_number_of_cells;
    using BASE::write_velocity;using BASE::write_levelset;using BASE::write_particles;using BASE::write_debug_data;using BASE::restart_data_write_rate;using BASE::turbulence;
    using BASE::write_removed_positive_particles;using BASE::write_removed_negative_particles;using BASE::write_strain;using BASE::write_ghost_values;using BASE::solid_affects_fluid;using BASE::fluid_affects_solid;
    using BASE::maximum_tree_depth;using BASE::Initialize_Density_And_Temperature;using BASE::write_restart_data;using BASE::modify_wall_tangential_velocities;using BASE::store_particle_ids;
    using BASE::use_separation_inside_water;using BASE::collision_bodies_affecting_fluid;using BASE::boundary_mac_slip;using BASE::collidable_contour_value;using BASE::collidable_phi_replacement_value;
    using BASE::rho_top;using BASE::rho_bottom;using BASE::density_buoyancy_constant;using BASE::density_buoyancy_threshold;

    ARRAY<VECTOR<T,2> > density_refinement_thresholds;
    ARRAY<VECTOR<T,2> > temperature_refinement_thresholds;

    T_PARTICLE_LEVELSET_EVOLUTION particle_levelset_evolution;
    T_INCOMPRESSIBLE incompressible;

    FLUIDS_PARAMETERS_DYADIC(const typename FLUIDS_PARAMETERS<T_GRID>::TYPE type);
    virtual ~FLUIDS_PARAMETERS_DYADIC();
        
//#####################################################################
    void Initialize_Grids();
    void Use_Fluid_Coupling_Defaults() PHYSBAM_OVERRIDE;
    void Use_Fluid_Slip_Coupling_Defaults() PHYSBAM_OVERRIDE;
    void Use_No_Fluid_Coupling_Defaults() PHYSBAM_OVERRIDE;
    void Get_Levelset_Velocity(const T_GRID& grid,T_LEVELSET& levelset,ARRAY<T>& face_velocities_levelset,const T time) const;
    void Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time);
    void Delete_Particles_Inside_Objects(const T time);
    template<class T_PARTICLES> void Delete_Particles_Inside_Objects(ARRAY<T_PARTICLES*>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time);
    void Delete_Particles_Inside_Objects(const BLOCK_DYADIC<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time);
    void Update_Fluid_Parameters(const T dt,const T time);
    void Get_Body_Force(ARRAY<TV>& force,const T dt,const T time);
    void Get_Neumann_And_Dirichlet_Boundary_Conditions(const T dt,const T time);
    void Set_Domain_Boundary_Conditions(const T time);
    void Blend_In_External_Velocity(const T dt,const T time);
    void Move_Grid(const T time);
    void Specify_Refinement_For_Cell_With_Thresholds_Helper(const T_CELL& cell,ARRAY<typename T_CELL::REFINE_ACTION>& refine_action,const ARRAY<VECTOR<T,2> >& refinement_thresholds,
        const ARRAY<T>& refinement_variable,const T dt,const T time);
    void Specify_Refinement_For_Cell_With_Levelset_Helper(const T_CELL& cell,ARRAY<typename T_CELL::REFINE_ACTION>& refine_action,const ARRAY<T>& phi,const T band_width_times_small_dx,
        const T dt,const T time);
    template<class OBJECT> void Specify_Refinement_For_Cell_With_Object_Helper(const T_CELL& cell,ARRAY<typename T_CELL::REFINE_ACTION>& refine_action,const OBJECT& object,
        const TRANSFORMATION_MATRIX& world_to_source,const T band_width_times_small_dx,const T dt,const T time);
    void Specify_Refinement(ARRAY<typename T_CELL::REFINE_ACTION>& refine_action,const ARRAY<T_CELL*>* cells_to_specify_for,const T dt,const T time);
    void Setup_Initial_Refinement(const T initial_time);
    template<class T_PARTICLES> int Total_Number_Of_Particles(const ARRAY<T_PARTICLES*>& particles) const;
    template<class T_PARTICLES> void Write_Particles(const STREAM_TYPE stream_type,const ARRAY<T_PARTICLES*>& particles,const std::string& output_directory,const std::string& prefix,
        const int frame) const;
    template<class T_PARTICLES> void Read_Particles(const STREAM_TYPE stream_type,ARRAY<T_PARTICLES*>& particles,const std::string& output_directory,const std::string& prefix,const int frame);
    void Read_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame);
    void Write_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int first_frame,const int frame) const;
//#####################################################################
};      
}
#endif
#endif

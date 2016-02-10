//#####################################################################
// Copyright 2004-2008, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Avi Robinson-Mosher, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Fourier_Transforms_Calculations/TURBULENCE_1D.h>
#include <PhysBAM_Tools/Fourier_Transforms_Calculations/TURBULENCE_2D.h>
#include <PhysBAM_Tools/Fourier_Transforms_Calculations/TURBULENCE_3D.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Dyadic_Advection/ADVECTION_SEMI_LAGRANGIAN_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Tools/Grids_RLE_Advection/ADVECTION_SEMI_LAGRANGIAN_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Boundaries/BOUNDARY_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/AVERAGING_RLE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_HAMILTON_JACOBI_WENO.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_REFLECTION_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Vectors/VECTOR_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Collisions/GRID_BASED_COLLISION_GEOMETRY_DYADIC.h>
#include <PhysBAM_Geometry/Grids_RLE_Collisions/GRID_BASED_COLLISION_GEOMETRY_RLE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Boundaries/BOUNDARY_MAC_GRID_SOLID_WALL_SLIP.h>
#include <PhysBAM_Dynamics/Boundaries/BOUNDARY_PHI_WATER.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_CALLBACKS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> FLUIDS_PARAMETERS<T_GRID>::
FLUIDS_PARAMETERS(const TYPE type,INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& incompressible_fluid_container)
    :smoke(type==SMOKE),fire(type==FIRE),water(type==WATER),sph(type==SPH),compressible(type==COMPRESSIBLE),quadtree(false),octree(false),
    number_of_ghost_cells(3),cfl((T).9),gravity((T)9.8),gravity_direction(-TV::Axis_Vector(TV::m==1?1:2)),grid(&incompressible_fluid_container.rigid_grid.grid),
    need_destroy_grid(false),maximum_tree_depth(1),
    levelset_refinement_bandwidth((T)6),minimal_air_bandwidth(false),
    phi_boundary_reflection(*new T_BOUNDARY_REFLECTION(VECTOR_UTILITIES::Complement(domain_walls))),phi_boundary_water(*new T_BOUNDARY_PHI_WATER),
    particle_half_bandwidth((T)3),reseeding_frame_rate(20),reinitialize_geometry_frame_rate(1),bias_towards_negative_particles(false),
    use_particle_levelset(true),number_particles_per_cell(64),
    use_removed_positive_particles(false),use_removed_negative_particles(false),
    store_particle_ids(false),reincorporate_removed_particle_velocity(false),removed_particle_mass_scaling((T)1),
    use_sph_for_removed_negative_particles(false),
    normal_flame_speed((T).5),curvature_flame_speed(0),
    fluid_boundary_water(*new T_BOUNDARY_SCALAR),
    boundary_mac_slip(*new T_BOUNDARY_MAC_GRID_SOLID_WALL_SLIP(VECTOR_UTILITIES::Complement(domain_walls))),
    incompressible_tolerance((T)1e-8),incompressible_iterations(20),show_residual(false),cg_restart_iterations(0),
    use_body_force(false),
    density((T)1e3),outside_density(0),density_fuel(1),
    surface_tension(0),variable_surface_tension(false),
    viscosity(0),viscosity_fuel(0),
    variable_viscosity(false),implicit_viscosity(true),use_explicit_part_of_implicit_viscosity(false),implicit_viscosity_iterations(40),implicit_viscosity_tolerance((T)1e-25),
    use_coupled_implicit_viscosity(false),
    use_vorticity_confinement(false),use_vorticity_confinement_fuel(true),use_variable_vorticity_confinement(false),
    confinement_parameter((T).05),confinement_parameter_fuel((T)1.2),
    use_non_zero_divergence(false),
    second_order_cut_cell_method(false),enforce_divergence_free_extrapolation(false),solve_neumann_regions(true),solve_single_cell_neumann_regions(false),evolution_solver_type(krylov_solver_cg),
    kolmogorov((T)0),turbulence_lowest_angular_frequency((T)16),turbulence_update_frame_rate(20),turbulence(*new T_TURBULENCE),
    use_external_velocity(false),
    use_soot(false),use_density(true),use_temperature(true),
    use_fixed_soot_boundary(true),use_fixed_density_boundary(true),use_fixed_temperature_boundary(true),
    soot_advection_order(2),ambient_soot(0),ambient_density(0),ambient_temperature((T)298),
    soot_container(*grid),soot_fuel_container(*grid),density_container(incompressible_fluid_container.density_container),temperature_container(incompressible_fluid_container.temperature_container),
    use_soot_fuel_combustion(false),burn_temperature_threshold((T)2000),burn_rate((T).2),
    soot_fuel_calorific_value((T)54000),
    soot_boundary(0),density_boundary(0),temperature_boundary(0),
    temperature_fuel((T)298),temperature_products((T)3000),
    density_buoyancy_threshold((T).5),density_buoyancy_constant((T).001),temperature_buoyancy_constant((T).001),removed_positive_particle_buoyancy_constant((T).3),
    rho_bottom((T).65),rho_top((T)1),
    use_strain(false),elastic_modulus((T)1000),plasticity_alpha(0),plasticity_gamma(0),adhesion_coefficient(0),adhesion_normal_strain(0),adhesion_half_bandwidth((T)1.5),
    strain_boundary_default(*new T_BOUNDARY_SYMMETRIC_MATRIX),object_friction(0),delete_fluid_inside_objects(false),
    move_grid(false),move_grid_explicitly(false),moving_grid_number_of_cells(number_of_ghost_cells),
    monitor_mass(true),
    write_levelset(true),write_particles(false),write_removed_positive_particles(false),write_removed_negative_particles(false),write_flattened_particles(false),
    write_velocity(true),write_strain(true),write_debug_data(true),write_restart_data(true),write_ghost_values(false),restart_data_write_rate(1),
    semi_lagrangian(*new T_ADVECTION_SEMI_LAGRANGIAN_SCALAR),
    hamilton_jacobi_weno(*new T_ADVECTION_HAMILTON_JACOBI_WENO_SCALAR),simulate(true),
    min_collision_distance_factor((T).1),max_collision_distance_factor((T)1),
    solid_affects_fluid(false),fluid_affects_solid(false),
    thin_shells_refine_near_objects(true),
    use_separation_inside_water(false),separation_velocity_tolerance((T)9.8/24), // dt*gravity where dt=1/24 is based on the length of a frame
    callbacks(0),refine_fmm_initialization_with_iterative_solver(true),modify_wall_tangential_velocities(true),
    collision_bodies_affecting_fluid(new T_GRID_BASED_COLLISION_GEOMETRY(*grid)),need_destroy_collision_bodies_affecting_fluid(true),
    collidable_contour_value(0),collidable_phi_replacement_value((T)1e-5),flood_fill_for_bubbles(false),
    use_maccormack_semi_lagrangian_advection(false),use_maccormack_compute_mask(true),use_maccormack_for_level_set(true),use_maccormack_for_incompressible(true),
    bandwidth_without_maccormack_near_interface(0),mass_conservation(false),mass_conservation_minimum_refinement_depth(1),mass_conservation_maximum_refinement_depth(0),analytic_test(false),
    compressible_boundary(0),compressible_pressure_boundary(0),compressible_eos(0),compressible_conservation_method(0),
    compressible_set_max_time_step(false),compressible_max_time_step((T)1.e8),compressible_spatial_order(3),compressible_rungekutta_order(3),compressible_tolerance((T)1e-8),
    compressible_iterations(0),compressible_timesplit(false),compressible_apply_isobaric_fix(false),compressible_monitor_conservation_error(false),
    compressible_perform_rungekutta_for_implicit_part(false),compressible_use_sound_speed_for_cfl(false),
    compressible_use_sound_speed_based_dt_multiple_for_cfl(false),compressible_multiplication_factor_for_sound_speed_based_dt((T)1.),
    compressible_apply_cavitation_correction(false),compressible_adaptive_time_step(false),compressible_log_extremas(true),
    use_poisson(false),solids_override_source_velocities(false),use_slip(false),
    use_slip_constraints_across_non_occluded_faces(false),use_preconditioner_for_slip_system(true),stokes_flow(false),
    thread_queue(0),number_of_threads(1),use_trapezoid_rule(false)
{
    phi_boundary=&phi_boundary_reflection;
    fluid_boundary=&boundary_mac_slip;
    strain_boundary=&strain_boundary_default;
    if(smoke){
        density=1;use_vorticity_confinement=true;temperature_container.Set_Cooling_Constant(0);}
    else if(fire){
        use_body_force=true;density=(T).2;use_vorticity_confinement=true;confinement_parameter=1;}
    else if(water){
        use_density=false;use_temperature=false;Set_Default_Number_Particles_Per_Cell(TV());
        phi_boundary=&phi_boundary_water;} // override default
    else if(sph){
        use_density=false;use_temperature=false;}
    else if(compressible){
        use_density=false;use_temperature=false;}
    turbulence_grid=Default_Turbulence_Grid();

    domain_walls.Fill(VECTOR<bool,2>::Constant_Vector(true));
    domain_walls(TV::m==1?1:2)(2)=false;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> FLUIDS_PARAMETERS<T_GRID>::
~FLUIDS_PARAMETERS()
{
    if(need_destroy_grid) delete grid;
    if(need_destroy_collision_bodies_affecting_fluid) delete collision_bodies_affecting_fluid;
    delete &phi_boundary_reflection;delete &phi_boundary_water;delete &fluid_boundary_water;delete &boundary_mac_slip;delete &turbulence;
    delete &strain_boundary_default;delete &semi_lagrangian;delete &hamilton_jacobi_weno;
}
//#####################################################################
// Function Initialize_Turbulence
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS<T_GRID>::
Initialize_Turbulence(const T time,const T frame_rate)
{
    turbulence.Set_Lowest_Angular_Frequency(turbulence_lowest_angular_frequency);
    turbulence.Initialize_Grid(turbulence_grid);
    turbulence.Generate_Initial_Turbulence(time,time+(T)turbulence_update_frame_rate/frame_rate);
}
//#####################################################################
// Function Initialize_Domain_Boundary_Conditions
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS<T_GRID>::
Initialize_Domain_Boundary_Conditions()
{
    VECTOR<VECTOR<bool,2>,TV::dimension> domain_open_boundaries=VECTOR_UTILITIES::Complement(domain_walls);
    if(phi_boundary) phi_boundary->Set_Constant_Extrapolation(domain_open_boundaries);
    if(fluid_boundary) fluid_boundary->Set_Constant_Extrapolation(domain_open_boundaries);

    if(soot_boundary) soot_container.Set_Custom_Boundary(*soot_boundary);
    soot_container.Set_Ambient_Density(ambient_soot);
    soot_container.boundary->Set_Fixed_Boundary(use_fixed_soot_boundary,ambient_soot);
    soot_container.Initialize_Domain_Boundary_Conditions(domain_walls);

    if(soot_boundary) soot_fuel_container.Set_Custom_Boundary(*soot_boundary);
    soot_fuel_container.Set_Ambient_Density(ambient_soot);
    soot_fuel_container.boundary->Set_Fixed_Boundary(use_fixed_soot_boundary,ambient_soot);
    soot_fuel_container.Initialize_Domain_Boundary_Conditions(domain_walls);

    if(density_boundary) density_container.Set_Custom_Boundary(*density_boundary);
    density_container.Set_Ambient_Density(ambient_density);
    density_container.boundary->Set_Fixed_Boundary(use_fixed_density_boundary,ambient_density);
    density_container.Initialize_Domain_Boundary_Conditions(domain_walls);

    if(temperature_boundary) temperature_container.Set_Custom_Boundary(*temperature_boundary);
    temperature_container.Set_Ambient_Temperature(ambient_temperature);
    temperature_container.boundary->Set_Fixed_Boundary(use_fixed_temperature_boundary,ambient_temperature);
    temperature_container.Initialize_Domain_Boundary_Conditions(domain_walls);

    if(strain_boundary) strain_boundary->Set_Constant_Extrapolation(domain_open_boundaries);
}
//#####################################################################
// Function Initialize_Soot
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS<T_GRID>::
Initialize_Soot(const T time)
{
    assert(use_soot);
    soot_container.Set_To_Ambient_Density();
    soot_fuel_container.Set_To_Ambient_Density();
}
//#####################################################################
// Function Evolve_Soot
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS<T_GRID>::
Evolve_Soot(const T dt,const T time)
{
    assert(use_soot);
    callbacks->Adjust_Soot_With_Sources(time);
    soot_container.Euler_Step(dt,time,number_of_ghost_cells);
    soot_fuel_container.Euler_Step(dt,time,number_of_ghost_cells);
}
//#####################################################################
// Function Initialize_Density_And_Temperature
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS<T_GRID>::
Initialize_Density_And_Temperature(const T time)
{
    if(use_density) density_container.Set_To_Ambient_Density();
    if(use_temperature) temperature_container.Set_To_Ambient_Temperature();
}
//#####################################################################
// Function Evolve_Density_And_Temperature
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS<T_GRID>::
Evolve_Density_And_Temperature(const T dt,const T time)
{
    if(use_density || use_temperature) callbacks->Adjust_Density_And_Temperature_With_Sources(time);
    if(use_density) density_container.Euler_Step(dt,time,number_of_ghost_cells);
    if(use_temperature) temperature_container.Euler_Step(dt,time,number_of_ghost_cells);
}
//#####################################################################
// Function Log_Parameters 
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS<T_GRID>::
Log_Parameters() const
{
    LOG::SCOPE scope("FLUIDS_PARAMETERS parameters");
    std::stringstream ss;
    ss<<"cfl="<<cfl<<std::endl;
    ss<<"solid_affects_fluid="<<solid_affects_fluid<<std::endl;
    ss<<"fluid_affects_solid="<<fluid_affects_solid<<std::endl;
    // compressible parameters
    ss<<"compressible_set_max_time_step="<<compressible_set_max_time_step<<std::endl;
    ss<<"compressible_max_time_step="<<compressible_max_time_step<<std::endl;
    ss<<"compressible_spatial_order="<<compressible_spatial_order<<std::endl;
    ss<<"compressible_rungekutta_order="<<compressible_rungekutta_order<<std::endl;
    ss<<"compressible_timesplit="<<compressible_timesplit<<std::endl;
    ss<<"compressible_monitor_conservation_error="<<compressible_monitor_conservation_error<<std::endl;
    ss<<"compressible_perform_rungekutta_for_implicit_part="<<compressible_perform_rungekutta_for_implicit_part<<std::endl;
    ss<<"compressible_use_sound_speed_for_cfl="<<compressible_use_sound_speed_for_cfl<<std::endl;
    ss<<"compressible_use_sound_speed_based_dt_multiple_for_cfl="<<compressible_use_sound_speed_based_dt_multiple_for_cfl<<std::endl;
    ss<<"compressible_multiplication_factor_for_sound_speed_based_dt="<<compressible_multiplication_factor_for_sound_speed_based_dt<<std::endl;
    ss<<"compressible_log_extremas="<<compressible_log_extremas<<std::endl;
    ss<<"compressible_apply_isobaric_fix="<<compressible_apply_isobaric_fix<<std::endl;
    LOG::filecout(ss.str());
}
//#####################################################################
template class FLUIDS_PARAMETERS<GRID<VECTOR<float,1> > >;
template class FLUIDS_PARAMETERS<GRID<VECTOR<float,2> > >;
template class FLUIDS_PARAMETERS<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template class FLUIDS_PARAMETERS<QUADTREE_GRID<float> >;
template class FLUIDS_PARAMETERS<OCTREE_GRID<float> >;
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template class FLUIDS_PARAMETERS<RLE_GRID_2D<float> >;
template class FLUIDS_PARAMETERS<RLE_GRID_3D<float> >;
#endif
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FLUIDS_PARAMETERS<GRID<VECTOR<double,1> > >;
template class FLUIDS_PARAMETERS<GRID<VECTOR<double,2> > >;
template class FLUIDS_PARAMETERS<GRID<VECTOR<double,3> > >;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template class FLUIDS_PARAMETERS<QUADTREE_GRID<double> >;
template class FLUIDS_PARAMETERS<OCTREE_GRID<double> >;
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template class FLUIDS_PARAMETERS<RLE_GRID_2D<double> >;
template class FLUIDS_PARAMETERS<RLE_GRID_3D<double> >;
#endif
#endif

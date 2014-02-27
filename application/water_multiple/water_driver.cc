//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parallel_Computation/DOMAIN_ITERATOR_THREADED.h>
#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_THREADED.h>
#include <PhysBAM_Tools/Vectors/VECTOR_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM.h>
#include <PhysBAM_Dynamics/Boundaries/BOUNDARY_PHI_WATER.h>

#include "application/water_multiple/water_driver.h"
#include "application/water_multiple/water_example.h"
#include "application/water_multiple/projection/laplace_solver_wrapper.h"
#include "application/water_multiple/projection/nimbus_pcg_sparse_mpi.h"
#include "application/water_multiple/projection/projection_helper.h"
#include "shared/dbg_modes.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"
#include "stdio.h"
#include "string.h"

using namespace PhysBAM;
namespace{
    template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
    {
        ((WATER_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
    }
};
//#####################################################################
// Initialize
//#####################################################################
template<class TV> WATER_DRIVER<TV>::
    WATER_DRIVER(WATER_EXAMPLE<TV>& example)
:example(example)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> WATER_DRIVER<TV>::
~WATER_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}
//#####################################################################
// Initialize
//#####################################################################
// TODO(quhang): Any memory allocation in this initialization function should be
// controlled by Nimbus.
template<class TV> void WATER_DRIVER<TV>::
Initialize(const nimbus::Job *job,
           const nimbus::DataArray &da,
           const bool set_boundary_conditions)
{
  typedef application::DataConfig DataConfig;
  DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);

  output_number=current_frame;
  if (init_phase) {
    time=example.Time_At_Frame(current_frame);
  }

  for (int i = 1; i <= TV::dimension; i++) {
    example.domain_boundary(i)(1) = true;
    example.domain_boundary(i)(2) = true;
  }

  example.domain_boundary(2)(2)=false;
  example.phi_boundary_water.Set_Velocity_Pointer(example.face_velocities);
  VECTOR<VECTOR<bool,2>,TV::dimension> domain_open_boundaries=VECTOR_UTILITIES::Complement(example.domain_boundary);
  example.phi_boundary=&example.phi_boundary_water;
  example.phi_boundary->Set_Constant_Extrapolation(domain_open_boundaries);
  example.boundary=&example.boundary_scalar;
  example.boundary->Set_Constant_Extrapolation(domain_open_boundaries);

  // Allocates array for levelset/particles/removed particles.  --quhang
  InitializeParticleLevelsetEvolutionHelper(
      example.data_config,
      example.mac_grid,
      &example.particle_levelset_evolution);

  example.particle_levelset_evolution.particle_levelset.Set_Band_Width(6);

  // Allocates array for valid mask and data structures used in projeciton.
  // --quhang
  InitializeIncompressibleProjectionHelper(
      example.data_config,
      example.mac_grid,
      &example.incompressible,
      &example.projection);
  // This initialization function initializes valid mask, which we
  // don't know. And this one also calls projection.Initialize_Grid, which is
  // duplicated.
  // example.incompressible.Initialize_Grids(example.mac_grid);
  // I think the only way to debug whether psi_D and psi_N is
  // passed around is to tune the initialization here. I will break this
  // function and make it possible to tune whether initialze psi or not. I
  // asked, but didn't receive a sure answer. Have to figure out ourselves.
  // example.projection.Initialize_Grid(example.mac_grid);

  example.collision_bodies_affecting_fluid.Initialize_Grids();
  if (example.data_config.GetFlag(DataConfig::VELOCITY)) {
      LOG::Time("Velocity memory allocated.\n");
      example.face_velocities.Resize(example.mac_grid);
  }
  // Initialize the face_velocities_ghoas here, it may not be the best place
  // when we get to water_multiple though. -omidm
  if (example.data_config.GetFlag(DataConfig::VELOCITY_GHOST)) {
      LOG::Time("Ghost Velocity memory allocated.\n");
      example.face_velocities_ghost.Resize(example.incompressible.grid,
                                           example.number_of_ghost_cells, false);
  }

  example.particle_levelset_evolution.Set_Time(time);
  example.particle_levelset_evolution.Set_CFL_Number((T).9);

  example.incompressible.Set_Custom_Advection(example.advection_scalar);
  example.particle_levelset_evolution.Levelset_Advection(1).Set_Custom_Advection(example.advection_scalar);

  example.particle_levelset_evolution.Set_Number_Particles_Per_Cell(16);
  example.particle_levelset_evolution.Set_Levelset_Callbacks(example);
  example.particle_levelset_evolution.Initialize_FMM_Initialization_Iterative_Solver(true);

  example.particle_levelset_evolution.particle_levelset.levelset.Set_Custom_Boundary(*example.phi_boundary);
  example.particle_levelset_evolution.Bias_Towards_Negative_Particles(false);
  example.particle_levelset_evolution.particle_levelset.Use_Removed_Positive_Particles();
  example.particle_levelset_evolution.particle_levelset.Use_Removed_Negative_Particles();
  example.particle_levelset_evolution.particle_levelset.Store_Unique_Particle_Id();
  example.particle_levelset_evolution.Use_Particle_Levelset(true);
  example.particle_levelset_evolution.particle_levelset.levelset.Set_Collision_Body_List(example.collision_bodies_affecting_fluid);
  example.particle_levelset_evolution.particle_levelset.levelset.Set_Face_Velocities_Valid_Mask(&example.incompressible.valid_mask);
  example.particle_levelset_evolution.particle_levelset.Set_Collision_Distance_Factors(.1,1);

  example.incompressible.Set_Custom_Boundary(*example.boundary);
  example.incompressible.projection.elliptic_solver->Set_Relative_Tolerance(1e-8);
  example.incompressible.projection.elliptic_solver->pcg.Set_Maximum_Iterations(40);
  example.incompressible.projection.elliptic_solver->pcg.evolution_solver_type=krylov_solver_cg;
  example.incompressible.projection.elliptic_solver->pcg.cg_restart_iterations=0;
  example.incompressible.projection.elliptic_solver->pcg.Show_Results();
  example.incompressible.projection.collidable_solver->Use_External_Level_Set(example.particle_levelset_evolution.particle_levelset.levelset);
  LAPLACE_COLLIDABLE_UNIFORM<T_GRID>* laplace_solver =
      dynamic_cast<LAPLACE_COLLIDABLE_UNIFORM<T_GRID>* >(
          example.projection.elliptic_solver);
  example.laplace_solver_wrapper.BindLaplaceAndInitialize(laplace_solver);

  if (init_phase) {
    example.collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false);
    example.collision_bodies_affecting_fluid.Rasterize_Objects();
    example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*example.mac_grid.Minimum_Edge_Length(),5);
    example.Initialize_Phi();
    example.Adjust_Phi_With_Sources(time);
    example.particle_levelset_evolution.Make_Signed_Distance();
  }
  else {
    // physbam init
    example.Load_From_Nimbus(job, da, current_frame);
    example.collision_bodies_affecting_fluid.Rasterize_Objects();
    example.collision_bodies_affecting_fluid.
        Compute_Occupied_Blocks(false, (T)2*example.mac_grid.Minimum_Edge_Length(),5);
  }

  example.collision_bodies_affecting_fluid.Compute_Grid_Visibility();
  example.particle_levelset_evolution.Set_Seed(2606);
  if (init_phase) {
    example.particle_levelset_evolution.Seed_Particles(time);
    // Comment: seems that particle should be not updated if loaded from
    // Nimbus.  -quhang
    example.particle_levelset_evolution.Delete_Particles_Outside_Grid();
  }

  //add forces
  example.incompressible.Set_Gravity();
  example.incompressible.Set_Body_Force(true);
  example.incompressible.projection.Use_Non_Zero_Divergence(false);
  // We don't want to deal with the additional burden caused by newmann regions.
  // Just set it to false.  --quhang
  example.incompressible.projection.elliptic_solver->Solve_Neumann_Regions(false);
  example.incompressible.projection.elliptic_solver->solve_single_cell_neumann_regions=false;
  example.incompressible.Use_Explicit_Part_Of_Implicit_Viscosity(false);
  example.incompressible.Set_Maximum_Implicit_Viscosity_Iterations(40);
  example.incompressible.Use_Variable_Vorticity_Confinement(false);
  example.incompressible.Set_Surface_Tension(0);
  example.incompressible.Set_Variable_Surface_Tension(false);
  example.incompressible.Set_Viscosity(0);
  example.incompressible.Set_Variable_Viscosity(false);
  example.incompressible.projection.Set_Density(1e3);

  if (set_boundary_conditions) {
    // TODO(quhang): Needs a better understanding what this block is doing. This
    // one is certainly doing something we haven't taken care of.
    ARRAY<T,TV_INT> exchanged_phi_ghost(example.mac_grid.Domain_Indices(8));
    example.particle_levelset_evolution.particle_levelset.levelset.boundary->Fill_Ghost_Cells(example.mac_grid,example.particle_levelset_evolution.phi,exchanged_phi_ghost,0,time,8);
    example.incompressible.Extrapolate_Velocity_Across_Interface(example.face_velocities,exchanged_phi_ghost,false,3,0,TV());
    example.Set_Boundary_Conditions(time); // get so CFL is correct
  }

  if (init_phase) {
    example.Save_To_Nimbus(job, da, current_frame);
    Write_Output_Files(example.first_frame);
  }
  // Comments by quhang:
  // The collision body should not matter, and the last two parameters should
  // not matter. So just add them in the initialization part.
  example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(true,0,0);
  // Don't know why this statement should be here.
  example.particle_levelset_evolution.Set_Number_Particles_Per_Cell(16);

}

template<class TV> bool WATER_DRIVER<TV>::
InitializeIncompressibleProjectionHelper(
    const application::DataConfig& data_config,
    const GRID<TV>& grid_input,
    INCOMPRESSIBLE_UNIFORM<GRID<TV> >* incompressible,
    PROJECTION_DYNAMICS_UNIFORM<GRID<TV> >* projection) {
  typedef application::DataConfig DataConfig;

  // T_FACE_ARRAYS_BOOL.
  if (data_config.GetFlag(DataConfig::VALID_MASK)) {
    incompressible->valid_mask.Resize(
        grid_input.Domain_Indices(3), true, true, true);
  }
  incompressible->grid = grid_input.Get_MAC_Grid();
  // Strain is not considered.
  assert(incompressible->strain == NULL);
  assert(grid_input.Is_MAC_Grid());
  projection->p_grid = grid_input;
  // Laplace solver is used.
  assert(projection->poisson == NULL);
  assert(projection->laplace != NULL);
  assert(grid_input.DX()==TV() || grid_input.Is_MAC_Grid());
  // projection->laplace->Initialize_Grid(grid_input);
  {
    LAPLACE_COLLIDABLE_UNIFORM<GRID<TV> >* laplace =
        dynamic_cast<LAPLACE_COLLIDABLE_UNIFORM<GRID<TV> >*>(
            projection->laplace);
    laplace->grid = grid_input;
    laplace->second_order_cut_cell_method = true;
    if (data_config.GetFlag(DataConfig::U_INTERFACE)) {
      laplace->u_interface.Resize(grid_input);
    } else {
      laplace->u_interface.Clean_Memory();
    }
    // T_ARRAYS_SCALAR.
    if (data_config.GetFlag(DataConfig::DIVERGENCE)) {
      laplace->f.Resize(grid_input.Domain_Indices(1));
    }
    // T_FACE_ARRAYS_BOOL.
    if (data_config.GetFlag(DataConfig::PSI_N)) {
      laplace->psi_N.Resize(grid_input, 1);
    }
    // T_ARRAYS_BOOL.
    if (data_config.GetFlag(DataConfig::PSI_D)) {
      laplace->psi_D.Resize(grid_input.Domain_Indices(1));
    }
    // T_ARRAYS_INT.
    if (data_config.GetFlag(DataConfig::REGION_COLORS)) {
        laplace->filled_region_colors.Resize(
            grid_input.Domain_Indices(1));
    }
    // assert(laplace->levelset != laplace->levelset_default);
    // Assume uniform region coloring.
    laplace->number_of_regions = 1;
    laplace->filled_region_touches_dirichlet.Resize(1);
    laplace->filled_region_touches_dirichlet(1) = true;
  }
  // Flag use_non_zero_divergence is expected to be false.
  assert(!projection->use_non_zero_divergence);
  projection->divergence.Clean_Memory();
  // T_ARRAYS_SCALAR.
  if (data_config.GetFlag(DataConfig::PRESSURE)) {
    projection->p.Resize(grid_input.Domain_Indices(1));
  }
  // T_ARRAYS_SCALAR.
  if (data_config.GetFlag(DataConfig::PRESSURE_SAVE)) {
    projection->p_save_for_projection.Resize(grid_input.Domain_Indices(1));
  }
  // T_FACE_ARRAYS_SCALAR.
  if (data_config.GetFlag(DataConfig::VELOCITY_SAVE)) {
    projection->face_velocities_save_for_projection.Resize(grid_input);
  }
  // dsd is not considered.
  assert(projection->dsd == NULL);
  return true;
}

template<class TV> bool WATER_DRIVER<TV>::
InitializeParticleLevelsetEvolutionHelper(
    const application::DataConfig& data_config,
    const GRID<TV>& grid_input,
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<TV> >*
    particle_levelset_evolution) {
  typedef application::DataConfig DataConfig;
  PARTICLE_LEVELSET_UNIFORM<GRID<TV> >* particle_levelset =
      &particle_levelset_evolution->particle_levelset;
  assert(grid_input.Is_MAC_Grid());
  particle_levelset_evolution->grid = grid_input;
  // Resizes phi here.
  if (data_config.GetFlag(DataConfig::LEVELSET)) {
    particle_levelset_evolution->phi.Resize(
        grid_input.Domain_Indices(particle_levelset->number_of_ghost_cells));
  }
  // Resizes particles.
  if (data_config.GetFlag(DataConfig::POSITIVE_PARTICLE)) {
    particle_levelset->positive_particles.Resize(
        particle_levelset->levelset.grid.Block_Indices(
            particle_levelset->number_of_ghost_cells));
  }
  if (data_config.GetFlag(DataConfig::NEGATIVE_PARTICLE)) {
    particle_levelset->negative_particles.Resize(
        particle_levelset->levelset.grid.Block_Indices(
            particle_levelset->number_of_ghost_cells));
  }
  particle_levelset->use_removed_positive_particles=true;
  particle_levelset->use_removed_negative_particles=true;
  // Resizes removed particles.
  if (data_config.GetFlag(DataConfig::REMOVED_POSITIVE_PARTICLE)) {
    particle_levelset->removed_positive_particles.Resize(
        particle_levelset->levelset.grid.Block_Indices(
            particle_levelset->number_of_ghost_cells));
  }
  if (data_config.GetFlag(DataConfig::REMOVED_NEGATIVE_PARTICLE)) {
    particle_levelset->removed_negative_particles.Resize(
        particle_levelset->levelset.grid.Block_Indices(
            particle_levelset->number_of_ghost_cells));
  }

  particle_levelset->Set_Minimum_Particle_Radius(
      (T).1*particle_levelset->levelset.grid.Minimum_Edge_Length());
  particle_levelset->Set_Maximum_Particle_Radius(
      (T).5*particle_levelset->levelset.grid.Minimum_Edge_Length());
  if (particle_levelset->half_band_width &&
      particle_levelset->levelset.grid.Minimum_Edge_Length()) {
   particle_levelset->Set_Band_Width(particle_levelset->half_band_width /
                   ((T).5*particle_levelset->levelset.grid.Minimum_Edge_Length()));
  } else {
    particle_levelset->Set_Band_Width();
  }
  particle_levelset->levelset.Initialize_Levelset_Grid_Values();
  if (particle_levelset_evolution->
      levelset_advection.semi_lagrangian_collidable) {
    particle_levelset->levelset.Initialize_Valid_Masks(grid_input);
  }
  return true;
}
//#####################################################################
// Run
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Run(RANGE<TV_INT>& domain,const T dt,const T time)
{
    T_FACE_ARRAYS_SCALAR face_velocities_ghost;face_velocities_ghost.Resize(example.incompressible.grid,3,false);
    example.incompressible.boundary->Fill_Ghost_Cells_Face(example.mac_grid,example.face_velocities,face_velocities_ghost,time+dt,example.number_of_ghost_cells);
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,TV> interpolation;
    PARTICLE_LEVELSET_UNIFORM<GRID<TV> >& pls=example.particle_levelset_evolution.particle_levelset;
    if(pls.use_removed_positive_particles) for(typename GRID<TV>::NODE_ITERATOR iterator(example.mac_grid,domain);iterator.Valid();iterator.Next()) if(pls.removed_positive_particles(iterator.Node_Index())){
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_positive_particles(iterator.Node_Index());
        for(int p=1;p<=particles.array_collection->Size();p++){
            TV X=particles.X(p),V=interpolation.Clamped_To_Array_Face(example.mac_grid,face_velocities_ghost,X);
            if(-pls.levelset.Phi(X)>1.5*particles.radius(p)) V-=-TV::Axis_Vector(2)*.3; // buoyancy
            particles.V(p)=V;}}
    if(pls.use_removed_negative_particles) for(typename GRID<TV>::NODE_ITERATOR iterator(example.mac_grid,domain);iterator.Valid();iterator.Next()) if(pls.removed_negative_particles(iterator.Node_Index())){
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_negative_particles(iterator.Node_Index());
        for(int p=1;p<=particles.array_collection->Size();p++) particles.V(p)+=-TV::Axis_Vector(2)*dt*9.8; // ballistic
        for(int p=1;p<=particles.array_collection->Size();p++) particles.V(p)+=dt*interpolation.Clamped_To_Array_Face(example.mac_grid,example.incompressible.force,particles.X(p));} // external forces
}



// Substep without reseeding and writing to frame.
// Operation on time should be solved carefully. --quhang
template<class TV> void WATER_DRIVER<TV>::
CalculateFrameImpl(const nimbus::Job *job,
                   const nimbus::DataArray &da,
                   const bool set_boundary_conditions,
                   const T dt) {
  // LOG::Time("Calculate Dt");
  // example.particle_levelset_evolution.Set_Number_Particles_Per_Cell(16);
  // T dt=example.cfl*example.incompressible.CFL(example.face_velocities);dt=min(dt,example.particle_levelset_evolution.CFL(false,false));
  // if(time+dt>=target_time){dt=target_time-time;done=true;}
  // else if(time+2*dt>=target_time){dt=.5*(target_time-time);}

  //LOG::Time("Compute Occupied Blocks");
  // T maximum_fluid_speed=example.face_velocities.Maxabs().Max();
  // T max_particle_collision_distance=example.particle_levelset_evolution.particle_levelset.max_collision_distance_factor*example.mac_grid.dX.Max();
  // example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(true,dt*maximum_fluid_speed+2*max_particle_collision_distance+(T).5*example.mac_grid.dX.Max(),10);
  //example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(true,dt*maximum_fluid_speed+2*max_particle_collision_distance+(T).5*example.mac_grid.dX.Max(),10);

  LOG::Time("Adjust Phi With Objects");
  T_FACE_ARRAYS_SCALAR face_velocities_ghost;face_velocities_ghost.Resize(example.incompressible.grid,example.number_of_ghost_cells,false);
  example.incompressible.boundary->Fill_Ghost_Cells_Face(example.mac_grid,example.face_velocities,face_velocities_ghost,time+dt,example.number_of_ghost_cells);

  //Advect Phi 3.6% (Parallelized)
  LOG::Time("Advect Phi");
  example.phi_boundary_water.Use_Extrapolation_Mode(false);
  assert(example.particle_levelset_evolution.runge_kutta_order_levelset == 1);
  example.particle_levelset_evolution.levelset_advection.Euler_Step(
      example.face_velocities,
      dt, time,
      example.particle_levelset_evolution.particle_levelset.
      number_of_ghost_cells);

  example.particle_levelset_evolution.time += dt;
  example.phi_boundary_water.Use_Extrapolation_Mode(true);

  //Advect Particles 12.1% (Parallelized)
  LOG::Time("Step Particles");
  example.particle_levelset_evolution.particle_levelset.Euler_Step_Particles(face_velocities_ghost,dt,time,true,true,false,false);

  //Advect removed particles (Parallelized)
  LOG::Time("Advect Removed Particles");
  RANGE<TV_INT> domain(example.mac_grid.Domain_Indices());domain.max_corner+=TV_INT::All_Ones_Vector();
  DOMAIN_ITERATOR_THREADED_ALPHA<WATER_DRIVER<TV>,TV>(domain,0).template Run<T,T>(*this,&WATER_DRIVER<TV>::Run,dt,time);

  //Advect Velocities 26% (Parallelized)
  LOG::Time("Advect V");
  example.incompressible.advection->Update_Advection_Equation_Face(example.mac_grid,example.face_velocities,face_velocities_ghost,face_velocities_ghost,*example.incompressible.boundary,dt,time);

  //Add Forces 0%
  LOG::Time("Forces");
  example.incompressible.Advance_One_Time_Step_Forces(example.face_velocities,dt,time,true,0,example.number_of_ghost_cells);

  //Modify Levelset with Particles 15% (Parallelizedish)
  LOG::Time("Modify Levelset");
  example.particle_levelset_evolution.particle_levelset.Exchange_Overlap_Particles();
  example.particle_levelset_evolution.Modify_Levelset_And_Particles(&face_velocities_ghost);

  //Adjust Phi 0%
  LOG::Time("Adjust Phi");
  example.Adjust_Phi_With_Sources(time+dt);

  //Delete Particles 12.5 (Parallelized)
  LOG::Time("Delete Particles");
  example.particle_levelset_evolution.Delete_Particles_Outside_Grid();                                                            //0.1%
  example.particle_levelset_evolution.particle_levelset.Delete_Particles_In_Local_Maximum_Phi_Cells(1);                           //4.9%
  example.particle_levelset_evolution.particle_levelset.Delete_Particles_Far_From_Interface(); // uses visibility                 //7.6%
  example.particle_levelset_evolution.particle_levelset.Identify_And_Remove_Escaped_Particles(face_velocities_ghost,1.5,time+dt); //2.4%

  //Reincorporate Particles 0% (Parallelized)
  LOG::Time("Reincorporate Particles");
  if(example.particle_levelset_evolution.particle_levelset.use_removed_positive_particles || example.particle_levelset_evolution.particle_levelset.use_removed_negative_particles)
    example.particle_levelset_evolution.particle_levelset.Reincorporate_Removed_Particles(1,1,0,true);

  //Project 7% (Parallelizedish)
  LOG::SCOPE *scope=0;
  scope=new LOG::SCOPE("Project");
  if (set_boundary_conditions)
    example.Set_Boundary_Conditions(time);
  example.incompressible.Set_Dirichlet_Boundary_Conditions(&example.particle_levelset_evolution.phi,0);
  example.projection.p*=dt;
  // example.projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method();
  example.incompressible.Advance_One_Time_Step_Implicit_Part(example.face_velocities,dt,time,true);
  example.projection.p*=(1/dt);
  example.incompressible.boundary->Apply_Boundary_Condition_Face(example.incompressible.grid,example.face_velocities,time+dt);
  // example.projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method(false);
  delete scope;

  //Extrapolate Velocity 7%
  LOG::Time("Extrapolate Velocity");
  T_ARRAYS_SCALAR exchanged_phi_ghost(example.mac_grid.Domain_Indices(8));
  example.particle_levelset_evolution.particle_levelset.levelset.boundary->Fill_Ghost_Cells(example.mac_grid,example.particle_levelset_evolution.phi,exchanged_phi_ghost,0,time+dt,8);
  example.incompressible.Extrapolate_Velocity_Across_Interface(example.face_velocities,exchanged_phi_ghost,false,3,0,TV());

  // TODO(quhang) Take care of this!
  // time+=dt;

  // Save State.
  example.Save_To_Nimbus(job, da, current_frame+1);
}

// Substep with reseeding and writing to frame.
// Operation on time should be solved carefully. --quhang
template<class TV> void WATER_DRIVER<TV>::
WriteFrameImpl(const nimbus::Job *job,
               const nimbus::DataArray &da,
               const bool set_boundary_conditions,
               const T dt) {
  // Comments(quhang): Notice time has already been increased here.
  // Not sure if the Set_Number_Particles_Per_Cell function should go to
  // initalization.
  // example.particle_levelset_evolution.Set_Number_Particles_Per_Cell(16);

  //Reseed
  LOG::Time("Reseed");
  example.particle_levelset_evolution.Reseed_Particles(time);
  example.particle_levelset_evolution.Delete_Particles_Outside_Grid();

  // I changed the order. --quhang
  Write_Output_Files(++output_number);

  //Save State
  example.Save_To_Nimbus(job, da, current_frame+1);
}

template<class TV> bool WATER_DRIVER<TV>::
ProjectionCalculateBoundaryConditionImpl (
    const nimbus::Job *job,
    const nimbus::DataArray &da,
    T dt) {
  INCOMPRESSIBLE_UNIFORM<GRID<TV> >& incompressible = example.incompressible;
  PROJECTION_DYNAMICS_UNIFORM<GRID<TV> >& projection = example.projection;
  LAPLACE_COLLIDABLE_UNIFORM<T_GRID>& laplace_solver =
      *dynamic_cast<LAPLACE_COLLIDABLE_UNIFORM<T_GRID>* >(
          projection.elliptic_solver);

  // Sets boundary conditions.
  // Local.
  // Read velocity and pressure. Write velocity, pressure, psi_D, and psi_N.
  example.Set_Boundary_Conditions(time);

  // Sets dirichlet boundary conditions in the air.
  // Remote: exchange psi_D and pressure afterwards.
  // Read levelset. Write psi_D and pressure.
  incompressible.Set_Dirichlet_Boundary_Conditions(
      &example.particle_levelset_evolution.phi, 0);

  // Scales pressure.
  // Read/Write pressure.
  projection.p *= dt;

  typedef typename INTERPOLATION_POLICY<GRID<TV> >::FACE_LOOKUP
      T_FACE_LOOKUP;

  // Computes divergence.
  // Local.
  // Read velocity. Write divergence(solver->f).
  projection.Compute_Divergence(
      T_FACE_LOOKUP(example.face_velocities),
      &laplace_solver);

  // Coloring.
  // Local.
  // Read psi_D, psi_N.
  // Write filled_region_colors.
  FillUniformRegionColor(
      laplace_solver.grid,
      laplace_solver.psi_D, laplace_solver.psi_N,
      false, &laplace_solver.filled_region_colors);

  return true;
}

// TAG_PROJECTION
template<class TV> bool WATER_DRIVER<TV>::
ProjectionConstructMatrixImpl (
    const nimbus::Job *job,
    const nimbus::DataArray &da,
    T dt) {
  // Read psi_N, psi_D, filled_region_colors, divergence, pressure.
  // Write A, b, x, tolerance, indexing.
  example.laplace_solver_wrapper.PrepareProjectionInput();
  return true;
}

/*
template<class TV> bool WATER_DRIVER<TV>::
ProjectionCoreImpl(
    const nimbus::Job *job,
    const nimbus::DataArray &da,
    T dt) {
  // MPI reference version:
  // laplace_mpi->Solve(A, x, b, q, s, r, k, z, tolerance, color);
  // color only used for MPI version.
  // laplace->pcg.Solve(A, x, b, q, s, r, k, z, laplace->tolerance);
  PCG_SPARSE<T> pcg_temp;
  pcg_temp.Set_Maximum_Iterations(40);
  pcg_temp.evolution_solver_type=krylov_solver_cg;
  pcg_temp.cg_restart_iterations=0;
  pcg_temp.Show_Results();

  NIMBUS_PCG_SPARSE_MPI pcg_mpi(pcg_temp);
  pcg_mpi.projection_data.matrix_index_to_cell_index =
      &laplace_solver_wrapper.matrix_index_to_cell_index_array(1);
  pcg_mpi.projection_data.cell_index_to_matrix_index =
      &laplace_solver_wrapper.cell_index_to_matrix_index;
  pcg_mpi.projection_data.matrix_a = &laplace_solver_wrapper.A_array(1);
  pcg_mpi.projection_data.vector_b = &laplace_solver_wrapper.b_array(1);
  pcg_mpi.projection_data.vector_x = &laplace_solver_wrapper.x;
  pcg_mpi.projection_data.local_tolerance = laplace_solver.tolerance;
  pcg_mpi.Initialize();
  pcg_mpi.CommunicateConfig();
  pcg_mpi.Parallel_Solve();

  return true;
}
*/

template<class TV> bool WATER_DRIVER<TV>::
ProjectionWrapupImpl(
    const nimbus::Job *job,
    const nimbus::DataArray &da,
    T dt) {
  // Read matrix_index_to_cell_index and x. Write u.
  example.laplace_solver_wrapper.TransformResult();

  // Applies pressure.
  // Local.
  // Read pressure(u/p), levelset, psi_D, psi_N, u_interface, velocity.
  // Write velocity.
  example.projection.Apply_Pressure(example.face_velocities, dt, time);

  // Scales pressure.
  // Read/Write pressure.
  example.projection.p *= (1/dt);
  return true;
}

/*
template<class TV> bool WATER_DRIVER<TV>::
ProjectionImpl(
    const nimbus::Job *job,
    const nimbus::DataArray &da,
    T dt) {
  ProjectionCalculateBoundaryConditionImpl(job, da, dt);
  ProjectionCoreImpl(job, da, dt);
  ProjectionWrapupImpl(job, da, dt);
  return true;
}
*/

/*
template<class TV> bool WATER_DRIVER<TV>::
ProjectionImpl(const nimbus::Job *job,
               const nimbus::DataArray &da,
               T dt) {
  INCOMPRESSIBLE_UNIFORM<GRID<TV> >& incompressible = example.incompressible;
  PROJECTION_DYNAMICS_UNIFORM<GRID<TV> >& projection = example.projection;
  LAPLACE_COLLIDABLE_UNIFORM<T_GRID>& laplace_solver =
      *dynamic_cast<LAPLACE_COLLIDABLE_UNIFORM<T_GRID>* >(
          projection.elliptic_solver);
  LaplaceSolverWrapper laplace_solver_wrapper(&laplace_solver);

  // [START] Code before entering INCOMPRESSIBLE.

  // Sets boundary conditions.
  // Local.
  // Read velocity and pressure. Write velocity, pressure, psi_D, and psi_N.
  example.Set_Boundary_Conditions(time);

  // Sets dirichlet boundary conditions in the air.
  // Remote: exchange psi_D and pressure afterwards.
  // Read levelset. Write psi_D and pressure.
  incompressible.Set_Dirichlet_Boundary_Conditions(
      &example.particle_levelset_evolution.phi, 0);
  // Scales pressure.
  // Read/Write pressure.
  projection.p *= dt;

  // [END] Code before entering INCOMPRESSIBLE.


  // [START] Code in INCOMPRESSIBLE:
  //     example.incompressible.Advance_One_Time_Step_Implicit_Part(
  //       example.face_velocities, dt, time, true);

  // Averaging common face of face_velocities.
  // incompressible.boundary->Apply_Boundary_Condition_Face(
  //    incompressible.projection.p_grid,
  //    example.face_velocities,
  //    time + dt);

  // [END] Code in INCOMPRESSIBLE.


  // [START] Code in PROJECTION:
  //    projection.Make_Divergence_Free(face_velocities, dt, time);

  typedef typename INTERPOLATION_POLICY<GRID<TV> >::FACE_LOOKUP
      T_FACE_LOOKUP;
  // Computes divergence.
  // Local.
  // Read velocity. Write divergence(solver->f).
  projection.Compute_Divergence(
      T_FACE_LOOKUP(example.face_velocities),
      &laplace_solver);

  // Coloring.
  // Local.
  // Read psi_D, psi_N.
  // Write filled_region_colors.
  FillUniformRegionColor(
      laplace_solver.grid,
      laplace_solver.psi_D, laplace_solver.psi_N,
      false, &laplace_solver.filled_region_colors);

  // projection.elliptic_solver->Solve(time,true);
  // Read psi_N, psi_D, filled_region_colors, divergence, pressure.
  // Write A, b, x, tolerance, indexing.
  laplace_solver_wrapper.PrepareProjectionInput();

  // MPI reference version:
  // laplace_mpi->Solve(A, x, b, q, s, r, k, z, tolerance, color);
  // color only used for MPI version.
  // laplace->pcg.Solve(A, x, b, q, s, r, k, z, laplace->tolerance);
  NIMBUS_PCG_SPARSE_MPI pcg_mpi(laplace_solver.pcg);
  pcg_mpi.projection_data.matrix_index_to_cell_index =
      &laplace_solver_wrapper.matrix_index_to_cell_index_array(1);
  pcg_mpi.projection_data.cell_index_to_matrix_index =
      &laplace_solver_wrapper.cell_index_to_matrix_index;
  pcg_mpi.projection_data.matrix_a = &laplace_solver_wrapper.A_array(1);
  pcg_mpi.projection_data.vector_b = &laplace_solver_wrapper.b_array(1);
  pcg_mpi.projection_data.vector_x = &laplace_solver_wrapper.x;
  pcg_mpi.projection_data.local_tolerance = laplace_solver.tolerance;
  pcg_mpi.Initialize();
  pcg_mpi.CommunicateConfig();
  pcg_mpi.Parallel_Solve();

  // Read matrix_index_to_cell_index and x. Write u.
  laplace_solver_wrapper.TransformResult();

  // Applies pressure.
  // Not clear.
  // pressure, psi_D, psi_N, u_interface is needed.
  projection.Apply_Pressure(example.face_velocities, dt, time);

  // [END] Code in PROJECTION.

  // Scales pressure.
  // Read/Write pressure.
  projection.p *= (1/dt);

  // incompressible.boundary->Apply_Boundary_Condition_Face(
  //    incompressible.grid,
  //    example.face_velocities,
  //    time + dt);

  return true;
}
*/

template<class TV> bool WATER_DRIVER<TV>::
ExtrapolationImpl (const nimbus::Job *job,
                 const nimbus::DataArray &da,
                 T dt) {
  T_ARRAYS_SCALAR exchanged_phi_ghost(example.mac_grid.Domain_Indices(8));
  example.particle_levelset_evolution.particle_levelset.levelset.boundary->Fill_Ghost_Cells(
      example.mac_grid,
      example.particle_levelset_evolution.phi,
      exchanged_phi_ghost,
      0, time+dt, 8);
  example.incompressible.Extrapolate_Velocity_Across_Interface(
      example.face_velocities,
      exchanged_phi_ghost,
      false, 3, 0, TV());

  return true;
}

template<class TV> bool WATER_DRIVER<TV>::
AdjustPhiWithObjectsImpl (const nimbus::Job *job,
                          const nimbus::DataArray &da,
                          T dt) {
  LOG::Time("Adjust Phi With Objects");
  example.incompressible.boundary->Fill_Ghost_Cells_Face(
      example.mac_grid, example.face_velocities, example.face_velocities_ghost,
      time + dt, example.number_of_ghost_cells);

  // Save State.
  example.Save_To_Nimbus(job, da, current_frame + 1);

  return true;
}

template<class TV> bool WATER_DRIVER<TV>::
ExtrapolatePhiImpl(const nimbus::Job *job,
                   const nimbus::DataArray &da,
                   T dt) {
  example.phi_boundary_water.Use_Extrapolation_Mode(false);
  assert(example.particle_levelset_evolution.runge_kutta_order_levelset == 1);
  {
    typedef typename LEVELSET_ADVECTION_POLICY<GRID<TV> >
        ::FAST_LEVELSET_ADVECTION_T T_FAST_LEVELSET_ADVECTION;
    typedef typename LEVELSET_POLICY<GRID<TV> >
        ::FAST_LEVELSET_T T_FAST_LEVELSET;
    T_FAST_LEVELSET_ADVECTION& levelset_advection =
        example.particle_levelset_evolution.levelset_advection;
    GRID<TV>& grid = ((T_FAST_LEVELSET*)levelset_advection.levelset)->grid;
    T_ARRAYS_SCALAR& phi = ((T_FAST_LEVELSET*)levelset_advection.levelset)->phi;
    assert(grid.Is_MAC_Grid() && levelset_advection.advection);
    T_ARRAYS_SCALAR phi_ghost(grid.Domain_Indices(
            example.particle_levelset_evolution.particle_levelset.
            number_of_ghost_cells));
    ((T_FAST_LEVELSET*)levelset_advection.levelset)->boundary->Fill_Ghost_Cells(
        grid, phi, phi_ghost, dt, time,
        example.particle_levelset_evolution.particle_levelset.number_of_ghost_cells);
    T_ARRAYS_SCALAR::Copy(phi_ghost, phi);
  }
  example.phi_boundary_water.Use_Extrapolation_Mode(true);

  // Save State.
  example.Save_To_Nimbus(job, da, current_frame + 1);

  return true;
}

template<class TV> bool WATER_DRIVER<TV>::
AdvectPhiImpl(const nimbus::Job *job,
              const nimbus::DataArray &da,
              T dt) {
  //Advect Phi 3.6% (Parallelized)
  LOG::Time("Advect Phi");
  example.phi_boundary_water.Use_Extrapolation_Mode(false);
  assert(example.particle_levelset_evolution.runge_kutta_order_levelset == 1);
  // I wrote and tested the following code, which broke levelset advection
  // into function calls. Because extrapolation and MPI calls are used
  // implicitly in this function call, I think we cannot get rid of it.
  {
    typedef typename LEVELSET_ADVECTION_POLICY<GRID<TV> >
        ::FAST_LEVELSET_ADVECTION_T T_FAST_LEVELSET_ADVECTION;
    typedef typename LEVELSET_POLICY<GRID<TV> >
        ::FAST_LEVELSET_T T_FAST_LEVELSET;
    T_FAST_LEVELSET_ADVECTION& levelset_advection =
        example.particle_levelset_evolution.levelset_advection;
    GRID<TV>& grid = ((T_FAST_LEVELSET*)levelset_advection.levelset)->grid;
    T_ARRAYS_SCALAR& phi = ((T_FAST_LEVELSET*)levelset_advection.levelset)->phi;
    assert(grid.Is_MAC_Grid() && levelset_advection.advection);
    T_ARRAYS_SCALAR phi_ghost(grid.Domain_Indices(
            example.particle_levelset_evolution.particle_levelset.
            number_of_ghost_cells));
    T_ARRAYS_SCALAR::Copy(phi, phi_ghost);
    levelset_advection.advection->Update_Advection_Equation_Cell(
        grid, phi, phi_ghost, example.face_velocities,
        *((T_FAST_LEVELSET*)levelset_advection.levelset)->boundary, dt, time);
    ((T_FAST_LEVELSET*)levelset_advection.levelset)->boundary->Apply_Boundary_Condition(
        grid, phi, time + dt);
  }
  /*
  example.particle_levelset_evolution.levelset_advection.Euler_Step(
      example.face_velocities,
      dt, time,
      example.particle_levelset_evolution.particle_levelset.
      number_of_ghost_cells);
  example.phi_boundary_water.Use_Extrapolation_Mode(true);
  */

  // Save State.
  example.Save_To_Nimbus(job, da, current_frame + 1);

  return true;
}

template<class TV> bool WATER_DRIVER<TV>::
StepParticlesImpl(const nimbus::Job *job,
                  const nimbus::DataArray &da,
                  T dt) {
  //Advect Particles 12.1% (Parallelized)
  LOG::Time("Step Particles");
  example.particle_levelset_evolution.particle_levelset.Euler_Step_Particles(
      example.face_velocities_ghost, dt, time, true, true, false, false);

  // Save State.
  example.Save_To_Nimbus(job, da, current_frame + 1);

  return true;
}

template<class TV> bool WATER_DRIVER<TV>::
AdvectRemovedParticlesImpl(const nimbus::Job *job,
                           const nimbus::DataArray &da,
                           T dt) {
  //Advect removed particles (Parallelized)
  LOG::Time("Advect Removed Particles");
  RANGE<TV_INT> domain(example.mac_grid.Domain_Indices());
  domain.max_corner += TV_INT::All_Ones_Vector();
  DOMAIN_ITERATOR_THREADED_ALPHA<WATER_DRIVER<TV>,TV>(domain,0).template Run<T,T>(
      *this, &WATER_DRIVER<TV>::Run, dt, time);

  // Save State.
  example.Save_To_Nimbus(job, da, current_frame + 1);

  return true;
}

template<class TV> bool WATER_DRIVER<TV>::
AdvectVImpl(const nimbus::Job *job,
            const nimbus::DataArray &da,
            T dt) {
  //Advect Velocities 26% (Parallelized)
  LOG::Time("Advect V");
  example.incompressible.advection->Update_Advection_Equation_Face(
      example.mac_grid, example.face_velocities, example.face_velocities_ghost,
      example.face_velocities_ghost, *example.incompressible.boundary, dt, time);

  // Save State.
  example.Save_To_Nimbus(job, da, current_frame + 1);

  return true;
}

template<class TV> bool WATER_DRIVER<TV>::
ApplyForcesImpl(const nimbus::Job *job,
           const nimbus::DataArray &da,
           T dt) {
  //Add Forces 0%
  LOG::Time("Forces");
  example.incompressible.Advance_One_Time_Step_Forces(
      example.face_velocities, dt, time, true, 0, example.number_of_ghost_cells);

  // Save State.
  example.Save_To_Nimbus(job, da, current_frame + 1);

  return true;
}

template<class TV> bool WATER_DRIVER<TV>::
ModifyLevelSetImpl(const nimbus::Job *job,
                   const nimbus::DataArray &da,
                   T dt) {
    LOG::Time("Modify Levelset ...\n");

    // modify levelset
    example.particle_levelset_evolution.particle_levelset.
        Exchange_Overlap_Particles();
    example.particle_levelset_evolution.
        Modify_Levelset_And_Particles(&example.face_velocities_ghost);

    // save state
    example.Save_To_Nimbus(job, da, current_frame+1);

    return true;
}

template<class TV> bool WATER_DRIVER<TV>::
AdjustPhiImpl(const nimbus::Job *job,
        const nimbus::DataArray &da,
        T dt) {
    LOG::Time("Adjust Phi ...\n");

    // adjust phi with sources
    example.Adjust_Phi_With_Sources(time+dt);

    // Save State.
    example.Save_To_Nimbus(job, da, current_frame + 1);

    return true;
}

template<class TV> bool WATER_DRIVER<TV>::
DeleteParticlesImpl(const nimbus::Job *job,
                    const nimbus::DataArray &da,
                    T dt) {
    LOG::Time("Delete Particles ...\n");

    // delete particles
    example.particle_levelset_evolution.Delete_Particles_Outside_Grid();
    example.particle_levelset_evolution.particle_levelset.
        Delete_Particles_In_Local_Maximum_Phi_Cells(1);
    example.particle_levelset_evolution.particle_levelset.
        Delete_Particles_Far_From_Interface(); // uses visibility
    example.particle_levelset_evolution.particle_levelset.
        Identify_And_Remove_Escaped_Particles(example.face_velocities_ghost,
                1.5,
                time + dt);


    // save state
    example.Save_To_Nimbus(job, da, current_frame+1);

    return true;
}

template<class TV> bool WATER_DRIVER<TV>::
ReincorporateParticlesImpl(const nimbus::Job *job,
                           const nimbus::DataArray &da,
                           T dt) {
    LOG::Time("Reincorporate Removed Particles ...\n");

    // reincorporate removed particles
    if (example.particle_levelset_evolution.particle_levelset.
            use_removed_positive_particles ||
            example.particle_levelset_evolution.particle_levelset.
            use_removed_negative_particles)
        example.particle_levelset_evolution.particle_levelset.
            Reincorporate_Removed_Particles(1, 1, 0, true);

    // save state
    example.Save_To_Nimbus(job, da, current_frame+1);

    return true;
}

//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Write_Substep(const std::string& title,const int substep,const int level)
{
    if(level<=example.write_substeps_level){
        example.frame_title=title;
        std::stringstream ss;ss<<"Writing substep ["<<example.frame_title<<"]: output_number="<<output_number+1<<", time="<<time<<", frame="<<current_frame<<", substep="<<substep<<std::endl;LOG::filecout(ss.str());
        Write_Output_Files(++output_number);example.frame_title="";}
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Write_Output_Files(const int frame)
{
    FILE_UTILITIES::Create_Directory(example.output_directory);
    FILE_UTILITIES::Create_Directory(example.output_directory+STRING_UTILITIES::string_sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+STRING_UTILITIES::string_sprintf("/%d/frame_title",frame),example.frame_title);
    if(frame==example.first_frame) 
        FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/first_frame",frame,"\n");
    example.Write_Output_Files(frame);
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/last_frame",frame,"\n");
}
//#####################################################################
template class WATER_DRIVER<VECTOR<float,3> >;

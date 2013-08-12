//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include "myinclude.h"
#include "water-driver.h"

#include "water-example.h"
#include "advection-velocity.h"
#include "job-context.h"
#include "global-repo.h"
#include <sched.h>
#include <unistd.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

using namespace PhysBAM;
namespace {
void Write_Substep_Helper(void* writer, const std::string& title, int substep,
    int level) {
  ((::PhysBAM::WATER_DRIVER*) writer)->Write_Substep(title, substep, level);
}
}

//#####################################################################
// Initialize
//#####################################################################
// ?? What is substep writer for.
WATER_DRIVER::WATER_DRIVER(WATER_EXAMPLE& example) :
    example(example), kinematic_evolution(example.rigid_geometry_collection,
        true), thread_queue(example.thread_queue) {
  DEBUG_SUBSTEPS::Set_Substep_Writer((void*) this, &Write_Substep_Helper);
}
//#####################################################################
// Initialize
//#####################################################################
WATER_DRIVER::~WATER_DRIVER() {
  DEBUG_SUBSTEPS::Clear_Substep_Writer((void*) this);
}

//#####################################################################
// Initialize
//#####################################################################
// Execute Main Program calls Initialize() and Simulate_To_Frame.
void WATER_DRIVER::Execute_Main_Program() {
  Initialize();
  Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
void WATER_DRIVER::Initialize() {
  DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);

  // setup time
  if (example.restart)
    current_frame = example.restart;
  else
    current_frame = example.first_frame;
  output_number = current_frame;
  time = example.Time_At_Frame(current_frame);

  // initialize collision objects
  kinematic_evolution.Get_Current_Kinematic_Keyframes(1 / example.frame_rate,
      time);
  kinematic_evolution.Set_External_Positions(
      example.rigid_geometry_collection.particles.X,
      example.rigid_geometry_collection.particles.rotation, time);
  kinematic_evolution.Set_External_Velocities(
      example.rigid_geometry_collection.particles.V,
      example.rigid_geometry_collection.particles.angular_velocity, time, time);

  example.phi_boundary_water.Set_Velocity_Pointer(example.face_velocities);

  {
    example.particle_levelset_evolution.Initialize_Domain(example.mac_grid);
    example.particle_levelset_evolution.particle_levelset.Set_Band_Width(6);
    example.incompressible.Initialize_Grids(example.mac_grid);
    example.projection.Initialize_Grid(example.mac_grid);
    example.collision_bodies_affecting_fluid.Initialize_Grids();
  }
  example.face_velocities.Resize(example.mac_grid);

  example.particle_levelset_evolution.Set_Time(time);
  example.particle_levelset_evolution.Set_CFL_Number((T) .9);

  if (example.mpi_grid)
    example.mpi_grid->Initialize(example.domain_boundary);
  example.incompressible.mpi_grid = example.mpi_grid;
  example.projection.elliptic_solver->mpi_grid = example.mpi_grid;
  example.particle_levelset_evolution.particle_levelset.mpi_grid =
      example.mpi_grid;
  if (example.mpi_grid) {
    example.boundary = new BOUNDARY_MPI<GRID<TV> >(example.mpi_grid,
        example.boundary_scalar);
    example.phi_boundary = new BOUNDARY_MPI<GRID<TV> >(example.mpi_grid,
        example.phi_boundary_water);
    example.particle_levelset_evolution.particle_levelset.last_unique_particle_id =
        example.mpi_grid->rank * 30000000;
  } else {
    example.boundary = &example.boundary_scalar;
    example.phi_boundary = &example.phi_boundary_water;
  }

  if (PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM < GRID<TV> > *refine =
      dynamic_cast<PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<GRID<TV> >*>(&example.projection)) {
    refine->boundary = example.boundary;
    refine->phi_boundary = example.phi_boundary;
  }
  example.rigid_geometry_collection.Update_Kinematic_Particles();

  VECTOR<VECTOR<bool, 2>, TV::dimension> domain_open_boundaries =
      VECTOR_UTILITIES::Complement(example.domain_boundary);
  example.phi_boundary->Set_Constant_Extrapolation(domain_open_boundaries);
  example.boundary->Set_Constant_Extrapolation(domain_open_boundaries);
  if (thread_queue) {
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<GRID<TV>, T>* threaded_advection_scalar =
        new ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<GRID<TV>, T>(thread_queue);
    example.particle_levelset_evolution.Levelset_Advection(1).Set_Custom_Advection(
        *threaded_advection_scalar);
    example.incompressible.Set_Custom_Advection(*threaded_advection_scalar);
    example.particle_levelset_evolution.particle_levelset.Set_Thread_Queue(
        thread_queue);
    example.particle_levelset_evolution.particle_levelset.levelset.thread_queue =
        thread_queue;
    if (PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM < GRID<TV> > *refinement =
        dynamic_cast<PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<GRID<TV> >*>(&example.projection))
      refinement->thread_queue = thread_queue;
  } else {
    example.particle_levelset_evolution.Levelset_Advection(1).Set_Custom_Advection(
        example.advection_scalar);
    example.incompressible.Set_Custom_Advection(example.advection_scalar);
  }

  example.particle_levelset_evolution.Set_Number_Particles_Per_Cell(16);
  example.particle_levelset_evolution.Set_Levelset_Callbacks(example);
  example.particle_levelset_evolution.Initialize_FMM_Initialization_Iterative_Solver(
      true);

  example.particle_levelset_evolution.particle_levelset.levelset.Set_Custom_Boundary(
      *example.phi_boundary);
  example.particle_levelset_evolution.Bias_Towards_Negative_Particles(false);
  example.particle_levelset_evolution.particle_levelset.Use_Removed_Positive_Particles();
  example.particle_levelset_evolution.particle_levelset.Use_Removed_Negative_Particles();
  example.particle_levelset_evolution.particle_levelset.Store_Unique_Particle_Id();
  example.particle_levelset_evolution.Use_Particle_Levelset(true);
  example.particle_levelset_evolution.particle_levelset.levelset.Set_Collision_Body_List(
      example.collision_bodies_affecting_fluid);
  example.particle_levelset_evolution.particle_levelset.levelset.Set_Face_Velocities_Valid_Mask(
      &example.incompressible.valid_mask);
  example.particle_levelset_evolution.particle_levelset.Set_Collision_Distance_Factors(
      .1, 1);

  example.incompressible.Set_Custom_Boundary(*example.boundary);
  example.incompressible.projection.elliptic_solver->Set_Relative_Tolerance(
      1e-8);
  example.incompressible.projection.elliptic_solver->pcg.Set_Maximum_Iterations(
      40);
  example.incompressible.projection.elliptic_solver->pcg.evolution_solver_type =
      krylov_solver_cg;
  example.incompressible.projection.elliptic_solver->pcg.cg_restart_iterations =
      0;
  example.incompressible.projection.elliptic_solver->pcg.Show_Results();
  example.incompressible.projection.collidable_solver->Use_External_Level_Set(
      example.particle_levelset_evolution.particle_levelset.levelset);

  if (example.restart) {
    example.Read_Output_Files(example.restart);
    example.collision_bodies_affecting_fluid.Rasterize_Objects();
    example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,
        (T) 2 * example.mac_grid.Minimum_Edge_Length(), 5);
  } // compute grid visibility (for advection later)
  else {
    example.collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(
        false);
    example.collision_bodies_affecting_fluid.Rasterize_Objects();
    example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,
        (T) 2 * example.mac_grid.Minimum_Edge_Length(), 5);
    example.Initialize_Phi();
    example.Adjust_Phi_With_Sources(time);
    example.particle_levelset_evolution.Make_Signed_Distance();
    example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time);
  }

  example.collision_bodies_affecting_fluid.Compute_Grid_Visibility();
  example.particle_levelset_evolution.Set_Seed(2606);
  if (!example.restart)
    example.particle_levelset_evolution.Seed_Particles(time);
  example.particle_levelset_evolution.Delete_Particles_Outside_Grid();

  //add forces
  example.incompressible.Set_Gravity();
  example.incompressible.Set_Body_Force(true);
  example.incompressible.projection.Use_Non_Zero_Divergence(false);
  example.incompressible.projection.elliptic_solver->Solve_Neumann_Regions(
      true);
  example.incompressible.projection.elliptic_solver->solve_single_cell_neumann_regions =
      false;
  example.incompressible.Use_Explicit_Part_Of_Implicit_Viscosity(false);
  example.incompressible.Set_Maximum_Implicit_Viscosity_Iterations(40);
  example.incompressible.Use_Variable_Vorticity_Confinement(false);
  example.incompressible.Set_Surface_Tension(0);
  example.incompressible.Set_Variable_Surface_Tension(false);
  example.incompressible.Set_Viscosity(0);
  example.incompressible.Set_Variable_Viscosity(false);
  example.incompressible.projection.Set_Density(1e3);

  ARRAY<T, TV_INT> exchanged_phi_ghost(example.mac_grid.Domain_Indices(8));
  example.particle_levelset_evolution.particle_levelset.levelset.boundary->Fill_Ghost_Cells(
      example.mac_grid, example.particle_levelset_evolution.phi,
      exchanged_phi_ghost, 0, time, 8);
  example.incompressible.Extrapolate_Velocity_Across_Interface(
      example.face_velocities, exchanged_phi_ghost, false, 3, 0, TV());

  example.Set_Boundary_Conditions(time); // get so CFL is correct
  if (!example.restart)
    Write_Output_Files(example.first_frame);
}
//#####################################################################
// Run
//#####################################################################
void WATER_DRIVER::Run(RANGE<TV_INT>& domain, const T dt, const T time) {
  T_FACE_ARRAYS_SCALAR face_velocities_ghost;
  face_velocities_ghost.Resize(example.incompressible.grid, 3, false);
  example.incompressible.boundary->Fill_Ghost_Cells_Face(example.mac_grid,
      example.face_velocities, face_velocities_ghost, time + dt,
      example.number_of_ghost_cells);
  LINEAR_INTERPOLATION_UNIFORM<GRID<TV>, TV> interpolation;
  PARTICLE_LEVELSET_UNIFORM < GRID<TV> > &pls =
      example.particle_levelset_evolution.particle_levelset;
  if (pls.use_removed_positive_particles)
    for (typename GRID<TV>::NODE_ITERATOR iterator(example.mac_grid, domain);
        iterator.Valid(); iterator.Next())
      if (pls.removed_positive_particles(iterator.Node_Index())) {
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> &particles =
            *pls.removed_positive_particles(iterator.Node_Index());
        for (int p = 1; p <= particles.array_collection->Size(); p++) {
          TV X = particles.X(p), V = interpolation.Clamped_To_Array_Face(
              example.mac_grid, face_velocities_ghost, X);
          if (-pls.levelset.Phi(X) > 1.5 * particles.radius(p))
            V -= -TV::Axis_Vector(2) * .3; // buoyancy
          particles.V(p) = V;
        }
      }
  if (pls.use_removed_negative_particles)
    for (typename GRID<TV>::NODE_ITERATOR iterator(example.mac_grid, domain);
        iterator.Valid(); iterator.Next())
      if (pls.removed_negative_particles(iterator.Node_Index())) {
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> &particles =
            *pls.removed_negative_particles(iterator.Node_Index());
        for (int p = 1; p <= particles.array_collection->Size(); p++)
          particles.V(p) += -TV::Axis_Vector(2) * dt * 9.8; // ballistic
        for (int p = 1; p <= particles.array_collection->Size(); p++)
          particles.V(p) += dt
              * interpolation.Clamped_To_Array_Face(example.mac_grid,
                  example.incompressible.force, particles.X(p));
      } // external forces
}
//#####################################################################
// Advance_To_Target_Time
//#####################################################################
void WATER_DRIVER::Advance_To_Target_Time(const T target_time) {
  struct JobContext *temp_job_context = new struct JobContext;
  temp_job_context->sub_step = 1;
  temp_job_context->dt = -1;
  temp_job_context->target_time = target_time;
  temp_job_context->done_after_this_substep = false;
  // Initiate the jobs. 
}

void WATER_DRIVER::CalculateDtAndSpawnJob(const struct JobContext &job_context) {
  struct JobContext *temp_job_context = job_context.clone();
  float &dt = temp_job_context->dt;
  const float target_time = temp_job_context->target_time;
  bool &done = temp_job_context->done_after_this_substep;
  LOG::Time("Calculate Dt");
  example.particle_levelset_evolution.Set_Number_Particles_Per_Cell(16);
  // Calculate time step.
  dt = example.cfl * example.incompressible.CFL(example.face_velocities);
  dt = min(dt, example.particle_levelset_evolution.CFL(false, false));
  if (time + dt >= target_time) {
    dt = target_time - time;
    done = true;
  } else if (time + 2 * dt >= target_time) {
    dt = .5 * (target_time - time);
  }
  // [TODO] Spawn all the jobs. 
}

void WATER_DRIVER::BeforeAdvectionJob(const struct JobContext &job_context) {
  struct JobContext *temp_job_context = job_context.clone();
  const float dt = temp_job_context->dt;

  LOG::Time("Compute Occupied Blocks");
  // What is the interaction? Why cannot delete?
  T maximum_fluid_speed = example.face_velocities.Maxabs().Max();
  T max_particle_collision_distance =
      example.particle_levelset_evolution.particle_levelset.max_collision_distance_factor
      * example.mac_grid.dX.Max();
  example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(true,
      dt * maximum_fluid_speed + 2 * max_particle_collision_distance
      + (T) .5 * example.mac_grid.dX.Max(), 10);

  // [TODO] Dangerous grammer! Extra care!
  T_FACE_ARRAYS_SCALAR &face_velocities_ghost = *(new T_FACE_ARRAYS_SCALAR);
  g_global_repo->face_velocities_ghost = (void*) &face_velocities_ghost;

  face_velocities_ghost.Resize(example.incompressible.grid,
      example.number_of_ghost_cells, false);
  example.incompressible.boundary->Fill_Ghost_Cells_Face(example.mac_grid,
      example.face_velocities, face_velocities_ghost, time + dt,
      example.number_of_ghost_cells);


  //Advect Phi 3.6% (Parallelized)
  LOG::Time("Advect Phi");
  // Advect the level set. What is extrapolation mode?
  example.phi_boundary_water.Use_Extrapolation_Mode(false);
  example.particle_levelset_evolution.Advance_Levelset(dt);
  example.phi_boundary_water.Use_Extrapolation_Mode(true);

  //Advect Particles 12.1% (Parallelized)
  LOG::Time("Step Particles");
  // Advect particle.
  example.particle_levelset_evolution.particle_levelset.Euler_Step_Particles(
      face_velocities_ghost, dt, time, true, true, false, false);

  //Advect removed particles (Parallelized)
  LOG::Time("Advect Removed Particles");
  example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time);
    RANGE<TV_INT> domain(example.mac_grid.Domain_Indices());
    domain.max_corner += TV_INT::All_Ones_Vector();
    DOMAIN_ITERATOR_THREADED_ALPHA<WATER_DRIVER, TV>(domain, 0).Run<T, T>(*this,
        &WATER_DRIVER::Run, dt, time);

  // Just finish... And pass back context.
}

void WATER_DRIVER::AdvectionJob(const struct JobContext &job_context) {
  struct JobContext *temp_job_context = job_context.clone();
  const float dt = temp_job_context->dt;
  ADVECT_VELOCITY_WORKER.face_velocities_ghost = 
      (T_FACE_ARRAYS_SCALAR*)g_global_repo->face_velocities_ghost;

    LOG::Time("Advect V");

    typename WATER_DRIVER::ADVECT_VELOCITY_WORKER_T &INFO =
      g_global_repo->water_driver->ADVECT_VELOCITY_WORKER;
    ADVECT_VELOCITY_NS::AVERAGING_TYPE averaging;
    ADVECT_VELOCITY_NS::INTERPOLATION_TYPE interpolation;
    INFO.dt = dt;
    INFO.range_all = g_global_repo->water_example->mac_grid.counts;
    INFO.range_x = INFO.range_y = INFO.range_z = INFO.range_re = INFO.range_all;
    INFO.range_re += TV_INT(2, 2, 2);
    INFO.range_x += TV_INT(2, 1, 1);
    INFO.range_y += TV_INT(1, 2, 1);
    INFO.range_z += TV_INT(1, 1, 2);
    TV_INT segment_start(1, 1, 1);
    INFO.fetcher_stop = false;
    while (!INFO.fetcher_stop) {
      ADVECT_VELOCITY_NS::run(*(g_global_repo->water_driver), *(g_global_repo->water_example), INFO, averaging, interpolation, segment_start);
      segment_start(3) += INFO.segment_len;
      if (segment_start(3) >= INFO.range_re(3)) {
        segment_start(3) = 1;
        segment_start(2) += INFO.segment_len;
        if (segment_start(2) >= INFO.range_re(2)) {
          segment_start(2) = 1;
          segment_start(1) += INFO.segment_len;
          if (segment_start(1) >= INFO.range_re(1)) {
            segment_start(1) = 1;
            INFO.fetcher_stop = true;
          }
        }
      }
    } // End while
  // Context not needed...
}


void WATER_DRIVER::AfterAdvectionJob(const struct JobContext &job_context) {
  struct JobContext *temp_job_context = job_context.clone();
  const float dt = temp_job_context->dt;
  T_FACE_ARRAYS_SCALAR &face_velocities_ghost = 
      *((T_FACE_ARRAYS_SCALAR*)g_global_repo->face_velocities_ghost);
    //Add Forces 0%
    LOG::Time("Forces");
    // Add the force for velocity.
    example.incompressible.Advance_One_Time_Step_Forces(example.face_velocities,
        dt, time, true, 0, example.number_of_ghost_cells);

    //Modify Levelset with Particles 15% (Parallelizedish)
    LOG::Time("Modify Levelset");
    example.particle_levelset_evolution.particle_levelset.Exchange_Overlap_Particles();
    example.particle_levelset_evolution.Modify_Levelset_And_Particles(
        &face_velocities_ghost);
    //example.particle_levelset_evolution.Make_Signed_Distance(); //TODO(mlentine) Figure out why this was needed

    //Adjust Phi 0%
    LOG::Time("Adjust Phi");
    example.Adjust_Phi_With_Sources(time + dt);
    example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time + dt);

    //Delete Particles 12.5 (Parallelized)
    LOG::Time("Delete Particles");
    example.particle_levelset_evolution.Delete_Particles_Outside_Grid(); //0.1%
    example.particle_levelset_evolution.particle_levelset.Delete_Particles_In_Local_Maximum_Phi_Cells(
        1);                           //4.9%
    example.particle_levelset_evolution.particle_levelset.Delete_Particles_Far_From_Interface(); // uses visibility                 //7.6%
    example.particle_levelset_evolution.particle_levelset.Identify_And_Remove_Escaped_Particles(
        face_velocities_ghost, 1.5, time + dt); //2.4%

    //Reincorporate Particles 0% (Parallelized)
    LOG::Time("Reincorporate Particles");
    if (example.particle_levelset_evolution.particle_levelset.use_removed_positive_particles
        || example.particle_levelset_evolution.particle_levelset.use_removed_negative_particles)
      example.particle_levelset_evolution.particle_levelset.Reincorporate_Removed_Particles(
          1, 1, 0, true);
    example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time + dt);






    //Project 7% (Parallelizedish)
    //LOG::Time("Project");
    {
      LOG::SCOPE *scope = 0;
      if (!thread_queue)
        scope = new LOG::SCOPE("Project");
      LOG::Time("Boundary Conditions");
      example.Set_Boundary_Conditions(time);
      example.incompressible.Set_Dirichlet_Boundary_Conditions(
          &example.particle_levelset_evolution.phi, 0);
      example.projection.p *= dt;
      {
        LOG::SCOPE *scope = 0;
        if (!thread_queue)
          scope = new LOG::SCOPE("Implicit Part");
        example.projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method();
        example.incompressible.Advance_One_Time_Step_Implicit_Part(
            example.face_velocities, dt, time, true);
        delete scope;
      }
      LOG::Time("Boundary Condition Face");
      example.projection.p *= (1 / dt);
      example.incompressible.boundary->Apply_Boundary_Condition_Face(
          example.incompressible.grid, example.face_velocities, time + dt);
      example.projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method(
          false);
      delete scope;
    }

    //Extrapolate Velocity 7%
    LOG::Time("Extrapolate Velocity");
    T_ARRAYS_SCALAR exchanged_phi_ghost(example.mac_grid.Domain_Indices(8));
    example.particle_levelset_evolution.particle_levelset.levelset.boundary->Fill_Ghost_Cells(
        example.mac_grid, example.particle_levelset_evolution.phi,
        exchanged_phi_ghost, 0, time + dt, 8);
    example.incompressible.Extrapolate_Velocity_Across_Interface(
        example.face_velocities, exchanged_phi_ghost, false, 3, 0, TV());

    time += dt;
    delete &face_velocities_ghost;
}


//#####################################################################
// Simulate_To_Frame
//#####################################################################
void WATER_DRIVER::Simulate_To_Frame(const int frame, const int tid) {
  while (current_frame < frame) {
    LOG::SCOPE scope("FRAME", "Frame %d", current_frame + 1, tid);
    /*kinematic_evolution.Get_Current_Kinematic_Keyframes(
     example.Time_At_Frame(current_frame + 1) - time, time);*/
    Advance_To_Target_Time(example.Time_At_Frame(current_frame + 1));
    LOG::Time("Reseed");
    if ((current_frame - example.first_frame) % 1 == 0) {
      example.particle_levelset_evolution.Reseed_Particles(time);
      example.particle_levelset_evolution.Delete_Particles_Outside_Grid();
    }
    Write_Output_Files(++output_number);
    current_frame++;
  }
}
//#####################################################################
// Function Write_Substep
//#####################################################################
void WATER_DRIVER::Write_Substep(const std::string& title, const int substep,
    const int level) {
  if (level <= example.write_substeps_level) {
    example.frame_title = title;
    std::stringstream ss;
    ss << "Writing substep [" << example.frame_title << "]: output_number="
        << output_number + 1 << ", time=" << time << ", frame=" << current_frame
        << ", substep=" << substep << std::endl;
    LOG::filecout(ss.str());
    Write_Output_Files(++output_number);
    example.frame_title = "";
  }
}
//#####################################################################
// Write_Output_Files
//#####################################################################
void WATER_DRIVER::Write_Output_Files(const int frame) {
  FILE_UTILITIES::Create_Directory(example.output_directory);
  FILE_UTILITIES::Create_Directory(
      example.output_directory
          + STRING_UTILITIES::string_sprintf("/%d", frame));
  FILE_UTILITIES::Create_Directory(example.output_directory + "/common");
  FILE_UTILITIES::Write_To_Text_File(
      example.output_directory
          + STRING_UTILITIES::string_sprintf("/%d/frame_title", frame),
      example.frame_title);
  if (frame == example.first_frame)
    FILE_UTILITIES::Write_To_Text_File(
        example.output_directory + "/common/first_frame", frame, "\n");
  example.Write_Output_Files(frame);
  FILE_UTILITIES::Write_To_Text_File(
      example.output_directory + "/common/last_frame", frame, "\n");
}
//#####################################################################
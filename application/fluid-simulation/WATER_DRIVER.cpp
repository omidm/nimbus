//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "myinclude.h"
#include "WATER_DRIVER.h"
#include "WATER_EXAMPLE.h"
#include "advection-velocity.h"
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
  ((WATER_DRIVER*) writer)->Write_Substep(title, substep, level);
}
static void* advect_velocity_worker_wrapper(void *arg) {
  return ADVECT_VELOCITY_NS::advect_velocity_worker(arg);
}
static void* advect_velocity_fetcher_wrapper(void *arg) {
  //return ADVECT_VELOCITY_NS::advect_velocity_fetcher_network(arg);
  return ADVECT_VELOCITY_NS::advect_velocity_fetcher(arg);
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

  // Assuming at most 64 cores.
  typename ADVECT_VELOCITY_WORKER_T::ThreadInfo tinfo[64];
  typename ADVECT_VELOCITY_WORKER_T::ThreadInfo fetcher_info;

  // Set the affinity of the scheduler to the last core.
  cpu_set_t temp_set;
  CPU_ZERO(&temp_set);
  CPU_SET(ADVECT_VELOCITY_WORKER.worker_num - 1, &temp_set);
  sched_setaffinity(0, sizeof(temp_set), &temp_set);

  pthread_mutex_init(&ADVECT_VELOCITY_WORKER.mutex_buffer, NULL);
  pthread_cond_init(&ADVECT_VELOCITY_WORKER.cond_buffer_any, NULL);
  pthread_cond_init(&ADVECT_VELOCITY_WORKER.cond_buffer_clear, NULL);
  pthread_cond_init(&ADVECT_VELOCITY_WORKER.cond_finish, NULL);
  pthread_mutex_init(&ADVECT_VELOCITY_WORKER.mutex_fetcher, NULL);
  pthread_cond_init(&ADVECT_VELOCITY_WORKER.cond_fetcher_ready, NULL);
  pthread_cond_init(&ADVECT_VELOCITY_WORKER.cond_fetcher_go, NULL);
  ADVECT_VELOCITY_WORKER.task_exec_buffer = new ADVECT_VELOCITY_WORKER_T::TaskList;
  ADVECT_VELOCITY_WORKER.task_exec_buffer->top = 0;
  ADVECT_VELOCITY_WORKER.task_recv_buffer = new ADVECT_VELOCITY_WORKER_T::TaskList;
  ADVECT_VELOCITY_WORKER.task_recv_buffer->top = 0;
  ADVECT_VELOCITY_WORKER.ongoing_worker_num = 0;

  ADVECT_VELOCITY_WORKER.fetcher_refresh = true;
  ADVECT_VELOCITY_WORKER.fetcher_stop = true;
  fetcher_info.assigned_core_num = ADVECT_VELOCITY_WORKER.worker_num;
  fetcher_info.driver = this;
  int t;
  t = pthread_create(&fetcher_info.thread_id, NULL,
      &(advect_velocity_fetcher_wrapper), &fetcher_info);
  if (t != 0) {
    printf("Cannot create threads!!!!\n");
    // [TODO] Handle the error.
  }

  for (int tnum = 0; tnum < ADVECT_VELOCITY_WORKER.worker_num - 1; tnum++) {
    tinfo[tnum].assigned_core_num = tnum;
    tinfo[tnum].driver = this;
    int s;
    s = pthread_create(&tinfo[tnum].thread_id, NULL,
        &(advect_velocity_worker_wrapper), &tinfo[tnum]);
    if (s != 0) {
      printf("Cannot create threads!!!!\n");
      // [TODO] Handle the error.
    }
  }
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
  bool done = false;
  for (int substep = 1; !done; substep++) {

    LOG::Time("Calculate Dt");
// Configure how many particles per cell?
    example.particle_levelset_evolution.Set_Number_Particles_Per_Cell(16);
// Calculate time step.
    T dt = example.cfl * example.incompressible.CFL(example.face_velocities);
    dt = min(dt, example.particle_levelset_evolution.CFL(false, false));
    /*if (example.mpi_grid)
     example.mpi_grid->Synchronize_Dt(dt);*/
    if (time + dt >= target_time) {
      dt = target_time - time;
      done = true;
    } else if (time + 2 * dt >= target_time) {
      dt = .5 * (target_time - time);
    }
    // Added by quhang for advection parallel.
    ADVECT_VELOCITY_WORKER.my_dt = dt;

    LOG::Time("Compute Occupied Blocks");
// What is the interaction? Why cannot delete?
    T maximum_fluid_speed = example.face_velocities.Maxabs().Max();
    T max_particle_collision_distance =
        example.particle_levelset_evolution.particle_levelset.max_collision_distance_factor
            * example.mac_grid.dX.Max();
    example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(true,
        dt * maximum_fluid_speed + 2 * max_particle_collision_distance
            + (T) .5 * example.mac_grid.dX.Max(), 10);

    //LOG::Time("Adjust Phi With Objects");
// Copy the velocity array for further computation? Complete copy?
    T_FACE_ARRAYS_SCALAR face_velocities_ghost;
    face_velocities_ghost.Resize(example.incompressible.grid,
        example.number_of_ghost_cells, false);
    // What is ghost cell????
    example.incompressible.boundary->Fill_Ghost_Cells_Face(example.mac_grid,
        example.face_velocities, face_velocities_ghost, time + dt,
        example.number_of_ghost_cells);
    // Added by quhang for advection parallel.
    ADVECT_VELOCITY_WORKER.my_face_velocities_ghost = &face_velocities_ghost;

    //example.Adjust_Phi_With_Objects(time);

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
    // What is removed particles???
    example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time);
    RANGE<TV_INT> domain(example.mac_grid.Domain_Indices());
    domain.max_corner += TV_INT::All_Ones_Vector();
    DOMAIN_ITERATOR_THREADED_ALPHA<WATER_DRIVER, TV>(domain, 0).Run<T, T>(*this,
        &WATER_DRIVER::Run, dt, time);

    pthread_mutex_lock(&ADVECT_VELOCITY_WORKER.mutex_fetcher);
    ADVECT_VELOCITY_WORKER.fetcher_stop = false;
    ADVECT_VELOCITY_WORKER.fetcher_refresh = true;
    pthread_cond_signal(&ADVECT_VELOCITY_WORKER.cond_fetcher_go);
    pthread_mutex_unlock(&ADVECT_VELOCITY_WORKER.mutex_fetcher);
    sleep(1);
    LOG::Time("Advect V by Hang Qu");

    ADVECT_VELOCITY_WORKER.ongoing_worker_num = 0;
    pthread_mutex_lock(&ADVECT_VELOCITY_WORKER.mutex_buffer);
    while (true) {
      while (ADVECT_VELOCITY_WORKER.task_exec_buffer->top > 0) {
        // Wait for workers to clear the task buffer.
        pthread_cond_wait(&ADVECT_VELOCITY_WORKER.cond_buffer_clear,
            &ADVECT_VELOCITY_WORKER.mutex_buffer);
      }

      pthread_mutex_lock(&ADVECT_VELOCITY_WORKER.mutex_fetcher);
      while (!ADVECT_VELOCITY_WORKER.fetcher_stop
	      && (ADVECT_VELOCITY_WORKER.task_recv_buffer->top == 0)) {
        pthread_cond_wait(&ADVECT_VELOCITY_WORKER.cond_fetcher_ready,
            &ADVECT_VELOCITY_WORKER.mutex_fetcher);
      }
      if (ADVECT_VELOCITY_WORKER.fetcher_stop
	  && ADVECT_VELOCITY_WORKER.task_recv_buffer->top == 0) {
        pthread_mutex_unlock(&ADVECT_VELOCITY_WORKER.mutex_fetcher);
	break;
      }
      std::swap(ADVECT_VELOCITY_WORKER.task_exec_buffer, ADVECT_VELOCITY_WORKER.task_recv_buffer);
      pthread_cond_signal(&ADVECT_VELOCITY_WORKER.cond_fetcher_go);
      pthread_mutex_unlock(&ADVECT_VELOCITY_WORKER.mutex_fetcher);

      pthread_cond_broadcast(&ADVECT_VELOCITY_WORKER.cond_buffer_any);
    }
    // Wait for the last worker to finish.
    while ((ADVECT_VELOCITY_WORKER.ongoing_worker_num != 0)
        || (ADVECT_VELOCITY_WORKER.task_exec_buffer->top != 0)) {
      pthread_cond_wait(&ADVECT_VELOCITY_WORKER.cond_finish,
          &ADVECT_VELOCITY_WORKER.mutex_buffer);
    }
    pthread_mutex_unlock(&ADVECT_VELOCITY_WORKER.mutex_buffer);

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
// What ghost cell is filled?
    example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time + dt);

    //Delete Particles 12.5 (Parallelized)
    LOG::Time("Delete Particles");
// Delete particles. How is the particle stored.
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
  }
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

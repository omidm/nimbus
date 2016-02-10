/*
 * Copyright 2013 Stanford University.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the
 *   distribution.
 *
 * - Neither the name of the copyright holders nor the names of
 *   its contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
 * THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */
/*
 * Author: Hang Qu <quhang@stanford.edu>
 */
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_THREADED.h>
#include "WATER_DRIVER.h"  // NOLINT
#include "WATER_EXAMPLE.h"  // NOLINT

namespace PhysBAM {

// Helper function to print out a grid.
template<typename GRID_T>
void PrintGrid(const GRID_T& grid) {
  printf("Counts %d %d\n", grid.counts(1), grid.counts(2));
  printf("Min %1.0f %1.0f\n", static_cast<float>(grid.domain.min_corner(1)),
      static_cast<float>(grid.domain.min_corner(2)));
  printf("Max %1.0f %1.0f\n", static_cast<float>(grid.domain.max_corner(1)),
      static_cast<float>(grid.domain.max_corner(2)));
}

// Main execution.
template<class TV> void WATER_DRIVER<TV>::Execute_Main_Program() {
  Initialize();
  Simulate_To_Frame(example.last_frame);
}

// Initialization.
template<class TV> void WATER_DRIVER<TV>::Initialize() {
  // Setup time.
  if (example.restart) {
    current_frame = example.restart;
  } else {
    current_frame = example.first_frame;
  }
  time = example.Time_At_Frame(current_frame);

  // MPI setup.
  if (example.mpi_grid) {
    example.mpi_grid->Initialize(example.domain_boundary);
  }
  example.projection.elliptic_solver->mpi_grid = example.mpi_grid;
  if (example.mpi_grid) {
    example.boundary = new BOUNDARY_MPI<GRID<TV> >(
        example.mpi_grid, example.boundary_scalar);
  } else if (example.thread_queue) {
      example.boundary = new BOUNDARY_THREADED<GRID<TV> >(
          *example.thread_queue, example.boundary_scalar);
  } else {
    example.boundary = &example.boundary_scalar;
  }

  // Threading setup.
  example.projection.elliptic_solver->thread_queue = example.thread_queue;

  // Setup grids and velocities.
  example.projection.Initialize_Grid(example.mac_grid);
  example.face_velocities.Resize(example.mac_grid);
  example.levelset.phi.Resize(example.mac_grid.Domain_Indices(3));
  example.Initialize_Fields();

  // Setup projection.
  example.projection.elliptic_solver->Set_Relative_Tolerance(1e-9);
  example.projection.elliptic_solver->pcg.Set_Maximum_Iterations(1000);
  example.projection.elliptic_solver->pcg.evolution_solver_type =
      krylov_solver_cg;
  example.projection.elliptic_solver->pcg.cg_restart_iterations = 40;

  if (example.restart) {
    // Not sure if it will work.
    assert(false);
    example.Read_Output_Files(example.restart);
  }

  // Setup domain boundaries.
  VECTOR<VECTOR<bool, 2>, TV::dimension> constant_extrapolation;
  constant_extrapolation.Fill(VECTOR<bool, 2>::Constant_Vector(true));
  example.boundary->Set_Constant_Extrapolation(constant_extrapolation);
  example.Set_Boundary_Conditions(time);

  // Setup particle_levelset.
  example.particle_levelset.mpi_grid = example.mpi_grid;
  example.particle_levelset.levelset.Set_Levelset_Callbacks(example);
  example.particle_levelset.levelset.Set_Custom_Boundary(*example.boundary);
  example.particle_levelset.Initialize_Particle_Levelset_Grid_Values();
  example.particle_levelset.random.Set_Seed(2008);
  example.particle_levelset.Seed_Particles(0, false);

  if (!example.restart) {
    Write_Output_Files(example.first_frame);
  }
  output_number = example.first_frame;
}

// Advect levelset.
template<class TV> void WATER_DRIVER<TV>::Scalar_Advance(const T dt,
    const T time) {
  // Add water source. Not threaded.
  example.Get_Scalar_Field_Sources(time);
  // Fill phi ghost. MPI.
  ARRAY<T, TV_INT> phi_ghost(example.mac_grid.Domain_Indices(3));
  example.boundary->Set_Fixed_Boundary(true, 0);
  example.boundary->Fill_Ghost_Cells(example.mac_grid, example.levelset.phi,
      phi_ghost, dt, time, 3);  // Communication.
  // Advect V. Not threaded.
  example.advection_scalar.Update_Advection_Equation_Cell(example.mac_grid,
      example.levelset.phi, phi_ghost, example.face_velocities,
      *example.boundary, dt, time);
  example.boundary->Set_Fixed_Boundary(false);
  // Threaded. Communication inside.
  example.levelset.Fast_Marching_Method();
}

// Advect velocity and particle-related functionality.
template<class TV> void WATER_DRIVER<TV>::Convect(const T dt, const T time) {
  ARRAY<T, FACE_INDEX<TV::dimension> > face_velocities_ghost(example.mac_grid,
      3, false);
  example.boundary->Fill_Ghost_Cells_Face(example.mac_grid,
      example.face_velocities, face_velocities_ghost, time, 3);

  example.particle_levelset.Reseed_Particles(time);
  example.particle_levelset.Delete_Particles_Outside_Grid();

  example.particle_levelset.Euler_Step_Particles(face_velocities_ghost, dt,
      time, true, true, false, false);
  example.particle_levelset.Exchange_Overlap_Particles();
  example.particle_levelset.Modify_Levelset_Using_Escaped_Particles(
      &face_velocities_ghost);
  example.particle_levelset.Delete_Particles_Outside_Grid();
  example.particle_levelset.Delete_Particles_In_Local_Maximum_Phi_Cells(1);
  example.particle_levelset.Delete_Particles_Far_From_Interface();
  example.particle_levelset.Identify_And_Remove_Escaped_Particles(
      face_velocities_ghost, 1.5, time + dt);

  example.advection_scalar.Update_Advection_Equation_Face(example.mac_grid,
      example.face_velocities, face_velocities_ghost, face_velocities_ghost,
      *example.boundary, dt, time);
}

// Projection.
template<class TV> void WATER_DRIVER<TV>::Project(const T dt, const T time) {
  example.Set_Boundary_Conditions(time + dt);
  example.projection.p *= dt;
  example.boundary->Apply_Boundary_Condition_Face(example.mac_grid,
      example.face_velocities, time + dt);
  example.projection.Make_Divergence_Free(example.face_velocities, dt, time);
  example.projection.p *= (1 / dt);

  const int ghost_cells = 7;
  T delta = 3 * example.mac_grid.dX.Max();
  ARRAY<T, TV_INT> phi_ghost(example.mac_grid.Domain_Indices(3));
  example.boundary->Fill_Ghost_Cells(example.mac_grid, example.levelset.phi,
      phi_ghost, dt, time, 3);
  for (int axis = 1; axis <= GRID<TV>::dimension; axis++) {
    GRID<TV> face_grid = example.mac_grid.Get_Face_Grid(axis);
    ARRAY<T, TV_INT> phi_face(face_grid.Domain_Indices(), false);
    T_ARRAYS_BASE& face_velocity = example.face_velocities.Component(axis);
    ARRAY<bool, TV_INT> fixed_face(face_grid.Domain_Indices());
    for (typename GRID<TV>::FACE_ITERATOR iterator(example.mac_grid, 0,
        GRID<TV>::WHOLE_REGION, 0, axis); iterator.Valid(); iterator.Next()) {
      TV_INT index = iterator.Face_Index();
      phi_face(index) = (T) .5
          * (phi_ghost(iterator.First_Cell_Index())
              + phi_ghost(iterator.Second_Cell_Index()));
      if (phi_face(index) <= 0)
        fixed_face(index) = true;
      if (phi_face(index) >= delta && !fixed_face(index))
        face_velocity(index) = (T) 0;
    }
    T_EXTRAPOLATION_SCALAR extrapolate(face_grid, phi_face, face_velocity,
        ghost_cells);
    extrapolate.Set_Band_Width(3);
    extrapolate.Set_Custom_Seed_Done(&fixed_face);
    extrapolate.Extrapolate();
  }
}

template<class TV> void WATER_DRIVER<TV>::Advance_To_Target_Time(
    const T target_time) {
  bool done = false;
  for (int substep = 1; !done; substep++) {
    // Cacluate dt. Threaded/MPI.
    T dt = example.cfl * example.CFL(example.face_velocities);
    if (example.mpi_grid)
      example.mpi_grid->Synchronize_Dt(dt);
    if (time + dt >= target_time) {
      dt = target_time - time;
      done = true;
    } else if (time + 2 * dt >= target_time) {
      dt = .5 * (target_time - time);
    }
    Scalar_Advance(dt, time);
    Convect(dt, time);
    for (typename GRID<TV>::FACE_ITERATOR iterator(example.mac_grid);
        iterator.Valid(); iterator.Next()) {
      int axis = iterator.Axis();
      if (axis != 2)
        continue;
      example.face_velocities.Component(axis)(iterator.Face_Index()) -=
          dt * 9.8;
    }
    Project(dt, time);
    time += dt;
  }
}

template<class TV> void WATER_DRIVER<TV>::Simulate_To_Frame(const int frame) {
  while (current_frame < frame) {
    Advance_To_Target_Time(example.Time_At_Frame(current_frame + 1));
    Write_Output_Files(++output_number);
    current_frame++;
  }
}

// Write data.
template<class TV> void WATER_DRIVER<TV>::Write_Output_Files(const int frame) {
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

template class WATER_DRIVER<VECTOR<float, 1> >;
template class WATER_DRIVER<VECTOR<float, 2> >;
template class WATER_DRIVER<VECTOR<float, 3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class WATER_DRIVER<VECTOR<double, 1> >;
template class WATER_DRIVER<VECTOR<double, 2> >;
template class WATER_DRIVER<VECTOR<double, 3> >;
#endif

}  // namespace PhysBAM

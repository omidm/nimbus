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
#include <PhysBAM_Tools/Parallel_Computation/DOMAIN_ITERATOR_THREADED.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_1D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_3D.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Uniform/READ_WRITE_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Dynamics/Read_Write/Particles/READ_WRITE_PARTICLES.h>
#include <pthread.h>
#include <algorithm>
#include <string>
#include "WATER_EXAMPLE.h"  // NOLINT
namespace PhysBAM {

template<class TV> WATER_EXAMPLE<TV>::WATER_EXAMPLE(
    const STREAM_TYPE stream_type_input, int number_of_threads) :
    stream_type(stream_type_input),
    initial_time(0),
    first_frame(0),
    last_frame(100),
    frame_rate(24),
    restart(0),
    write_debug_data(false),
    output_directory("output"),
    cfl(.9),
    mac_grid(TV_INT(), RANGE<TV>::Unit_Box(), true),
    mpi_grid(0),
    thread_queue(number_of_threads > 1 ?
                 new THREAD_QUEUE(number_of_threads) : 0),
    projection(mac_grid, false, false, thread_queue),
    advection_scalar(thread_queue),
    boundary(0),
    levelset(mac_grid, *new ARRAY<T, TV_INT>()),  // Memory leakage here.
    particle_levelset(mac_grid, levelset.phi, 3) {
  for (int i = 1; i <= TV::dimension; i++) {
    domain_boundary(i)(1) = true;
    domain_boundary(i)(2) = true;
  }
  pthread_mutex_init(&lock, 0);
  Initialize_Particles();
  // Initialize_Read_Write_General_Structures();
}

template<class TV> WATER_EXAMPLE<TV>::~WATER_EXAMPLE() {
  if (mpi_grid)
    delete boundary;
}

template<class TV> typename TV::SCALAR WATER_EXAMPLE<TV>::CFL(
    ARRAY<T, FACE_INDEX<TV::dimension> >& face_velocities) {
  T dt = FLT_MAX;
  DOMAIN_ITERATOR_THREADED_ALPHA<WATER_EXAMPLE<TV>, TV>(
      mac_grid.Domain_Indices(), thread_queue).template Run<
      ARRAY<T, FACE_INDEX<TV::dimension> >&, T&>(*this,
      &WATER_EXAMPLE::CFL_Threaded, face_velocities, dt);
  return dt;
}

template<class TV> void WATER_EXAMPLE<TV>::CFL_Threaded(
    RANGE<TV_INT>& domain,  // NOLINT
    ARRAY<T, FACE_INDEX<TV::dimension> >& face_velocities, T& dt) {
  T dt_convection = 0;
  for (typename GRID<TV>::CELL_ITERATOR iterator(mac_grid, domain);
      iterator.Valid(); iterator.Next()) {
    TV_INT cell = iterator.Cell_Index();
    T local_V_norm = 0;
    for (int axis = 1; axis <= GRID<TV>::dimension; axis++)
      local_V_norm += mac_grid.one_over_dX[axis]
          * maxabs(
              face_velocities(axis,
                  mac_grid.First_Face_Index_In_Cell(axis, cell)),
              face_velocities(axis,
                  mac_grid.Second_Face_Index_In_Cell(axis, cell)));
    dt_convection = max(dt_convection, local_V_norm);
  }
  pthread_mutex_lock(&lock);
  dt = min(dt, (T) 1.0 / dt_convection);
  pthread_mutex_unlock(&lock);
}

template<class TV> void WATER_EXAMPLE<TV>::Set_Boundary_Conditions(
    const T time) {
  projection.elliptic_solver->psi_D.Fill(false);
  projection.elliptic_solver->psi_N.Fill(false);
  for (int axis = 1; axis <= TV::dimension; axis++)
    for (int axis_side = 1; axis_side <= 2; axis_side++) {
      int side = 2 * (axis - 1) + axis_side;
      TV_INT interior_cell_offset =
          axis_side == 1 ? TV_INT() : -TV_INT::Axis_Vector(axis);
      TV_INT exterior_cell_offset =
          axis_side == 1 ? -TV_INT::Axis_Vector(axis) : TV_INT();
      TV_INT boundary_face_offset =
          axis_side == 1 ?
              TV_INT::Axis_Vector(axis) : -TV_INT::Axis_Vector(axis);
      if (domain_boundary(axis)(axis_side)) {
        for (typename GRID<TV>::FACE_ITERATOR iterator(mac_grid, 1,
            GRID<TV>::BOUNDARY_REGION, side); iterator.Valid();
            iterator.Next()) {
          TV_INT face = iterator.Face_Index() + boundary_face_offset;
          if (levelset.phi(face + interior_cell_offset) <= 0) {
            if (face_velocities.Component(axis).Valid_Index(face)) {
              projection.elliptic_solver->psi_N.Component(axis)(face) = true;
              face_velocities.Component(axis)(face) = 0;
            }
          } else {
            TV_INT cell = face + exterior_cell_offset;
            projection.elliptic_solver->psi_D(cell) = true;
            projection.p(cell) = 0;
          }
        }
      } else {
        for (typename GRID<TV>::FACE_ITERATOR iterator(mac_grid, 1,
            GRID<TV>::BOUNDARY_REGION, side); iterator.Valid();
            iterator.Next()) {
          TV_INT cell = iterator.Face_Index() + interior_cell_offset;
          projection.elliptic_solver->psi_D(cell) = true;
          projection.p(cell) = 0;
        }
      }
    }
  for (typename GRID<TV>::FACE_ITERATOR iterator(mac_grid); iterator.Valid();
      iterator.Next()) {
    if (time <= 3 && source.Lazy_Inside(iterator.Location())) {
      projection.elliptic_solver->psi_N(iterator.Full_Index()) = true;
      if ((TV::dimension == 2 && iterator.Axis() == 1)
          || (TV::dimension == 3 && iterator.Axis() == 3))
        face_velocities(iterator.Full_Index()) = 1;
      else
        face_velocities(iterator.Full_Index()) = 0;
    }
  }
  for (typename GRID<TV>::CELL_ITERATOR iterator(projection.p_grid);
      iterator.Valid(); iterator.Next())
    if (levelset.phi(iterator.Cell_Index()) > 0) {
      projection.elliptic_solver->psi_D(iterator.Cell_Index()) = true;
      projection.p(iterator.Cell_Index()) = 0;
    }
  if (projection.elliptic_solver->mpi_grid) {
    projection.elliptic_solver->mpi_grid->Exchange_Boundary_Cell_Data(
        projection.elliptic_solver->psi_D, 1, false);
    projection.elliptic_solver->mpi_grid->Exchange_Boundary_Cell_Data(
        projection.p, 1, false);
  }
}

template<class TV> void WATER_EXAMPLE<TV>::Write_Output_Files(const int frame) {
  std::string f = STRING_UTILITIES::string_sprintf("%d", frame);
  FILE_UTILITIES::Write_To_File(stream_type,
      output_directory + "/" + f + "/mac_velocities", face_velocities);
  if (mpi_grid)
    FILE_UTILITIES::Write_To_File(stream_type,
        output_directory + "/common/global_grid", mpi_grid->global_grid);
  FILE_UTILITIES::Write_To_File(stream_type,
      output_directory + "/" + f + "/grid", mac_grid);
  FILE_UTILITIES::Write_To_File(stream_type, output_directory + "/common/grid",
      mac_grid);
  FILE_UTILITIES::Write_To_File(stream_type,
      output_directory + "/" + f + "/levelset", levelset);
  FILE_UTILITIES::Write_To_File(stream_type,
      output_directory + "/" + f + "/positive_particles",
      particle_levelset.positive_particles);
  FILE_UTILITIES::Write_To_File(stream_type,
      output_directory + "/" + f + "/negative_particles",
      particle_levelset.negative_particles);
  if (write_debug_data) {
    FILE_UTILITIES::Write_To_File(stream_type,
        output_directory + "/" + f + "/pressure", projection.p);
    FILE_UTILITIES::Write_To_File(stream_type,
        output_directory + "/" + f + "/psi_N",
        projection.elliptic_solver->psi_N);
    FILE_UTILITIES::Write_To_File(stream_type,
        output_directory + "/" + f + "/psi_D",
        projection.elliptic_solver->psi_D);
  }
}
template<class TV> void WATER_EXAMPLE<TV>::Read_Output_Files(const int frame) {
  // TODO(quhang) Particle IO not added.
  std::string f = STRING_UTILITIES::string_sprintf("%d", frame);
  FILE_UTILITIES::Read_From_File(stream_type,
      output_directory + "/" + f + "/levelset", levelset);
  std::string filename;
  filename = output_directory + "/" + f + "/mac_velocities";
  if (FILE_UTILITIES::File_Exists(filename)) {
    std::stringstream ss;
    ss << "Reading mac_velocities " << filename << std::endl;
    FILE_UTILITIES::Read_From_File(stream_type, filename, face_velocities);
    LOG::filecout(ss.str());
  }
  filename = output_directory + "/" + f + "/pressure";
  if (FILE_UTILITIES::File_Exists(filename)) {
    std::stringstream ss;
    ss << "Reading pressure " << filename << std::endl;
    FILE_UTILITIES::Read_From_File(stream_type, filename, projection.p);
    LOG::filecout(ss.str());
  }
}

template class WATER_EXAMPLE<VECTOR<float, 1> >;
template class WATER_EXAMPLE<VECTOR<float, 2> >;
template class WATER_EXAMPLE<VECTOR<float, 3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class WATER_EXAMPLE<VECTOR<double, 1> >;
template class WATER_EXAMPLE<VECTOR<double, 2> >;
template class WATER_EXAMPLE<VECTOR<double, 3> >;
#endif

}  // namespace PhysBAM

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
#ifndef __WATER_EXAMPLE__  // NOLINT
#define __WATER_EXAMPLE__  // NOLINT
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/PROJECTION_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_UNIFORM.h>
#include <algorithm>
#include <string>
namespace PhysBAM {

template<class TV>
class WATER_EXAMPLE: public LEVELSET_CALLBACKS<GRID<TV> > {
 public:
  typedef typename TV::SCALAR T;
  typedef typename TV::template REBIND<int>::TYPE TV_INT;
  typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;
  enum workaround1 {
    d = TV::m
  };

  STREAM_TYPE stream_type;
  T initial_time;
  int first_frame, last_frame;
  T frame_rate;
  int restart;
  std::string frame_title;
  int write_substeps_level;
  bool write_debug_data;
  std::string output_directory;

  T cfl;

  GRID<TV> mac_grid;
  MPI_UNIFORM_GRID<GRID<TV> > *mpi_grid;
  THREAD_QUEUE* thread_queue;
  PROJECTION_UNIFORM<GRID<TV> > projection;
  ARRAY<T, FACE_INDEX<TV::dimension> > face_velocities;
  ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<GRID<TV>, T,
      AVERAGING_UNIFORM<GRID<TV>, FACE_LOOKUP_UNIFORM<GRID<TV> > >,
      LINEAR_INTERPOLATION_UNIFORM<GRID<TV>, T, FACE_LOOKUP_UNIFORM<GRID<TV> > > > advection_scalar;
  BOUNDARY_UNIFORM<GRID<TV>, T> boundary_scalar;
  BOUNDARY_UNIFORM<GRID<TV>, T> *boundary;
  T_LEVELSET levelset;
  VECTOR<VECTOR<bool, 2>, TV::dimension> domain_boundary;
  RANGE<TV> source;
  pthread_mutex_t lock;

  PARTICLE_LEVELSET_UNIFORM<GRID<TV> > particle_levelset;

  WATER_EXAMPLE(const STREAM_TYPE stream_type_input, int refine = 0);
  virtual ~WATER_EXAMPLE();

  T CFL(ARRAY<T, FACE_INDEX<TV::dimension> >& face_velocities);  // NOLINT
  void CFL_Threaded(RANGE<TV_INT>& domain,  // NOLINT
      ARRAY<T, FACE_INDEX<TV::dimension> >& face_velocities, T& dt);

  T Time_At_Frame(const int frame) const {
    return initial_time + (frame - first_frame) / frame_rate;
  }

  void Initialize_Grid(TV_INT counts, RANGE<TV> domain) {
    mac_grid.Initialize(counts, domain, true);
  }

  void Initialize_Fields() {
    for (typename GRID<TV>::FACE_ITERATOR iterator(mac_grid); iterator.Valid();
        iterator.Next())
      face_velocities(iterator.Full_Index()) = 0;
    for (typename GRID<TV>::CELL_ITERATOR iterator(mac_grid); iterator.Valid();
        iterator.Next())
      levelset.phi(iterator.Cell_Index()) = iterator.Location()(2)
          - mac_grid.dX(2) * 0;
  }

  void Get_Scalar_Field_Sources(const T time) {
    if (time > 3)
      return;
    for (typename GRID<TV>::CELL_ITERATOR iterator(mac_grid); iterator.Valid();
        iterator.Next()) {
      TV_INT index = iterator.Cell_Index();
      T distance = abs((iterator.Location() - source.min_corner).Min());
      distance = min(distance,
          abs((iterator.Location() - source.max_corner).Min()));
      T phi = 0;
      if (source.Lazy_Inside(iterator.Location()))
        phi = -1 * distance;
      else
        phi = distance;
      levelset.phi(index) = min(levelset.phi(index), phi);
    }
  }

  virtual void Write_Output_Files(const int frame);
  virtual void Read_Output_Files(const int frame);
  virtual void Set_Boundary_Conditions(const T time);
};
}  // namespace PhysBAM
#endif  // NOLINT

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
 * This class serves as the bridge between WATER_DRIVER and PROJECTION_DRIVER.
 * It is a wrapper around LAPLACE_COLLIDABLE_UNIFORM so that it can extract
 * necessary data from WATER_DRIVER/WATER_EXAMPLE and calculates all the data
 * needed for PROJECTION_DRIVER. After this point, PROJECTION_DRIVER gets all
 * its data to do projection and WATER_DRIVER is no more needed.
 *
 * Application first binds LAPLACE_COLLIDABLE_UNIFORM of WATER_EXAMPLE to this
 * class, and then calls PrePareProjectionInput. After that, all the data needed
 * for projection is stored in this class.
 *
 * This class is written based on original PhysBAM implementation but is fully
 * rewritten.
 *
 * Author: Hang Qu <quhang@stanford.edu>
 */

#ifndef NIMBUS_APPLICATION_WATER_MULTIPLE_PROJECTION_LAPLACE_SOLVER_WRAPPER_H_
#define NIMBUS_APPLICATION_WATER_MULTIPLE_PROJECTION_LAPLACE_SOLVER_WRAPPER_H_

#include "applications/physbam/water//app_utils.h"

namespace PhysBAM {

typedef application::T T;
typedef application::TV TV;
typedef application::TV_INT TV_INT;

typedef GRID<TV> T_GRID;

class LaplaceSolverWrapper {
 public:
  LaplaceSolverWrapper() {}
  void BindLaplaceAndInitialize(
      LAPLACE_COLLIDABLE_UNIFORM<GRID<TV> >* laplace_input);
  ~LaplaceSolverWrapper() {}

  int local_n, interior_n;
  // region_id -> matrix_id -> (dim_t, dim_t, dim_t)
  ARRAY<ARRAY<TV_INT> > matrix_index_to_cell_index_array;
  // (dim_t, dim_t, dim_t) -> matrix_id
  ARRAY<int, TV_INT> cell_index_to_matrix_index;
  // region_id -> matrix
  ARRAY<SPARSE_MATRIX_FLAT_NXN<T> > A_array;
  // region_id -> vector
  ARRAY<VECTOR_ND<T> > b_array;
  VECTOR_ND<T> x;

  // Input refers to A, x, b, indexing, tolerance.
  void PrepareProjectionInput();
 private:
  LAPLACE_COLLIDABLE_UNIFORM<GRID<TV> >* laplace;
};

}  // namespace PhysBAM

#endif  // NIMBUS_APPLICATION_WATER_MULTIPLE_PROJECTION_LAPLACE_SOLVER_WRAPPER_H_

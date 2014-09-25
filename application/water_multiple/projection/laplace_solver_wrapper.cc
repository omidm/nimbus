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
 * Author: Hang Qu<quhang@stanford.edu>
 */

#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>

#include "application/water_multiple/app_utils.h"

#include "application/water_multiple/projection/data_sparse_matrix.h"
#include "application/water_multiple/projection/data_raw_array_m2c.h"
#include "application/water_multiple/projection/data_raw_grid_array.h"
#include "application/water_multiple/projection/data_raw_vector_nd.h"
#include "application/water_multiple/projection/projection_helper.h"

#include "application/water_multiple/projection/laplace_solver_wrapper.h"

namespace PhysBAM {

void LaplaceSolverWrapper::BindLaplaceAndInitialize(
    LAPLACE_COLLIDABLE_UNIFORM<GRID<TV> >* laplace_input) {
  laplace = laplace_input;
  const int number_of_regions = 1;
  matrix_index_to_cell_index_array.Resize(number_of_regions);
  cell_index_to_matrix_index.Resize(laplace->grid.Domain_Indices(1));
  // Intialize array since it may contain old data/ garbage when read from
  // cache (added by Chinmayee)
  ARRAYS_COMPUTATIONS::Fill(cell_index_to_matrix_index.array, 0);
  cell_index_to_matrix_index.hash_code = 0;
  A_array.Resize(number_of_regions);
  b_array.Resize(number_of_regions);
}

void LaplaceSolverWrapper::PrepareProjectionInput() {
  const int number_of_regions = laplace->number_of_regions;
  assert(number_of_regions == 1);
  // region_id -> sum
  ARRAY<int, VECTOR<int, 1> > filled_region_cell_count(-1, number_of_regions);

  // Count the cells in each region.
  for (typename T_GRID::CELL_ITERATOR iterator(laplace->grid, 1);
       iterator.Valid();
       iterator.Next()) {
    /*
    // TODO(quhang) I am not sure whether the areas specified in this if
    // statement are handled well.
    if ((iterator.Cell_Index().x <= 0)
       + (iterator.Cell_Index().x > laplace->grid.counts.x)
       + (iterator.Cell_Index().y <= 0)
       + (iterator.Cell_Index().y > laplace->grid.counts.y)
       + (iterator.Cell_Index().z <= 0)
       + (iterator.Cell_Index().z > laplace->grid.counts.z) >= 2) {
      laplace->filled_region_colors(iterator.Cell_Index()) = -1;
    }
    */
    filled_region_cell_count(
        laplace->filled_region_colors(iterator.Cell_Index()))++;
  }

  // Assume only one color.
  const int color = 1;

  matrix_index_to_cell_index_array(color).Resize(
      filled_region_cell_count(color));

  // int temp1 = filled_region_cell_count(-1);
  // int temp2 = filled_region_cell_count(1);

  // Reusing this array in order to make the indirection arrays.
  filled_region_cell_count.Fill(0);

  // MPI reference version.
  // laplace_mpi->Find_Matrix_Indices(filled_region_cell_count,
  //                                  cell_index_to_matrix_index,
  //                                  matrix_index_to_cell_index_array);

  FindMatrixIndices(
      laplace->grid,
      laplace->filled_region_colors,
      &filled_region_cell_count,
      &cell_index_to_matrix_index,
      &matrix_index_to_cell_index_array(1),
      &local_n,
      &interior_n);

  RANGE<TV_INT> domain = laplace->grid.Domain_Indices(1);

  if (interior_n != 0) {
    // Construct both A and b.
    laplace->Find_A(
        domain, A_array, b_array,
        filled_region_cell_count, cell_index_to_matrix_index);

    laplace->pcg.Enforce_Compatibility(
        !laplace->filled_region_touches_dirichlet(color) &&
        laplace->enforce_compatibility);

    SPARSE_MATRIX_FLAT_NXN<T>& A = A_array(color);
    VECTOR_ND<T>& b = b_array(color);
    A.Negate();
    b *= (T) -1;

    laplace->Find_Tolerance(b);
    A.Initialize_Diagonal_Index();

  } else {
    laplace->tolerance = 0;
    dbg(APP_LOG, "No water in this region!\n");
  }
}

}  // namespace PhysBAM

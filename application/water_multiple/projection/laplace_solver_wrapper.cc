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

#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/projection/nimbus_pcg_sparse_mpi.h"
#include "application/water_multiple/projection/laplace_solver_wrapper.h"

namespace PhysBAM {

void LaplaceSolverWrapper::PrepareProjectionInput() {
  const int number_of_regions = laplace->number_of_regions;
  assert(number_of_regions == 1);
  // region_id -> sum
  ARRAY<int, VECTOR<int, 1> > filled_region_cell_count(-1, number_of_regions);

  // Count the cells in each region.
  for (typename T_GRID::CELL_ITERATOR iterator(laplace->grid, 1);
       iterator.Valid();
       iterator.Next()) {
    filled_region_cell_count(
        laplace->filled_region_colors(iterator.Cell_Index()))++;
  }

  // Assume only one color.
  const int color = 1;

  matrix_index_to_cell_index_array(color).Resize(
      filled_region_cell_count(color));

  // Reusing this array in order to make the indirection arrays.
  filled_region_cell_count.Fill(0);

  // MPI reference version.
  // laplace_mpi->Find_Matrix_Indices(filled_region_cell_count,
  //                                  cell_index_to_matrix_index,
  //                                  matrix_index_to_cell_index_array);
  laplace->Compute_Matrix_Indices(
      laplace->grid.Domain_Indices(1),
      filled_region_cell_count,
      matrix_index_to_cell_index_array,
      cell_index_to_matrix_index);

  RANGE<TV_INT> domain = laplace->grid.Domain_Indices(1);
  // Construct both A and b.
  laplace->Find_A(
      domain, A_array, b_array,
      filled_region_cell_count, cell_index_to_matrix_index);

  laplace->pcg.Enforce_Compatibility(
     !laplace->filled_region_touches_dirichlet(color) &&
     laplace->enforce_compatibility);

  ARRAY<TV_INT>& matrix_index_to_cell_index =
      matrix_index_to_cell_index_array(color);
  SPARSE_MATRIX_FLAT_NXN<T>& A = A_array(color);
  VECTOR_ND<T>& b = b_array(color);

  int number_of_unknowns = matrix_index_to_cell_index.m;
  A.Negate();
  b *= (T) -1;
  x.Resize(number_of_unknowns);
  VECTOR_ND<T> x(number_of_unknowns);
  for (int i = 1; i <= number_of_unknowns; i++) {
    x(i) = laplace->u(matrix_index_to_cell_index(i));
  }

  laplace->Find_Tolerance(b);
}

void LaplaceSolverWrapper::TransformResult() {
  // Assume only one color.
  const int color = 1;
  ARRAY<TV_INT>& matrix_index_to_cell_index =
      matrix_index_to_cell_index_array(color);
  int number_of_unknowns = matrix_index_to_cell_index.m;
  for (int i = 1; i <= number_of_unknowns; i++) {
    TV_INT cell_index = matrix_index_to_cell_index(i);
    laplace->u(cell_index) = x(i);
  }
  // Set some velocity to zero.
  // for (typename T_GRID::CELL_ITERATOR iterator(laplace->grid, 1);
  //     iterator.Valid();
  //     iterator.Next()) {
  //  int filled_region_color =
  //      laplace->filled_region_colors(iterator.Cell_Index());
  //  if (filled_region_color > 0 &&
  //      !laplace->filled_region_touches_dirichlet(filled_region_color))
  //    laplace->u(iterator.Cell_Index()) = 0;
  // }
}

}  // namespace PhysBAM

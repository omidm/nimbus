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

#include "applications/physbam/water//app_utils.h"
#include "applications/physbam/water//projection/projection_helper.h"

namespace PhysBAM {

bool All_Cell_Faces_Neumann(const TV_INT& cell_index, const T_PSI_N& psi_N);

void FillUniformRegionColor(
    const T_GRID& grid,
    const T_PSI_D& psi_D,
    const T_PSI_N& psi_N,
    const bool solve_single_cell_neumann_regions,
    T_COLOR* filled_region_colors) {
  for (typename T_GRID::CELL_ITERATOR iterator(grid, 1, T_GRID::GHOST_REGION);
       iterator.Valid();
       iterator.Next()) {
    (*filled_region_colors)(iterator.Cell_Index())=-1;
  }
  for (typename T_GRID::CELL_ITERATOR iterator(grid);
       iterator.Valid();
       iterator.Next()) {
    if (psi_D(iterator.Cell_Index()) ||
        (!solve_single_cell_neumann_regions &&
         All_Cell_Faces_Neumann(iterator.Cell_Index(), psi_N))) {
      (*filled_region_colors)(iterator.Cell_Index()) = -1;
    } else {
      (*filled_region_colors)(iterator.Cell_Index()) = 1;
    }
  }
}

bool All_Cell_Faces_Neumann(
    const TV_INT& cell_index, const T_PSI_N& psi_N) {
  for (int axis = 1; axis <= T_GRID::dimension; axis++)
    if (!psi_N.Component(axis)(cell_index) ||
        !psi_N.Component(axis)(cell_index+TV_INT::Axis_Vector(axis))) {
      return false;
    }
  return true;
}

void Find_Matrix_Indices_In_Region(
    const T_GRID& local_grid,
    const int region_index,
    const RANGE<TV_INT>& region,
    const ARRAY<int, TV_INT>& filled_region_colors,
    ARRAY<int, VECTOR<int,1> >* filled_region_cell_count,
    ARRAY<int, TV_INT>* cell_index_to_matrix_index,
    ARRAY<TV_INT>* matrix_index_to_cell_index,
    int* local_n = NULL,
    int* interior_n = NULL) {
  typedef UNIFORM_GRID_ITERATOR_CELL<TV> CELL_ITERATOR;
  for (CELL_ITERATOR iterator(local_grid, region);
       iterator.Valid();
       iterator.Next()) {
    TV_INT c = iterator.Cell_Index();
    int color = filled_region_colors(c);
    if (color < 1) {
      continue;
    }
    assert(color == 1);
    int new_index = ++(*filled_region_cell_count)(color);
    (*cell_index_to_matrix_index)(c) = new_index;
    (*matrix_index_to_cell_index)(new_index) = c;
  }
  if (region_index == 0) {
    assert(interior_n != NULL);
    *interior_n = (*filled_region_cell_count)(1);
  } else if (region_index == 6) {
    assert(local_n != NULL);
    *local_n = (*filled_region_cell_count)(1);
  }
}

void FindMatrixIndices(
    const T_GRID& local_grid,
    const ARRAY<int, TV_INT>& filled_region_colors,
    ARRAY<int, VECTOR<int, 1> >* filled_region_cell_count,
    ARRAY<int, TV_INT>* cell_index_to_matrix_index,
    ARRAY<TV_INT >* matrix_index_to_cell_index,
    int* local_n,
    int* interior_n) {
  assert(local_grid.Is_MAC_Grid());
  int m = local_grid.counts.x;
  int n = local_grid.counts.y;
  int mn = local_grid.counts.z;
  dbg(APP_LOG, "Local grid in finding matrix indices: %d,%d,%d.\n",
      m, n, mn);
  Find_Matrix_Indices_In_Region(
      local_grid,
      0, RANGE<VECTOR<int,3> >(1,m,1,n,1,mn),
      filled_region_colors,
      filled_region_cell_count,
      cell_index_to_matrix_index,
      matrix_index_to_cell_index,
      NULL,
      interior_n);
  Find_Matrix_Indices_In_Region(
      local_grid,
      1, RANGE<VECTOR<int,3> >(0,0,0,n+1,0,mn+1),
      filled_region_colors,
      filled_region_cell_count,
      cell_index_to_matrix_index,
      matrix_index_to_cell_index);
  Find_Matrix_Indices_In_Region(
      local_grid,
      2, RANGE<VECTOR<int,3> >(m+1,m+1,0,n+1,0,mn+1),
      filled_region_colors,
      filled_region_cell_count,
      cell_index_to_matrix_index,
      matrix_index_to_cell_index);
  Find_Matrix_Indices_In_Region(
      local_grid,
      3, RANGE<VECTOR<int,3> >(1,m,0,0,0,mn+1),
      filled_region_colors,
      filled_region_cell_count,
      cell_index_to_matrix_index,
      matrix_index_to_cell_index);
  Find_Matrix_Indices_In_Region(
      local_grid,
      4, RANGE<VECTOR<int,3> >(1,m,n+1,n+1,0,mn+1),
      filled_region_colors,
      filled_region_cell_count,
      cell_index_to_matrix_index,
      matrix_index_to_cell_index);
  Find_Matrix_Indices_In_Region(
      local_grid,
      5, RANGE<VECTOR<int,3> >(1,m,1,n,0,0),
      filled_region_colors,
      filled_region_cell_count,
      cell_index_to_matrix_index,
      matrix_index_to_cell_index);
  Find_Matrix_Indices_In_Region(
      local_grid,
      6, RANGE<VECTOR<int,3> >(1,m,1,n,mn+1,mn+1),
      filled_region_colors,
      filled_region_cell_count,
      cell_index_to_matrix_index,
      matrix_index_to_cell_index,
      local_n,
      NULL);
}

}  // namespace PhysBAM

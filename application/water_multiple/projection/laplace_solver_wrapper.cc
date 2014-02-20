#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/projection/laplace_solver_wrapper.h"

namespace PhysBAM {

void LaplaceSolverWrapper::Solve() {
  const int number_of_regions = laplace->number_of_regions;
  // int -> (dim_t, dim_t, dim_t)
  ARRAY<ARRAY<TV_INT> > matrix_index_to_cell_index_array(number_of_regions);

  // (dim_t, dim_t, dim_t) -> int
  ARRAY<int, TV_INT> cell_index_to_matrix_index(
      laplace->grid.Domain_Indices(1));

  // color_id -> int
  ARRAY<int, VECTOR<int, 1> > filled_region_cell_count(-1, number_of_regions);

  ARRAY<SPARSE_MATRIX_FLAT_NXN<T> > A_array(number_of_regions);

  ARRAY<VECTOR_ND<T> > b_array(number_of_regions);

  // Count the cells in each region.
  for (typename T_GRID::CELL_ITERATOR iterator(laplace->grid, 1);
       iterator.Valid();
       iterator.Next()) {
    filled_region_cell_count(
        laplace->filled_region_colors(iterator.Cell_Index()))++;
  }

  for (int color = 1; color <= number_of_regions; color++)
    if (laplace->filled_region_touches_dirichlet(color) ||
        laplace->solve_neumann_regions) {
      matrix_index_to_cell_index_array(color).Resize(
          filled_region_cell_count(color));
    }

  // reusing this array in order to make the indirection arrays
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
  laplace->Find_A(
      domain, A_array, b_array,
      filled_region_cell_count, cell_index_to_matrix_index);
  for (int color = 1; color <= number_of_regions; color++)
    if (filled_region_cell_count(color) > 0 &&
        (laplace->filled_region_touches_dirichlet(color) ||
         laplace->solve_neumann_regions)) {
      laplace->pcg.Enforce_Compatibility(
          !laplace->filled_region_touches_dirichlet(color) &&
          laplace->enforce_compatibility);
      SolveSubregion(
          matrix_index_to_cell_index_array(color),
          A_array(color),
          b_array(color),
          color);
    }

  // Set some velocity to zero.
  if (!laplace->solve_neumann_regions)
    for (typename T_GRID::CELL_ITERATOR iterator(laplace->grid, 1);
         iterator.Valid();
         iterator.Next()) {
      int filled_region_color =
          laplace->filled_region_colors(iterator.Cell_Index());
      if (filled_region_color > 0 &&
          !laplace->filled_region_touches_dirichlet(filled_region_color))
        laplace->u(iterator.Cell_Index()) = 0;
    }
}

void LaplaceSolverWrapper::SolveSubregion(
    ARRAY<TV_INT>& matrix_index_to_cell_index,
    SPARSE_MATRIX_FLAT_NXN<T>& A,
    VECTOR_ND<T>& b,
    const int color) {
  // "m" means Size.
  int number_of_unknowns = matrix_index_to_cell_index.m;
  A.Negate();
  b *= (T) -1;
  VECTOR_ND<T> x(number_of_unknowns), q, s, r, k, z;
  for (int i = 1; i <= number_of_unknowns; i++)
    x(i) = laplace->u(matrix_index_to_cell_index(i));

  laplace->Find_Tolerance(b); // needs to happen after b is completely set up

  // MPI reference version:
  // laplace_mpi->Solve(A, x, b, q, s, r, k, z, tolerance, color);
  // color only used for MPI version.
  laplace->pcg.Solve(A, x, b, q, s, r, k, z, laplace->tolerance);

  for (int i = 1; i <= number_of_unknowns; i++) {
    TV_INT cell_index = matrix_index_to_cell_index(i);
    laplace->u(cell_index) = x(i);
  }
}

}  // namespace PhysBAM

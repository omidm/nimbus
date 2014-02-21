#include "application/water_multiple/app_utils.h"

#ifndef NIMBUS_APPLICATION_WATER_MULTIPLE_PROJECTION_LAPLACE_SOLVER_WRAPPER_H_
#define NIMBUS_APPLICATION_WATER_MULTIPLE_PROJECTION_LAPLACE_SOLVER_WRAPPER_H_

namespace PhysBAM {

typedef application::T T;
typedef application::TV TV;
typedef application::TV_INT TV_INT;

typedef GRID<TV> T_GRID;

class LaplaceSolverWrapper {
 public:
  LaplaceSolverWrapper(
      LAPLACE_COLLIDABLE_UNIFORM<GRID<TV> >* laplace_input) {
    laplace = laplace_input;
  }
  void Solve();
 private:
  LAPLACE_COLLIDABLE_UNIFORM<GRID<TV> >* laplace;
  void SolveSubregion(
      ARRAY<TV_INT>& matrix_index_to_cell_index,
      SPARSE_MATRIX_FLAT_NXN<T>& A,
      VECTOR_ND<T>& b,
      const int color);
};

}  // namespace PhysBAM

#endif  // NIMBUS_APPLICATION_WATER_MULTIPLE_PROJECTION_LAPLACE_SOLVER_WRAPPER_HELPER_H_

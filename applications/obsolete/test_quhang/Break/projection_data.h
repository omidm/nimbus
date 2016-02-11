#ifndef __PROJECTION_DATA__
#define __PROJECTION_DATA__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
namespace PhysBAM{
template<class TV>
class ProjectionData {
  typedef typename TV::SCALAR T;
  typedef typename TV::template REBIND<int>::TYPE TV_INT;
 public:
  ARRAY< ARRAY<TV_INT> >* matrix_index_to_cell_index_array;

  ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >* A_array;

  ARRAY<VECTOR_ND<T> >* b_array;

  ARRAY<int, TV_INT>* domain_index;

  VECTOR_ND<T>* x;

  T tolerance;
};


template<class TV>
class ProjectionInternalData {
  typedef typename TV::SCALAR T;
  typedef typename TV::template REBIND<int>::TYPE TV_INT;
 public:
  VECTOR_ND<T> *temp, *p, *z_interior, *x_interior, *b_interior, *p_interior, *temp_interior;
  int global_n;
  T global_tolerance;
  double rho, rho_old;
  T beta, alpha;
  int iteration;
  T residual;
  bool move_on;
};

}  // namespace PhysBAM

#endif

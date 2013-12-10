/*
 * The internal data structures of PhysBAM projection.
 * These data structures should not be here in future.
 * It should either be passed by Nimbus, or be config constantants.
 *
 * Author: Hang Qu <quhang@stanford.edu>
 */

#ifndef __PROJECTION_DATA__
#define __PROJECTION_DATA__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>

namespace PhysBAM{
template<class TV>
class ProjectionData {
 private:
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

// All data needed to run PCG.
template<class TV>
class ProjectionInternalData {
 private:
  typedef typename TV::SCALAR T;
  typedef typename TV::template REBIND<int>::TYPE TV_INT;
 public:
  VECTOR_ND<T> *temp, *p, *p_interior;
};
}  // namespace PhysBAM

#endif


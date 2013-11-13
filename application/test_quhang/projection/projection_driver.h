#ifndef __PROJECTION_DRIVER__
#define __PROJECTION_DRIVER__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/LAPLACE_UNIFORM.h>
#include "projection_data.h"
#include "nimbus_pcg_sparse_mpi.h"

namespace PhysBAM{

template<class TV> class PROJECTION_EXAMPLE;
template<class T_GRID> class NIMBUS_PCG_SPARSE_MPI;

template<class TV>
class PROJECTION_DRIVER {
 public:
  typedef typename TV::SCALAR T;
  typedef typename TV::template REBIND<int>::TYPE TV_INT;

  T time;
  int current_frame;
  int output_number;
  ProjectionData<TV> *projection_data;
  ProjectionInternalData<TV> *projection_internal_data;
  NIMBUS_PCG_SPARSE_MPI< GRID<TV> >* pcg_mpi;

  PROJECTION_EXAMPLE<TV>& example;
  PROJECTION_DRIVER(PROJECTION_EXAMPLE<TV>& example);
  virtual ~PROJECTION_DRIVER();

  void PrepareForProjection();
  void PrepareForOneRegion();
  void WindUpForOneRegion();
  void ApplyPressureAndFinish();

  void Execute_Main_Program();

  void Initialize();
  void Write_Output_Files(const int frame);
  void Write_Substep(const std::string& title,const int substep,const int level=0);
};

}
#endif

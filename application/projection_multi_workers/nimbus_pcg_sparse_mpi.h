/* 
 * This file should be changed a lot in the future. The interface is not good.
 *
 * Author: Hang Qu <quhang@stanford.edu>
 */
#ifndef __NIMBUS_PCG_SPARSE_MPI__
#define __NIMBUS_PCG_SPARSE_MPI__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include "projection_data.h"

namespace PhysBAM{

class SPARSE_MATRIX_PARTITION;

template<class T_GRID>
class NIMBUS_PCG_SPARSE_MPI : public NONCOPYABLE
{
 public:
  typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
  PCG_SPARSE<T>& pcg;
  MPI::Intracomm& comm;
  SPARSE_MATRIX_PARTITION& partition;
  ARRAY<MPI::Datatype> boundary_datatypes,ghost_datatypes;
  ARRAY<ARRAY<int> > columns_to_send;
  ARRAY<ARRAY<int> > columns_to_receive;

  NIMBUS_PCG_SPARSE_MPI(PCG_SPARSE<T>& pcg_input,MPI::Intracomm& comm_input,SPARSE_MATRIX_PARTITION& partition_input)
      : pcg(pcg_input),comm(comm_input),partition(partition_input) {}

  virtual ~NIMBUS_PCG_SPARSE_MPI() {
    MPI_UTILITIES::Free_Elements_And_Clean_Memory(boundary_datatypes);MPI_UTILITIES::Free_Elements_And_Clean_Memory(ghost_datatypes);
  }

  // [TODO] Global sum should not use MPI.
  template<class TYPE> TYPE Global_Sum(const TYPE& input) {
    TYPE output;
    MPI_UTILITIES::Reduce(input,output,MPI::SUM,comm);
    return output;
  }

  // [TODO] Global max should not use MPI.
  template<class TYPE> TYPE Global_Max(const TYPE& input) {
     TYPE output;MPI_UTILITIES::Reduce(input,output,MPI::MAX,comm);
    return output;
  }

  // [TODO] Ghost cell transmission should not use MPI.
  virtual void Fill_Ghost_Cells(VECTOR_ND<T>& v) {
    ARRAY<MPI::Request> requests;
    requests.Preallocate(2*partition.number_of_sides);
    for(int s=1;s<=partition.number_of_sides;s++)
      if(boundary_datatypes(s)!=MPI::DATATYPE_NULL)
        requests.Append(comm.Isend(v.x-1,1,boundary_datatypes(s),partition.neighbor_ranks(s),s));
    for(int s=1;s<=partition.number_of_sides;s++)
      if(ghost_datatypes(s)!=MPI::DATATYPE_NULL)
        requests.Append(comm.Irecv(v.x-1,1,ghost_datatypes(s),partition.neighbor_ranks(s),((s-1)^1)+1));
    MPI_UTILITIES::Wait_All(requests);
  }
  
  virtual void Initialize_Datatypes();

  // Functions called during projection.
  // Currently all the needed data is passed by "projection_internal_data" and "projection_data".
  // [TODO] Change data abstraction.

  // [TODO] Data is initialized here, to be moved.
  void Initialize(
      ProjectionInternalData<TV>* projection_internal_data,
      ProjectionData<TV>* projection_data);
  // [TODO] Data is initialized here, to be moved.
  // [TODO] Communcation involved to exchange config parameters.
  void CommunicateConfig( 
      ProjectionInternalData<TV>* projection_internal_data,
      ProjectionData<TV>* projection_data);
  void ExchangePressure(
      ProjectionInternalData<TV>* projection_internal_data,
      ProjectionData<TV>* projection_data);
  void InitializeResidual(
      ProjectionInternalData<TV>* projection_internal_data,
      ProjectionData<TV>* projection_data);
  void SpawnFirstIteration(
      ProjectionInternalData<TV>* projection_internal_data,
      ProjectionData<TV>* projection_data);
  void DoPrecondition(
      ProjectionInternalData<TV>* projection_internal_data,
      ProjectionData<TV>* projection_data);
  void CalculateBeta(
      ProjectionInternalData<TV>* projection_internal_data,
      ProjectionData<TV>* projection_data);
  void UpdateSearchVector(
      ProjectionInternalData<TV>* projection_internal_data,
      ProjectionData<TV>* projection_data);
  void ExchangeSearchVector(
      ProjectionInternalData<TV>* projection_internal_data,
      ProjectionData<TV>* projection_data);
  void UpdateTempVector(
      ProjectionInternalData<TV>* projection_internal_data,
      ProjectionData<TV>* projection_data);
  void CalculateAlpha(
      ProjectionInternalData<TV>* projection_internal_data,
      ProjectionData<TV>* projection_data);
  void UpdateOtherVectors(
      ProjectionInternalData<TV>* projection_internal_data,
      ProjectionData<TV>* projection_data);
  void CalculateResidual(
      ProjectionInternalData<TV>* projection_internal_data,
      ProjectionData<TV>* projection_data);
  void DecideToSpawnNextIteration(
      ProjectionInternalData<TV>* projection_internal_data,
      ProjectionData<TV>* projection_data);
};
}  // namespace PhysBAM
#endif

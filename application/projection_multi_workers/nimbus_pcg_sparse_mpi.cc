#ifdef USE_MPI

#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <PhysBAM_Tools/Vectors/SPARSE_VECTOR_ND.h>
#include "nimbus_pcg_sparse_mpi.h"

using namespace PhysBAM;

template<class T_GRID> void NIMBUS_PCG_SPARSE_MPI<T_GRID>::Initialize(
    ProjectionInternalData<TV>* projection_internal_data,
    ProjectionData<TV>* projection_data) {
  int color = 1;
  Initialize_Datatypes();
  int local_n=(*projection_data->A_array)(color).n;
  int interior_n=partition.interior_indices.Size()+1;

  projection_internal_data->temp = new VECTOR_ND<T>(local_n,false);
  projection_internal_data->p = new VECTOR_ND<T>(local_n,false);
  projection_internal_data->z_interior = new VECTOR_ND<T>(interior_n,false);

  projection_internal_data->x_interior = new VECTOR_ND<T>;
  projection_internal_data->x_interior->Set_Subvector_View(*projection_data->x,partition.interior_indices);
  projection_internal_data->b_interior = new VECTOR_ND<T>;
  projection_internal_data->b_interior->Set_Subvector_View((*projection_data->b_array)(color),partition.interior_indices);
  projection_internal_data->p_interior = new VECTOR_ND<T>;
  projection_internal_data->p_interior->Set_Subvector_View(*projection_internal_data->p,partition.interior_indices);
  projection_internal_data->temp_interior = new VECTOR_ND<T>;
  projection_internal_data->temp_interior->Set_Subvector_View(*projection_internal_data->temp,partition.interior_indices);
  projection_internal_data->rho=0;
  projection_internal_data->rho_old=0;
  projection_internal_data->alpha=0;
  projection_internal_data->beta=0;
}

template<class T_GRID> void NIMBUS_PCG_SPARSE_MPI<T_GRID>::
Initialize_Datatypes() {
    MPI_UTILITIES::Free_Elements_And_Clean_Memory(boundary_datatypes);MPI_UTILITIES::Free_Elements_And_Clean_Memory(ghost_datatypes);
    boundary_datatypes.Resize(partition.number_of_sides);ghost_datatypes.Resize(partition.number_of_sides);
    for(int s=1;s<=partition.number_of_sides;s++) if(partition.neighbor_ranks(s)!=MPI::PROC_NULL){
        if(partition.boundary_indices(s).m){
            const ARRAY<int>& displacements=partition.boundary_indices(s);
            ARRAY<int> block_lengths(displacements.m,false);ARRAYS_COMPUTATIONS::Fill(block_lengths,1);
            boundary_datatypes(s)=MPI_UTILITIES::Datatype<T>().Create_indexed(displacements.m,&block_lengths(1),&displacements(1));
            boundary_datatypes(s).Commit();}
        int ghost_indices_length=partition.ghost_indices(s).Size()+1;
        if(ghost_indices_length){
            ghost_datatypes(s)=MPI_UTILITIES::Datatype<T>().Create_indexed(1,&ghost_indices_length,&partition.ghost_indices(s).min_corner);
            ghost_datatypes(s).Commit();}}
}

template<class T_GRID> void NIMBUS_PCG_SPARSE_MPI<T_GRID>::
CommunicateConfig(
    ProjectionInternalData<TV>* projection_internal_data,
    ProjectionData<TV>* projection_data) {
  // [TODO] Communication!
  int interior_n=partition.interior_indices.Size()+1;
  int global_n=Global_Sum(interior_n);
  // [TODO] Communication!
  T global_tolerance=Global_Max(projection_data->tolerance);
  projection_internal_data->global_n = global_n;
  projection_internal_data->global_tolerance = global_tolerance;
  int desired_iterations=global_n;
  if(pcg.maximum_iterations) desired_iterations=min(desired_iterations,pcg.maximum_iterations);
}

// Projection is broken to "smallest" code piece to allow future changes.
// Data pointers are first initialized,
// and then computation functions are invoked.

template<class T_GRID> void NIMBUS_PCG_SPARSE_MPI<T_GRID>::
ExchangePressure(
    ProjectionInternalData<TV>* projection_internal_data,
    ProjectionData<TV>* projection_data) {
  VECTOR_ND<T>& x = (*projection_data->x); 
  Fill_Ghost_Cells(x);
}

template<class T_GRID> void NIMBUS_PCG_SPARSE_MPI<T_GRID>::
InitializeResidual(
    ProjectionInternalData<TV>* projection_internal_data,
    ProjectionData<TV>* projection_data) {
  int color = 1;
  SPARSE_MATRIX_FLAT_NXN<T>& A = (*projection_data->A_array)(color);
  VECTOR_ND<T>& x = (*projection_data->x); 
  VECTOR_ND<T>& temp = (*projection_internal_data->temp);
  VECTOR_ND<T>& b_interior = (*projection_internal_data->b_interior);
  VECTOR_ND<T>& temp_interior = (*projection_internal_data->temp_interior);
  A.Times(x,temp);b_interior-=temp_interior;
}

template<class T_GRID> void NIMBUS_PCG_SPARSE_MPI<T_GRID>::
SpawnFirstIteration(
    ProjectionInternalData<TV>* projection_internal_data,
    ProjectionData<TV>* projection_data) {
  int color = 1;
  SPARSE_MATRIX_FLAT_NXN<T>& A = (*projection_data->A_array)(color);
  const bool recompute_preconditioner = true;
  VECTOR_ND<T>& b_interior = (*projection_internal_data->b_interior);
  // T global_tolerance = projection_internal_data->global_tolerance;
  projection_internal_data->partial_norm = b_interior.Max_Abs();
  // TODO fine-grained control flow for this.
  //if(Global_Max(b_interior.Max_Abs())<=global_tolerance) {
  //  projection_internal_data->move_on = false;
  //  return;
  //}
  // find an incomplete cholesky preconditioner - actually an LU that saves square roots, and an inverted diagonal to save on divides
  if(pcg.incomplete_cholesky && (recompute_preconditioner || !A.C)){
      delete A.C;A.C=A.Create_Submatrix(partition.interior_indices);
      A.C->In_Place_Incomplete_Cholesky_Factorization(pcg.modified_incomplete_cholesky,pcg.modified_incomplete_cholesky_coefficient,
          pcg.preconditioner_zero_tolerance,pcg.preconditioner_zero_replacement);}

  projection_internal_data->rho=0,projection_internal_data->rho_old=0;
  projection_internal_data->move_on = true;
}

template<class T_GRID> void NIMBUS_PCG_SPARSE_MPI<T_GRID>::
DoPrecondition(
    ProjectionInternalData<TV>* projection_internal_data,
    ProjectionData<TV>* projection_data) {
  int color = 1;
  SPARSE_MATRIX_FLAT_NXN<T>& A = (*projection_data->A_array)(color);
  VECTOR_ND<T>& z_interior = (*projection_internal_data->z_interior);
  VECTOR_ND<T>& b_interior = (*projection_internal_data->b_interior);
  VECTOR_ND<T>& temp_interior = (*projection_internal_data->temp_interior);
  A.C->Solve_Forward_Substitution(b_interior,temp_interior,true); // diagonal should be treated as the identity
  A.C->Solve_Backward_Substitution(temp_interior,z_interior,false,true); // diagonal is inverted to save on divides
}

template<class T_GRID> void NIMBUS_PCG_SPARSE_MPI<T_GRID>::
CalculateBeta(
    ProjectionInternalData<TV>* projection_internal_data,
    ProjectionData<TV>* projection_data) {
  VECTOR_ND<T>& z_interior = (*projection_internal_data->z_interior);
  VECTOR_ND<T>& b_interior = (*projection_internal_data->b_interior);
  projection_internal_data->rho_old=projection_internal_data->rho;
  projection_internal_data->rho=Global_Sum(VECTOR_ND<T>::Dot_Product_Double_Precision(z_interior,b_interior));
  projection_internal_data->beta=(T)(projection_internal_data->rho/projection_internal_data->rho_old);
}

template<class T_GRID> void NIMBUS_PCG_SPARSE_MPI<T_GRID>::
UpdateSearchVector(
    ProjectionInternalData<TV>* projection_internal_data,
    ProjectionData<TV>* projection_data) {
  int interior_n=partition.interior_indices.Size()+1;
  VECTOR_ND<T>& z_interior = (*projection_internal_data->z_interior);
  VECTOR_ND<T>& p_interior = (*projection_internal_data->p_interior);
  if(projection_internal_data->iteration==1)
    p_interior=z_interior;
  else {
    for(int i=1;i<=interior_n;i++) p_interior(i)=z_interior(i)+projection_internal_data->beta*p_interior(i);
  }
}

template<class T_GRID> void NIMBUS_PCG_SPARSE_MPI<T_GRID>::
ExchangeSearchVector(
    ProjectionInternalData<TV>* projection_internal_data,
    ProjectionData<TV>* projection_data) {
  VECTOR_ND<T>& p = (*projection_internal_data->p);
  Fill_Ghost_Cells(p);
}

template<class T_GRID> void NIMBUS_PCG_SPARSE_MPI<T_GRID>::
UpdateTempVector(
    ProjectionInternalData<TV>* projection_internal_data,
    ProjectionData<TV>* projection_data) {
  int color = 1;
  SPARSE_MATRIX_FLAT_NXN<T>& A = (*projection_data->A_array)(color);
  VECTOR_ND<T>& temp = (*projection_internal_data->temp);
  VECTOR_ND<T>& p = (*projection_internal_data->p);
  A.Times(p,temp);
}

template<class T_GRID> void NIMBUS_PCG_SPARSE_MPI<T_GRID>::
CalculateAlpha(
    ProjectionInternalData<TV>* projection_internal_data,
    ProjectionData<TV>* projection_data) {
  VECTOR_ND<T>& p_interior = (*projection_internal_data->p_interior);
  VECTOR_ND<T>& temp_interior = (*projection_internal_data->temp_interior);
  projection_internal_data->alpha=(T)(projection_internal_data->rho/Global_Sum(VECTOR_ND<T>::Dot_Product_Double_Precision(p_interior,temp_interior)));
}

template<class T_GRID> void NIMBUS_PCG_SPARSE_MPI<T_GRID>::
UpdateOtherVectors(
    ProjectionInternalData<TV>* projection_internal_data,
    ProjectionData<TV>* projection_data) {
  int interior_n=partition.interior_indices.Size()+1;
  VECTOR_ND<T>& x_interior = (*projection_internal_data->x_interior);
  VECTOR_ND<T>& p_interior = (*projection_internal_data->p_interior);
  VECTOR_ND<T>& b_interior = (*projection_internal_data->b_interior);
  VECTOR_ND<T>& temp_interior = (*projection_internal_data->temp_interior);
  for(int i=1;i<=interior_n;i++) {
    x_interior(i)+=projection_internal_data->alpha*p_interior(i);
    b_interior(i)-=projection_internal_data->alpha*temp_interior(i);
  }
}

template<class T_GRID> void NIMBUS_PCG_SPARSE_MPI<T_GRID>::
CalculateResidual(
    ProjectionInternalData<TV>* projection_internal_data,
    ProjectionData<TV>* projection_data) {
  VECTOR_ND<T>& b_interior = (*projection_internal_data->b_interior);
  // projection_internal_data->residual=Global_Max(b_interior.Max_Abs());
  projection_internal_data->partial_norm=b_interior.Max_Abs();
}

template<class T_GRID> void NIMBUS_PCG_SPARSE_MPI<T_GRID>::
DecideToSpawnNextIteration(
    ProjectionInternalData<TV>* projection_internal_data,
    ProjectionData<TV>* projection_data) {
  int desired_iterations=projection_internal_data->global_n;
  if(projection_internal_data->residual<=projection_internal_data->global_tolerance) {
    projection_internal_data->move_on = false;
    return;
  }
  if(projection_internal_data->iteration==desired_iterations) {
    projection_internal_data->move_on = false;
    return;
  }
  projection_internal_data->move_on = true;
}

template class NIMBUS_PCG_SPARSE_MPI<GRID<VECTOR<float,2> > >;
template class NIMBUS_PCG_SPARSE_MPI<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class NIMBUS_PCG_SPARSE_MPI<GRID<VECTOR<double,2> > >;
template class NIMBUS_PCG_SPARSE_MPI<GRID<VECTOR<double,3> > >;
#endif

#endif

/*
   // The original code for projection.
    // find initial residual, r=b-Ax - reusing b for the residual
    // Communication!
    Fill_Ghost_Cells(x);
    A.Times(x,temp);b_interior-=temp_interior;
    // Communication!
    if(Global_Max(b_interior.Max_Abs())<=global_tolerance) {
      return;
    }

    // find an incomplete cholesky preconditioner - actually an LU that saves square roots, and an inverted diagonal to save on divides
    if(pcg.incomplete_cholesky && (recompute_preconditioner || !A.C)){
        delete A.C;A.C=A.Create_Submatrix(partition.interior_indices);
        A.C->In_Place_Incomplete_Cholesky_Factorization(pcg.modified_incomplete_cholesky,pcg.modified_incomplete_cholesky_coefficient,
            pcg.preconditioner_zero_tolerance,pcg.preconditioner_zero_replacement);}

    double rho=0,rho_old=0;
    for(int iteration=1;;iteration++){
        if(pcg.incomplete_cholesky){
            // solve Mz=r
            A.C->Solve_Forward_Substitution(b_interior,temp_interior,true); // diagonal should be treated as the identity
            A.C->Solve_Backward_Substitution(temp_interior,z_interior,false,true);} // diagonal is inverted to save on divides
        else z_interior=b_interior; // set z=r when there is no preconditioner

        // update search direction
	// Communication.
        rho_old=rho;rho=Global_Sum(VECTOR_ND<T>::Dot_Product_Double_Precision(z_interior,b_interior));
        T beta=0;
	if(iteration==1)
	  p_interior=z_interior;
	else{beta=(T)(rho/rho_old);for(int i=1;i<=interior_n;i++) p_interior(i)=z_interior(i)+beta*p_interior(i);} // when iteration=1, beta=0

        // update solution and residual
        Fill_Ghost_Cells(p);
        A.Times(p,temp);
        T alpha=(T)(rho/Global_Sum(VECTOR_ND<T>::Dot_Product_Double_Precision(p_interior,temp_interior)));
        for(int i=1;i<=interior_n;i++){x_interior(i)+=alpha*p_interior(i);b_interior(i)-=alpha*temp_interior(i);}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        T residual=Global_Max(b_interior.Max_Abs());

        // check for convergence
        std::stringstream ss;
        if(pcg.show_residual) ss<<residual<<std::endl;
        if(residual<=global_tolerance){if(pcg.show_results) ss<<"NUMBER OF ITERATIONS = "<<iteration<<std::endl;break;}
        if(iteration==desired_iterations){if(pcg.show_results) ss<<"DID NOT CONVERGE IN "<<iteration<<" ITERATIONS"<<std::endl;break;}
        LOG::filecout(ss.str());
#endif
    }
  
    Fill_Ghost_Cells(x);
*/

#ifdef USE_MPI

#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <PhysBAM_Tools/Vectors/SPARSE_VECTOR_ND.h>
#include "nimbus_pcg_sparse_mpi.h"

using namespace PhysBAM;

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

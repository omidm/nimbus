//#####################################################################
// Copyright 2008, Avi Robinson-Mosher, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONJUGATE_RESIDUAL_SPARSE_MPI
//#####################################################################
#ifndef __CONJUGATE_RESIDUAL_SPARSE_MPI__
#define __CONJUGATE_RESIDUAL_SPARSE_MPI__

#ifdef USE_MPI

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
namespace PhysBAM{

class SPARSE_MATRIX_PARTITION;
template<class TV> class GENERALIZED_VELOCITY;
template<class TV> class FLUID_SYSTEM_MPI;
template<class TV> class SOLID_SYSTEM_MPI;

template<class TV>
class CONJUGATE_RESIDUAL_SPARSE_MPI:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    MPI::Intracomm& global_comm;
    ARRAY<MPI::Intracomm>* fluid_comm;
    ARRAY<SPARSE_MATRIX_PARTITION>* partitions;
    ARRAY<ARRAY<MPI::Datatype> > boundary_datatypes_array,ghost_datatypes_array;
    bool print_diagnostics,print_residuals;
    T nullspace_tolerance; // don't attempt to invert eigenvalues approximately less than nullspace_tolerance*max_eigenvalue
    int* iterations_used;
    T residual_magnitude_squared,nullspace_measure; // extra convergence information
    int restart_iterations;
    int maximum_iterations;


    bool show_residual,show_results;
    bool incomplete_cholesky; // true when using this preconditioner or the modified version below
    bool modified_incomplete_cholesky; // true when using this preconditioner
    T modified_incomplete_cholesky_coefficient;
    T preconditioner_zero_tolerance,preconditioner_zero_replacement;

    CONJUGATE_RESIDUAL_SPARSE_MPI(MPI::Intracomm& global_comm_input,ARRAY<MPI::Intracomm>* fluid_comm_input,ARRAY<SPARSE_MATRIX_PARTITION>* partitions_input)
        :global_comm(global_comm_input),fluid_comm(fluid_comm_input),partitions(partitions_input),print_diagnostics(true),print_residuals(true),nullspace_tolerance((T)1e-5),iterations_used(0),
        restart_iterations(0),maximum_iterations(200),show_residual(false),show_results(false),preconditioner_zero_tolerance((T)1e-8),preconditioner_zero_replacement((T)1e-8)
    {
        Use_Modified_Incomplete_Cholesky();
    }

    ~CONJUGATE_RESIDUAL_SPARSE_MPI()
    {for(int i=1;i<=boundary_datatypes_array.m;i++) MPI_UTILITIES::Free_Elements_And_Clean_Memory(boundary_datatypes_array(i));
    for(int i=1;i<=ghost_datatypes_array.m;i++) MPI_UTILITIES::Free_Elements_And_Clean_Memory(ghost_datatypes_array(i));}

    // use .97 for octrees and .99 for uniform grids!
    void Use_Modified_Incomplete_Cholesky(const T modified_incomplete_cholesky_coefficient_input=(T).97)
    {incomplete_cholesky=true;modified_incomplete_cholesky=true;modified_incomplete_cholesky_coefficient=modified_incomplete_cholesky_coefficient_input;} // note that both are true

    template<class TV2> TV2 Global_Sum(const TV2& input,MPI::Intracomm& comm)
    {TV2 output;MPI_UTILITIES::Reduce(input,output,MPI::SUM,comm);return output;}

    template<class TV2> TV2 Global_Max(const TV2& input,MPI::Intracomm& comm)
    {TV2 output;MPI_UTILITIES::Reduce(input,output,MPI::MAX,comm);return output;}

    void Fill_Ghost_Cells(GENERALIZED_VELOCITY<TV>& V) {} // stub for solids

//#####################################################################
    void Parallel_Solve_Fluid_Part(FLUID_SYSTEM_MPI<TV>& fluid_system,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& x_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& b_array,
        KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& p_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& ap_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& ar_array,
        KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& r_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& z_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& zaq_array,
        const int min_iterations,const int max_iterations,const T tolerance,const bool recompute_preconditioner);
    void Parallel_Solve_Solid_Part(SOLID_SYSTEM_MPI<TV>& solid_system,GENERALIZED_VELOCITY<TV>& x_array,GENERALIZED_VELOCITY<TV>& b_array,GENERALIZED_VELOCITY<TV>& p_array,
        GENERALIZED_VELOCITY<TV>& ap_array,GENERALIZED_VELOCITY<TV>& ar_array,GENERALIZED_VELOCITY<TV>& r_array,GENERALIZED_VELOCITY<TV>& z_array,GENERALIZED_VELOCITY<TV>& zaq_array,
        const int min_iterations,const int max_iterations,const T tolerance);
    template<class T_SYSTEM,class TV2> bool Parallel_Solve(T_SYSTEM& system,TV2& x_array,const TV2& b_array,TV2& p_array,TV2& ap_array,TV2& ar_array,TV2& r_array,TV2& z_array,
        TV2& zaq_array,const int min_iterations,const int maximum_iterations,const T tolerance=1e-7);
private:
    void Fill_Ghost_Cells(KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& v_array);
    void Initialize_Datatypes();
//#####################################################################
};
}
#endif
#endif

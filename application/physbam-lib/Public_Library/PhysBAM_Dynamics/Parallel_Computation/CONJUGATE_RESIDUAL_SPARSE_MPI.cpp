//#####################################################################
// Copyright 2008, Avi Robinson-Mosher, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifdef USE_MPI
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <PhysBAM_Tools/Vectors/SPARSE_VECTOR_ND.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
#include <PhysBAM_Dynamics/Parallel_Computation/CONJUGATE_RESIDUAL_SPARSE_MPI.h>
#include <PhysBAM_Dynamics/Parallel_Computation/FLUID_SYSTEM_MPI.h>
#include <PhysBAM_Dynamics/Parallel_Computation/SOLID_SYSTEM_MPI.h>
#include <limits>
using namespace PhysBAM;

//#####################################################################
// Function Parallel_Solve_Fluid_Part
//#####################################################################
template<class TV> void CONJUGATE_RESIDUAL_SPARSE_MPI<TV>::
Parallel_Solve_Fluid_Part(FLUID_SYSTEM_MPI<TV>& fluid_system,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& x_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& b_array,
    KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& p_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& ap_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& ar_array,
    KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& r_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& z_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& zaq_array,
    const int min_iterations,const int max_iterations,const T tolerance,const bool recompute_preconditioner)
{
    Initialize_Datatypes();

    ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array=fluid_system.A_array;
    // find an incomplete cholesky preconditioner - actually an LU that saves square roots, and an inverted diagonal to save on divides
    for(int color=1;color<=A_array.m;color++){
        SPARSE_MATRIX_FLAT_NXN<T>& A=A_array(color);
        if(incomplete_cholesky && (recompute_preconditioner || !A.C)){
            if(color<=partitions->m){
                delete A.C;A.C=A.Create_Submatrix((*partitions)(color).interior_indices);
                A.C->In_Place_Incomplete_Cholesky_Factorization(modified_incomplete_cholesky,modified_incomplete_cholesky_coefficient,
                    preconditioner_zero_tolerance,preconditioner_zero_replacement);}
            else A.Construct_Incomplete_Cholesky_Factorization(modified_incomplete_cholesky,modified_incomplete_cholesky_coefficient,
                preconditioner_zero_tolerance,preconditioner_zero_replacement);}}

    Parallel_Solve(fluid_system,x_array,b_array,p_array,ap_array,ar_array,r_array,z_array,zaq_array,min_iterations,max_iterations,tolerance);
}
//#####################################################################
// Function Parallel_Solve
//#####################################################################
template<class TV> void CONJUGATE_RESIDUAL_SPARSE_MPI<TV>::
Parallel_Solve_Solid_Part(SOLID_SYSTEM_MPI<TV>& solid_system,GENERALIZED_VELOCITY<TV>& x_array,GENERALIZED_VELOCITY<TV>& b_array,GENERALIZED_VELOCITY<TV>& p_array,
    GENERALIZED_VELOCITY<TV>& ap_array,GENERALIZED_VELOCITY<TV>& ar_array,GENERALIZED_VELOCITY<TV>& r_array,GENERALIZED_VELOCITY<TV>& z_array,GENERALIZED_VELOCITY<TV>& zaq_array,
    const int min_iterations,const int max_iterations,const T tolerance)
{
    Parallel_Solve(solid_system,x_array,b_array,p_array,ap_array,ar_array,r_array,z_array,zaq_array,min_iterations,max_iterations,tolerance);
}
//#####################################################################
// Function Parallel_Solve
//#####################################################################
template<class TV> template<class T_SYSTEM,class TV2> bool CONJUGATE_RESIDUAL_SPARSE_MPI<TV>::
Parallel_Solve(T_SYSTEM& system,TV2& x_array,const TV2& b_array,TV2& p_array,TV2& ap_array,TV2& ar_array,TV2& r_array,TV2& z_array,
    TV2& zaq_array,const int min_iterations,const int max_iterations,const T tolerance)
{
    static const T small_number=std::numeric_limits<T>::epsilon();
    system.Set_Boundary_Conditions(x_array);

    double rho_old=0;T convergence_norm=0;
    int iteration=0;
    for(;;iteration++){
        bool restart=!iteration || (restart_iterations && iteration%restart_iterations==0);
        if(restart){
            r_array=b_array;
            Fill_Ghost_Cells(x_array);
            system.Multiply(x_array,p_array);
            r_array-=p_array;
            system.Project(r_array);
            if(show_residual) LOG::filecout("restarting conjugate_residual\n");
        }
        // stopping conditions
        convergence_norm=Global_Max(system.Convergence_Norm(r_array),global_comm);
        if(print_residuals) {std::stringstream ss;ss<<convergence_norm<<std::endl;LOG::filecout(ss.str());}
        residual_magnitude_squared=(T)Global_Sum(system.Inner_Product(r_array,r_array),global_comm); // reduce
        nullspace_measure=residual_magnitude_squared?(T)abs(rho_old/residual_magnitude_squared):0;
        if((convergence_norm<=tolerance || (iteration && nullspace_measure<=nullspace_tolerance)) &&
            (iteration>=min_iterations || convergence_norm<small_number)){ // TODO: get the stopping criterion right
            if(print_diagnostics) LOG::Stat("conjugate_residual iterations",iteration);if(iterations_used) *iterations_used=iteration;
            Fill_Ghost_Cells(x_array);return true;}
        if(iteration==max_iterations) break;

        // TODO: fix me for non-preconditioned
        system.Precondition(r_array,z_array);
        system.Project(z_array);system.Project_Nullspace(z_array);
        //TV2& mr_array=system.Precondition(r_array,z_array);
        TV2& mr_array=z_array;

        Fill_Ghost_Cells(mr_array);
        system.Multiply(mr_array,ar_array);
        system.Project(ar_array);

        T rho=(T)Global_Sum(system.Inner_Product(mr_array,ar_array),global_comm);

        if(restart){p_array=mr_array;ap_array=ar_array;}
        else{T beta=rho/(T)rho_old;T_SYSTEM::Copy(beta,p_array,mr_array,p_array);T_SYSTEM::Copy(beta,ap_array,ar_array,ap_array);}
        const TV2& map_array=debug_cast<const TV2&>(system.Precondition(ap_array,zaq_array));
        T alpha=rho/(T)Global_Sum(system.Inner_Product(map_array,ap_array),global_comm);
        T_SYSTEM::Copy(alpha,p_array,x_array,x_array);
        T_SYSTEM::Copy(-alpha,ap_array,r_array,r_array);
        rho_old=rho;}

    if(print_diagnostics && iteration==max_iterations){
        LOG::Stat("conjugate_residual iterations",iteration);
        std::stringstream ss;ss<<"convergence norm after maximum number of conjugate_residual iterations = "<<convergence_norm<<std::endl;LOG::filecout(ss.str());}

    Fill_Ghost_Cells(x_array);
    return false;
}
//#####################################################################
// Function Initialize_Datatypes
//#####################################################################
template<class TV> void CONJUGATE_RESIDUAL_SPARSE_MPI<TV>::
Fill_Ghost_Cells(KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& v_array)
{
    for(int p=1;p<=partitions->m;p++){
        SPARSE_MATRIX_PARTITION& partition=(*partitions)(p);
        VECTOR_ND<T>& v=v_array.v(p);
        ARRAY<MPI::Datatype>& boundary_datatypes=boundary_datatypes_array(p);
        ARRAY<MPI::Datatype>& ghost_datatypes=ghost_datatypes_array(p);
        ARRAY<MPI::Request> requests;requests.Preallocate(2*partition.number_of_sides);
        for(int s=1;s<=partition.number_of_sides;s++)if(boundary_datatypes(s)!=MPI::DATATYPE_NULL) requests.Append((*fluid_comm)(p).Isend(v.x-1,1,boundary_datatypes(s),partition.neighbor_ranks(s),s));
        for(int s=1;s<=partition.number_of_sides;s++)if(ghost_datatypes(s)!=MPI::DATATYPE_NULL) requests.Append((*fluid_comm)(p).Irecv(v.x-1,1,ghost_datatypes(s),partition.neighbor_ranks(s),((s-1)^1)+1));
        MPI_UTILITIES::Wait_All(requests);}
}
//#####################################################################
// Function Initialize_Datatypes
//#####################################################################
template<class TV> void CONJUGATE_RESIDUAL_SPARSE_MPI<TV>::
Initialize_Datatypes()
{
    boundary_datatypes_array.Resize(partitions->m);
    ghost_datatypes_array.Resize(partitions->m);
    for(int p=1;p<=partitions->m;p++){
        SPARSE_MATRIX_PARTITION& partition=(*partitions)(p);
        ARRAY<MPI::Datatype>& boundary_datatypes=boundary_datatypes_array(p);
        ARRAY<MPI::Datatype>& ghost_datatypes=ghost_datatypes_array(p);
        MPI_UTILITIES::Free_Elements_And_Clean_Memory(boundary_datatypes);MPI_UTILITIES::Free_Elements_And_Clean_Memory(ghost_datatypes);
        boundary_datatypes.Resize(partition.number_of_sides);ghost_datatypes.Resize(partition.number_of_sides);
        for(int s=1;s<=partition.number_of_sides;s++) if(partition.neighbor_ranks(s)!=MPI::PROC_NULL){
            if(partition.boundary_indices(s).m){
                ARRAY<int>& displacements=partition.boundary_indices(s);
                ARRAY<int> block_lengths(displacements.m,false);ARRAYS_COMPUTATIONS::Fill(block_lengths,1);
                boundary_datatypes(s)=MPI_UTILITIES::Datatype<T>().Create_indexed(displacements.m,&block_lengths(1),&displacements(1)); // TODO: collapse consecutive elements into blocks
                boundary_datatypes(s).Commit();}
            int ghost_indices_length=partition.ghost_indices(s).Size()+1;
            if(ghost_indices_length){
                ghost_datatypes(s)=MPI_UTILITIES::Datatype<T>().Create_indexed(1,&ghost_indices_length,&partition.ghost_indices(s).min_corner);
                ghost_datatypes(s).Commit();}}}
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template class CONJUGATE_RESIDUAL_SPARSE_MPI<VECTOR<T,d> >; \
    template bool CONJUGATE_RESIDUAL_SPARSE_MPI<VECTOR<T,d> >::Parallel_Solve(SOLID_SYSTEM_MPI<VECTOR<T,d> >& system,GENERALIZED_VELOCITY<VECTOR<T,d> >& x_array, \
        const GENERALIZED_VELOCITY<VECTOR<T,d> >& b_array,GENERALIZED_VELOCITY<VECTOR<T,d> >& p_array,GENERALIZED_VELOCITY<VECTOR<T,d> >& ap_array, \
        GENERALIZED_VELOCITY<VECTOR<T,d> >& ar_array,GENERALIZED_VELOCITY<VECTOR<T,d> >& r_array,GENERALIZED_VELOCITY<VECTOR<T,d> >& z_array, \
        GENERALIZED_VELOCITY<VECTOR<T,d> >& zaq_array,const int min_iterations,const int max_iterations,const T tolerance); \
    template bool CONJUGATE_RESIDUAL_SPARSE_MPI<VECTOR<T,d> >::Parallel_Solve(FLUID_SYSTEM_MPI<VECTOR<T,d> >& system,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& x_array, \
        const KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& b_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& p_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& ap_array, \
        KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& ar_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& r_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& z_array, \
        KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& zaq_array,const int min_iterations,const int max_iterations,const T tolerance);

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif
#endif

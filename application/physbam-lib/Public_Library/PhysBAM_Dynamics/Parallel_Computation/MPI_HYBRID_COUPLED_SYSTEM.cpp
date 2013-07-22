//#####################################################################
// Copyright 2012, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_HYBRID_COUPLED_SYSTEM
//#####################################################################
#ifdef USE_MPI
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#endif
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_HYBRID_COUPLED_SYSTEM.h>
using namespace PhysBAM;

#ifdef USE_MPI
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPI_HYBRID_COUPLED_SYSTEM<TV>::
MPI_HYBRID_COUPLED_SYSTEM(const SPARSE_MATRIX_FLAT_NXN<T>& A_input,const ARRAY<TV>& positions,const ARRAY<int>& global_ids,MPI_GRID<GRID<TV> >& mpi_grid_input,bool use_diagonal_preconditioner)
    :KRYLOV_SYSTEM_BASE<T>(use_diagonal_preconditioner,false),A(A_input),mpi_grid(mpi_grid_input)
{
    // Note: global id must be consistent across processes
    LOG::SCOPE scope("Initialize MPI system");
    neighbor_ranks=mpi_grid.all_neighbor_ranks;
    neighbor_directions=mpi_grid.all_neighbor_directions;
    send_indices.Resize(neighbor_ranks.m);
    recv_indices.Resize(neighbor_ranks.m);
    HASHTABLE<TV_INT,int> direction_to_neighbor;
    for(int i=1;i<=neighbor_ranks.m;i++){
        direction_to_neighbor.Set(neighbor_directions(i),i);}

    // set up recv
    for(int i=1;i<=positions.m;i++){
        TV_INT direction=Get_Direction(positions(i));
        if(direction==TV_INT()) interior.Append(i);}
    ARRAY<bool> touched_by_local_matrix(positions.m);ARRAYS_COMPUTATIONS::Fill(touched_by_local_matrix,false);
    for(int i=1;i<=interior.m;i++){
        int row=interior(i);
        int start=A.offsets(row),end=A.offsets(row+1);
        for(int j=start;j<end;j++) touched_by_local_matrix(A.A(j).j)=true;}
    for(int i=1;i<=positions.m;i++) if(touched_by_local_matrix(i)) {
        TV_INT direction=Get_Direction(positions(i));
        if(direction!=TV_INT()){
            int* neighbor_index=direction_to_neighbor.Get_Pointer(direction);
            PHYSBAM_ASSERT(neighbor_index);
            recv_indices(*neighbor_index).Append(i);}}
    touched_by_local_matrix.Clean_Memory();
    
    // set up send
    ARRAY<int> send_counts(neighbor_ranks.m),recv_counts(neighbor_ranks.m);
    ARRAY<MPI::Request> requests;
    for(int i=1;i<=recv_indices.m;i++){
        recv_counts(i)=recv_indices(i).m;
        requests.Append(mpi_grid.comm->Isend(&recv_counts(i),1,MPI::INT,neighbor_ranks(i),mpi_grid.Get_Send_Tag(neighbor_directions(i))));
        requests.Append(mpi_grid.comm->Irecv(&send_counts(i),1,MPI::INT,neighbor_ranks(i),mpi_grid.Get_Recv_Tag(neighbor_directions(i))));}
    MPI_UTILITIES::Wait_All(requests);
    ARRAY<MPI_PACKAGE> packages;
    requests.Resize(0);
    ARRAY<ARRAY<int> > recv_global_ids(recv_indices.m);
    for(int i=1;i<=send_indices.m;i++){
        send_indices(i).Resize(send_counts(i));
        recv_global_ids(i).Resize(recv_counts(i));
        for(int j=1;j<=recv_counts(i);j++) recv_global_ids(i)(j)=global_ids(recv_indices(i)(j));
        MPI_PACKAGE send_package(recv_global_ids(i));
        packages.Append(send_package);requests.Append(send_package.Isend(*mpi_grid.comm,neighbor_ranks(i),mpi_grid.Get_Send_Tag(neighbor_directions(i))));        
        MPI_PACKAGE recv_package(send_indices(i));
        packages.Append(recv_package);requests.Append(recv_package.Irecv(*mpi_grid.comm,neighbor_ranks(i),mpi_grid.Get_Recv_Tag(neighbor_directions(i))));}
    MPI_UTILITIES::Wait_All(requests);MPI_PACKAGE::Free_All(packages);
    HASHTABLE<int,int> global_to_local;
    for(int i=1;i<=global_ids.m;i++) global_to_local.Set(global_ids(i),i);
    for(int i=1;i<=send_indices.m;i++) for(int j=1;j<=send_indices(i).m;j++){
        int* local_index_pt=global_to_local.Get_Pointer(send_indices(i)(j));
        PHYSBAM_ASSERT(local_index_pt);
        send_indices(i)(j)=*local_index_pt;}
}
//#####################################################################
// Function Fill_Ghost_Values
//#####################################################################
template<class TV> void MPI_HYBRID_COUPLED_SYSTEM<TV>::
Fill_Ghost_Values(VECTOR_ND<T>& v)const
{
    ARRAY<MPI_PACKAGE> packages;ARRAY<MPI::Request> requests;
    // send
    ARRAY<INDIRECT_ARRAY<ARRAY_VIEW<T> >* > send_values(neighbor_ranks.m);
    for(int n=1;n<=neighbor_ranks.m;n++)if(neighbor_ranks(n)!=MPI::PROC_NULL && send_indices(n).m>0){
        send_values(n)=new INDIRECT_ARRAY<ARRAY_VIEW<T> >(v,send_indices(n));
        MPI_PACKAGE package(*send_values(n));
        packages.Append(package);requests.Append(package.Isend(*mpi_grid.comm,neighbor_ranks(n),mpi_grid.Get_Send_Tag(neighbor_directions(n))));}
    // receive
    ARRAY<INDIRECT_ARRAY<ARRAY_VIEW<T> >* > recv_values(neighbor_ranks.m);
    for(int n=1;n<=neighbor_ranks.m;n++)if(neighbor_ranks(n)!=MPI::PROC_NULL && recv_indices(n).m>0){
        recv_values(n)=new INDIRECT_ARRAY<ARRAY_VIEW<T> >(v,recv_indices(n));
        MPI_PACKAGE package(*recv_values(n));
        packages.Append(package);requests.Append(package.Irecv(*mpi_grid.comm,neighbor_ranks(n),mpi_grid.Get_Recv_Tag(neighbor_directions(n))));}
    // finish
    MPI_UTILITIES::Wait_All(requests);MPI_PACKAGE::Free_All(packages);
    send_values.Delete_Pointers_And_Clean_Memory();
    recv_values.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void MPI_HYBRID_COUPLED_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& bp,KRYLOV_VECTOR_BASE<T>& bresult) const
{
    const KRYLOV_VECTOR_T& p=debug_cast<const KRYLOV_VECTOR_T&>(bp);KRYLOV_VECTOR_T& result=debug_cast<KRYLOV_VECTOR_T&>(bresult);
    const KRYLOV_VECTOR_T* use_p=&p;
    Fill_Ghost_Values(use_p->v);
    A.Times(use_p->v,result.v);
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double MPI_HYBRID_COUPLED_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& bx,const KRYLOV_VECTOR_BASE<T>& by) const
{
    const KRYLOV_VECTOR_T& x=debug_cast<const KRYLOV_VECTOR_T&>(bx),&y=debug_cast<const KRYLOV_VECTOR_T&>(by);
    double local_sum=0,global_sum;
    for(int i=1;i<=interior.m;i++){int j=interior(i);
        local_sum+=(double)x.v(j)*(double)y.v(j);}
    MPI_UTILITIES::Reduce(local_sum,global_sum,MPI::SUM,*mpi_grid.comm);
    return global_sum;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR MPI_HYBRID_COUPLED_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& bx) const
{
    const KRYLOV_VECTOR_T& x=debug_cast<const KRYLOV_VECTOR_T&>(bx);
    T local_max=0,global_max;
    for(int i=1;i<=interior.m;i++){int j=interior(i);
        local_max=PhysBAM::max(local_max,abs(x.v(j)));}
    MPI_UTILITIES::Reduce(local_max,global_max,MPI::MAX,*mpi_grid.comm);
    return global_max;
}
//#####################################################################
template class MPI_HYBRID_COUPLED_SYSTEM<VECTOR<float,1> >;
template class MPI_HYBRID_COUPLED_SYSTEM<VECTOR<float,2> >;
template class MPI_HYBRID_COUPLED_SYSTEM<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MPI_HYBRID_COUPLED_SYSTEM<VECTOR<double,1> >;
template class MPI_HYBRID_COUPLED_SYSTEM<VECTOR<double,2> >;
template class MPI_HYBRID_COUPLED_SYSTEM<VECTOR<double,3> >;
#endif
#endif

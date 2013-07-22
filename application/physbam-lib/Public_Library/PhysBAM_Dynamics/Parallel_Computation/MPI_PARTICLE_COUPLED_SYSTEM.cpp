//#####################################################################
// Copyright 2011, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_PARTICLE_COUPLED_SYSTEM
//#####################################################################
#ifdef USE_MPI
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#endif
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_PARTICLE_COUPLED_SYSTEM.h>
using namespace PhysBAM;

#ifdef USE_MPI
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPI_PARTICLE_COUPLED_SYSTEM<TV>::
MPI_PARTICLE_COUPLED_SYSTEM(const SPARSE_MATRIX_FLAT_MXN<T>& K_input,const VECTOR_ND<T>& Mi_input,VECTOR_ND<T>& temp_v_input,const VECTOR<int,2>& diagonal_range_input,const ARRAY<int>& filter_input,const ARRAY<int>& mapping_input,const ARRAY<int>& back_mapping_input,MPI_PARTICLES<GRID<TV> >& mpi_particles_input)
    :KRYLOV_SYSTEM_BASE<T>(false,false),K(K_input),Mi(Mi_input),temp_v(temp_v_input),diagonal_range(diagonal_range_input),filter(filter_input),mapping(mapping_input),back_mapping(back_mapping_input),mpi_particles(mpi_particles_input)
{
    int pn=mapping.m;
    for(int n=diagonal_range(1);n<=diagonal_range(2);n++){// surface tension segments
        int j1=K.A(K.offsets(n)).j,j2=K.A(K.offsets(n)+1).j;
        int k1=(j1-1)/pn+1,k2=(j2-1)/pn+1;
        int i1=mapping(k1),i2=mapping(k2);
        if(filter(i1)!=filter(i2)) cross_over_indices.Append(n);
        else if(filter(i1)==0 && filter(i2)==0) interior_indices.Append(n);}
    for(int n=1;n<=pn;n++){int i=mapping(n);
        if(filter(i)==0) interior_indices.Append(n+diagonal_range(2));}
    // recv
    recv_indices.Resize(mpi_particles.send_regions.m);
    recv_indices_scalar.Resize(mpi_particles.send_regions.m);
    backward_filter.Resize(back_mapping.m);ARRAYS_COMPUTATIONS::Fill(backward_filter,0);
    for(int n=1;n<=mpi_particles.send_regions.m;n++)if(mpi_particles.neighbor_ranks(n)!=MPI::PROC_NULL && !mpi_particles.recv_particle_range(n).Empty()){
        for(int i=mpi_particles.recv_particle_range(n).min_corner(1);i<=mpi_particles.recv_particle_range(n).max_corner(1);i++){
            if(filter(i)!=1) continue;// receive locally ghost values
            int k=back_mapping(i);
            if(k==0) continue;// must be an unknown
            backward_filter(i)=1;
            recv_indices_scalar(n).Append(k);}
        for(int axis=1;axis<=TV::m;axis++) for(int i=1;i<=recv_indices_scalar(n).m;i++) recv_indices(n).Append((axis-1)*pn+recv_indices_scalar(n)(i));
        for(int i=1;i<=recv_indices_scalar(n).m;i++) recv_indices_scalar(n)(i)+=diagonal_range(2);}
    // send
    ARRAY_VIEW<int> backward_filter_view(backward_filter);
    ARRAY<ARRAY<int> > backward_filters(mpi_particles.send_regions.m);
    ARRAY<ARRAY_VIEW<int>*> backward_filter_views(mpi_particles.send_regions.m);
    for(int n=1;n<=mpi_particles.send_regions.m;n++) if(mpi_particles.neighbor_ranks(n)!=MPI::PROC_NULL && mpi_particles.send_particles(n).m>0){
        backward_filters(n).Resize(backward_filter.m);ARRAYS_COMPUTATIONS::Fill(backward_filters(n),0);
        backward_filter_views(n)=new ARRAY_VIEW<int>(backward_filters(n));}
    else backward_filter_views(n)=0;
    mpi_particles.Backwards_Exchange_Boundary_Particle_Values(backward_filter_view,backward_filter_views);
    send_indices.Resize(mpi_particles.send_regions.m);
    send_indices_scalar.Resize(mpi_particles.send_regions.m);
    for(int n=1;n<=mpi_particles.send_regions.m;n++)if(mpi_particles.neighbor_ranks(n)!=MPI::PROC_NULL && mpi_particles.send_particles(n).m>0){
        for(int j=1;j<=mpi_particles.send_particles(n).m;j++){
            int i=mpi_particles.send_particles(n)(j);
            int k=back_mapping(i);
            if(!(k>0 && filter(i)==0 && backward_filters(n)(i)==1)) continue;// send locally interior and remotely ghost values
            send_indices_scalar(n).Append(k);}
        for(int axis=1;axis<=TV::m;axis++) for(int i=1;i<=send_indices_scalar(n).m;i++) send_indices(n).Append((axis-1)*pn+send_indices_scalar(n)(i));
        for(int i=1;i<=send_indices_scalar(n).m;i++) send_indices_scalar(n)(i)+=diagonal_range(2);}
    backward_filter_views.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Fill_Scalar_Ghost_Values
//#####################################################################
template<class TV> void MPI_PARTICLE_COUPLED_SYSTEM<TV>::
Fill_Scalar_Ghost_Values(VECTOR_ND<T>& v)const
{
    ARRAY<MPI_PACKAGE> packages;ARRAY<MPI::Request> requests;
    // send
    for(int n=1;n<=mpi_particles.send_regions.m;n++)if(mpi_particles.neighbor_ranks(n)!=MPI::PROC_NULL && send_indices_scalar(n).m>0){
        MPI_PACKAGE package(v,send_indices_scalar(n));
        packages.Append(package);requests.Append(package.Isend(*mpi_particles.mpi_grid.comm,mpi_particles.neighbor_ranks(n),mpi_particles.mpi_grid.Get_Send_Tag(mpi_particles.neighbor_directions(n))));}
    // receive
    for(int n=1;n<=mpi_particles.send_regions.m;n++)if(mpi_particles.neighbor_ranks(n)!=MPI::PROC_NULL && recv_indices_scalar(n).m>0){
        MPI_PACKAGE package(v,recv_indices_scalar(n));
        packages.Append(package);requests.Append(package.Irecv(*mpi_particles.mpi_grid.comm,mpi_particles.neighbor_ranks(n),mpi_particles.mpi_grid.Get_Recv_Tag(mpi_particles.neighbor_directions(n))));}
    // finish
    MPI_UTILITIES::Wait_All(requests);MPI_PACKAGE::Free_All(packages);
}
//#####################################################################
// Function Fill_Vector_Ghost_Values
//#####################################################################
template<class TV> void MPI_PARTICLE_COUPLED_SYSTEM<TV>::
Fill_Vector_Ghost_Values(VECTOR_ND<T>& v)const
{
    ARRAY<MPI_PACKAGE> packages;ARRAY<MPI::Request> requests;
    // send
    for(int n=1;n<=mpi_particles.send_regions.m;n++)if(mpi_particles.neighbor_ranks(n)!=MPI::PROC_NULL && send_indices(n).m>0){
        MPI_PACKAGE package(v,send_indices(n));
        packages.Append(package);requests.Append(package.Isend(*mpi_particles.mpi_grid.comm,mpi_particles.neighbor_ranks(n),mpi_particles.mpi_grid.Get_Send_Tag(mpi_particles.neighbor_directions(n))));}
    // receive
    for(int n=1;n<=mpi_particles.send_regions.m;n++)if(mpi_particles.neighbor_ranks(n)!=MPI::PROC_NULL && recv_indices(n).m>0){
        MPI_PACKAGE package(v,recv_indices(n));
        packages.Append(package);requests.Append(package.Irecv(*mpi_particles.mpi_grid.comm,mpi_particles.neighbor_ranks(n),mpi_particles.mpi_grid.Get_Recv_Tag(mpi_particles.neighbor_directions(n))));}
    // finish
    MPI_UTILITIES::Wait_All(requests);MPI_PACKAGE::Free_All(packages);
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void MPI_PARTICLE_COUPLED_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& bp,KRYLOV_VECTOR_BASE<T>& bresult) const
{
    const KRYLOV_VECTOR_T& p=debug_cast<const KRYLOV_VECTOR_T&>(bp);KRYLOV_VECTOR_T& result=debug_cast<KRYLOV_VECTOR_T&>(bresult);
    const KRYLOV_VECTOR_T* use_p=&p;
    Fill_Scalar_Ghost_Values(use_p->v);
    K.Transpose_Times(use_p->v,temp_v);
    temp_v*=Mi;
    Fill_Vector_Ghost_Values(temp_v);
    K.Times(temp_v,result.v);
    for(int i=diagonal_range(1);i<=diagonal_range(2);i++) result.v(i)+=use_p->v(i);
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double MPI_PARTICLE_COUPLED_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& bx,const KRYLOV_VECTOR_BASE<T>& by) const
{
    const KRYLOV_VECTOR_T& x=debug_cast<const KRYLOV_VECTOR_T&>(bx),&y=debug_cast<const KRYLOV_VECTOR_T&>(by);
    double local_sum=0,global_sum;
    for(int i=1;i<=cross_over_indices.m;i++){int j=cross_over_indices(i);
        local_sum+=(double)0.5*(double)x.v(j)*(double)y.v(j);}
    for(int i=1;i<=interior_indices.m;i++){int j=interior_indices(i);
        local_sum+=(double)x.v(j)*(double)y.v(j);}
    MPI_UTILITIES::Reduce(local_sum,global_sum,MPI::SUM,*mpi_particles.mpi_grid.comm);
    return global_sum;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR MPI_PARTICLE_COUPLED_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& bx) const
{
    const KRYLOV_VECTOR_T& x=debug_cast<const KRYLOV_VECTOR_T&>(bx);
    T local_max=0,global_max;
    for(int i=1;i<=cross_over_indices.m;i++){int j=cross_over_indices(i);
        local_max=PhysBAM::max(local_max,abs(x.v(j)));}
    for(int i=1;i<=interior_indices.m;i++){int j=interior_indices(i);
        local_max=PhysBAM::max(local_max,abs(x.v(j)));}
    MPI_UTILITIES::Reduce(local_max,global_max,MPI::MAX,*mpi_particles.mpi_grid.comm);
    return global_max;
}
//#####################################################################
template class MPI_PARTICLE_COUPLED_SYSTEM<VECTOR<float,1> >;
template class MPI_PARTICLE_COUPLED_SYSTEM<VECTOR<float,2> >;
template class MPI_PARTICLE_COUPLED_SYSTEM<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MPI_PARTICLE_COUPLED_SYSTEM<VECTOR<double,1> >;
template class MPI_PARTICLE_COUPLED_SYSTEM<VECTOR<double,2> >;
template class MPI_PARTICLE_COUPLED_SYSTEM<VECTOR<double,3> >;
#endif
#endif

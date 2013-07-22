//#####################################################################
// Copyright 2008, Avi Robinson-Mosher, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_SOLID_FLUID
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <PhysBAM_Dynamics/Parallel_Computation/CONJUGATE_RESIDUAL_SPARSE_MPI.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#ifdef USE_MPI
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#endif
namespace PhysBAM{

#ifdef USE_MPI

//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPI_SOLID_FLUID<TV>::
MPI_SOLID_FLUID()
    :number_of_solid_processes(1),solid_node(0),comm(0),group(0),solid_group(0),fluid_group(0),fluid_ranks(0),current_tag(0)
{
    LOG::SCOPE scope("MPI INITIALIZE","Initializing MPI_SOLID_FLUID");
    number_of_processes=MPI::COMM_WORLD.Get_size();
    std::stringstream ss1;
    ss1<<"number of processes = "<<number_of_processes<<std::endl;
    LOG::filecout(ss1.str());

    comm=new MPI::Intracomm;
    *comm=MPI::COMM_WORLD.Dup();
    rank=comm->Get_rank();
    std::stringstream ss;
    ss<<"global rank = "<<rank<<std::endl;
    LOG::filecout(ss.str());

    group=new MPI::Group(comm->Get_group());
    solid_ranks.Resize(number_of_solid_processes);
    for(int i=1;i<=number_of_solid_processes;i++)solid_ranks(i)=i-1;
    solid_group=new MPI::Group(group->Incl(solid_ranks.n,&solid_ranks(1)));
    fluid_ranks.Resize(number_of_processes-number_of_solid_processes);
    for(int i=1;i<=fluid_ranks.n;i++)fluid_ranks(i)=i+number_of_solid_processes-1;
    fluid_group=new MPI::Group(group->Incl(fluid_ranks.n,&fluid_ranks(1)));
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPI_SOLID_FLUID<TV>::
~MPI_SOLID_FLUID()
{
    if(comm){comm->Free();delete comm;}
    if(group){group->Free();delete group;}
    if(solid_group){solid_group->Free();delete solid_group;}
}
//#####################################################################
// Function Exchange_Solid_Positions_And_Velocities
//#####################################################################
namespace {
template<class T> void Exchange_Solid_Positions_And_Velocities_Helper(const MPI_SOLID_FLUID<VECTOR<T,1> >& mpi,SOLID_BODY_COLLECTION<VECTOR<T,1> >& solid_body_collection)
{
    int tag=mpi.Get_Unique_Tag();
    PARTICLES<VECTOR<T,1> >& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_PARTICLES<VECTOR<T,1> >& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particle;
    if(mpi.Solid_Node()){
        ARRAY<ARRAY<char> > send_buffers(mpi.fluid_ranks.n);ARRAY<MPI::Request> requests;
        for(int i=1;i<=mpi.fluid_ranks.n;i++){
            int buffer_size=MPI_UTILITIES::Pack_Size(particles.X,rigid_body_particles.X,particles.V,rigid_body_particles.angular_velocity,*mpi.comm)+1;
            send_buffers(i).Resize(buffer_size);int position=0;
            MPI_UTILITIES::Pack(particles.X,rigid_body_particles.X,particles.V,rigid_body_particles.angular_velocity,send_buffers(i),position,*mpi.comm);
            requests.Append(mpi.comm->Isend(&(send_buffers(i)(1)),position,MPI::PACKED,mpi.fluid_ranks(i),tag));}
        MPI_UTILITIES::Wait_All(requests);}
    else{
        int buffer_size=MPI_UTILITIES::Pack_Size(particles.X,rigid_body_particles.X,particles.V,rigid_body_particles.angular_velocity,*mpi.comm)+1;
        ARRAY<char> buffer(buffer_size);int position=0;
        mpi.comm->Recv(&buffer(1),buffer_size,MPI::PACKED,mpi.solid_node,tag);
        MPI_UTILITIES::Unpack(particles.X,rigid_body_particles.X,particles.V,rigid_body_particles.angular_velocity,buffer,position,*mpi.comm);}
}
};
template<> void MPI_SOLID_FLUID<VECTOR<float,1> >::Exchange_Solid_Positions_And_Velocities(SOLID_BODY_COLLECTION<VECTOR<float,1> >& solid_body_collection) const {Exchange_Solid_Positions_And_Velocities_Helper(*this,solid_body_collection);}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template<> void MPI_SOLID_FLUID<VECTOR<double,1> >::Exchange_Solid_Positions_And_Velocities(SOLID_BODY_COLLECTION<VECTOR<double,1> >& solid_body_collection) const {Exchange_Solid_Positions_And_Velocities_Helper(*this,solid_body_collection);}
#endif
template<class TV> void MPI_SOLID_FLUID<TV>::
Exchange_Solid_Positions_And_Velocities(SOLID_BODY_COLLECTION<TV>& solid_body_collection) const
{
    int tag=Get_Unique_Tag();
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particle;
    if(Solid_Node()){
        ARRAY<ARRAY<char> > send_buffers(fluid_ranks.n);ARRAY<MPI::Request> requests;
        for(int i=1;i<=fluid_ranks.n;i++){
            int buffer_size=MPI_UTILITIES::Pack_Size(particles.X,rigid_body_particles.X,rigid_body_particles.rotation,particles.V,rigid_body_particles.angular_velocity,
                rigid_body_particles.angular_momentum,*comm)+1;
            send_buffers(i).Resize(buffer_size);int position=0;
            MPI_UTILITIES::Pack(particles.X,rigid_body_particles.X,rigid_body_particles.rotation,particles.V,rigid_body_particles.angular_velocity,rigid_body_particles.angular_momentum,
                send_buffers(i),position,*comm);
            requests.Append(comm->Isend(&(send_buffers(i)(1)),position,MPI::PACKED,fluid_ranks(i),tag));}
        MPI_UTILITIES::Wait_All(requests);}
    else{
        int buffer_size=MPI_UTILITIES::Pack_Size(particles.X,rigid_body_particles.X,rigid_body_particles.rotation,particles.V,rigid_body_particles.angular_velocity,
            rigid_body_particles.angular_momentum,*comm)+1;
        ARRAY<char> buffer(buffer_size);int position=0;
        comm->Recv(&buffer(1),buffer_size,MPI::PACKED,solid_node,tag);
        MPI_UTILITIES::Unpack(particles.X,rigid_body_particles.X,rigid_body_particles.rotation,particles.V,rigid_body_particles.angular_velocity,rigid_body_particles.angular_momentum,
            buffer,position,*comm);}
}
//#####################################################################
// Function Fluid_Node
//#####################################################################
template<class TV> bool MPI_SOLID_FLUID<TV>::
Fluid_Node() const
{
    return fluid_group->Get_rank()!=MPI::UNDEFINED;
}
//#####################################################################
// Function Solid_Node
//#####################################################################
template<class TV> bool MPI_SOLID_FLUID<TV>::
Solid_Node() const
{
    return solid_group->Get_rank()!=MPI::UNDEFINED;
}
//#####################################################################
// Function Create_Fluid_Comm_For_Solid_Nodes
//#####################################################################
template<class TV> void MPI_SOLID_FLUID<TV>::
Create_Fluid_Comm_For_Solid_Nodes() const
{
    MPI::COMM_WORLD.Create(*fluid_group);
}
//#####################################################################
// Function Reduce_Add
//#####################################################################
template<class TV> template<class T2> void MPI_SOLID_FLUID<TV>::
Reduce_Add(const T2& input,T2& output) const
{
    MPI_UTILITIES::Reduce(input,output,MPI::SUM,*comm);
}
//#####################################################################
// Function Reduce_Min
//#####################################################################
template<class TV> typename TV::SCALAR MPI_SOLID_FLUID<TV>::
Reduce_Min(const T local_value) const
{
    T global_value;
    MPI_UTILITIES::Reduce(local_value,global_value,MPI::MIN,*comm);
    return global_value;
}
//#####################################################################
// Function Reduce_Max
//#####################################################################
template<class TV> typename TV::SCALAR MPI_SOLID_FLUID<TV>::
Reduce_Max(const T local_value) const
{
    T global_value;
    MPI_UTILITIES::Reduce(local_value,global_value,MPI::MAX,*comm);
    return global_value;
}
//#####################################################################
// Function Parallel_Solve_Fluid_Part
//#####################################################################
template<class TV> void MPI_SOLID_FLUID<TV>::
Parallel_Solve_Fluid_Part(FLUID_SYSTEM_MPI<TV>& fluid_system,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& x_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& b_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& p_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& ap_array,
    KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& ar_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& r_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& z_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& zaq_array,const int min_iterations,const int max_iterations,
    const T tolerance,const bool recompute_preconditioner,ARRAY<MPI::Intracomm>* fluid_comm,ARRAY<SPARSE_MATRIX_PARTITION>* partitions)
{
    CONJUGATE_RESIDUAL_SPARSE_MPI<TV> cr_mpi(*comm,fluid_comm,partitions);
    cr_mpi.nullspace_tolerance=(T)0;
    cr_mpi.restart_iterations=100;
    //cr_mpi.maximum_iterations=solids_parameters.cg_iterations;
    cr_mpi.Parallel_Solve_Fluid_Part(fluid_system,x_array,b_array,p_array,ap_array,ar_array,r_array,z_array,zaq_array,min_iterations,max_iterations,tolerance,recompute_preconditioner);
}
//#####################################################################
// Function Parallel_Solve_Solid_Part
//#####################################################################
template<class TV> void MPI_SOLID_FLUID<TV>::
Parallel_Solve_Solid_Part(SOLID_SYSTEM_MPI<TV>& solid_system,GENERALIZED_VELOCITY<TV>& x_array,GENERALIZED_VELOCITY<TV>& b_array,GENERALIZED_VELOCITY<TV>& p_array,GENERALIZED_VELOCITY<TV>& ap_array,
    GENERALIZED_VELOCITY<TV>& ar_array,GENERALIZED_VELOCITY<TV>& r_array,GENERALIZED_VELOCITY<TV>& z_array,GENERALIZED_VELOCITY<TV>& zaq_array,const int min_iterations,const int max_iterations,const T tolerance)
{
    CONJUGATE_RESIDUAL_SPARSE_MPI<TV> cr_mpi(*comm,0,0);
    cr_mpi.nullspace_tolerance=(T)0;
    cr_mpi.restart_iterations=100;
    //cr_mpi.maximum_iterations=solids_parameters.cg_iterations;
    cr_mpi.Parallel_Solve_Solid_Part(solid_system,x_array,b_array,p_array,ap_array,ar_array,r_array,z_array,zaq_array,min_iterations,max_iterations,tolerance);
}
//#####################################################################
// Function Exchange_Coupled_Deformable_Particle_List
//#####################################################################
template<class TV> void MPI_SOLID_FLUID<TV>::
Exchange_Coupled_Deformable_Particle_List(ARRAY<int>* fluid_list,ARRAY<ARRAY<int> >* results)
{
    int tag=Get_Unique_Tag();
    if(Solid_Node()){
        for(int i=1;i<=fluid_ranks.n;i++){
            MPI::Status status;
            comm->Probe(MPI::ANY_SOURCE,tag,status);
            int source=status.Get_source();
            ARRAY<char> buffer(status.Get_count(MPI::PACKED));int position=0;
            comm->Recv(&buffer(1),buffer.m,MPI::PACKED,source,tag);
            MPI_UTILITIES::Unpack((*results)(source),buffer,position,*comm);}}
    else{
        int buffer_size=MPI_UTILITIES::Pack_Size(*fluid_list,*comm)+1;
        ARRAY<char> buffer(buffer_size);int position=0;
        MPI_UTILITIES::Pack(*fluid_list,buffer,position,*comm);
        comm->Send(&buffer(1),buffer_size,MPI::PACKED,solid_node,tag);}
}
//#####################################################################
// Function Distribute_Lists_From_Solid_Node
//#####################################################################
template<class TV> void MPI_SOLID_FLUID<TV>::
Distribute_Lists_From_Solid_Node(GENERALIZED_VELOCITY<TV>& F) const
{
    int tag=Get_Unique_Tag();
    if(Solid_Node()){
        ARRAY<ARRAY<char> > send_buffers(fluid_ranks.n);ARRAY<MPI::Request> requests;
        for(int i=1;i<=fluid_ranks.n;i++){
            int buffer_size=MPI_UTILITIES::Pack_Size(F.V,F.rigid_V,*comm)+1;
            send_buffers(i).Resize(buffer_size);int position=0;
            MPI_UTILITIES::Pack(F.V,F.rigid_V,send_buffers(i),position,*comm);
            requests.Append(comm->Isend(&(send_buffers(i)(1)),position,MPI::PACKED,fluid_ranks(i),tag));}
        MPI_UTILITIES::Wait_All(requests);}
    else{
        MPI::Status status;
        comm->Probe(MPI::ANY_SOURCE,tag,status);
        ARRAY<char> buffer(status.Get_count(MPI::PACKED));int position=0;
        comm->Recv(&buffer(1),buffer.m,MPI::PACKED,solid_node,tag);
        MPI_UTILITIES::Unpack(F.V,F.rigid_V,buffer,position,*comm);}
}
//#####################################################################
// Function Aggregate_Lists_To_Solid_Node
//#####################################################################
template<class TV> void MPI_SOLID_FLUID<TV>::
Aggregate_Lists_To_Solid_Node(GENERALIZED_VELOCITY<TV>& F)
{
    int F_tag=Get_Unique_Tag(),rigid_F_tag=Get_Unique_Tag();
    ARRAY<TV> F_buffer;F_buffer.Resize(F.V.array.m);
    ARRAY<TWIST<TV> > rigid_F_buffer;rigid_F_buffer.Resize(F.rigid_V.array.m);
    if(Solid_Node()){
        for(int i=1;i<=fluid_ranks.n;++i){
            { // Receive F_buffer
                MPI::Status status;
                comm->Probe(MPI::ANY_SOURCE,F_tag,status);
                int source=status.Get_source();
                ARRAY<char> buffer(status.Get_count(MPI::PACKED));int position=0;
                comm->Recv(&buffer(1),buffer.m,MPI::PACKED,source,F_tag);
                MPI_UTILITIES::Unpack(F_buffer,buffer,position,*comm);}
            { // Receive rigid_F_buffer
                MPI::Status status;
                comm->Probe(MPI::ANY_SOURCE,rigid_F_tag,status);
                int source=status.Get_source();
                ARRAY<char> buffer(status.Get_count(MPI::PACKED));int position=0;
                comm->Recv(&buffer(1),buffer.m,MPI::PACKED,source,rigid_F_tag);
                MPI_UTILITIES::Unpack(rigid_F_buffer,buffer,position,*comm);}
            for(int j=1;j<=F.V.array.m;++j) F.V.array(j)+=F_buffer(j);
            for(int j=1;j<=F.rigid_V.array.m;++j) F.rigid_V.array(j)+=rigid_F_buffer(j);}}
    else{
        { // Send F_buffer
            int buffer_size=MPI_UTILITIES::Pack_Size(F.V.array,*comm)+1;
            ARRAY<char> buffer(buffer_size);int position=0;
            MPI_UTILITIES::Pack(F.V.array,buffer,position,*comm);
            comm->Send(&buffer(1),buffer_size,MPI::PACKED,solid_node,F_tag);}
        { // Send rigid_F_buffer
            int buffer_size=MPI_UTILITIES::Pack_Size(F.rigid_V.array,*comm)+1;
            ARRAY<char> buffer(buffer_size);int position=0;
            MPI_UTILITIES::Pack(F.rigid_V.array,buffer,position,*comm);
            comm->Send(&buffer(1),buffer_size,MPI::PACKED,solid_node,rigid_F_tag);}}
}
//#####################################################################
#else
//#####################################################################
template<class TV> MPI_SOLID_FLUID<TV>::MPI_SOLID_FLUID(){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> MPI_SOLID_FLUID<TV>::~MPI_SOLID_FLUID(){}
template<class TV> void MPI_SOLID_FLUID<TV>::Exchange_Solid_Positions_And_Velocities(SOLID_BODY_COLLECTION<TV>& solid_body_collection) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool MPI_SOLID_FLUID<TV>::Fluid_Node() const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool MPI_SOLID_FLUID<TV>::Solid_Node() const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> typename TV::SCALAR MPI_SOLID_FLUID<TV>::Reduce_Max(const T local_value) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_SOLID_FLUID<TV>::Create_Fluid_Comm_For_Solid_Nodes() const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> template<class T2> void MPI_SOLID_FLUID<TV>::Reduce_Add(const T2& input,T2& output) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> typename TV::SCALAR MPI_SOLID_FLUID<TV>::Reduce_Min(const T local_value) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_SOLID_FLUID<TV>::Parallel_Solve_Fluid_Part(FLUID_SYSTEM_MPI<TV>& fluid_system,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& x_array, \
    KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& b_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& p_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& ap_array, \
    KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& ar_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& r_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& z_array, \
    KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& zaq_array,const int min_iterations,const int max_iterations, const T tolerance,const bool recompute_preconditioner, \
    ARRAY<MPI::Intracomm>* fluid_comm,ARRAY<SPARSE_MATRIX_PARTITION>* partitions){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_SOLID_FLUID<TV>::Parallel_Solve_Solid_Part(SOLID_SYSTEM_MPI<TV>& solid_system,GENERALIZED_VELOCITY<TV>& x_array,GENERALIZED_VELOCITY<TV>& b_array, \
    GENERALIZED_VELOCITY<TV>& p_array,GENERALIZED_VELOCITY<TV>& ap_array,GENERALIZED_VELOCITY<TV>& ar_array,GENERALIZED_VELOCITY<TV>& r_array,GENERALIZED_VELOCITY<TV>& z_array, \
    GENERALIZED_VELOCITY<TV>& zaq_array,const int min_iterations,const int max_iterations,const T tolerance){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_SOLID_FLUID<TV>::Distribute_Lists_From_Solid_Node(GENERALIZED_VELOCITY<TV>& F) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_SOLID_FLUID<TV>::Exchange_Coupled_Deformable_Particle_List(ARRAY<int>* fluid_list,ARRAY<ARRAY<int> >* results){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_SOLID_FLUID<TV>::Aggregate_Lists_To_Solid_Node(GENERALIZED_VELOCITY<TV>& F){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
#endif
template void MPI_SOLID_FLUID<VECTOR<float,1> >::Reduce_Add(const double&,double&) const;
template void MPI_SOLID_FLUID<VECTOR<float,2> >::Reduce_Add(const double&,double&) const;
template void MPI_SOLID_FLUID<VECTOR<float,3> >::Reduce_Add(const double&,double&) const;

#define INSTANTIATION_HELPER_UNIFORM(T,d) \
    template class MPI_SOLID_FLUID<VECTOR<T,d> >; \
    template void MPI_SOLID_FLUID<VECTOR<T,d> >::Reduce_Add(const T&,T&) const; \
    template void MPI_SOLID_FLUID<VECTOR<T,d> >::Reduce_Add(const ARRAY<T>&,ARRAY<T>&) const;
#define INSTANTIATION_HELPER(T) \
    INSTANTIATION_HELPER_UNIFORM(T,1);INSTANTIATION_HELPER_UNIFORM(T,2);INSTANTIATION_HELPER_UNIFORM(T,3);
INSTANTIATION_HELPER(float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double);
#endif
}

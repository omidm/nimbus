//#####################################################################
// Copyright 2008-2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_SYSTEM_MPI_SLIP
//#####################################################################
#include <PhysBAM_Tools/Arrays/IDENTITY_ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Vectors/Dot_Product.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID_SLIP.h>
#include <PhysBAM_Dynamics/Parallel_Computation/SOLID_SYSTEM_MPI_SLIP.h>
#ifdef USE_MPI
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#endif
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLID_SYSTEM_MPI_SLIP<TV>::
SOLID_SYSTEM_MPI_SLIP(const bool use_preconditioner_input,BACKWARD_EULER_SYSTEM<TV>& solid_system_input,MPI_SOLID_FLUID_SLIP<TV>* mpi_solid_fluid_input,
    ARRAY<ARRAY<int> >& coupled_deformable_particle_indices_input,NEWMARK_EVOLUTION<TV>& newmark_evolution_input,const int rigid_V_size)
    :BASE(use_preconditioner_input,!use_preconditioner_input),solid_system(solid_system_input),mpi_solid_fluid(mpi_solid_fluid_input),
    coupled_deformable_particle_indices(coupled_deformable_particle_indices_input),solid_mass(0),newmark_evolution(newmark_evolution_input)
{
    recv_fluid_V_boundary_arrays.Resize(coupled_deformable_particle_indices.m);
    recv_fluid_rigid_V_boundary_arrays.Resize(coupled_deformable_particle_indices.m);
    solid_mass=&solid_system.projection_data.mass;
    for(int i=1;i<=coupled_deformable_particle_indices.m;i++){
        recv_fluid_V_boundary_arrays(i).Resize(coupled_deformable_particle_indices(i).m);
        recv_fluid_rigid_V_boundary_arrays(i).Resize(rigid_V_size);}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLID_SYSTEM_MPI_SLIP<TV>::
~SOLID_SYSTEM_MPI_SLIP()
{
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void SOLID_SYSTEM_MPI_SLIP<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& BV) const
{
    VECTOR_T& V=debug_cast<VECTOR_T&>(BV);
    solid_system.Set_Global_Boundary_Conditions(V,newmark_evolution.X_save,newmark_evolution.rigid_X_save,newmark_evolution.rigid_rotation_save,newmark_evolution.rigid_velocity_save,
        newmark_evolution.rigid_angular_momentum_save,newmark_evolution.V_save,newmark_evolution.solids_parameters.implicit_solve_parameters.test_system,
        newmark_evolution.solids_parameters.implicit_solve_parameters.print_matrix);
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void SOLID_SYSTEM_MPI_SLIP<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF) const
{
    const VECTOR_T& V=debug_cast<const VECTOR_T&>(BV);VECTOR_T& F=debug_cast<VECTOR_T&>(BF);
    solid_system.Force(V,F);
    Send_Generalized_Velocity_To_Fluid(V); // MPI
    for(int i=1;i<=V.V.Size();i++)
        F.V(i)=F.V(i)*solid_system.dt-(*solid_mass).mass(i)*V.V(i);
    for(int i=1;i<=V.rigid_V.Size();i++)
        F.rigid_V(i)=F.rigid_V(i)*solid_system.dt-(*solid_mass).world_space_rigid_mass(i)*V.rigid_V(i);

    Get_Generalized_Velocity_From_Fluid(F);
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double SOLID_SYSTEM_MPI_SLIP<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& BV1,const KRYLOV_VECTOR_BASE<T>& BV2) const
{
    const VECTOR_T& V1=debug_cast<const VECTOR_T&>(BV1);const VECTOR_T& V2=debug_cast<const VECTOR_T&>(BV2);
    double inner_product_solid_velocities=0;
    for(int i=1;i<=V1.V.Size();i++) inner_product_solid_velocities+=Dot_Product(V1.V(i),V2.V(i));
    for(int i=1;i<=V1.rigid_V.Size();i++){
        inner_product_solid_velocities+=Dot_Product(V1.rigid_V(i).linear,V2.rigid_V(i).linear);
        inner_product_solid_velocities+=Dot_Product(V1.rigid_V(i).angular,V2.rigid_V(i).angular);}
    return inner_product_solid_velocities;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR SOLID_SYSTEM_MPI_SLIP<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const
{
    const VECTOR_T& R=debug_cast<const VECTOR_T&>(BR);
    packed_solid_velocities.Resize(R.Raw_Size());
    R.Pack(packed_solid_velocities);
    T convergence_norm=packed_solid_velocities.Maximum_Magnitude();
    return convergence_norm;
}
#ifdef USE_MPI
//#####################################################################
// Function void Send_Generalized_Velocity_To_Solid
//#####################################################################
template<class TV> void SOLID_SYSTEM_MPI_SLIP<TV>::
Send_Generalized_Velocity_To_Fluid(const GENERALIZED_VELOCITY<TV>& V) const
{
    int tag_solid_to_fluid=mpi_solid_fluid->Get_Unique_Tag();
    ARRAY<ARRAY<char> > send_buffers(mpi_solid_fluid->fluid_ranks.n);ARRAY<MPI::Request> requests;
    for(int i=1;i<=mpi_solid_fluid->fluid_ranks.n;i++){
        const INDIRECT_ARRAY<ARRAY_VIEW<TV> > V_boundary(V.V.array,coupled_deformable_particle_indices(i));
        int buffer_size=MPI_UTILITIES::Pack_Size(V_boundary,V.rigid_V,*mpi_solid_fluid->comm)+1;
        send_buffers(i).Resize(buffer_size);int position=0;
        MPI_UTILITIES::Pack(V_boundary,V.rigid_V,send_buffers(i),position,*mpi_solid_fluid->comm);
        requests.Append(mpi_solid_fluid->comm->Isend(&(send_buffers(i)(1)),position,MPI::PACKED,mpi_solid_fluid->fluid_ranks(i),tag_solid_to_fluid));}
    MPI_UTILITIES::Wait_All(requests);
}
//#####################################################################
// Function Get_Generalized_Velocity_From_Fluid
//#####################################################################
template<class TV> void SOLID_SYSTEM_MPI_SLIP<TV>::
Get_Generalized_Velocity_From_Fluid(GENERALIZED_VELOCITY<TV>& V) const
{
    int tag_fluid_to_solid=mpi_solid_fluid->Get_Unique_Tag();
    ARRAY<ARRAY<char> > receive_buffers(mpi_solid_fluid->fluid_ranks.n);ARRAY<MPI::Request> requests;
    for(int i=1;i<=mpi_solid_fluid->fluid_ranks.n;i++){
        INDIRECT_ARRAY<ARRAY<TV>,IDENTITY_ARRAY<int> > fluid_V_boundary(recv_fluid_V_boundary_arrays(i),IDENTITY_ARRAY<int>(recv_fluid_V_boundary_arrays(i).Size()));
        INDIRECT_ARRAY<ARRAY<TWIST<TV> >,IDENTITY_ARRAY<int> > fluid_rigid_V_boundary(recv_fluid_rigid_V_boundary_arrays(i),IDENTITY_ARRAY<int>(recv_fluid_rigid_V_boundary_arrays(i).Size()));
        int buffer_size=MPI_UTILITIES::Pack_Size(fluid_V_boundary,fluid_rigid_V_boundary,*mpi_solid_fluid->comm)+1;
        receive_buffers(i).Resize(buffer_size);
        requests.Append(mpi_solid_fluid->comm->Irecv(&(receive_buffers(i)(1)),buffer_size,MPI::PACKED,i,tag_fluid_to_solid));}
    MPI_UTILITIES::Wait_All(requests);
    for(int i=1;i<=mpi_solid_fluid->fluid_ranks.n;i++){
        INDIRECT_ARRAY<ARRAY_VIEW<TV> > fluid_V_boundary(V.V.array,coupled_deformable_particle_indices(i));
        INDIRECT_ARRAY<ARRAY<TV>,IDENTITY_ARRAY<int> > recv_fluid_V_boundary(recv_fluid_V_boundary_arrays(i),IDENTITY_ARRAY<int>(recv_fluid_V_boundary_arrays(i).Size()));
        INDIRECT_ARRAY<ARRAY<TWIST<TV> >,IDENTITY_ARRAY<int> > recv_fluid_rigid_V_boundary(recv_fluid_rigid_V_boundary_arrays(i),IDENTITY_ARRAY<int>(recv_fluid_rigid_V_boundary_arrays(i).Size()));
        int position=0;
        MPI_UTILITIES::Unpack(recv_fluid_V_boundary,recv_fluid_rigid_V_boundary,receive_buffers(i),position,*mpi_solid_fluid->comm);
        fluid_V_boundary+=recv_fluid_V_boundary;V.rigid_V+=recv_fluid_rigid_V_boundary;}
}
#else
template<class TV> void SOLID_SYSTEM_MPI_SLIP<TV>::Send_Generalized_Velocity_To_Fluid(const GENERALIZED_VELOCITY<TV>& V) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void SOLID_SYSTEM_MPI_SLIP<TV>::Get_Generalized_Velocity_From_Fluid(GENERALIZED_VELOCITY<TV>& V) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
#endif
#define INSTANTIATION_HELPER(T,d) \
    template class SOLID_SYSTEM_MPI_SLIP<VECTOR<T,d> >;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif


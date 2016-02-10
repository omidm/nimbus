//#####################################################################
// Copyright 2008, Avi Robinson-Mosher, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_SYSTEM_MPI
//#####################################################################
#include <PhysBAM_Tools/Arrays/IDENTITY_ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/INNER_PRODUCT.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX_0X0.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <PhysBAM_Dynamics/Parallel_Computation/SOLID_SYSTEM_MPI.h>
#ifdef USE_MPI
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#endif
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLID_SYSTEM_MPI<TV>::
SOLID_SYSTEM_MPI(BACKWARD_EULER_SYSTEM<TV>& solid_system_input,ARRAY<T_DIAGONAL_MATRIX>& fluid_mass_input,ARRAY<T_DIAGONAL_MATRIX>& rigid_body_fluid_mass_input,
    ARRAY<T_INERTIA_TENSOR>& modified_world_space_rigid_inertia_tensor_input,MPI_SOLID_FLUID<TV>* mpi_solid_fluid_input,
    ARRAY<ARRAY<int> >& coupled_deformable_particle_indices_input,NEWMARK_EVOLUTION<TV>& newmark_evolution_input,const int rigid_V_size,bool precondition)
    :KRYLOV_SYSTEM_BASE<typename TV::SCALAR>(precondition,!precondition),
    solid_system(solid_system_input),fluid_mass(fluid_mass_input),rigid_body_fluid_mass(rigid_body_fluid_mass_input),
    modified_world_space_rigid_inertia_tensor(modified_world_space_rigid_inertia_tensor_input),mpi_solid_fluid(mpi_solid_fluid_input),
    coupled_deformable_particle_indices(coupled_deformable_particle_indices_input),newmark_evolution(newmark_evolution_input)
{
    modified_world_space_rigid_inertia_tensor_inverse.Resize(modified_world_space_rigid_inertia_tensor.m);
    for(int i=1;i<=modified_world_space_rigid_inertia_tensor.m;i++) modified_world_space_rigid_inertia_tensor_inverse(i)=modified_world_space_rigid_inertia_tensor(i).Inverse();
    modified_world_space_rigid_mass.Resize(rigid_body_fluid_mass.m);modified_world_space_rigid_mass_inverse.Resize(rigid_body_fluid_mass.m);
    for(int i=1;i<=modified_world_space_rigid_mass_inverse.m;i++){modified_world_space_rigid_mass(i)=rigid_body_fluid_mass(i)+solid_system.projection_data.mass.world_space_rigid_mass(i).mass;
        modified_world_space_rigid_mass_inverse(i)=modified_world_space_rigid_mass(i).Inverse();}
    modified_mass.Resize(fluid_mass.m);one_over_modified_mass.Resize(fluid_mass.m);
    for(int i=1;i<=one_over_modified_mass.m;i++){modified_mass(i)=fluid_mass(i)+solid_system.projection_data.mass.mass(i);
        one_over_modified_mass(i)=modified_mass(i).Inverse();}
    recv_fluid_V_boundary_arrays.Resize(coupled_deformable_particle_indices.m);
    recv_fluid_rigid_V_boundary_arrays.Resize(coupled_deformable_particle_indices.m);
    for(int i=1;i<=coupled_deformable_particle_indices.m;i++){
        recv_fluid_V_boundary_arrays(i).Resize(coupled_deformable_particle_indices(i).m);
        recv_fluid_rigid_V_boundary_arrays(i).Resize(rigid_V_size);}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLID_SYSTEM_MPI<TV>::
~SOLID_SYSTEM_MPI()
{
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void SOLID_SYSTEM_MPI<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& BV) const
{
    VECTOR_T& V=dynamic_cast<VECTOR_T&>(BV);
    solid_system.Set_Global_Boundary_Conditions(V,newmark_evolution.X_save,newmark_evolution.rigid_X_save,newmark_evolution.rigid_rotation_save,newmark_evolution.rigid_velocity_save,
        newmark_evolution.rigid_angular_momentum_save,newmark_evolution.V_save,newmark_evolution.solids_parameters.implicit_solve_parameters.test_system,
        newmark_evolution.solids_parameters.implicit_solve_parameters.print_matrix);
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void SOLID_SYSTEM_MPI<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF) const
{
    const VECTOR_T& V=debug_cast<const VECTOR_T&>(BV);VECTOR_T& F=debug_cast<VECTOR_T&>(BF);
    Send_Generalized_Velocity_To_Fluid(V); // MPI
    solid_system.Force(V,F);
    for(int i=1;i<=V.V.Size();i++) F.V(i)=Solid_Sign()*(modified_mass(i)*V.V(i)-solid_system.dt*F.V(i));
    for(int i=1;i<=V.rigid_V.Size();i++){
        F.rigid_V(i).linear=Solid_Sign()*(modified_world_space_rigid_mass(i)*V.rigid_V(i).linear-solid_system.dt*F.rigid_V(i).linear);
        F.rigid_V(i).angular=Solid_Sign()*(modified_world_space_rigid_inertia_tensor(i)*V.rigid_V(i).angular-solid_system.dt*F.rigid_V(i).angular);}

    Get_Generalized_Velocity_From_Fluid(F);
    for(int i=1;i<=V.V.Size();i++) F.V(i)=one_over_modified_mass(i)*F.V(i);
    for(int i=1;i<=V.rigid_V.Size();i++){
        F.rigid_V(i).linear=modified_world_space_rigid_mass_inverse(i)*F.rigid_V(i).linear;
        F.rigid_V(i).angular=modified_world_space_rigid_inertia_tensor_inverse(i)*F.rigid_V(i).angular;}
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double SOLID_SYSTEM_MPI<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& BV1,const KRYLOV_VECTOR_BASE<T>& BV2) const
{
    const VECTOR_T& V1=debug_cast<const VECTOR_T&>(BV1),&V2=debug_cast<const VECTOR_T&>(BV2);
    double inner_product=ARRAYS_COMPUTATIONS::Inner_Product_Double_Precision(modified_mass,V1.V,V2.V);
    for(int i=1;i<=V1.rigid_V.Size();i++){
        inner_product+=Dot_Product(V1.rigid_V(i).linear,modified_world_space_rigid_mass(i)*V2.rigid_V(i).linear);
        inner_product+=Dot_Product(V1.rigid_V(i).angular,modified_world_space_rigid_inertia_tensor(i)*V2.rigid_V(i).angular);}
    return inner_product;
}
#ifdef USE_MPI
//#####################################################################
// Function void Send_Generalized_Velocity_To_Solid
//#####################################################################
template<class TV> void SOLID_SYSTEM_MPI<TV>::
Send_Generalized_Velocity_To_Fluid(const GENERALIZED_VELOCITY<TV>& BV) const
{
    const VECTOR_T& V=debug_cast<const VECTOR_T&>(BV);
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
template<class TV> void SOLID_SYSTEM_MPI<TV>::
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
template<class TV> void SOLID_SYSTEM_MPI<TV>::Send_Generalized_Velocity_To_Fluid(const GENERALIZED_VELOCITY<TV>& V) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void SOLID_SYSTEM_MPI<TV>::Get_Generalized_Velocity_From_Fluid(GENERALIZED_VELOCITY<TV>& V) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
#endif
#define INSTANTIATION_HELPER(T,d) \
    template class SOLID_SYSTEM_MPI<VECTOR<T,d> >;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif

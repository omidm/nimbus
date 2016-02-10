//#####################################################################
// Copyright 2007-2008, Avi (Snarky) Robinson-Mosher, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_FLUID_SYSTEM
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Vectors/Dot_Product.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_SYSTEM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class T_MATRIX> SOLID_FLUID_SYSTEM<TV,T_MATRIX>::
SOLID_FLUID_SYSTEM(BACKWARD_EULER_SYSTEM<TV>& solid_system_input,const ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& J_deformable_array_input,
    const ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& J_rigid_array_input,const ARRAY<T_DIAGONAL_MATRIX>& fluid_mass_input,
    const ARRAY<T_DIAGONAL_MATRIX>& rigid_body_fluid_mass_input,const ARRAY<T_INERTIA_TENSOR>& modified_world_space_rigid_inertia_tensor_input,
    const T fluid_tolerance_input,const T solid_tolerance_input,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array_input)
    :BASE(false,true),solid_system(solid_system_input),J_deformable_array(J_deformable_array_input),
    J_rigid_array(J_rigid_array_input),fluid_mass(fluid_mass_input),rigid_body_fluid_mass(rigid_body_fluid_mass_input),
    modified_world_space_rigid_inertia_tensor(modified_world_space_rigid_inertia_tensor_input),fluid_tolerance(fluid_tolerance_input),solid_tolerance(solid_tolerance_input),
    A_array(A_array_input)
{
    temp_array.Resize(A_array.m);
    modified_world_space_rigid_inertia_tensor_inverse.Resize(modified_world_space_rigid_inertia_tensor.m);
    for(int i=1;i<=modified_world_space_rigid_inertia_tensor.m;i++) modified_world_space_rigid_inertia_tensor_inverse(i)=modified_world_space_rigid_inertia_tensor(i).Inverse();
    modified_world_space_rigid_mass.Resize(rigid_body_fluid_mass.m);modified_world_space_rigid_mass_inverse.Resize(rigid_body_fluid_mass.m);
    for(int i=1;i<=rigid_body_fluid_mass.m;i++){
        modified_world_space_rigid_mass(i)=rigid_body_fluid_mass(i)+solid_system.projection_data.mass.world_space_rigid_mass(i).mass;
        modified_world_space_rigid_mass_inverse(i)=modified_world_space_rigid_mass(i).Inverse();}
    modified_mass.Resize(fluid_mass.m);one_over_modified_mass.Resize(fluid_mass.m);
    for(int i=1;i<=fluid_mass.m;i++){
        modified_mass(i)=fluid_mass(i)+solid_system.solid_body_collection.deformable_body_collection.particles.mass(solid_system.solid_body_collection.deformable_body_collection.dynamic_particles(i));
        one_over_modified_mass(i)=modified_mass(i).Inverse();}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,class T_MATRIX> SOLID_FLUID_SYSTEM<TV,T_MATRIX>::
~SOLID_FLUID_SYSTEM()
{
}
#if 0
//#####################################################################
// Function Print_Matrix
//#####################################################################
template<class TV,class T_MATRIX> void SOLID_FLUID_SYSTEM<TV,T_MATRIX>::
Print_Matrix(VECTOR_T& V,VECTOR_T& F)
{
    std::stringstream ss;
    // treat vector as: p, deformable v, rigid v (linear, angular)
    for(int i=1;i<=fluid_system.A.n;i++){ // fluid rows
        for(int j=1;j<=fluid_system.A.n;j++) if(fluid_system.A.Element_Present(i,j)) ss<<fluid_system.A(i,j)<<" "; else ss<<"0\t";
        for(int j=1;j<=J.m;j++) if(J.Element_Present(j,i)) ss<<J(j,i)<<"\t"; else ss<<"0\t";
        for(int j=1;j<=J_rigid.m;j++) if(J_rigid.Element_Present(j,i)) ss<<J_rigid(j,i)<<"\t";else ss<<"0\t";
        ss<<std::endl;
    }
    for(int i=1;i<=J.m;i++){ // solid rows
        for(int j=1;j<=J.n;j++) if(J.Element_Present(i,j)) ss<<J(i,j)<<"\t";else ss<<"0\t";
        // need to get solid terms
        int our_dynamic_particle_index=(i-1)/TV::dimension+1;int our_axis=i-(our_dynamic_particle_index-1)*TV::dimension;

        for(int j=1;j<=TV::dimension*V.solid_velocity.V.Size();j++){
            int dynamic_particle_index=(j-1)/TV::dimension+1;int axis=j-(dynamic_particle_index-1)*TV::dimension;
            V.solid_velocity*=0;
            V.solid_velocity.V(dynamic_particle_index)(axis)=(T)1;
            solid_system.Force(V.solid_velocity,F.solid_velocity);
            for(int k=1;k<=V.solid_velocity.V.Size();k++) F.solid_velocity.V(k)=Solid_Sign()*(modified_mass(k)*V.solid_velocity.V(k)-solid_system.dt*F.solid_velocity.V(k));

            ss<<F.solid_velocity.V(our_dynamic_particle_index)(our_axis)<<"\t";}
        ss<<std::endl;
    }
    for(int i=1;i<=J_rigid.m;i++){ // solid rows
        for(int j=1;j<=J_rigid.n;j++) if(J_rigid.Element_Present(i,j)) ss<<J_rigid(i,j)<<"\t";else ss<<"0\t";
        // need to get solid terms
        int our_rigid_particle_index=(i-1)/rows_per_rigid_body+1;int our_component=i-(our_rigid_particle_index-1)*rows_per_rigid_body;

        for(int j=1;j<=rows_per_rigid_body*V.solid_velocity.rigid_V.Size();j++){
            int rigid_particle_index=(j-1)/rows_per_rigid_body+1;int component=j-(rigid_particle_index-1)*rows_per_rigid_body;

            V.solid_velocity*=0;
            if(component<=TV::dimension) V.solid_velocity.rigid_V(rigid_particle_index).linear(component)=(T)1;
            else V.solid_velocity.rigid_V(rigid_particle_index).angular(component-TV::dimension)=(T)1;

            solid_system.Force(V.solid_velocity,F.solid_velocity);
            for(int k=1;k<=F.solid_velocity.rigid_V.Size();k++){
                F.solid_velocity.rigid_V(k).linear=Solid_Sign()*(modified_world_space_rigid_mass(k)*V.solid_velocity.rigid_V(k).linear-solid_system.dt*F.solid_velocity.rigid_V(k).linear);
                F.solid_velocity.rigid_V(k).angular=Solid_Sign()*(modified_world_space_rigid_inertia_tensor(k)*V.solid_velocity.rigid_V(k).angular-solid_system.dt*F.solid_velocity.rigid_V(k).angular);}
            
            if(our_component<=TV::dimension) ss<<F.solid_velocity.rigid_V(our_rigid_particle_index).linear(our_component)<<"\t";
            else ss<<F.solid_velocity.rigid_V(our_rigid_particle_index).angular(our_component-TV::dimension)<<"\t";}
        ss<<std::endl;
    }
    LOG::filecout(ss.str());
}
#endif
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV,class T_MATRIX> void SOLID_FLUID_SYSTEM<TV,T_MATRIX>::
Multiply(const KRYLOV_VECTOR_BASE<T>& bV,KRYLOV_VECTOR_BASE<T>& bF) const
{
    const VECTOR_T& V=debug_cast<const VECTOR_T&>(bV);VECTOR_T& F=debug_cast<VECTOR_T&>(bF);
    solid_system.Force(V.solid_velocity,F.solid_velocity);
    // fluid_mass has components as vectors; need to scale velocities by each mass separately.
    for(int i=1;i<=V.solid_velocity.V.Size();i++) F.solid_velocity.V(i)=Solid_Sign()*(modified_mass(i)*V.solid_velocity.V(i)-solid_system.dt*F.solid_velocity.V(i));
    for(int i=1;i<=V.solid_velocity.rigid_V.Size();i++){
        F.solid_velocity.rigid_V(i).linear=Solid_Sign()*(modified_world_space_rigid_mass(i)*V.solid_velocity.rigid_V(i).linear-solid_system.dt*F.solid_velocity.rigid_V(i).linear);
        F.solid_velocity.rigid_V(i).angular=Solid_Sign()*(modified_world_space_rigid_inertia_tensor(i)*V.solid_velocity.rigid_V(i).angular-solid_system.dt*F.solid_velocity.rigid_V(i).angular);}

    for(int i=1;i<=A_array.m;i++){
        const SPARSE_MATRIX_FLAT_NXN<T>& A=A_array(i);
        const SPARSE_MATRIX_FLAT_MXN<T>& J_deformable=J_deformable_array(i);
        const SPARSE_MATRIX_FLAT_MXN<T>& J_rigid=J_rigid_array(i);
        F.pressure(i).Resize(A.n);A.Times(V.pressure(i),F.pressure(i));

        Add_J_Rigid_Transpose_Times_Velocity(J_rigid,V.solid_velocity,F.pressure(i));
        Add_J_Deformable_Transpose_Times_Velocity(J_deformable,V.solid_velocity,F.pressure(i));

        Add_J_Rigid_Times_Pressure(J_rigid,V.pressure(i),F.solid_velocity);
        Add_J_Deformable_Times_Pressure(J_deformable,V.pressure(i),F.solid_velocity);}

    for(int i=1;i<=V.solid_velocity.V.Size();i++) F.solid_velocity.V(i)=one_over_modified_mass(i)*F.solid_velocity.V(i);
    for(int i=1;i<=V.solid_velocity.rigid_V.Size();i++){
        F.solid_velocity.rigid_V(i).linear=modified_world_space_rigid_mass_inverse(i)*F.solid_velocity.rigid_V(i).linear;
        F.solid_velocity.rigid_V(i).angular=modified_world_space_rigid_inertia_tensor_inverse(i)*F.solid_velocity.rigid_V(i).angular;}
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV,class T_MATRIX> void SOLID_FLUID_SYSTEM<TV,T_MATRIX>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& bV) const
{
    // fluid_system.Set_Boundary_Conditions(V.pressure); // TODO: nullspace stuff for fluids
    //solid_system.Set_Global_Boundary_Conditions(V.solid_velocity);
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV,class T_MATRIX> void SOLID_FLUID_SYSTEM<TV,T_MATRIX>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& bV,KRYLOV_VECTOR_BASE<T>& bR) const  // solve MR=V
{
    const VECTOR_T& V=debug_cast<const VECTOR_T&>(bV);VECTOR_T& R=debug_cast<VECTOR_T&>(bR);
    for(int i=1;i<=A_array.m;i++){
        const SPARSE_MATRIX_FLAT_NXN<T>& A=A_array(i);
        temp_array(i).Resize(A.n);
        A.C->Solve_Forward_Substitution(V.pressure(i),temp_array(i),true); // diagonal should be treated as the identity
        A.C->Solve_Backward_Substitution(temp_array(i),R.pressure(i),false,true);} // diagonal is inverted to save on divides
    R.solid_velocity=V.solid_velocity;R.solid_velocity*=(T)-1;
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV,class T_MATRIX> void SOLID_FLUID_SYSTEM<TV,T_MATRIX>::
Project(KRYLOV_VECTOR_BASE<T>& bV) const
{
    // fluid_system.Project(V.pressure);
    solid_system.Project(debug_cast<VECTOR_T&>(bV).solid_velocity);
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV,class T_MATRIX> double SOLID_FLUID_SYSTEM<TV,T_MATRIX>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& bV1,const KRYLOV_VECTOR_BASE<T>& bV2) const
{
    const VECTOR_T& V1=debug_cast<const VECTOR_T&>(bV1),&V2=debug_cast<const VECTOR_T&>(bV2);
    double fluid_inner_product=0.0;for(int i=1;i<=V1.pressure.m;i++) fluid_inner_product+=Dot_Product_Double_Precision(V1.pressure(i),V2.pressure(i));
    double solid_inner_product=0.0;
    for(int i=1;i<=V1.solid_velocity.V.Size();i++) solid_inner_product+=Dot_Product(V1.solid_velocity.V(i),modified_mass(i)*V2.solid_velocity.V(i));
    for(int i=1;i<=V1.solid_velocity.rigid_V.Size();i++){
        solid_inner_product+=Dot_Product(V1.solid_velocity.rigid_V(i).linear,modified_world_space_rigid_mass(i)*V2.solid_velocity.rigid_V(i).linear);
        solid_inner_product+=Dot_Product(V1.solid_velocity.rigid_V(i).angular,modified_world_space_rigid_inertia_tensor(i)*V2.solid_velocity.rigid_V(i).angular);}
    //LOG::cout<<"Inner product: solid: "<<solid_inner_product<<" fluid: "<<fluid_inner_product<<std::endl;
    return fluid_inner_product+solid_inner_product;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV,class T_MATRIX> typename TV::SCALAR SOLID_FLUID_SYSTEM<TV,T_MATRIX>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& bR) const
{
    const VECTOR_T& R=debug_cast<const VECTOR_T&>(bR);
    T fluid_convergence_norm=(T)0;for(int i=1;i<=R.pressure.m;i++) fluid_convergence_norm=max(fluid_convergence_norm,R.pressure(i).Maximum_Magnitude());
    T solid_convergence_norm=solid_system.Convergence_Norm(R.solid_velocity);T scaled_fluid_convergence_norm=solid_tolerance/fluid_tolerance*fluid_convergence_norm;
    T convergence_norm=max(solid_convergence_norm,scaled_fluid_convergence_norm);
    //LOG::cout<<"residual norm: "<<Inner_Product(R,R)<<std::endl;
    /*LOG::cout<<"Convergence norm: "<<convergence_norm<<" solid: "<<solid_convergence_norm<<" fluid: "<<fluid_convergence_norm<<" weighted fluid: "<<scaled_fluid_convergence_norm<<"
     * fluid tolerance: "<<fluid_tolerance<<std::endl;*/
    return convergence_norm;
}
//#####################################################################
// Function Add_J_Deformable_Transpose_Times_Velocity
//#####################################################################
template<class TV,class T_MATRIX> void SOLID_FLUID_SYSTEM<TV,T_MATRIX>::
Add_J_Deformable_Transpose_Times_Velocity(const SPARSE_MATRIX_FLAT_MXN<T>& J_deformable,const GENERALIZED_VELOCITY<TV>& V,VECTOR_ND<T>& pressure)
{
    assert(pressure.n==J_deformable.n && J_deformable.m==V.V.indices.Size()*TV::dimension);
    // computes pressure+=J*V.V
    int index=J_deformable.offsets(1);
    for(int i=1;i<=J_deformable.m;i++){
        const int end=J_deformable.offsets(i+1);
        for(;index<end;index++){
            int dynamic_particle_index=(i-1)/TV::dimension+1;int axis=i-(dynamic_particle_index-1)*TV::dimension;
            pressure(J_deformable.A(index).j)+=J_deformable.A(index).a*V.V(dynamic_particle_index)(axis);}}
}
//#####################################################################
// Function Add_J_Rigid_Transpose_Times_Velocity
//#####################################################################
template<class TV,class T_MATRIX> void SOLID_FLUID_SYSTEM<TV,T_MATRIX>::
Add_J_Rigid_Transpose_Times_Velocity(const SPARSE_MATRIX_FLAT_MXN<T>& J_rigid,const GENERALIZED_VELOCITY<TV>& V,VECTOR_ND<T>& pressure)
{
    assert(pressure.n==J_rigid.n && J_rigid.m==V.rigid_V.indices.Size()*rows_per_rigid_body);
    // computes pressure+=J*V.V
    int index=J_rigid.offsets(1);
    for(int i=1;i<=J_rigid.m;i++){
        const int end=J_rigid.offsets(i+1);
        for(;index<end;index++){
            int rigid_particle_index=(i-1)/rows_per_rigid_body+1;int component=i-(rigid_particle_index-1)*rows_per_rigid_body;
            if(component<=TV::dimension) pressure(J_rigid.A(index).j)+=J_rigid.A(index).a*V.rigid_V(rigid_particle_index).linear(component);
            else pressure(J_rigid.A(index).j)+=J_rigid.A(index).a*V.rigid_V(rigid_particle_index).angular(component-TV::dimension);}}
}
//#####################################################################
// Function Add_J_Deformable_Times_Pressure
//#####################################################################
template<class TV,class T_MATRIX> void SOLID_FLUID_SYSTEM<TV,T_MATRIX>::
Add_J_Deformable_Times_Pressure(const SPARSE_MATRIX_FLAT_MXN<T>& J_deformable,const VECTOR_ND<T>& pressure,GENERALIZED_VELOCITY<TV>& V)
{
    assert(pressure.n==J_deformable.n && J_deformable.m==V.V.indices.Size()*TV::dimension);
    int index=J_deformable.offsets(1);
    for(int i=1;i<=J_deformable.m;i++){
        const int end=J_deformable.offsets(i+1);
        for(;index<end;index++){
            int dynamic_particle_index=(i-1)/TV::dimension+1;int axis=i-(dynamic_particle_index-1)*TV::dimension;
            V.V(dynamic_particle_index)(axis)+=pressure(J_deformable.A(index).j)*J_deformable.A(index).a;}}
}
//#####################################################################
// Function Add_J_Rigid_Times_Pressure
//#####################################################################
template<class TV,class T_MATRIX> void SOLID_FLUID_SYSTEM<TV,T_MATRIX>::
Add_J_Rigid_Times_Pressure(const SPARSE_MATRIX_FLAT_MXN<T>& J_rigid,const VECTOR_ND<T>& pressure,GENERALIZED_VELOCITY<TV>& V)
{
    assert(pressure.n==J_rigid.n && J_rigid.m==V.rigid_V.indices.Size()*rows_per_rigid_body);
    // computes pressure+=J*V.V
    int index=J_rigid.offsets(1);
    for(int i=1;i<=J_rigid.m;i++){
        const int end=J_rigid.offsets(i+1);
        for(;index<end;index++){
            int rigid_particle_index=(i-1)/rows_per_rigid_body+1;int component=i-(rigid_particle_index-1)*rows_per_rigid_body;
            if(component<=TV::dimension) V.rigid_V(rigid_particle_index).linear(component)+=J_rigid.A(index).a*pressure(J_rigid.A(index).j);
            else V.rigid_V(rigid_particle_index).angular(component-TV::dimension)+=J_rigid.A(index).a*pressure(J_rigid.A(index).j);}}
}
//#####################################################################
template class SOLID_FLUID_SYSTEM<VECTOR<float,1>,SPARSE_MATRIX_FLAT_NXN<float> >;
template class SOLID_FLUID_SYSTEM<VECTOR<float,2>,SPARSE_MATRIX_FLAT_NXN<float> >;
template class SOLID_FLUID_SYSTEM<VECTOR<float,3>,SPARSE_MATRIX_FLAT_NXN<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SOLID_FLUID_SYSTEM<VECTOR<double,1>,SPARSE_MATRIX_FLAT_NXN<double> >;
template class SOLID_FLUID_SYSTEM<VECTOR<double,2>,SPARSE_MATRIX_FLAT_NXN<double> >;
template class SOLID_FLUID_SYSTEM<VECTOR<double,3>,SPARSE_MATRIX_FLAT_NXN<double> >;
#endif

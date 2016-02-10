//#####################################################################
// Copyright 2008-2009, Elliot English, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Vectors/Dot_Product.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SLIP_SYSTEM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SLIP_SYSTEM<TV>::
SLIP_SYSTEM(const bool use_preconditioner_input,BACKWARD_EULER_SYSTEM<TV>* solid_system_input,SPARSE_MATRIX_FLAT_MXN<T>& W_input,SPARSE_MATRIX_FLAT_MXN<T>& N_input,
        SPARSE_MATRIX_FLAT_MXN<T>& div_input,VECTOR_ND<T>& M_inverse_input,SPARSE_MATRIX_FLAT_MXN<T>& J_input,SPARSE_MATRIX_FLAT_MXN<T>& P_input,SPARSE_MATRIX_FLAT_MXN<T> div_precondition,
        GENERALIZED_VELOCITY<TV>& temp_generalized_velocity_input,GENERALIZED_VELOCITY<TV>& second_temp_generalized_velocity_input,const T dt_input,const bool leakproof_solve_input,
        const INTERVAL<int> divergence_indices_input,const bool using_slip)
    :BASE(use_preconditioner_input,true),solid_system(solid_system_input),J(J_input),N(N_input),div(div_input),W(W_input),PP(P_input),M_inverse(M_inverse_input),solid_mass(0),dt(dt_input),
        leakproof_solve(leakproof_solve_input)
{
    N.Transpose(N_transpose);
    div.Transpose(div_transpose);
    W.Transpose(W_transpose);
    pressures_size_vector.Resize(div.m);

    if(using_slip)
        C_f=div+PP*N*W;
    else
        C_f=div+PP*W;

    C_f.Transpose(C_f_transpose);

    if(!leakproof_solve){
        J.Transpose(J_transpose);
        solid_mass=&solid_system->projection_data.mass;
        if(using_slip)
            C_s=PP*N*J;
        else
            C_s=PP*J;
        C_s*=-1;
#ifdef DEBUG_OUTPUT
        std::stringstream ss;
        ss<<"C_s:"<<std::endl;
        for(int i=1;i<=C_s.m;i++){
            ss<<"Row "<<i<<": ";
            for(int j=1;j<=C_s.n;j++){
                if(C_s.Element_Present(i,j))
                    ss<<"("<<j<<", "<<C_s(i,j)<<") ";
            }
            ss<<std::endl;}
        LOG::filecout(ss.str());
#endif

        C_s.Transpose(C_s_transpose);
        solid_velocities_size_vector.Resize(J.n);
    }
        
    if(use_preconditioner){
        SPARSE_MATRIX_FLAT_MXN<T> div_precondition_transpose;
        div_precondition.Transpose(div_precondition_transpose);
        div_M_inverse_div_transpose_precondition=(div_precondition*div_precondition_transpose.Scale_Rows(M_inverse)).Create_NXN_Matrix();
        //LOG::cout<<div_M_inverse_div_transpose_precondition<<std::endl;
        div_M_inverse_div_transpose_precondition.Construct_Incomplete_Cholesky_Factorization();
        preconditioned_pressures_size_vector.Resize(div_precondition.m);
    }
    fluid_velocities_size_vector.Resize(M_inverse.Size());
    lagrange_multipliers_size_vector.Resize(C_f.m);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SLIP_SYSTEM<TV>::
~SLIP_SYSTEM()
{
}
template<class TV> void SLIP_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& bV,KRYLOV_VECTOR_BASE<T>& bF) const
{
    const VECTOR_T& V=debug_cast<const VECTOR_T&>(bV);VECTOR_T& F=debug_cast<VECTOR_T&>(bF);
    C_f_transpose.Times(V.lagrange_multipliers,fluid_velocities_size_vector);
    fluid_velocities_size_vector*=M_inverse;
    C_f.Times(fluid_velocities_size_vector,F.lagrange_multipliers);
        
    if(!leakproof_solve){
        V.solid_velocity.Pack(solid_velocities_size_vector);
        C_s.Times(solid_velocities_size_vector,lagrange_multipliers_size_vector);
        F.lagrange_multipliers+=lagrange_multipliers_size_vector;

        solid_system->Force(V.solid_velocity,F.solid_velocity);

        for(int i=1;i<=V.solid_velocity.V.Size();i++)
            F.solid_velocity.V(i)=F.solid_velocity.V(i)*dt-solid_mass->mass(i)*V.solid_velocity.V(i);
        for(int i=1;i<=V.solid_velocity.rigid_V.Size();i++)
            F.solid_velocity.rigid_V(i)=F.solid_velocity.rigid_V(i)*dt-solid_mass->world_space_rigid_mass(i)*V.solid_velocity.rigid_V(i);

        C_s_transpose.Times(V.lagrange_multipliers,solid_velocities_size_vector);
        F.solid_velocity.Unpack_And_Add(solid_velocities_size_vector);
        
        for(int i=1;i<=V.solid_velocity.V.Size();i++)
            F.solid_velocity.V(i)=solid_mass->one_over_mass(i)*F.solid_velocity.V(i);
        for(int i=1;i<=V.solid_velocity.rigid_V.Size();i++)
            F.solid_velocity.rigid_V(i)=solid_mass->world_space_rigid_mass_inverse(i)*F.solid_velocity.rigid_V(i);
    }
}
template<class TV> void SLIP_SYSTEM<TV>::
Apply_Fluid(const VECTOR_T& V,VECTOR_ND<T>& result_dual_cells_size_vector) const
{
#ifdef DEBUG_OUTPUT
    std::stringstream ss;
    ss<<"Writing out matrix (fluid part) after solve"<<std::endl;
    SPARSE_MATRIX_FLAT_MXN<T> full_matrix_fluid_part=C_f_transpose.Scale_Rows(M_inverse);
    for(int i=1;i<=full_matrix_fluid_part.m;i++){
        ss<<"Row "<<i<<": ";
        for(int j=1;j<=full_matrix_fluid_part.n;j++){
            if(full_matrix_fluid_part.Element_Present(i,j))
                ss<<"("<<j<<", "<<full_matrix_fluid_part(i,j)<<") ";
        }
        ss<<std::endl;}
    LOG::filecout(ss.str());
#endif

    result_dual_cells_size_vector.Resize(div_transpose.m);
    C_f_transpose.Times(V.lagrange_multipliers,fluid_velocities_size_vector);
    fluid_velocities_size_vector *= M_inverse;
    result_dual_cells_size_vector=fluid_velocities_size_vector;
}
template<class TV> void SLIP_SYSTEM<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& bV,KRYLOV_VECTOR_BASE<T>& bR) const  // solve MR=V
{
    const VECTOR_T& V=debug_cast<const VECTOR_T&>(bV);VECTOR_T& R=debug_cast<VECTOR_T&>(bR);
    R.Copy(1,V);
    div_M_inverse_div_transpose_precondition.C->Solve_Forward_Substitution(R.lagrange_multipliers,preconditioned_pressures_size_vector,true);
    div_M_inverse_div_transpose_precondition.C->Solve_Backward_Substitution(preconditioned_pressures_size_vector,R.lagrange_multipliers,false,true);
    R.solid_velocity*=-1;
}
template<class TV> double SLIP_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& bV1,const KRYLOV_VECTOR_BASE<T>& bV2) const
{
    const VECTOR_T& V1=debug_cast<const VECTOR_T&>(bV1);const VECTOR_T& V2=debug_cast<const VECTOR_T&>(bV2);
    double inner_product_lagrange_multipliers=Dot_Product_Double_Precision(V1.lagrange_multipliers,V2.lagrange_multipliers);
    double inner_product_solid_velocities=leakproof_solve?0:solid_system->Inner_Product(V1.solid_velocity,V2.solid_velocity);
    return inner_product_lagrange_multipliers+inner_product_solid_velocities;
}
template<class TV> typename TV::SCALAR SLIP_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& bR) const
{
    const VECTOR_T& R=debug_cast<const VECTOR_T&>(bR);
    T convergence_norm_solid=leakproof_solve?0:solid_system->Convergence_Norm(R.solid_velocity);
    T convergence_norm_fluid=R.lagrange_multipliers.Maximum_Magnitude();
    T convergence_norm=max(convergence_norm_fluid,convergence_norm_solid);
    return convergence_norm;
}
//#####################################################################
template class SLIP_SYSTEM<VECTOR<float,1> >;
template class SLIP_SYSTEM<VECTOR<float,2> >;
template class SLIP_SYSTEM<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SLIP_SYSTEM<VECTOR<double,1> >;
template class SLIP_SYSTEM<VECTOR<double,2> >;
template class SLIP_SYSTEM<VECTOR<double,3> >;
#endif

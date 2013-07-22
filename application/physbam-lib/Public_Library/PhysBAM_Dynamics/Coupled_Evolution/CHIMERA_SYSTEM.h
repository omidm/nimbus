//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CHIMERA_SYSTEM
//#####################################################################
#ifndef __CHIMERA_SYSTEM__
#define __CHIMERA_SYSTEM__
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE_SYSTEM.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/CHIMERA_VECTOR.h>
#include <fstream>
namespace PhysBAM{
//#####################################################################
// Class CHIMERA_SYSTEM
//#####################################################################
template<class TV>
class CHIMERA_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef CHIMERA_VECTOR<TV> VECTOR_T;
    typedef KRYLOV_SYSTEM_BASE<T> BASE;
public:
    using BASE::use_preconditioner;using BASE::preconditioner_commutes_with_projection;

    bool least_squares;

    SPARSE_MATRIX_FLAT_MXN<T> inverse_mass_matrix;
    
    SPARSE_MATRIX_FLAT_MXN<T> divergence_matrix;
    SPARSE_MATRIX_FLAT_MXN<T> kinematic_matrix;
    
    SPARSE_MATRIX_FLAT_MXN<T> divergence_matrix_transpose; //gradient
    SPARSE_MATRIX_FLAT_MXN<T> kinematic_matrix_transpose; //kinematic velocity forcing matrix
    
    void Set_Matrices(SPARSE_MATRIX_FLAT_MXN<T>& inverse_mass_matrix_input,SPARSE_MATRIX_FLAT_MXN<T>& divergence_matrix_input,SPARSE_MATRIX_FLAT_MXN<T>& kinematic_matrix_input)
    {
        inverse_mass_matrix=inverse_mass_matrix_input;
        divergence_matrix=divergence_matrix_input;
        kinematic_matrix=kinematic_matrix_input;
        
        divergence_matrix.Transpose(divergence_matrix_transpose);
        kinematic_matrix.Transpose(kinematic_matrix_transpose);
    }
    
    CHIMERA_SYSTEM():
        BASE(false,false),least_squares(false)
    {
    }
    
    CHIMERA_SYSTEM(SPARSE_MATRIX_FLAT_MXN<T>& inverse_mass_matrix_input,SPARSE_MATRIX_FLAT_MXN<T>& divergence_matrix_input,SPARSE_MATRIX_FLAT_MXN<T>& kinematic_matrix_input):
        BASE(false,false),least_squares(false)
    {
        Set_Matrices(inverse_mass_matrix_input,divergence_matrix_input,kinematic_matrix_input);
        
        //LOG::cout << "preconditioning " << use_preconditioner << " " << preconditioner_commutes_with_projection << std::endl;

        /*LOG::filecout("chimera system\n");
        ARRAY<SPARSE_MATRIX_FLAT_MXN<T> > constraint_matrices;
        constraint_matrices.Append(divergence_matrix);
        constraint_matrices.Append(grid_coupling_matrix);
        SPARSE_MATRIX_FLAT_MXN<T> constraint_matrix=SPARSE_MATRIX_FLAT_MXN<T>::Concatenate_Matrices(constraint_matrices);
        //LOG::cout << "divergence_matrix"  << std::endl << divergence_matrix << std::endl;
        //LOG::cout << "grid_coupling_matrix" << std::endl << grid_coupling_matrix << std::endl;
        //LOG::cout << "constraint_matrix" << std::endl << constraint_matrix << std::endl;

        SPARSE_MATRIX_FLAT_MXN<T> constraint_matrix_transpose;
        constraint_matrix.Transpose(constraint_matrix_transpose);

        SPARSE_MATRIX_FLAT_MXN<T> A=constraint_matrix*inverse_mass_matrix*constraint_matrix_transpose;
        //LOG::cout << "schur complement " << std::endl << A << std::endl;*/

        //std::ofstream f("constraint_matrix.m");
        //f << "D=[" << divergence_matrix << "]\n";        
        //f << "K=[" << kinematic_matrix << "]\n";
        //f << "M=[" << inverse_mass_matrix << "]\n";
        //LOG::cout << "writing constraint matrix " << f.good() << std::endl;
        //f.close();

        /*LOG::cout << "divergence_matrix " << std::endl << divergence_matrix << std::endl;
        LOG::cout << "grid_coupling_matrix " << std::endl << grid_coupling_matrix << std::endl;
        LOG::cout << "kinematic_matrix " << std::endl << kinematic_matrix << std::endl;*/
    }
    ~CHIMERA_SYSTEM()
    {}

    void Multiply(const KRYLOV_VECTOR_BASE<T>& bV,KRYLOV_VECTOR_BASE<T>& bF) const PHYSBAM_OVERRIDE
    {
        Multiply_Single(bV,bF);
        if(least_squares)
            Multiply_Single(bF,bF);
    }

    void Multiply_Single(const KRYLOV_VECTOR_BASE<T>& bV,KRYLOV_VECTOR_BASE<T>& bF) const
    {
        const VECTOR_T& V=debug_cast<const VECTOR_T&>(bV);
        VECTOR_T& F=debug_cast<VECTOR_T&>(bF);
        
        VECTOR_ND<T> tmp_v1(divergence_matrix.n);
        VECTOR_ND<T> tmp_v2(divergence_matrix.n);
        
        divergence_matrix_transpose.Times(V.pressure,tmp_v1);
        kinematic_matrix_transpose.Times_Add(V.kinematic_force,tmp_v1);
        
        inverse_mass_matrix.Times(tmp_v1,tmp_v2);
        
        divergence_matrix.Times(tmp_v2,F.pressure);
        kinematic_matrix.Times(tmp_v2,F.kinematic_force);
    }

    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& bV) const PHYSBAM_OVERRIDE
    {}
    
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& bV,KRYLOV_VECTOR_BASE<T>& bR) const PHYSBAM_OVERRIDE
    {
        bR=bV;
    }

    void Project(KRYLOV_VECTOR_BASE<T>& bV) const PHYSBAM_OVERRIDE
    {}

    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& bV1,const KRYLOV_VECTOR_BASE<T>& bV2) const PHYSBAM_OVERRIDE
    {
        const VECTOR_T& V1=debug_cast<const VECTOR_T&>(bV1);
        const VECTOR_T& V2=debug_cast<const VECTOR_T&>(bV2);
        
        return VECTOR_ND<T>::Dot_Product_Double_Precision(V1.pressure,V2.pressure)+VECTOR_ND<T>::Dot_Product_Double_Precision(V1.kinematic_force,V2.kinematic_force);
    }
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& bR) const PHYSBAM_OVERRIDE
    {
        const VECTOR_T& R=debug_cast<const VECTOR_T&>(bR);
        T pressure_norm=R.pressure.Maximum_Magnitude();
        T kinematic_norm=R.kinematic_force.Maximum_Magnitude();
        T norm=max(pressure_norm,kinematic_norm);
        //LOG::cout << "convergence_norm " << pressure_norm << " " << kinematic_norm << std::endl;
        return norm;
    }

    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& bV) const PHYSBAM_OVERRIDE {} // TODO
//#####################################################################
};
}
#endif

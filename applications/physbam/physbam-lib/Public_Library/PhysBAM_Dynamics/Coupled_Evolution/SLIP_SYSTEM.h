//#####################################################################
// Copyright 2007, Avi Robinson-Mosher, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SLIP_SYSTEM
//#####################################################################
#ifndef __SLIP_SYSTEM__
#define __SLIP_SYSTEM__
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE_SYSTEM.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR.h>
namespace PhysBAM{
//#####################################################################
// Class SLIP_SYSTEM
//#####################################################################
template<class TV>
class SLIP_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR<TV> VECTOR_T;
    typedef KRYLOV_SYSTEM_BASE<T> BASE;
public:
    using BASE::use_preconditioner;

    BACKWARD_EULER_SYSTEM<TV>* solid_system;
    SPARSE_MATRIX_FLAT_MXN<T>& J;
    SPARSE_MATRIX_FLAT_MXN<T>& N;
    SPARSE_MATRIX_FLAT_MXN<T>& div;
    SPARSE_MATRIX_FLAT_MXN<T>& W;
    SPARSE_MATRIX_FLAT_MXN<T>& PP;
    VECTOR_ND<T>& M_inverse;
    GENERALIZED_MASS<TV>* solid_mass;
    const T dt;
    const bool leakproof_solve;

    SPARSE_MATRIX_FLAT_MXN<T> J_transpose;
    SPARSE_MATRIX_FLAT_MXN<T> N_transpose;
    SPARSE_MATRIX_FLAT_MXN<T> div_transpose;
    SPARSE_MATRIX_FLAT_MXN<T> W_transpose;
    SPARSE_MATRIX_FLAT_MXN<T> P_transpose;

    SPARSE_MATRIX_FLAT_MXN<T> C_s,C_f;
    SPARSE_MATRIX_FLAT_MXN<T> C_s_transpose,C_f_transpose;

    SPARSE_MATRIX_FLAT_NXN<T> div_M_inverse_div_transpose_precondition;

    mutable VECTOR_ND<T> fluid_velocities_size_vector;
    mutable VECTOR_ND<T> solid_velocities_size_vector;
    mutable VECTOR_ND<T> pressures_size_vector;
    mutable VECTOR_ND<T> preconditioned_pressures_size_vector;
    mutable VECTOR_ND<T> lagrange_multipliers_size_vector;
    
    SLIP_SYSTEM(
        const bool use_preconditioner_input,
        BACKWARD_EULER_SYSTEM<TV>* solid_system_input,
        SPARSE_MATRIX_FLAT_MXN<T>& W_input,
        SPARSE_MATRIX_FLAT_MXN<T>& N_input,
        SPARSE_MATRIX_FLAT_MXN<T>& div_input,
        VECTOR_ND<T>& M_inverse_input,
        SPARSE_MATRIX_FLAT_MXN<T>& J_input,
        SPARSE_MATRIX_FLAT_MXN<T>& P_input,
        SPARSE_MATRIX_FLAT_MXN<T> div_precondition,
        GENERALIZED_VELOCITY<TV>& temp_generalized_velocity_input,
        GENERALIZED_VELOCITY<TV>& second_temp_generalized_velocity_input,
        const T dt_input,
        const bool leakproof_solve_input,
        const INTERVAL<int> divergence_indices_input,
        const bool using_slip);
    virtual ~SLIP_SYSTEM();
    void Multiply(const KRYLOV_VECTOR_BASE<T>& bV,KRYLOV_VECTOR_BASE<T>& bF) const PHYSBAM_OVERRIDE;
    void Apply_Fluid(const VECTOR_T& V,VECTOR_ND<T>& result_dual_cells_size_vector) const;

    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& bV) const PHYSBAM_OVERRIDE
    {}

    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& bV,KRYLOV_VECTOR_BASE<T>& bR) const PHYSBAM_OVERRIDE;  // solve MR=V

    void Project(KRYLOV_VECTOR_BASE<T>& bV) const PHYSBAM_OVERRIDE
    {}

    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& bV1,const KRYLOV_VECTOR_BASE<T>& bV2) const PHYSBAM_OVERRIDE;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& bR) const PHYSBAM_OVERRIDE;

    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& bV) const PHYSBAM_OVERRIDE {} // TODO
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2007-2008, Avi (Snarky) Robinson-Mosher, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_FLUID_SYSTEM
//#####################################################################
#ifndef __SOLID_FLUID_SYSTEM__
#define __SOLID_FLUID_SYSTEM__
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_POLICY.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/PRESSURE_VELOCITY_VECTOR.h>
namespace PhysBAM{

template<class T> class SPARSE_MATRIX_FLAT_MXN;
template<class TV> class BACKWARD_EULER_SYSTEM;
template<class TV> class GRID;

//#####################################################################
// Class SOLID_FLUID_SYSTEM
//#####################################################################
template<class TV,class T_MATRIX>
class SOLID_FLUID_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef typename GRID<TV>::VECTOR_INT TV_INT;typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;
    typedef ARRAY<PAIR<int,T> > FACE_WEIGHT_ELEMENTS;typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS::template REBIND<FACE_WEIGHT_ELEMENTS*>::TYPE T_FACE_ARRAYS_FACE_WEIGHT_ELEMENTS;
    typedef typename TV::SPIN T_SPIN;
    typedef typename RIGID_BODY_POLICY<TV>::WORLD_SPACE_INERTIA_TENSOR T_INERTIA_TENSOR;typedef KRYLOV_SYSTEM_BASE<typename TV::SCALAR> BASE;
    typedef typename MATRIX_POLICY<TV>::DIAGONAL_MATRIX T_DIAGONAL_MATRIX;
public:
    typedef PRESSURE_VELOCITY_VECTOR<TV> VECTOR_T;
    static const int rows_per_rigid_body=TV::dimension+T_SPIN::dimension;
    using BASE::use_preconditioner;

    BACKWARD_EULER_SYSTEM<TV>& solid_system;
    const ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& J_deformable_array;
    const ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& J_rigid_array;
    const ARRAY<T_DIAGONAL_MATRIX>& fluid_mass;
    const ARRAY<T_DIAGONAL_MATRIX>& rigid_body_fluid_mass;
    ARRAY<T_DIAGONAL_MATRIX> modified_mass,one_over_modified_mass;
    ARRAY<T_DIAGONAL_MATRIX> modified_world_space_rigid_mass,modified_world_space_rigid_mass_inverse;
    const ARRAY<T_INERTIA_TENSOR>& modified_world_space_rigid_inertia_tensor;
    ARRAY<T_INERTIA_TENSOR> modified_world_space_rigid_inertia_tensor_inverse;
    const T fluid_tolerance,solid_tolerance;
    const ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array;
    mutable ARRAY<VECTOR_ND<T> > temp_array;

    SOLID_FLUID_SYSTEM(BACKWARD_EULER_SYSTEM<TV>& solid_system_input,const ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& J_deformable_array_input,
        const ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& J_rigid_array_input,const ARRAY<T_DIAGONAL_MATRIX>& fluid_mass_input,
        const ARRAY<T_DIAGONAL_MATRIX>& rigid_body_fluid_mass_input,const ARRAY<T_INERTIA_TENSOR>& modified_world_space_rigid_inertia_tensor_input,
        const T fluid_tolerance_input,const T solid_tolerance_input,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array_input);

    virtual ~SOLID_FLUID_SYSTEM();

    double Check_Positive_Definite(VECTOR_T& V,VECTOR_T& F,VECTOR_T& R) const
    {
        Multiply(F,R);return Inner_Product(V,R);
    }

    const T Solid_Sign() const
    {return (T)-1;}

    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE {} // TODO

//#####################################################################
    // void Print_Matrix(VECTOR_T& V,VECTOR_T& F);
    void Multiply(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& F) const PHYSBAM_OVERRIDE;
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const PHYSBAM_OVERRIDE;  // solve MR=V
    void Project(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& V1,const KRYLOV_VECTOR_BASE<T>& V2) const PHYSBAM_OVERRIDE;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& R) const PHYSBAM_OVERRIDE;
    static void Add_J_Deformable_Transpose_Times_Velocity(const SPARSE_MATRIX_FLAT_MXN<T>& J_deformable,const GENERALIZED_VELOCITY<TV>& V,VECTOR_ND<T>& pressure);
    static void Add_J_Rigid_Transpose_Times_Velocity(const SPARSE_MATRIX_FLAT_MXN<T>& J_rigid,const GENERALIZED_VELOCITY<TV>& V,VECTOR_ND<T>& pressure);
    static void Add_J_Deformable_Times_Pressure(const SPARSE_MATRIX_FLAT_MXN<T>& J_deformable,const VECTOR_ND<T>& pressure,GENERALIZED_VELOCITY<TV>& V);
    static void Add_J_Rigid_Times_Pressure(const SPARSE_MATRIX_FLAT_MXN<T>& J_rigid,const VECTOR_ND<T>& pressure,GENERALIZED_VELOCITY<TV>& V);
//#####################################################################
};
}
#endif

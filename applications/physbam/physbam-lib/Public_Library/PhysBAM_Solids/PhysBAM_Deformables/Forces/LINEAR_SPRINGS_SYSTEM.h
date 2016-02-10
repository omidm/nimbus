//#####################################################################
// Copyright 2010, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_SPRINGS_SYSTEM
//#####################################################################
#ifndef __LINEAR_SPRINGS_SYSTEM__
#define __LINEAR_SPRINGS_SYSTEM__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>

namespace PhysBAM{
template<class T> class SPARSE_MATRIX_FLAT_MXN;
//#####################################################################
// Class LINEAR_SPRINGS_SYSTEM_VECTOR
//#####################################################################
template<class TV>
class LINEAR_SPRINGS_SYSTEM_VECTOR:public KRYLOV_VECTOR_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef KRYLOV_VECTOR_BASE<T> BASE;
public:
    VECTOR_ND<T> X;

    LINEAR_SPRINGS_SYSTEM_VECTOR(const int size);
    ~LINEAR_SPRINGS_SYSTEM_VECTOR();

    LINEAR_SPRINGS_SYSTEM_VECTOR& operator=(const LINEAR_SPRINGS_SYSTEM_VECTOR& v);
    BASE& operator+=(const BASE& bv) PHYSBAM_OVERRIDE;
    BASE& operator-=(const BASE& bv) PHYSBAM_OVERRIDE;
    BASE& operator*=(const T a) PHYSBAM_OVERRIDE;
    void Copy(const T c,const BASE& bv) PHYSBAM_OVERRIDE;
    void Copy(const T c1,const BASE& bv1,const BASE& bv2) PHYSBAM_OVERRIDE;
    void Print() const;
    int Raw_Size() const PHYSBAM_OVERRIDE;
    T& Raw_Get(int i) PHYSBAM_OVERRIDE;
};
//#####################################################################
// Class LINEAR_SPRINGS_SYSTEM
//#####################################################################
template<class TV>
class LINEAR_SPRINGS_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
public:
    typedef LINEAR_SPRINGS_SYSTEM_VECTOR<TV> VECTOR_T;

    SPARSE_MATRIX_FLAT_MXN<T>& J;

    LINEAR_SPRINGS_SYSTEM(SPARSE_MATRIX_FLAT_MXN<T>& J_input);

    virtual ~LINEAR_SPRINGS_SYSTEM();

//#####################################################################
    void Force(const VECTOR_T& V,VECTOR_T& F) const;
    void Multiply(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& F) const PHYSBAM_OVERRIDE;
    void Project(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& V1,const KRYLOV_VECTOR_BASE<T>& V2) const PHYSBAM_OVERRIDE;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& R) const PHYSBAM_OVERRIDE;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const PHYSBAM_OVERRIDE;
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

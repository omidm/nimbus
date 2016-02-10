//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COUPLED_SYSTEM_VECTOR
//#####################################################################
#ifndef __COUPLED_SYSTEM_VECTOR__
#define __COUPLED_SYSTEM_VECTOR__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COUPLING_CONSTRAINT_ID.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/FORCE_AGGREGATE_ID.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/VISCOUS_FORCE_ID.h>
namespace PhysBAM{

template<class TV>
class COUPLED_SYSTEM_VECTOR:public KRYLOV_VECTOR_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef KRYLOV_VECTOR_BASE<T> BASE;
public:
    VECTOR_ND<T> pressure;
    ARRAY<T,COUPLING_CONSTRAINT_ID> lambda;
    ARRAY<T,FORCE_AGGREGATE_ID> force_coefficients;
    ARRAY<T,VISCOUS_FORCE_ID> viscous_force_coefficients;

    COUPLED_SYSTEM_VECTOR();
    ~COUPLED_SYSTEM_VECTOR();

    COUPLED_SYSTEM_VECTOR& operator=(const COUPLED_SYSTEM_VECTOR& v);
    BASE& operator+=(const BASE& bv) PHYSBAM_OVERRIDE;
    BASE& operator-=(const BASE& bv) PHYSBAM_OVERRIDE;
    BASE& operator*=(const T a) PHYSBAM_OVERRIDE;
    void Copy(const T c,const BASE& bv) PHYSBAM_OVERRIDE;
    void Copy(const T c1,const BASE& bv1,const BASE& bv2) PHYSBAM_OVERRIDE;
    void Print() const;
    int Raw_Size() const PHYSBAM_OVERRIDE;
    T& Raw_Get(int i) PHYSBAM_OVERRIDE;
};
}
#endif

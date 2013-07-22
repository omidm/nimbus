//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR
//#####################################################################
#ifndef __LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR__
#define __LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
namespace PhysBAM{
template<class T> class VECTOR_ND;
template<class TV> class GENERALIZED_VELOCITY;

template<class TV>
class LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR : public KRYLOV_VECTOR_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef KRYLOV_VECTOR_BASE<T> BASE;
public:
    VECTOR_ND<T>& lagrange_multipliers;
    GENERALIZED_VELOCITY<TV>& solid_velocity;

    LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR(VECTOR_ND<T>& lagrange_multipliers_input,GENERALIZED_VELOCITY<TV>& solid_velocity_input);

    LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR& operator=(const LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR& v);
    BASE& operator+=(const BASE& bv) PHYSBAM_OVERRIDE;
    BASE& operator-=(const BASE& bv) PHYSBAM_OVERRIDE;
    BASE& operator*=(const T a) PHYSBAM_OVERRIDE;
    void Copy(const T c,const BASE& bv) PHYSBAM_OVERRIDE;
    void Copy(const T c1,const BASE& bv1,const BASE& bv2) PHYSBAM_OVERRIDE;
};
}
#endif

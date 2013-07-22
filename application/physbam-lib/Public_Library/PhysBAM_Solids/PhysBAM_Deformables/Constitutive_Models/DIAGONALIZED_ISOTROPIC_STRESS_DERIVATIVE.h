//#####################################################################
// Copyright 2004-2007, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE
//#####################################################################
#ifndef __DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE__
#define __DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE__

#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T,int d> class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE;

template<class T>
class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>
{
public:
    T x1111,x2211,x2222; // 2x2 block
    T x2121,x2112; // 2x2 block

    MATRIX<T,2> Differential(const MATRIX<T,2>& dF) const
    {return MATRIX<T,2>(x1111*dF.x[0]+x2211*dF.x[3],x2121*dF.x[1]+x2112*dF.x[2],
        x2112*dF.x[1]+x2121*dF.x[2],x2211*dF.x[0]+x2222*dF.x[3]);}

//#####################################################################
    void Enforce_Definiteness();
//#####################################################################
};

template<class T>
class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>
{
public:
    T x1111,x2211,x3311,x2222,x3322,x3333; // 3x3 block
    T x2121,x2112; // 2x2 block
    T x3131,x3113; // 2x2 block
    T x3223,x3232; // 2x2 block

    MATRIX<T,3> Differential(const MATRIX<T,3>& dF) const
    {return MATRIX<T,3>(x1111*dF.x[0]+x2211*dF.x[4]+x3311*dF.x[8],x2121*dF.x[1]+x2112*dF.x[3],x3131*dF.x[2]+x3113*dF.x[6],
        x2112*dF.x[1]+x2121*dF.x[3],x2211*dF.x[0]+x2222*dF.x[4]+x3322*dF.x[8],x3232*dF.x[5]+x3223*dF.x[7],
        x3113*dF.x[2]+x3131*dF.x[6],x3223*dF.x[5]+x3232*dF.x[7],x3311*dF.x[0]+x3322*dF.x[4]+x3333*dF.x[8]);}

//#####################################################################
    void Enforce_Definiteness(const T eigenvalue_clamp_percentage=(T)0,const T epsilon=(T)1e-4);
//#####################################################################
};
}
#endif

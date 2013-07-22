//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANISOTROPIC_CONSTITUTIVE_MODEL
//##################################################################### 
#ifndef __ANISOTROPIC_CONSTITUTIVE_MODEL__
#define __ANISOTROPIC_CONSTITUTIVE_MODEL__

#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/CONSTITUTIVE_MODEL.h>
namespace PhysBAM{

template<class T,int d>
class ANISOTROPIC_CONSTITUTIVE_MODEL:public CONSTITUTIVE_MODEL<T,d>
{
public:
    bool use_isotropic_component_of_stress_derivative_only;

    ANISOTROPIC_CONSTITUTIVE_MODEL();
    virtual ~ANISOTROPIC_CONSTITUTIVE_MODEL();

//#####################################################################
    virtual MATRIX<T,d> dP_From_dF(const MATRIX<T,d>& dF,const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& V,const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dPi_dF,
        const T scale,const int simplex) const;
    virtual MATRIX<T,d> dP_From_dF(const MATRIX<T,d>& dF,const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& V,const DIAGONALIZED_STRESS_DERIVATIVE<T,d>& dP_dF,
        const T scale,const int tetrahedron) const;
    virtual MATRIX<T,d> P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& V,const T scale,const int simplex) const=0;
    virtual void Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& V,DIAGONALIZED_STRESS_DERIVATIVE<T,d>& dP_dF,const int simplex) const=0;
    virtual void Update_State_Dependent_Auxiliary_Variables(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& V,const int simplex);
//#####################################################################
};
}
#endif

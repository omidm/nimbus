//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ANISOTROPIC_CONSTITUTIVE_MODEL.h>
using namespace PhysBAM;
template<class T,int d> ANISOTROPIC_CONSTITUTIVE_MODEL<T,d>::
ANISOTROPIC_CONSTITUTIVE_MODEL()
    :use_isotropic_component_of_stress_derivative_only(false)
{
}
template<class T,int d> ANISOTROPIC_CONSTITUTIVE_MODEL<T,d>::
~ANISOTROPIC_CONSTITUTIVE_MODEL()
{
}
template<class T,int d> MATRIX<T,d> ANISOTROPIC_CONSTITUTIVE_MODEL<T,d>::
dP_From_dF(const MATRIX<T,d>& dF,const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& V,const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dPi_dF,const T scale,const int simplex) const
{
    return scale*dPi_dF.Differential(dF);
}
template<class T,int d> MATRIX<T,d> ANISOTROPIC_CONSTITUTIVE_MODEL<T,d>::
dP_From_dF(const MATRIX<T,d>& dF,const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& V,const DIAGONALIZED_STRESS_DERIVATIVE<T,d>& dP_dF,const T scale,const int tetrahedron) const
{
    return scale*dP_dF.Differential(dF);
}
template<class T,int d> void ANISOTROPIC_CONSTITUTIVE_MODEL<T,d>::
Update_State_Dependent_Auxiliary_Variables(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& V,const int simplex)
{
}

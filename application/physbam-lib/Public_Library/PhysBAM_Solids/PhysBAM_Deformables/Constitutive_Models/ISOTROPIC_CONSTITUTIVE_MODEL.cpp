//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
using namespace PhysBAM;
template<class T,int d> ISOTROPIC_CONSTITUTIVE_MODEL<T,d>::
ISOTROPIC_CONSTITUTIVE_MODEL()
{
}
template<class T,int d> ISOTROPIC_CONSTITUTIVE_MODEL<T,d>::
~ISOTROPIC_CONSTITUTIVE_MODEL()
{
}
template<class T,int d> MATRIX<T,d> ISOTROPIC_CONSTITUTIVE_MODEL<T,d>::
dP_From_dF(const MATRIX<T,d>& dF,const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dPi_dF,const T scale,const int simplex) const
{
    return scale*dPi_dF.Differential(dF);
}
template<class T,int d> void ISOTROPIC_CONSTITUTIVE_MODEL<T,d>::
Update_State_Dependent_Auxiliary_Variables(const DIAGONAL_MATRIX<T,d>& F,const int simplex)
{
}
template<class T,int d> T ISOTROPIC_CONSTITUTIVE_MODEL<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int simplex) const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
template class ISOTROPIC_CONSTITUTIVE_MODEL<float,2>;
template class ISOTROPIC_CONSTITUTIVE_MODEL<float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ISOTROPIC_CONSTITUTIVE_MODEL<double,2>;
template class ISOTROPIC_CONSTITUTIVE_MODEL<double,3>;
#endif

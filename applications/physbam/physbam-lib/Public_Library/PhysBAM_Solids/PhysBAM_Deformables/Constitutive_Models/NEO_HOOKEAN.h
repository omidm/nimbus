//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NEO_HOOKEAN
//#####################################################################
#ifndef __NEO_HOOKEAN__
#define __NEO_HOOKEAN__

#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
namespace PhysBAM{

using ::std::log;

template<class T,int d> class DIAGONAL_MATRIX;
template<class T,int d>
class NEO_HOOKEAN:public ISOTROPIC_CONSTITUTIVE_MODEL<T,d>
{
    typedef VECTOR<T,d> TV;
public:
    typedef ISOTROPIC_CONSTITUTIVE_MODEL<T,d> BASE;
    using BASE::enforce_definiteness;using BASE::constant_lambda;using BASE::constant_mu;using BASE::constant_alpha;using BASE::constant_beta;

    T youngs_modulus,poissons_ratio;
    T failure_threshold;
private:
    T dth_root_failure_threshold;
public:

    NEO_HOOKEAN(const T youngs_modulus_input=3e6,const T poissons_ratio_input=.475,const T Rayleigh_coefficient=.05,const T failure_threshold_input=.25);
    virtual ~NEO_HOOKEAN();

private:
    DIAGONAL_MATRIX<T,2> Clamp_To_Hyperbola(const DIAGONAL_MATRIX<T,2>& F) const;
    DIAGONAL_MATRIX<T,3> Clamp_To_Hyperbola(const DIAGONAL_MATRIX<T,3>& F) const;
public:

    // clamp to hyperbola to avoid indefiniteness "automatically"
    DIAGONAL_MATRIX<T,d> P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const T scale,const int simplex) const PHYSBAM_OVERRIDE;
    T Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int simplex) const PHYSBAM_OVERRIDE;

/*
    // alternate version using standard indefiniteness fix
    DIAGONAL_MATRIX<T,d> P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const T scale,const int simplex) const PHYSBAM_OVERRIDE
    {T scale_mu=scale*constant_mu,scale_lambda=scale*constant_lambda;
    DIAGONAL_MATRIX<T,d> F_threshold=F.Clamp_Min(failure_threshold),F_inverse=F_threshold.Inverse();
    T J=F_threshold.Determinant();
    T scale_mu_minus_lambda_log_J=scale_mu-scale_lambda*log(J);
    DIAGONAL_MATRIX<T,d> P_threshold=scale_mu*F_threshold-scale_mu_minus_lambda_log_J*F_inverse;
    if(F_threshold==F) return P_threshold;
    // otherwise, compute dP_dF, apply indefiniteness fix, and extrapolate
    if(J>0) scale_mu_minus_lambda_log_J=scale_mu;
    SYMMETRIC_MATRIX<T,d> dP_dF=scale_mu+scale_mu_minus_lambda_log_J*sqr(F_inverse)+scale_lambda*SYMMETRIC_MATRIX<T,d>::Outer_Product(F_inverse.To_Vector());
    //return P_threshold+DIAGONAL_MATRIX<T,d>(dP_dF.Positive_Definite_Part()*(F-F_threshold).To_Vector());}
*/

    MATRIX<T,d> P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const PHYSBAM_OVERRIDE;
    void Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,const int triangle) const;
    void Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF,const int tetrahedron) const;
    int P_From_Strain_Rate_Forces_Size() const PHYSBAM_OVERRIDE;
    void P_From_Strain_Rate_First_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<T> aggregate,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const PHYSBAM_OVERRIDE;
    MATRIX<T,d> P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,const ARRAY_VIEW<const T> aggregate,const T scale,const int simplex) const PHYSBAM_OVERRIDE;

//#####################################################################
};
}
#endif

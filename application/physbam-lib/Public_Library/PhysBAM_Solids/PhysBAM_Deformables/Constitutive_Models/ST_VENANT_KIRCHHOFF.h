//#####################################################################
// Copyright 2003-2005, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ST_VENANT_KIRCHHOFF
//##################################################################### 
#ifndef __ST_VENANT_KIRCHHOFF__
#define __ST_VENANT_KIRCHHOFF__

#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
namespace PhysBAM{

template<class T,int d>
class ST_VENANT_KIRCHHOFF:public ISOTROPIC_CONSTITUTIVE_MODEL<T,d>
{
    typedef VECTOR<T,d> TV;
public:
    typedef ISOTROPIC_CONSTITUTIVE_MODEL<T,d> BASE;
    using BASE::enforce_definiteness;using BASE::constant_lambda;using BASE::constant_mu;using BASE::constant_alpha;using BASE::constant_beta;

    T youngs_modulus,poissons_ratio;
    T failure_threshold;
    
    ST_VENANT_KIRCHHOFF(const T youngs_modulus_input=3e6,const T poissons_ratio_input=.475,const T Rayleigh_coefficient=.05)
        :youngs_modulus(youngs_modulus_input),poissons_ratio(poissons_ratio_input)
    {
        assert(-1<poissons_ratio && poissons_ratio<.5);
        constant_lambda=youngs_modulus*poissons_ratio/((1+poissons_ratio)*(1-2*poissons_ratio));
        constant_mu=youngs_modulus/(2*(1+poissons_ratio));
        constant_alpha=Rayleigh_coefficient*constant_lambda;
        constant_beta=Rayleigh_coefficient*constant_mu;
        if(d==2) failure_threshold=sqrt((constant_mu+constant_lambda)/(3*constant_mu+3*constant_lambda));
        else failure_threshold=sqrt((constant_mu+(T)1.5*constant_lambda)/(3*constant_mu+(T)4.5*constant_lambda));
    }

    DIAGONAL_MATRIX<T,d> P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const T scale,const int simplex) const PHYSBAM_OVERRIDE
    {DIAGONAL_MATRIX<T,d> F_threshold=F.Max(failure_threshold),twice_strain=F_threshold*F_threshold-1;
    return F_threshold*(scale*constant_mu*twice_strain+(T).5*scale*constant_lambda*twice_strain.Trace());}
    
    MATRIX<T,d> P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const PHYSBAM_OVERRIDE
    {DIAGONAL_MATRIX<T,d> F_threshold=F.Max(failure_threshold);
    SYMMETRIC_MATRIX<T,d> strain_rate=(F_threshold*F_dot).Symmetric_Part();
    return F_threshold*(2*scale*constant_beta*strain_rate+scale*constant_alpha*strain_rate.Trace());}

    void Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,const int triangle) const PHYSBAM_OVERRIDE
    {DIAGONAL_MATRIX<T,2> F_threshold=F.Max(failure_threshold);
    T lambda_tr_G_minus_mu=(T).5*constant_lambda*(F_threshold*F_threshold-1).Trace()-constant_mu,three_mu_plus_lambda=3*constant_mu+constant_lambda;
    SYMMETRIC_MATRIX<T,2> F_outer=SYMMETRIC_MATRIX<T,2>::Outer_Product(VECTOR<T,2>(F_threshold.x11,F_threshold.x22));
    dP_dF.x1111=lambda_tr_G_minus_mu+three_mu_plus_lambda*F_outer.x11;//alpha+beta+gamma
    dP_dF.x2222=lambda_tr_G_minus_mu+three_mu_plus_lambda*F_outer.x22;
    dP_dF.x2211=constant_lambda*F_outer.x21;//gamma
    dP_dF.x2121=lambda_tr_G_minus_mu+constant_mu*(F_outer.x22+F_outer.x11);//alpha
    dP_dF.x2112=constant_mu*F_outer.x21;//beta
    if(enforce_definiteness) dP_dF.Fix_Indefinite_Blocks();}

    void Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF,const int tetrahedron) const PHYSBAM_OVERRIDE
    {DIAGONAL_MATRIX<T,3> F_threshold=F.Max(failure_threshold);
    T lambda_tr_G_minus_mu=(T).5*constant_lambda*(F_threshold*F_threshold-1).Trace()-constant_mu,three_mu_plus_lambda=3*constant_mu+constant_lambda;
    SYMMETRIC_MATRIX<T,3> F_outer=SYMMETRIC_MATRIX<T,3>::Outer_Product(VECTOR<T,3>(F_threshold.x11,F_threshold.x22,F_threshold.x33));
    //alpha+beta+gamma
    dPi_dF.x1111=lambda_tr_G_minus_mu+three_mu_plus_lambda*F_outer.x11;
    dPi_dF.x2222=lambda_tr_G_minus_mu+three_mu_plus_lambda*F_outer.x22;
    dPi_dF.x3333=lambda_tr_G_minus_mu+three_mu_plus_lambda*F_outer.x33;
    //gamma
    dPi_dF.x2211=constant_lambda*F_outer.x21;
    dPi_dF.x3311=constant_lambda*F_outer.x31;
    dPi_dF.x3322=constant_lambda*F_outer.x32;
    //alpha
    dPi_dF.x2121=lambda_tr_G_minus_mu+constant_mu*(F_outer.x22+F_outer.x11);
    dPi_dF.x3131=lambda_tr_G_minus_mu+constant_mu*(F_outer.x33+F_outer.x11);
    dPi_dF.x3232=lambda_tr_G_minus_mu+constant_mu*(F_outer.x33+F_outer.x22);
    //beta
    dPi_dF.x2112=constant_mu*F_outer.x21;
    dPi_dF.x3113=constant_mu*F_outer.x31;
    dPi_dF.x3223=constant_mu*F_outer.x32;
    if(enforce_definiteness) dPi_dF.Enforce_Definiteness();}

//#####################################################################
};
}
#endif

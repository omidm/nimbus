//#####################################################################
// Copyright 2003-2006, Geoffrey Irving, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ROTATED_LINEAR
//#####################################################################
#ifndef __ROTATED_LINEAR__
#define __ROTATED_LINEAR__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
namespace PhysBAM{

template<class T,int d>
class ROTATED_LINEAR:public ISOTROPIC_CONSTITUTIVE_MODEL<T,d>
{
    typedef VECTOR<T,d> TV;
public:
    typedef ISOTROPIC_CONSTITUTIVE_MODEL<T,d> BASE;
    using BASE::lambda;using BASE::mu;using BASE::alpha;using BASE::beta;using BASE::constant_lambda;using BASE::constant_mu;using BASE::constant_alpha;using BASE::constant_beta;

    T constant_youngs_modulus,constant_poissons_ratio;
    ARRAY<T> youngs_modulus,poissons_ratio;

    template<class T_FIELD>
    ROTATED_LINEAR(const T_FIELD& youngs_modulus_input,const T_FIELD& poissons_ratio_input=(T).475,const T_FIELD& Rayleigh_coefficient=(T).05)
    {
        Set_Parameters(youngs_modulus_input,poissons_ratio_input,Rayleigh_coefficient);
    }

protected:
    ROTATED_LINEAR()
    {}
public:

    void Set_Parameters(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient)
    {youngs_modulus.Clean_Memory();poissons_ratio.Clean_Memory();lambda.Clean_Memory();mu.Clean_Memory();alpha.Clean_Memory();beta.Clean_Memory();
    constant_youngs_modulus=youngs_modulus_input;
    constant_poissons_ratio=poissons_ratio_input;
    assert(constant_poissons_ratio>-1 && constant_poissons_ratio<(T).5);
    constant_lambda=constant_youngs_modulus*constant_poissons_ratio/((1+constant_poissons_ratio)*(1-2*constant_poissons_ratio));
    constant_mu=constant_youngs_modulus/(2*(1+constant_poissons_ratio));
    constant_alpha=Rayleigh_coefficient*constant_lambda;
    constant_beta=Rayleigh_coefficient*constant_mu;}

    void Set_Parameters(ARRAY_VIEW<const T> youngs_modulus_input,ARRAY_VIEW<const T> poissons_ratio_input,ARRAY_VIEW<const T> Rayleigh_coefficient)
    {constant_youngs_modulus=constant_poissons_ratio=constant_lambda=constant_mu=constant_alpha=constant_beta=0;
    youngs_modulus=youngs_modulus_input;poissons_ratio=poissons_ratio_input;
    assert(youngs_modulus.m==poissons_ratio.m);
    lambda.Resize(youngs_modulus.m);mu.Resize(youngs_modulus.m);alpha.Resize(youngs_modulus.m);beta.Resize(youngs_modulus.m);
    for(int e=1;e<=youngs_modulus.m;e++){
        assert(poissons_ratio(e)>-1 && poissons_ratio(e)<(T).5);
        lambda(e)=youngs_modulus(e)*poissons_ratio(e)/((1+poissons_ratio(e))*(1-2*poissons_ratio(e)));
        mu(e)=youngs_modulus(e)/(2*(1+poissons_ratio(e)));
        alpha(e)=Rayleigh_coefficient(e)*lambda(e);
        beta(e)=Rayleigh_coefficient(e)*mu(e);}}

    DIAGONAL_MATRIX<T,d> P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const T scale,const int simplex) const PHYSBAM_OVERRIDE
    {DIAGONAL_MATRIX<T,d> strain=F-1;
    if(!mu.m) return 2*scale*constant_mu*strain+scale*constant_lambda*strain.Trace();
    else return 2*scale*mu(simplex)*strain+scale*lambda(simplex)*strain.Trace();}
    
    MATRIX<T,d> P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const PHYSBAM_OVERRIDE
    {SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part();
    if(!beta.m) return 2*scale*constant_beta*strain_rate+scale*constant_alpha*strain_rate.Trace();
    else return 2*scale*beta(simplex)*strain_rate+scale*alpha(simplex)*strain_rate.Trace();}

//#####################################################################
};
}
#endif

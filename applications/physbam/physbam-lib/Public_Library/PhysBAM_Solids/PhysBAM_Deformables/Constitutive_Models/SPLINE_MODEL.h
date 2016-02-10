//#####################################################################
// Copyright 2003-2008, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPLINE_MODEL
//#####################################################################
#ifndef __SPLINE_MODEL__
#define __SPLINE_MODEL__

#include <PhysBAM_Tools/Arrays/ARRAY_RATIO.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ROTATED_LINEAR.h>
namespace PhysBAM{

template<class T,int d>
class SPLINE_MODEL:public ROTATED_LINEAR<T,d>
{
    typedef VECTOR<T,d> TV;
public:
    typedef ROTATED_LINEAR<T,d> BASE;
    using BASE::constant_lambda;using BASE::constant_mu;using BASE::lambda;using BASE::mu;

private:
    ARRAY<T> hardening_deformation,hardening_strength;
    ARRAY<T> coefficient,base;
public:

    template<class T_FIELD1,class T_FIELD2>
    SPLINE_MODEL(const T_FIELD1& youngs_modulus_input=3e6,const T_FIELD1& poissons_ratio_input=.475,const T_FIELD2& hardening_deformation_input=.5,
        const T_FIELD2& hardening_strength_input=7,const T_FIELD1& Rayleigh_coefficient=.05)
    {
        Set_Parameters(youngs_modulus_input,poissons_ratio_input,hardening_deformation_input,hardening_strength_input,Rayleigh_coefficient);
    }

    template<class T_FIELD1,class T_FIELD2>
    void Set_Parameters(const T_FIELD1& youngs_modulus_input,const T_FIELD1& poissons_ratio_input,const T_FIELD2& hardening_deformation_input,const T_FIELD2& hardening_strength_input,
        const T_FIELD1& Rayleigh_coefficient)
    {ROTATED_LINEAR<T,d>::Set_Parameters(youngs_modulus_input,poissons_ratio_input,Rayleigh_coefficient);
    Set(hardening_deformation,hardening_deformation_input);
    Set(hardening_strength,hardening_strength_input);
    coefficient=(hardening_strength-1)/(3*(hardening_deformation*hardening_deformation));
    base=hardening_deformation*(hardening_strength-1)*(T)two_thirds;}

private:
    static void Set(ARRAY<T>& parameter,const T value)
    {parameter.Resize(1);parameter(1)=value;}

    static void Set(ARRAY<T>& parameter,ARRAY_VIEW<const T> value)
    {parameter=value;}
public:

    DIAGONAL_MATRIX<T,2> P_From_Strain(const DIAGONAL_MATRIX<T,2>& F,const T scale,const int simplex) const
    {DIAGONAL_MATRIX<T,2> strain=F-1,strain_abs=strain.Abs(),strain_sign=strain.Sign();
    DIAGONAL_MATRIX<T,2> D;
    int index=hardening_deformation.m>1?simplex:1;
    T hardening_deformation_=hardening_deformation(index),hardening_strength_=hardening_strength(index),coefficient_=coefficient(index),base_=base(index);
    D.x11=strain_abs.x11<hardening_deformation_?strain_abs.x11*(1+coefficient_*sqr(strain_abs.x11)):hardening_strength_*strain_abs.x11-base_;
    D.x22=strain_abs.x22<hardening_deformation_?strain_abs.x22*(1+coefficient_*sqr(strain_abs.x22)):hardening_strength_*strain_abs.x22-base_;
    if(!mu.m) return 2*scale*constant_mu*strain_sign*D+scale*constant_lambda*strain.Trace();
    else return 2*scale*mu(simplex)*strain_sign*D+scale*lambda(simplex)*strain.Trace();}

    DIAGONAL_MATRIX<T,3> P_From_Strain(const DIAGONAL_MATRIX<T,3>& F,const T scale,const int simplex) const
    {DIAGONAL_MATRIX<T,3> strain=F-1,strain_abs=strain.Abs(),strain_sign=strain.Sign();
    DIAGONAL_MATRIX<T,3> D;
    int index=hardening_deformation.m>1?simplex:1;
    T hardening_deformation_=hardening_deformation(index),hardening_strength_=hardening_strength(index),coefficient_=coefficient(index),base_=base(index);
    D.x11=strain_abs.x11<hardening_deformation_?strain_abs.x11*(1+coefficient_*sqr(strain_abs.x11)):hardening_strength_*strain_abs.x11-base_;
    D.x22=strain_abs.x22<hardening_deformation_?strain_abs.x22*(1+coefficient_*sqr(strain_abs.x22)):hardening_strength_*strain_abs.x22-base_;
    D.x33=strain_abs.x33<hardening_deformation_?strain_abs.x33*(1+coefficient_*sqr(strain_abs.x33)):hardening_strength_*strain_abs.x33-base_;
    if(!mu.m) return 2*scale*constant_mu*strain_sign*D+scale*constant_lambda*strain.Trace();
    else return 2*scale*mu(simplex)*strain_sign*D+scale*lambda(simplex)*strain.Trace();}

//#####################################################################
};
}
#endif

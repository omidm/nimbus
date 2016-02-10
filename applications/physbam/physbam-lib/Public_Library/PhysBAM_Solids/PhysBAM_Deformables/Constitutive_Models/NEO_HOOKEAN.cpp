//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Math_Tools/pow.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> NEO_HOOKEAN<T,d>::
NEO_HOOKEAN(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient,const T failure_threshold_input)
    :youngs_modulus(youngs_modulus_input),poissons_ratio(poissons_ratio_input),failure_threshold(failure_threshold_input)
{
    assert(poissons_ratio>-1&&poissons_ratio<.5);
    constant_lambda=youngs_modulus*poissons_ratio/((1+poissons_ratio)*(1-2*poissons_ratio));
    constant_mu=youngs_modulus/(2*(1+poissons_ratio));
    constant_alpha=Rayleigh_coefficient*constant_lambda;
    constant_beta=Rayleigh_coefficient*constant_mu;
    dth_root_failure_threshold=pow<1,d>(failure_threshold);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,int d> NEO_HOOKEAN<T,d>::
~NEO_HOOKEAN()
{
}
//#####################################################################
// Function Clamp_To_Hyperbola
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,2> NEO_HOOKEAN<T,d>::
Clamp_To_Hyperbola(const DIAGONAL_MATRIX<T,2>& F) const
{
    if(sqr(F.x11)>failure_threshold) return DIAGONAL_MATRIX<T,2>(F.x11,failure_threshold/F.x11);
    else return DIAGONAL_MATRIX<T,2>(dth_root_failure_threshold,dth_root_failure_threshold);
}
//#####################################################################
// Function Clamp_To_Hyperbola
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,3> NEO_HOOKEAN<T,d>::
Clamp_To_Hyperbola(const DIAGONAL_MATRIX<T,3>& F) const
{
    if(F.x11*F.x22*F.x22>failure_threshold) return DIAGONAL_MATRIX<T,3>(F.x11,F.x22,failure_threshold/(F.x11*F.x22));
    else if(cube(F.x11)>failure_threshold){
        T clamped=sqrt(failure_threshold/F.x11);
        return DIAGONAL_MATRIX<T,3>(F.x11,clamped,clamped);}
    else return DIAGONAL_MATRIX<T,3>(dth_root_failure_threshold,dth_root_failure_threshold,dth_root_failure_threshold);
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
// clamp to hyperbola to avoid indefiniteness "automatically"
template<class T,int d> DIAGONAL_MATRIX<T,d> NEO_HOOKEAN<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const T scale,const int simplex) const
{
    T scale_mu=scale*constant_mu,scale_lambda=scale*constant_lambda,J=F.Determinant();
    if(J>=failure_threshold) return scale_mu*F-(scale_mu-scale_lambda*log(J))*F.Inverse();
    DIAGONAL_MATRIX<T,d> F_clamp=Clamp_To_Hyperbola(F),dF=F-F_clamp,F_inverse=F_clamp.Inverse();
    T scale_mu_minus_lambda_log_J=scale_mu-scale_lambda*log(failure_threshold);
    return scale_mu*F+scale_mu_minus_lambda_log_J*(sqr(F_inverse)*dF-F_inverse)+scale_lambda*DIAGONAL_MATRIX<T,d>::Inner_Product(F_inverse,dF)*F_inverse;
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> NEO_HOOKEAN<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    return 2*scale*constant_beta*strain_rate+scale*constant_alpha*strain_rate.Trace();
}
//#####################################################################
// Function P_From_Strain_Rate_Forces_Size
//#####################################################################
template<class T,int d> int NEO_HOOKEAN<T,d>::
P_From_Strain_Rate_Forces_Size() const
{
    return sizeof(MATRIX<T,d>)/sizeof(T);
}
//#####################################################################
// Function P_From_Strain_Rate_First_Half
//#####################################################################
template<class T,int d> void NEO_HOOKEAN<T,d>::
P_From_Strain_Rate_First_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<T> aggregate,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=scale*F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*constant_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(constant_alpha/TV::dimension+dd*dd)-dd;
    SYMMETRIC_MATRIX<T,d> s=sb*strain_rate+sa*strain_rate.Trace();
    *(MATRIX<T,d>*)aggregate.Get_Array_Pointer()+=s;
}
//#####################################################################
// Function P_From_Strain_Rate_Second_Half
//#####################################################################
template<class T,int d> MATRIX<T,d> NEO_HOOKEAN<T,d>::
P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<const T> aggregate,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=scale*(*(const MATRIX<T,d>*)aggregate.Get_Array_Pointer()).Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*constant_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(constant_alpha/TV::dimension+dd*dd)-dd;
    return sb*strain_rate+sa*strain_rate.Trace();
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void NEO_HOOKEAN<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,const int triangle) const
{
    DIAGONAL_MATRIX<T,2> F_inverse=F.Clamp_Min(failure_threshold).Inverse();
    T mu_minus_lambda_logJ=constant_mu+constant_lambda*log(F_inverse.Determinant());
    SYMMETRIC_MATRIX<T,2> F_inverse_outer=SYMMETRIC_MATRIX<T,2>::Outer_Product(F_inverse.To_Vector());
    dP_dF.x1111=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x11;//alpha+beta+gamma
    dP_dF.x2222=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x22;
    dP_dF.x2211=constant_lambda*F_inverse_outer.x21;//gamma
    dP_dF.x2121=constant_mu;//alpha
    dP_dF.x2112=mu_minus_lambda_logJ*F_inverse_outer.x21;//beta
    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void NEO_HOOKEAN<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF,const int tetrahedron) const
{
    DIAGONAL_MATRIX<T,3> F_inverse=F.Clamp_Min(failure_threshold).Inverse();
    T mu_minus_lambda_logJ=constant_mu+constant_lambda*log(F_inverse.Determinant());
    SYMMETRIC_MATRIX<T,3> F_inverse_outer=SYMMETRIC_MATRIX<T,3>::Outer_Product(F_inverse.To_Vector());
    dPi_dF.x1111=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x11;
    dPi_dF.x2222=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x22;
    dPi_dF.x3333=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x33;
    dPi_dF.x2211=constant_lambda*F_inverse_outer.x21;
    dPi_dF.x3311=constant_lambda*F_inverse_outer.x31;
    dPi_dF.x3322=constant_lambda*F_inverse_outer.x32;
    dPi_dF.x2121=constant_mu;dPi_dF.x3131=constant_mu;dPi_dF.x3232=constant_mu;
    dPi_dF.x2112=mu_minus_lambda_logJ*F_inverse_outer.x21;
    dPi_dF.x3113=mu_minus_lambda_logJ*F_inverse_outer.x31;
    dPi_dF.x3223=mu_minus_lambda_logJ*F_inverse_outer.x32;
    if(enforce_definiteness) dPi_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class T,int d> T NEO_HOOKEAN<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int simplex) const
{
    T I1=(F*F.Transposed()).Trace(),J=F.Determinant();
    if(J<=0) return (T)1e16; // TODO: Do something smarter here.
    T log_J=log(J);
    return constant_mu*((T).5*(I1-TV::m)-log_J)+(T).5*constant_lambda*sqr(log_J);
}
template class NEO_HOOKEAN<float,2>;
template class NEO_HOOKEAN<float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class NEO_HOOKEAN<double,2>;
template class NEO_HOOKEAN<double,3>;
#endif

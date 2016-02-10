//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MOONEY_RIVLIN_3D
//##################################################################### 
#ifndef __MOONEY_RIVLIN_3D__
#define __MOONEY_RIVLIN_3D__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
namespace PhysBAM{

template<class T>
class MOONEY_RIVLIN_3D:public ISOTROPIC_CONSTITUTIVE_MODEL<T,3>
{
    typedef VECTOR<T,3> TV;
public:
    typedef ISOTROPIC_CONSTITUTIVE_MODEL<T,3> BASE;
    using BASE::enforce_definiteness;using BASE::constant_lambda;using BASE::constant_mu;using BASE::constant_alpha;using BASE::constant_beta;

    T mu_10,mu_01,kappa;
    T failure_threshold;

    MOONEY_RIVLIN_3D(const T mu_10_input=(T)6e4,const T mu_01_input=(T)2e4,const T kappa_input=(T)6e4,const T Rayleigh_coefficient=.05,const T failure_threshold_input=.25)
        :mu_10(mu_10_input),mu_01(mu_01_input),kappa(kappa_input),failure_threshold(failure_threshold_input)
    {
        constant_mu=2*(mu_10+mu_01);
        constant_lambda=(T)one_third*(8*mu_01-4*mu_10)+kappa;
        constant_alpha=Rayleigh_coefficient*constant_lambda;
        constant_beta=Rayleigh_coefficient*constant_mu;
    }

    DIAGONAL_MATRIX<T,3> P_From_Strain(const DIAGONAL_MATRIX<T,3>& F,const T scale,const int tetrahedron) const PHYSBAM_OVERRIDE
    {DIAGONAL_MATRIX<T,3> F_threshold=F.Clamp_Min(failure_threshold),C=F_threshold*F_threshold,F_cube=C*F_threshold,F_inverse=F_threshold.Inverse();
    T I_C=C.Trace(),II_C=(C*C).Trace(),J=F_threshold.Determinant(),Jcc=pow(J,-(T)two_thirds);
    return (scale*2*Jcc*(mu_10+Jcc*mu_01*I_C))*F_threshold-(scale*2*Jcc*Jcc*mu_01)*F_cube+(scale*(kappa*log(J)-(T)two_thirds*Jcc*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))))*F_inverse;}
    
    MATRIX<T,3> P_From_Strain_Rate(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& F_dot,const T scale,const int tetrahedron) const PHYSBAM_OVERRIDE
    {SYMMETRIC_MATRIX<T,3> strain_rate=F_dot.Symmetric_Part(); 
    return 2*scale*constant_beta*strain_rate+scale*constant_alpha*strain_rate.Trace();}

    void Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF,const int tetrahedron) const PHYSBAM_OVERRIDE
    {DIAGONAL_MATRIX<T,3> F_threshold=F.Clamp_Min(failure_threshold),C=F_threshold*F_threshold,F_cube=C*F_threshold,F_inverse=F_threshold.Inverse();
    T I_C=C.Trace(),II_C=(C*C).Trace(),J=F_threshold.Determinant(),Jcc=pow(J,-(T)two_thirds);
    SYMMETRIC_MATRIX<T,3> alpha;
    alpha.x11=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x11-C.x11));alpha.x21=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x22-C.x11));alpha.x31=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x33-C.x11));
    alpha.x22=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x22-C.x22));alpha.x32=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x33-C.x22));alpha.x33=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x33-C.x33));
    SYMMETRIC_MATRIX<T,3> beta;
    beta.x11=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x11*F_inverse.x11-2*Jcc*Jcc*mu_01*F_threshold.x11*F_threshold.x11;
    beta.x21=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x22*F_inverse.x11-2*Jcc*Jcc*mu_01*F_threshold.x22*F_threshold.x11;
    beta.x31=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x33*F_inverse.x11-2*Jcc*Jcc*mu_01*F_threshold.x33*F_threshold.x11;
    beta.x22=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x22*F_inverse.x22-2*Jcc*Jcc*mu_01*F_threshold.x22*F_threshold.x22;
    beta.x32=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x33*F_inverse.x22-2*Jcc*Jcc*mu_01*F_threshold.x33*F_threshold.x22;
    beta.x33=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x33*F_inverse.x33-2*Jcc*Jcc*mu_01*F_threshold.x33*F_threshold.x33;
    SYMMETRIC_MATRIX<T,3> eta;
    eta.x11=4*Jcc*Jcc*mu_01;
    eta.x31=-Jcc*(T)one_third*(4*mu_10+8*Jcc*mu_01*I_C);
    eta.x32=8*(T)one_third*Jcc*Jcc*mu_01;
    eta.x33=Jcc*(T)one_ninth*(4*mu_10*I_C+Jcc*8*mu_01*(I_C*I_C-II_C))+kappa;
    MATRIX<T,3> F_base(F_threshold.x11,F_threshold.x22,F_threshold.x33,F_cube.x11,F_cube.x22,F_cube.x33,F_inverse.x11,F_inverse.x22,F_inverse.x33);
    SYMMETRIC_MATRIX<T,3> gamma=SYMMETRIC_MATRIX<T,3>::Conjugate(F_base,eta);
    dPi_dF.x1111=alpha.x11+beta.x11+gamma.x11;dPi_dF.x2222=alpha.x22+beta.x22+gamma.x22;dPi_dF.x3333=alpha.x33+beta.x33+gamma.x33;
    dPi_dF.x2211=gamma.x21;dPi_dF.x3311=gamma.x31;dPi_dF.x3322=gamma.x32;
    dPi_dF.x2121=alpha.x21;dPi_dF.x3131=alpha.x31;dPi_dF.x3232=alpha.x32;
    dPi_dF.x2112=beta.x21;dPi_dF.x3113=beta.x31;dPi_dF.x3223=beta.x31;
    if(enforce_definiteness) dPi_dF.Enforce_Definiteness();}

//#####################################################################
};
}
#endif

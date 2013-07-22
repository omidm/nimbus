//#####################################################################
// Copyright 2003-2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SIMPLE_PLASTICITY
//##################################################################### 
#ifndef __SIMPLE_PLASTICITY__
#define __SIMPLE_PLASTICITY__

#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/PLASTICITY_MODEL.h>
namespace PhysBAM{

template<class T,int d>
class SIMPLE_PLASTICITY:public PLASTICITY_MODEL<T,d>
{
    typedef VECTOR<T,d> TV;
public:
    typedef PLASTICITY_MODEL<T,d> BASE;
    using BASE::Fp_inverse;

    T yield_ratio,plastic_clamp_ratio;
private:
    T sqr_log_yield_ratio,sqr_log_plastic_clamp_ratio;
public:

    SIMPLE_PLASTICITY(const int elements,const T yield_ratio_input,const T plastic_clamp_ratio_input)
        :PLASTICITY_MODEL<T,d>(elements),yield_ratio(yield_ratio_input),plastic_clamp_ratio(plastic_clamp_ratio_input)
    {
        assert(yield_ratio>0 && plastic_clamp_ratio>0);
        sqr_log_yield_ratio=sqr(log(yield_ratio));sqr_log_plastic_clamp_ratio=sqr(log(plastic_clamp_ratio));
    }
    
    bool Project_Fe(const DIAGONAL_MATRIX<T,d>& Fe_trial,DIAGONAL_MATRIX<T,d>& Fe_project) const PHYSBAM_OVERRIDE
    {DIAGONAL_MATRIX<T,d> Fe_log=log(Fe_trial.Clamp_Min((T)1e-4));T dilation=Fe_log.Dilational();
    DIAGONAL_MATRIX<T,d> Fe_deviatoric=Fe_log-dilation;T deviatoric_sqr_norm=Fe_deviatoric.Frobenius_Norm_Squared();
    if(deviatoric_sqr_norm<=sqr_log_yield_ratio)return false;Fe_deviatoric*=sqrt(sqr_log_yield_ratio/deviatoric_sqr_norm);
    Fe_project=exp(Fe_deviatoric+dilation);return true;}
    
    void Project_Fp(const int simplex,const MATRIX<T,d>& Fp_trial) PHYSBAM_OVERRIDE
    {DIAGONAL_MATRIX<T,d> Fp_hat;MATRIX<T,d> U,V;Fp_trial.Fast_Singular_Value_Decomposition(U,Fp_hat,V);
    DIAGONAL_MATRIX<T,d> Fp_log_hat=log(Fp_hat.Clamp_Min((T)1e-4));Fp_log_hat-=Fp_log_hat.Dilational();T sqr_norm=Fp_log_hat.Frobenius_Norm_Squared();
    if(sqr_norm>sqr_log_plastic_clamp_ratio)Fp_log_hat*=sqrt(sqr_log_plastic_clamp_ratio/sqr_norm);
    Fp_inverse(simplex)=SYMMETRIC_MATRIX<T,d>::Conjugate(V,exp(-Fp_log_hat));}

//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2010, Mridul Aanjaneya, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __EOS_CLAMPED_INTERNAL_ENERGY__
#define __EOS_CLAMPED_INTERNAL_ENERGY__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS.h>
#include <cmath>
namespace PhysBAM{

using ::std::sqrt;

template<class T1>
class EOS_CLAMPED_INTERNAL_ENERGY:public EOS<T1>
{
private:
    const EOS<T1>& base_eos;
    T1 e_min;
    T1 epsilon;
    T1 a1,a2,a3;

public:
    EOS_CLAMPED_INTERNAL_ENERGY(const EOS<T1>& base_eos_input,T1 e_min_input,T1 epsilon_input):
    base_eos(base_eos_input),e_min(e_min_input),epsilon(epsilon_input)
    {
        Compute_Coefficients();
    }

    virtual ~EOS_CLAMPED_INTERNAL_ENERGY()
    {}

//#####################################################################
    // sound speed  
    virtual T1 p_rho(const T1 rho,const T1 e) const {T1 e_clamped=Smooth_Clamp_Internal_Energy(e);return base_eos.p_rho(rho,e_clamped);}
    virtual T1 p_e(const T1 rho,const T1 e) const {T1 e_clamped=Smooth_Clamp_Internal_Energy(e);return base_eos.p_e(rho,e_clamped)*L_e(e);}
    // pressure
    virtual T1 p(const T1 rho,const T1 e) const {T1 e_clamped=Smooth_Clamp_Internal_Energy(e);return base_eos.p(rho,e_clamped);}        
    virtual T1 rho_From_p_And_e(const T1 p,const T1 e) const {T1 e_clamped=Smooth_Clamp_Internal_Energy(e);return base_eos.rho_From_p_And_e(p,e_clamped);}
    virtual T1 e_From_p_And_rho(const T1 p,const T1 rho) const {return L_Inverse(base_eos.e_From_p_And_rho(p,rho));} 
    // temperature
    virtual T1 T(const T1 rho,const T1 e) const {T1 e_clamped=Smooth_Clamp_Internal_Energy(e);return base_eos.T(rho,e_clamped);}       
    virtual T1 rho_From_T_And_e(const T1 T,const T1 e) const {T1 e_clamped=Smooth_Clamp_Internal_Energy(e);return base_eos.rho_From_T_And_e(T,e_clamped);}
    virtual T1 e_From_T_And_rho(const T1 T,const T1 rho) const {return L_Inverse(base_eos.e_From_T_And_rho(T,rho));}  
    // entropy
    virtual T1 S(const T1 rho,const T1 e) const {T1 e_clamped=Smooth_Clamp_Internal_Energy(e);return base_eos.S(rho,e_clamped);}     
    virtual T1 rho_From_S_And_e(const T1 S,const T1 e) const {T1 e_clamped=Smooth_Clamp_Internal_Energy(e);return base_eos.rho_From_S_And_e(S,e_clamped);}
    virtual T1 e_From_S_And_rho(const T1 S,const T1 rho) const {return L_Inverse(base_eos.e_From_S_And_rho(S,rho));} 
    // pressure and temperature
    virtual T1 rho_From_p_And_T(const T1 p,const T1 T) const {return base_eos.rho_From_p_And_T(p,T);}
    virtual T1 e_From_p_And_T(const T1 p,const T1 T) const {return L_Inverse(base_eos.e_From_p_And_T(p,T));}
    // pressure and entropy
    virtual T1 rho_From_p_And_S(const T1 p,const T1 S) const {return base_eos.rho_From_p_And_S(p,S);}
    virtual T1 e_From_p_And_S(const T1 p,const T1 S) const {return L_Inverse(base_eos.e_From_p_And_S(p,S));} 
    // sound speed
    virtual T1 c(const T1 rho,const T1 e) const {T1 e_clamped=Smooth_Clamp_Internal_Energy(e);return base_eos.c(rho,e_clamped);}
    virtual T1 one_over_c(const T1 rho,const T1 e) const {T1 e_clamped=Smooth_Clamp_Internal_Energy(e);return base_eos.one_over_c(rho,e_clamped);}

private:

    void Compute_Coefficients()
    {
        T1 one_over_epsilon = 1/epsilon;
        a1 = one_over_epsilon*.25;
        a2 = -(e_min - epsilon)*one_over_epsilon*.5;
        a3 = (e_min + epsilon)*(e_min + epsilon)*one_over_epsilon*.25;
    }

    T1 Smooth_Clamp_Internal_Energy(const T1 e) const
    {
        if(e > e_min + epsilon) return e;
        else if(e >= e_min - epsilon) return (a1*e*e + a2*e + a3);
        else return e_min;
    }

    T1 L_Inverse(T1 e) const
    {
        if(e>e_min+epsilon) return e;
        else PHYSBAM_FATAL_ERROR("WARNING: Inverse of internal energy not possible!!!"); 
    }

    T1 L_e(T1 e) const
    {
        if(e > e_min + epsilon) return (T1)1;
        else if(e >= e_min - epsilon) return (2*a1*e + a2);
        else return (T1)0;
    }
//#####################################################################
};   
}
#endif


//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EOS_MG  
//##################################################################### 
//
// Mie Gruneisen equation of state.
// Inherits the virtual base class EOS, overwriting its functions.
//
//#####################################################################
#ifndef __EOS_MG__
#define __EOS_MG__

#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS.h>
namespace PhysBAM{     

template<class T>
class EOS_MG:public EOS<T>
{
private:
    T gamma; 
    T s;       // scale factor
    T c_not;   // ambient sound speed???
    T rho_not; // ambient density 

public:
    EOS_MG()
    {
        Set_gamma();   // default is HMX
        Set_s();       // default is HMX
        Set_c_not();   // default is HMX
        Set_rho_not(); // default is HMX
    }          

    void Set_gamma(const T gamma_input = 1.7)
    {gamma=gamma_input;} 

    void Set_s(const T s_input = 1.79) // no units 
    {s=s_input;}
    
    void Set_c_not(const T c_not_input = 3.07e3) // units are m/s 
    {c_not=c_not_input;}

    void Set_rho_not(const T rho_not_input = 1891) // units are kg/m^3 
    {rho_not=rho_not_input;}

//#####################################################################
// Function f
//#####################################################################   
// p=(gamma-1)*rho*e+f(rho)
T f(const T rho)
    {
        return (1-(gamma-1)/2*(rho/rho_not-1))*rho_not*sqr(c_not)/s*
               (1-1/(1-s*(1-rho_not/rho)));
    } 
//#####################################################################
// Function f_prime
//#####################################################################   
// derivative of f with respect to rho
    T f_prime(const T rho)
    {
        return rho_not*sqr(c_not)/s*(-(gamma-1)/(2*rho_not)*(1-1/(1-s*(1-rho_not/rho)))
               -s*rho_not/pow(rho*(1-s*(1-rho_not/rho)),2)*(1-(gamma-1)/2*(rho/rho_not-1)));
    } 
//#####################################################################
// Function p_rho
//#####################################################################
// partial derivative of the pressure
    T p_rho(const T rho,const T e) const PHYSBAM_OVERRIDE
    {
        return (gamma-1)*e+f_prime(rho);
    }
//#####################################################################
// Function p_e
//#####################################################################
// partial derivative of the pressure - e is not needed
    T p_e(const T rho,const T e) const PHYSBAM_OVERRIDE
    {
        return (gamma-1)*rho;
    }
//#####################################################################
// Function p
//#####################################################################
// pressure
    T p(const T rho,const T e) const PHYSBAM_OVERRIDE
    {
        return (gamma-1)*rho*e+f(rho);
    }
//#####################################################################
// Function e_From_p_And_rho
//#####################################################################   
    T e_From_p_And_rho(const T p,const T rho) const PHYSBAM_OVERRIDE
    {
        return (p-f(rho))/((gamma-1)*rho);
    } 
//#####################################################################
};
}    
#endif


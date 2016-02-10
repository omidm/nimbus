//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EOS_JWL  
//##################################################################### 
//
// JWL gas.
// Inherits the virtual base class EOS, overwriting its functions.
//
//#####################################################################
#ifndef __EOS_JWL__
#define __EOS_JWL__

#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS.h>
namespace PhysBAM{

template<class T>
class EOS_JWL:public EOS<T>
{
private:
    T gamma,rho_not,A1,A2,R1,R2; 

public:
    EOS_JWL()
    {
        Set_gamma();  
        Set_rho_not(); 
        Set_A1();      
        Set_A2();      
        Set_R1();     
        Set_R2();     
    }

    void Set_gamma(const T gamma_input = 1.28) // no units  
    {gamma=gamma_input;}
    
    void Set_rho_not(const T rho_not_input = 1630) // units are kg/m^3  
    {rho_not=rho_not_input;}
    
    void Set_A1(const T A1_input = 5.484e11) // units are Pa  
    {A1=A1_input;}         

    void Set_A2(const T A2_input = 9.375e9) // units are Pa  
    {A2=A2_input;}
    
    void Set_R1(const T R1_input = 4.94) // no units  
    {R1=R1_input;}         

    void Set_R2(const T R2_input = 1.21) // no units  
    {R2=R2_input;} 
    
//#####################################################################
// Function p_rho
//#####################################################################
// partial derivative of the pressure 
    T p_rho(const T rho,const T e) const PHYSBAM_OVERRIDE
    {
        return (gamma-1)*e+A1*(R1*rho_not/sqr(rho)-(gamma-1)*(1/rho+1/(R1*rho_not)))*exp(-R1*rho_not/rho)+
                                        A2*(R2*rho_not/sqr(rho)-(gamma-1)*(1/rho+1/(R2*rho_not)))*exp(-R2*rho_not/rho);
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
        return (gamma-1)*rho*e+A1*(1-(gamma-1)*rho/(R1*rho_not))*exp(-R1*rho_not/rho)+
                                               A2*(1-(gamma-1)*rho/(R2*rho_not))*exp(-R2*rho_not/rho);
    }
//#####################################################################
// Function e_From_p_And_rho
//#####################################################################   
    T e_From_p_And_rho(const T p,const T rho) const PHYSBAM_OVERRIDE
    {
        return (p-A1*(1-(gamma-1)*rho/(R1*rho_not))*exp(-R1*rho_not/rho)
                      -A2*(1-(gamma-1)*rho/(R2*rho_not))*exp(-R2*rho_not/rho))/((gamma-1)*rho);
    } 
//#####################################################################
};   
}
#endif


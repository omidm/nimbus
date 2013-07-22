//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EOS_Tait  
//##################################################################### 
//
// Tait equation of state for water.
// Inherits the virtual base class EOS, overwriting its functions.
//
//#####################################################################
#ifndef __EOS_TAIT__
#define __EOS_TAIT__

#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS.h>
namespace PhysBAM{

using ::std::pow;

template<class T>
class EOS_TAIT:public EOS<T>
{
private:
    T gamma;
    T A;       // scale factor
    T rho_not; // ambient density
    T p_not;   // ambient pressure
public:
    T p_min;   // minimum cavitation pressure

public:
    EOS_TAIT()
    {
        Set_gamma();   // default
        Set_A();       // default
        Set_rho_not(); // default
        Set_p_not();   // default
        Set_p_min();   // default
    }         

    void Set_gamma(const T gamma_input = 7.15) // default is water
    {gamma=gamma_input;} 

    void Set_A(const T A_input = 3.31e8) // units are Pa
    {A=A_input;}
    
    void Set_rho_not(const T rho_not_input = 1e3) // units are kg/m^3
    {rho_not=rho_not_input;}
    
    void Set_p_not(const T p_not_input = 1e5) // units are Pa
    {p_not=p_not_input;}

    void Set_p_min(const T p_min_input = 22.02726) // units are Pa
    {p_min=p_min_input;}


//#####################################################################
// Function p_rho
//#####################################################################
// partial derivative of the pressure - e is not needed
    T p_rho(const T rho,const T e) const PHYSBAM_OVERRIDE
    {
        if(p(rho,e) == p_min) return 0; // pressure no longer depends on density
        else return A*gamma*pow(rho,gamma-1)/pow(rho_not,gamma);
    }
//#####################################################################
// Function p_e
//#####################################################################
// partial derivative of the pressure - neither rho nor p is needed
    T p_e(const T rho,const T e) const PHYSBAM_OVERRIDE
    {
        return 0; // pressure doesn't depend on internal energy
    }
//#####################################################################
// Function p
//#####################################################################
// pressure - e is not needed
    T p(const T rho,const T e) const PHYSBAM_OVERRIDE
    {
        return max(A*(pow(rho/rho_not,gamma)-1)+p_not,p_min);

    }
//#####################################################################
// Function rho_From_p_And_e
//#####################################################################   
// e is not needed
    T rho_From_p_And_e(const T p,const T e) const PHYSBAM_OVERRIDE
    {
        if(p <= p_min){
            std::stringstream ss;ss << "rho_From_p_And_e ERROR: WHEN p<=p_min, rho IS UNKNOWN!" << std::endl;LOG::filecout(ss.str()); 
            return rho_not*pow((p_min-p_not)/A+1,1/gamma);}
        else return rho_not*pow((p-p_not)/A+1,1/gamma);
    }
//#####################################################################
// Function e_From_p_And_rho
//#####################################################################
// p is not needed
    T e_From_p_And_rho(const T p,const T rho) const PHYSBAM_OVERRIDE
    {
        return A/pow(rho_not,gamma)*pow(rho,gamma-1)/(gamma-1)+(A-p_not)/rho;
    }  
//#####################################################################
};   
}
#endif


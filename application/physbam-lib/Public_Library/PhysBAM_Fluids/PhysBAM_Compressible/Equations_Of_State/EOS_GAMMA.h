//#####################################################################
// Copyright 2002-2005, Doug Enright, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EOS_GAMMA  
//##################################################################### 
//
// Gamma Law gas.
// Inherits the virtual base class EOS, overwriting its functions.
//
//#####################################################################
#ifndef __EOS_GAMMA__
#define __EOS_GAMMA__

#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS.h>
namespace PhysBAM{

using std::pow;

template<class T1>
class EOS_GAMMA:public EOS<T1>
{
public:
    typedef T1 SCALAR; 

    T1 gamma; // ratio of specific heats
    T1 Ru;    // universal gas constant - J/(mol K)
    T1 MW;    // molecular weight
    T1 R;     // specific gas constant
    T1 Cv;  // specific heat at constant volume 
    T1 e_not; // internal energy at 0K

public:
    EOS_GAMMA():EOS<T1>()
    {
        Ru=(T1)8.31451; 
        Cv=0;
        Set_Molecular_Weight(); 
        Set_Gamma();   
        Set_Internal_Energy_Of_Formation();
        Calculate_Cv();
    }

    void Set_Molecular_Weight(const T1 MW_input=29) // units are g/mol - default is air  
    {MW=MW_input/(T1)1000; // converts to  kg/mol
    R=Ru/MW; // J/(kg K)
    if(Cv) Calculate_Cv();} // update Cv if MW is changed

    void Set_Gamma(const T1 gamma_input=1.4) // default is air
    {gamma=gamma_input;
    if(Cv) Calculate_Cv();} // update Cv if gamma is changed

    void Set_Internal_Energy_Of_Formation(const T1 e_not_input=0)
    {e_not=e_not_input;}

    void Calculate_Cv()
    {Cv = R/(gamma-1);} // units are J/(kg K)
    
//#####################################################################
// Function p_rho
//#####################################################################
// partial derivative of the pressure - rho is not needed
    T1 p_rho(const T1 rho,const T1 e) const PHYSBAM_OVERRIDE
    {
        return (gamma-1)*(e-e_not);
    }
//#####################################################################
// Function p_e
//#####################################################################
// partial derivative of the pressure - e is not needed
    T1 p_e(const T1 rho,const T1 e) const PHYSBAM_OVERRIDE
    {
        return (gamma-1)*rho;
    }
//#####################################################################
// Function p
//#####################################################################
// pressure
    T1 p(const T1 rho,const T1 e) const PHYSBAM_OVERRIDE
    {
        return (gamma-1)*rho*(e-e_not);
    }
//#####################################################################
// Function rho_From_p_And_e
//#####################################################################   
    T1 rho_From_p_And_e(const T1 p,const T1 e) const PHYSBAM_OVERRIDE
    {
        return p/((gamma-1)*(e-e_not));
    }
//#####################################################################
// Function e_From_p_And_rho
//#####################################################################   
    T1 e_From_p_And_rho(const T1 p,const T1 rho) const PHYSBAM_OVERRIDE
    {
        return p/((gamma-1)*rho)+e_not;
    } 
//#####################################################################
// Function T
//#####################################################################   
// temperature - rho is not needed
    T1 T(const T1 rho,const T1 e) const PHYSBAM_OVERRIDE
    {
        return (e-e_not)/Cv;
    }       
//#####################################################################
// Function e_From_T_And_rho
//#####################################################################   
// rho is not needed
    T1 e_From_T_And_rho(const T1 T,const T1 rho) const PHYSBAM_OVERRIDE
    {
        return e_not+Cv*T;
    }  
//#####################################################################
// Function S
//#####################################################################   
// entropy
    T1 S(const T1 rho,const T1 e) const PHYSBAM_OVERRIDE
    {
        return p(rho,e)/pow(rho,gamma);
    }     
//#####################################################################
// Function rho_From_S_And_e
//#####################################################################
    T1 rho_From_S_And_e(const T1 S,const T1 e) const PHYSBAM_OVERRIDE
    {
        return pow((gamma-1)*e/S,1/(gamma-1));
    }
//#####################################################################
// Function e_From_S_And_rho
//#####################################################################
    T1 e_From_S_And_rho(const T1 S,const T1 rho) const PHYSBAM_OVERRIDE
    {
        return e_not+S*pow(rho,gamma-1)/(gamma-1);
    } 
//#####################################################################
// Function rho_From_p_And_T
//#####################################################################
    T1 rho_From_p_And_T(const T1 p,const T1 T) const PHYSBAM_OVERRIDE
    {
        return p/(R*T);
    }
//#####################################################################
// Function e_From_p_And_T
//#####################################################################
// p is not needed
    T1 e_From_p_And_T(const T1 p,const T1 T) const PHYSBAM_OVERRIDE
    {
        return e_not+Cv*T;
    }
//#####################################################################
// Function rho_From_p_And_S
//#####################################################################
    T1 rho_From_p_And_S(const T1 p,const T1 S) const PHYSBAM_OVERRIDE
    {
        return pow(p/S,1/gamma);
    }
//#####################################################################
// Function e_From_p_And_S
//#####################################################################
    T1 e_From_p_And_S(const T1 p,const T1 S) const PHYSBAM_OVERRIDE
    {
        return e_From_S_And_rho(S,rho_From_p_And_S(p,S));
    }  
//#####################################################################
};   
}
#endif


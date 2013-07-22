//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REACTIVE_EOS_GAMMA  
//##################################################################### 
//
// Gamma Law gas.
// Inherits the virtual base class EOS, overwriting its functions.
// assume Y is the mass fraction of gas #1 
//
//#####################################################################
#ifndef __REACTIVE_EOS_GAMMA__
#define __REACTIVE_EOS_GAMMA__

#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/REACTIVE_EOS.h>
namespace PhysBAM{

template<class TS>
class REACTIVE_EOS_GAMMA:public REACTIVE_EOS<TS>
{
    using REACTIVE_EOS<TS>::e_not1;using REACTIVE_EOS<TS>::e_not2;
public:
    TS Ru;    // universal gas constant - J/(mol K)
private:
    TS MW1,MW2; // molecular weights
    TS R1,R2; // specific gas constants
    TS gamma1,gamma2; // ratio of specific heats
    TS Cv1,Cv2;  // specific heat at constant volume 

public:
    REACTIVE_EOS_GAMMA()
    {
        Ru=8.31451; 
        Cv1=0;Cv2=0;
        Set_Molecular_Weight(); 
        Set_Gamma();   
        Set_Internal_Energy_Of_Formation();
        Calculate_Cv();
    }

    void Set_Molecular_Weight(const TS MW1_input=29,const TS MW2_input=29) // units are g/mol - default is air  
    {MW1=MW1_input*1000;MW2=MW2_input*1000;  // converts to  kg/mol
     R1=Ru/MW1;R2=Ru/MW2;  // J/(kg K)
     if(Cv1 || Cv2) Calculate_Cv();} // update Cv if MW is changed          

    void Set_Gamma(const TS gamma1_input=1.4,const TS gamma2_input=1.4) // default is air
    {gamma1=gamma1_input;gamma2=gamma2_input;
     if(Cv1 || Cv2) Calculate_Cv();} // update Cv if gamma is changed

    void Set_Internal_Energy_Of_Formation(const TS e_not1_input=0,const TS e_not2_input=0)
    {e_not1=e_not1_input;e_not2=e_not2_input;}

    void Calculate_Cv()
    {Cv1=R1/(gamma1-1);Cv2=R2/(gamma2-1);} // units are J/(kg K)
    
    TS MW(const TS Y)
    {return 1/(Y/MW1+(1-Y)/MW2);}

    TS R(const TS Y)
    {return Y*R1+(1-Y)*R2;}
    
    TS gamma(const TS Y)
    {TS Cv_total=Y*Cv1+(1-Y)*Cv2,R_total=Y*R1+(1-Y)*R2;return (Cv_total+R_total)/Cv_total;}
    
    TS Cv(const TS Y)
    {return Y*Cv1+(1-Y)*Cv2;}
    
    TS e_not(const TS Y)
    {return Y*e_not1+(1-Y)*e_not2;}

//#####################################################################
// Function p_rho
//#####################################################################
// partial derivative of the pressure - rho is not needed
    TS p_rho(const TS rho,const TS e,const TS Y) PHYSBAM_OVERRIDE
    {
        return (gamma(Y)-1)*(e-e_not(Y));
    }
//#####################################################################
// Function p_e
//#####################################################################
// partial derivative of the pressure - e is not needed
    TS p_e(const TS rho,const TS e,const TS Y) PHYSBAM_OVERRIDE
    {
        return (gamma(Y)-1)*rho;
    }
//#####################################################################
// Function p
//#####################################################################
// pressure
    TS p(const TS rho,const TS e,const TS Y) PHYSBAM_OVERRIDE
    {
        return (gamma(Y)-1)*rho*(e-e_not(Y));
    }
//#####################################################################
// Function rho_From_p_And_e
//#####################################################################   
    TS rho_From_p_And_e(const TS p,const TS e,const TS Y) PHYSBAM_OVERRIDE
    {
        return p/((gamma(Y)-1)*(e-e_not(Y)));
    }
//#####################################################################
// Function e_From_p_And_rho
//#####################################################################   
    TS e_From_p_And_rho(const TS p,const TS rho,const TS Y) PHYSBAM_OVERRIDE
    {
        return p/((gamma(Y)-1)*rho)+e_not(Y);
    } 
//#####################################################################
// Function T
//#####################################################################   
// temperature - rho is not needed
    TS T(const TS rho,const TS e,const TS Y) PHYSBAM_OVERRIDE
    {
        return (e-e_not(Y))/Cv(Y);
    }       
//#####################################################################
// Function e_From_T_And_rho
//#####################################################################   
// rho is not needed
    TS e_From_T_And_rho(const TS T,const TS rho,const TS Y) PHYSBAM_OVERRIDE
    {
        return e_not(Y)+Cv(Y)*T;
    }  
//#####################################################################
// Function rho_From_p_And_T
//#####################################################################
    TS rho_From_p_And_T(const TS p,const TS T,const TS Y) PHYSBAM_OVERRIDE
    {
        return p/(R(Y)*T);
    }
//#####################################################################
// Function e_From_p_And_T
//#####################################################################
// p is not needed
    TS e_From_p_And_T(const TS p,const TS T,const TS Y) PHYSBAM_OVERRIDE
    {
        return e_not(Y)+Cv(Y)*T;
    }
//#####################################################################
};   
}
#endif


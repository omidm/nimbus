//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REACTIVE_EOS
//##################################################################### 
//
// Virtual base class.
// All functions are overwritten by the equations of state classes that 
// inherit this virtual base class.
//
//#####################################################################
#ifndef __REACTIVE_EOS__
#define __REACTIVE_EOS__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <cmath>
namespace PhysBAM{

template<class TS>
class REACTIVE_EOS
{
public:
    TS e_not1,e_not2; // internal energy at 0K
private:
    bool check_sound_speed; // (1) check for c<=0, (0) don't check for c<=0

protected:    
    REACTIVE_EOS()
    {      
        Check_Sound_Speed();
    }

public:
    virtual ~REACTIVE_EOS()
    {}

    void Check_Sound_Speed(bool check_input=true)
    {check_sound_speed=check_input;}

//#####################################################################
    // sound speed  
    virtual TS p_rho(const TS rho,const TS e,const TS Y){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual TS p_e(const TS rho,const TS e,const TS Y){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    // pressure
    virtual TS p(const TS rho,const TS e,const TS Y){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}        
    virtual TS rho_From_p_And_e(const TS p,const TS e,const TS Y){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual TS e_From_p_And_rho(const TS p,const TS rho,const TS Y){PHYSBAM_FUNCTION_IS_NOT_DEFINED();} 
    // temperature
    virtual TS T(const TS rho,const TS e,const TS Y){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}       
    virtual TS rho_From_T_And_e(const TS T,const TS e,const TS Y){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual TS e_From_T_And_rho(const TS T,const TS rho,const TS Y){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}  
    // pressure and temperature
    virtual TS rho_From_p_And_T(const TS p,const TS T,const TS Y){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual TS e_From_p_And_T(const TS p,const TS T,const TS Y){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//
//#####################################################################
// Function c
//#####################################################################
// sound speed
    virtual TS c(const TS rho,const TS e,const TS Y)
    {
        TS c_squared=p_rho(rho,e,Y)+p(rho,e,Y)/sqr(rho)*p_e(rho,e,Y);
        if(c_squared > 0) return sqrt(c_squared);
        else{ // sound speed not hyperbolic
            if(check_sound_speed){std::stringstream ss;ss << "IMAGINARY SOUND SPEED: c^2 = " << c_squared << std::endl;LOG::filecout(ss.str());}
            return 0;} 
    }
//#####################################################################
};   
}
#endif


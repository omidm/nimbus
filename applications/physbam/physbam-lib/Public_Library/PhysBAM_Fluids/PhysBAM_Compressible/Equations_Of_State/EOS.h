//#####################################################################
// Copyright 2002-2003, Doug Enright, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EOS  
//##################################################################### 
//
// Virtual base class.
// All functions are overwritten by the equations of state classes that 
// inherit this virtual base class.
//
//#####################################################################
#ifndef __EOS__
#define __EOS__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <cmath>
namespace PhysBAM{

using ::std::sqrt;

template<class T1>
class EOS
{
private:
    bool check_sound_speed; // (true) check for c<=0, (false) don't check for c<=0

protected:    
    EOS()
    {      
        Check_Sound_Speed();
    }

public:
    virtual ~EOS()
    {}

    void Check_Sound_Speed(bool check_input=true)
    {check_sound_speed=check_input;}

//#####################################################################
    // sound speed  
    virtual T1 p_rho(const T1 rho,const T1 e) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T1 p_e(const T1 rho,const T1 e) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    // pressure
    virtual T1 p(const T1 rho,const T1 e) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}        
    virtual T1 rho_From_p_And_e(const T1 p,const T1 e) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T1 e_From_p_And_rho(const T1 p,const T1 rho) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();} 
    // temperature
    virtual T1 T(const T1 rho,const T1 e) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}       
    virtual T1 rho_From_T_And_e(const T1 T,const T1 e) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T1 e_From_T_And_rho(const T1 T,const T1 rho) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}  
    // entropy
    virtual T1 S(const T1 rho,const T1 e) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}     
    virtual T1 rho_From_S_And_e(const T1 S,const T1 e) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T1 e_From_S_And_rho(const T1 S,const T1 rho) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();} 
    // pressure and temperature
    virtual T1 rho_From_p_And_T(const T1 p,const T1 T) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T1 e_From_p_And_T(const T1 p,const T1 T) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    // pressure and entropy
    virtual T1 rho_From_p_And_S(const T1 p,const T1 S) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T1 e_From_p_And_S(const T1 p,const T1 S) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();} 
//#####################################################################
// Function c
//#####################################################################
// sound speed
    virtual T1 c(const T1 rho,const T1 e) const
    {
        T1 c_squared=p_rho(rho,e)+p(rho,e)/sqr(rho)*p_e(rho,e);
        assert(rho>0);assert(e>0);
        if(c_squared > 0) return sqrt(c_squared);
        else{ // sound speed not hyperbolic
            if(check_sound_speed){std::stringstream ss;ss << "IMAGINARY SOUND SPEED: c^2 = " << c_squared << std::endl;LOG::filecout(ss.str());}
            return (T1)0;} 
    }
    virtual T1 c_squared(const T1 rho,const T1 e) const
    {
        T1 c_squared=p_rho(rho,e)+p(rho,e)/sqr(rho)*p_e(rho,e);
        assert(rho>0);assert(e>0);
        if(c_squared > 0) return c_squared;
        else{ // sound speed not hyperbolic
            if(check_sound_speed){std::stringstream ss;ss << "IMAGINARY SOUND SPEED: c^2 = " << c_squared << std::endl;LOG::filecout(ss.str());}
            return (T1)0;} 
    }
    virtual T1 one_over_c(const T1 rho,const T1 e) const
    {
        return (T1)1/c(rho,e);
    }
//#####################################################################
};   
}
#endif


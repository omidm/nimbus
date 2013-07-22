//#####################################################################
// Copyright 2009, Nipun Kwatra, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EOS_SMOOTH_TRANSITION_INCOMPRESSIBLE  
//##################################################################### 
#ifndef __EOS_SMOOTH_TRANSITION_INCOMPRESSIBLE__
#define __EOS_SMOOTH_TRANSITION_INCOMPRESSIBLE__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS.h>
namespace PhysBAM{

template<class T_EOS>
class EOS_SMOOTH_TRANSITION_INCOMPRESSIBLE:public T_EOS
{
    typedef typename T_EOS::SCALAR T;

public:
    T t_current,t_start_transition,t_end_transition;
    T one_over_c_incompressible;
    bool use_inverse_polynomial;
    T exponent;

public:
    EOS_SMOOTH_TRANSITION_INCOMPRESSIBLE(const T t_start_transition_input,const T t_end_transition_input,
        const T one_over_c_incompressible_input,const bool use_inverse_polynomial_input=false,const T exponent_input=(T)1)
        :T_EOS(),t_current(0),t_start_transition(t_start_transition_input),t_end_transition(t_end_transition_input),
        one_over_c_incompressible(one_over_c_incompressible_input),use_inverse_polynomial(use_inverse_polynomial_input),
        exponent(exponent_input)
    {if(use_inverse_polynomial && one_over_c_incompressible!=0) PHYSBAM_FATAL_ERROR("one_over_c_incompressible should be 0 for inverse polynomial interpolation");}

    void Set_Current_Time(const T time_current_input)
    {t_current=time_current_input;}

    T Inverse_Polynomial(const T a,const T alpha) const
    {return a*pow((T)1-alpha,exponent);}
  
    virtual T one_over_c(const T rho,const T e) const PHYSBAM_OVERRIDE
    {T one_over_c_real=T_EOS::one_over_c(rho,e);
    if(t_current<t_start_transition) return one_over_c_real;
    else if(t_current>t_end_transition) return one_over_c_incompressible;
    else{T alpha=(t_current-t_start_transition)/(t_end_transition-t_start_transition);
        if(use_inverse_polynomial){return Inverse_Polynomial(one_over_c_real,alpha);}
        else return one_over_c_real+alpha*(one_over_c_incompressible-one_over_c_real);}}
//#####################################################################
};
}
#endif

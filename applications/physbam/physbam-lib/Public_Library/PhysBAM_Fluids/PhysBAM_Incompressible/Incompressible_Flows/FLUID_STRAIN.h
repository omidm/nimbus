//#####################################################################
// Copyright 2004, Doug Enright, Geoffrey Irving, Andy Lutimirski, Paul-James White.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_STRAIN  
//#####################################################################
#ifndef __FLUID_STRAIN__
#define __FLUID_STRAIN__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

template<class T>
class FLUID_STRAIN:public NONCOPYABLE
{
public:
    T viscosity_index,strainrate_time;
    T elastic_modulus;
    T plasticity_alpha,plasticity_gamma;

    FLUID_STRAIN()
    {
        Set_Elastic_Modulus();
        Set_Plasticity_Components();
        Set_Viscosity_Index();
        Set_Strainrate_Time();
    }
    
    ~FLUID_STRAIN()
    {}

    void Set_Elastic_Modulus(const T elastic_modulus_input=0)
    {elastic_modulus=elastic_modulus_input;}

    void Set_Plasticity_Components(const T plasticity_alpha_input=0,const T plasticity_gamma_input=0)
    {plasticity_alpha=plasticity_alpha_input;plasticity_gamma=plasticity_gamma_input;}

    void Set_Viscosity_Index(const T viscosity_index_input=1)
    {viscosity_index=viscosity_index_input;}

    void Set_Strainrate_Time(const T strainrate_time_input=0)
    {strainrate_time=strainrate_time_input;}

//#####################################################################
};
}
#endif


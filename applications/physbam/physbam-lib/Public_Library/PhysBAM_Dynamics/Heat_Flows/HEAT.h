//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HEAT
//#####################################################################
#ifndef __HEAT__
#define __HEAT__

#include <cfloat>
namespace PhysBAM{

template<class T>
class HEAT
{
protected:               
    T max_time_step;
    T density;
    T specific_heat;
    T thermal_conductivity; 
    T kappa;

    HEAT()
    {
        Set_Max_Time_Step();    
        Set_Density();
        Set_Specific_Heat();
        Set_Thermal_Conductivity();
    }
    
public: 
    void Set_Max_Time_Step(const T max_time_step_input=FLT_MAX)
    {max_time_step=max_time_step_input;}

    void Set_Density(const T density_input=1)
    {density=density_input;Calculate_Kappa();}

    void Set_Specific_Heat(const T specific_heat_input=1)
    {specific_heat=specific_heat_input;Calculate_Kappa();}

    void Set_Thermal_Conductivity(const T thermal_conductivity_input=1)
    {thermal_conductivity=thermal_conductivity_input;Calculate_Kappa();}

    void Calculate_Kappa()
    {if(density && specific_heat) kappa=thermal_conductivity/(density*specific_heat);}

//#####################################################################
};
}
#endif


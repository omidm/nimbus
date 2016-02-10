//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_LAGRANGE 
//#####################################################################
//
// Inherited by GRID_LAGRANGE_1D, GRID_LAGRANGE_2D, & GRID_LAGRANGE_3D. 
//
//#####################################################################
#ifndef __GRID_LAGRANGE__
#define __GRID_LAGRANGE__

#include <cfloat>
namespace PhysBAM{

template<class T>
class GRID_LAGRANGE
{
private:
    T max_convection_time_step; // default=FLT_MAX

public:
    GRID_LAGRANGE()
    {
        Set_Max_Convection_Time_Step();
    }

    void Set_Max_Convection_Time_Step(const T max_convection_time_step_input=FLT_MAX)
    {max_convection_time_step=max_convection_time_step_input;}

//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HAMILTONIAN_NORMAL_ADVECTION_1D 
//##################################################################### 
#ifndef __HAMILTONIAN_NORMAL_ADVECTION_1D__
#define __HAMILTONIAN_NORMAL_ADVECTION_1D__   

#include <PhysBAM_Dynamics/Level_Sets/HAMILTONIAN_1D.h>
namespace PhysBAM{

template<class T>
class HAMILTONIAN_NORMAL_ADVECTION_1D:public HAMILTONIAN_1D<T>
{
    typedef VECTOR<T,1> TV;
public:
    T speed;

    HAMILTONIAN_NORMAL_ADVECTION_1D(GRID<TV>& grid_input)
        :HAMILTONIAN_1D<T>(grid_input)
    {
        Set_Speed(1);
    }

    void Set_Speed(const T speed_input)
    {speed=speed_input;}

//#####################################################################
    T H(const T phi_x,const int i=0,const T t=0) PHYSBAM_OVERRIDE;
    T Maxabs_H1(const T phi_x_1,const T phi_x_2,const int i=0,const T t=0) PHYSBAM_OVERRIDE;
//#####################################################################
};   
}
#endif


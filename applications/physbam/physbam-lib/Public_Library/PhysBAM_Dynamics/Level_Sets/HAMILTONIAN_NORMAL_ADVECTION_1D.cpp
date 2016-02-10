//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Dynamics/Level_Sets/HAMILTONIAN_NORMAL_ADVECTION_1D.h>
#include <cmath>
using ::std::abs;
using namespace PhysBAM;
//#####################################################################
// Function H
//#####################################################################
template<class T> T HAMILTONIAN_NORMAL_ADVECTION_1D<T>::
H(const T phi_x,const int i,const T t)
{       
    return speed*abs(phi_x);
}
//#####################################################################
// Function Maxabs_H1
//#####################################################################
template<class T> T HAMILTONIAN_NORMAL_ADVECTION_1D<T>::
Maxabs_H1(const T phi_x_1,const T phi_x_2,const int i,const T t)
{       
    return abs(speed);
}
//#####################################################################
template class HAMILTONIAN_NORMAL_ADVECTION_1D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class HAMILTONIAN_NORMAL_ADVECTION_1D<double>;
#endif

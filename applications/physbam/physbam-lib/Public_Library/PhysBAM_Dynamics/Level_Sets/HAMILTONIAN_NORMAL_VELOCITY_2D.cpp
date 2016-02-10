//#####################################################################
// Copyright 2002, Ronald Fedkiw, Frederic Gibou.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/maxabs.h>
#include <PhysBAM_Tools/Math_Tools/minabs.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Dynamics/Level_Sets/HAMILTONIAN_NORMAL_VELOCITY_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function H
//#####################################################################
template<class T> T HAMILTONIAN_NORMAL_VELOCITY_2D<T>::
H(const T phi_x,const T phi_y,const int i,const int j,const T t)
{       
    return speed(i,j)*sqrt(sqr(phi_x)+sqr(phi_y));
}
//#####################################################################
// Function Maxabs_H1
//#####################################################################
template<class T> T HAMILTONIAN_NORMAL_VELOCITY_2D<T>::
Maxabs_H1(const T phi_x_1,const T phi_x_2,const T phi_y_1,const T phi_y_2,const int i,const int j,const T t)
{       
    if(phi_y_1*phi_y_2 <= 0) return abs(speed(i,j)); // use phi_y=0
    T phi_x=maxabs(phi_x_1,phi_x_2),phi_y=minabs(phi_y_1,phi_y_2);
    return abs(speed(i,j))*phi_x/sqrt(sqr(phi_x)+sqr(phi_y));
}
//#####################################################################
// Function Maxabs_H2
//#####################################################################
template<class T> T HAMILTONIAN_NORMAL_VELOCITY_2D<T>::
Maxabs_H2(const T phi_x_1,const T phi_x_2,const T phi_y_1,const T phi_y_2,const int i,const int j,const T t)
{       
    if(phi_x_1*phi_x_2 <= 0) return abs(speed(i,j)); // use phi_x=0
    T phi_x=minabs(phi_x_1,phi_x_2),phi_y=maxabs(phi_y_1,phi_y_2);
    return abs(speed(i,j))*phi_y/sqrt(sqr(phi_x)+sqr(phi_y));
}
//#####################################################################
template class HAMILTONIAN_NORMAL_VELOCITY_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class HAMILTONIAN_NORMAL_VELOCITY_2D<double>;
#endif

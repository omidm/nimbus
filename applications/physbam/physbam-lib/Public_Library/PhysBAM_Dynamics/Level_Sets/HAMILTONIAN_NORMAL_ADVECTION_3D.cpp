//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/maxabs.h>
#include <PhysBAM_Tools/Math_Tools/minabs.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Dynamics/Level_Sets/HAMILTONIAN_NORMAL_ADVECTION_3D.h>
using namespace PhysBAM;
//#####################################################################
// Function H
//#####################################################################
template<class T> T HAMILTONIAN_NORMAL_ADVECTION_3D<T>::
H(const T phi_x,const T phi_y,const T phi_z,const int i,const int j,const int k,const T t)
{       
    return speed*sqrt(sqr(phi_x)+sqr(phi_y)+sqr(phi_z));
}
//#####################################################################
// Function Maxabs_H1
//#####################################################################
template<class T> T HAMILTONIAN_NORMAL_ADVECTION_3D<T>::
Maxabs_H1(const T phi_x_1,const T phi_x_2,const T phi_y_1,const T phi_y_2,const T phi_z_1,const T phi_z_2,const int i,const int j,const int k,const T t)
{       
    T phi_y=0,phi_z=0;
    if(phi_y_1*phi_y_2 <= 0){
        if(phi_z_1*phi_z_2 <= 0) return abs(speed); // use phi_y=0 and phi_z=0
        else phi_z=minabs(phi_z_1,phi_z_2);}
    else{
        phi_y=minabs(phi_y_1,phi_y_2);
        if(phi_z_1*phi_z_2 > 0) phi_z=minabs(phi_z_1,phi_z_2);}
    T phi_x=maxabs(phi_x_1,phi_x_2);
    return abs(speed)*phi_x/sqrt(sqr(phi_x)+sqr(phi_y)+sqr(phi_z));
}
//#####################################################################
// Function Maxabs_H2
//#####################################################################
template<class T> T HAMILTONIAN_NORMAL_ADVECTION_3D<T>::
Maxabs_H2(const T phi_x_1,const T phi_x_2,const T phi_y_1,const T phi_y_2,const T phi_z_1,const T phi_z_2,const int i,const int j,const int k,const T t)
{       
    T phi_x=0,phi_z=0;
    if(phi_x_1*phi_x_2 <= 0){
        if(phi_z_1*phi_z_2 <= 0) return abs(speed); // use phi_x=0 and phi_z=0
        else phi_z=minabs(phi_z_1,phi_z_2);}
    else{
        phi_x=minabs(phi_x_1,phi_x_2);
        if(phi_z_1*phi_z_2 > 0) phi_z=minabs(phi_z_1,phi_z_2);}
    T phi_y=maxabs(phi_y_1,phi_y_2);
    return abs(speed)*phi_y/sqrt(sqr(phi_x)+sqr(phi_y)+sqr(phi_z));
}
//#####################################################################
// Function Maxabs_H3
//#####################################################################
template<class T> T HAMILTONIAN_NORMAL_ADVECTION_3D<T>::
Maxabs_H3(const T phi_x_1,const T phi_x_2,const T phi_y_1,const T phi_y_2,const T phi_z_1,const T phi_z_2,const int i,const int j,const int k,const T t)
{       
    T phi_x=0,phi_y=0;
    if(phi_x_1*phi_x_2 <= 0){
        if(phi_y_1*phi_y_2 <= 0) return abs(speed); // use phi_x=0 and phi_y=0
        else phi_y=minabs(phi_y_1,phi_y_2);}
    else{
        phi_x=minabs(phi_x_1,phi_x_2);
        if(phi_y_1*phi_y_2 > 0) phi_y=minabs(phi_y_1,phi_y_2);}
    T phi_z=maxabs(phi_z_1,phi_z_2);
    return abs(speed)*phi_z/sqrt(sqr(phi_x)+sqr(phi_y)+sqr(phi_z));
}
//#####################################################################
template class HAMILTONIAN_NORMAL_ADVECTION_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class HAMILTONIAN_NORMAL_ADVECTION_3D<double>;
#endif

//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/GRID_LAGRANGE_1D.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T> void GRID_LAGRANGE_1D<T>::
Euler_Step(const ARRAY<T,VECTOR<int,1> >& u,const T dt)
{       
    for(int i=1;i<=m;i++) x(i)+=dt*u(i); 
}
//#####################################################################
// Function Get_Lengths
//#####################################################################
// size (1,m-1)
template<class T> void GRID_LAGRANGE_1D<T>::
Get_Lengths(ARRAY<T,VECTOR<int,1> >& L)
{       
    for(int i=1;i<=m-1;i++) L(i)=abs(x(i+1)-x(i)); 
}
//#####################################################################
// Function Get_Midpoints
//#####################################################################
// size (1,m-1)
template<class T> void GRID_LAGRANGE_1D<T>::
Get_Midpoints(ARRAY<T,VECTOR<int,1> >& M)
{       
    for(int i=1;i<=m-1;i++) M(i)=(x(i)+x(i+1))/2; 
}
//#####################################################################
template class GRID_LAGRANGE_1D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_LAGRANGE_1D<double>;
#endif

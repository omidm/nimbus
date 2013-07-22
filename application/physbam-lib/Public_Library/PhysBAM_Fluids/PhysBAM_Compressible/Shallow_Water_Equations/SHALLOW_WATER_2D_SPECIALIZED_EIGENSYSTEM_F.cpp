//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Eran Guendelman, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_F  
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Shallow_Water_Equations/SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_F.h>
using namespace PhysBAM;
//#####################################################################
// Function Flux
//#####################################################################
// F(U) for i in (-2,m+3) - 3 ghost cells
template<class T> void SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_F<T>::
Flux(const int m,const ARRAY<VECTOR<T,2> ,VECTOR<int,1> >& U,ARRAY<VECTOR<T,2> ,VECTOR<int,1> >& F,ARRAY<VECTOR<T,2> ,VECTOR<int,1> >* U_clamped)       
{
    for(int i=-2;i<=m+3;i++){
        F(i)(1)=U(i)(1)*U(i)(2); // h*u
        F(i)(2)=(T).5*sqr(U(i)(2))+gravity*eta_ghost(i,slice_index.y);}  // .5*u^2+g*eta
}
//#####################################################################
// Function Eigenvalues
//#####################################################################
// eigenvalues for F(U) at flux i and and at points i and i+1
template<class T> bool SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_F<T>::
Eigenvalues(const ARRAY<VECTOR<T,2> ,VECTOR<int,1> >& U,const int i,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right)
{
    // eigenvalues on the left - at point i
    T u=U(i)(2);
    assert(U(i)(1)>=0); // assume height never becomes negative
    T celerity=sqrt(gravity*U(i)(1));
    lambda_left(1)=u-celerity;
    lambda_left(2)=u+celerity;
        
    // eigenvalues on the right - at point i+1
    u=U(i+1)(2);
    assert(U(i+1)(1)>=0); // assume height never becomes negative
    celerity=sqrt(gravity*U(i+1)(1));
    lambda_right(1)=u-celerity;
    lambda_right(2)=u+celerity;
        
    // eigenvalues in the center - at flux i
    T h=(T).5*(U(i)(1)+U(i+1)(1));
    u=(T).5*(U(i)(2)+U(i+1)(2));
    assert(h>=0); // assume height never becomes negative
    celerity=sqrt(gravity*h);
    lambda(1)=u-celerity;
    lambda(2)=u+celerity;

    return true; // eigensystem is well defined
}  
//#####################################################################
// Function Eigenvectors
//#####################################################################
// eigenvectors for F(U) at flux i between points i and i+1
template<class T> void SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_F<T>::
Eigenvectors(const ARRAY<VECTOR<T,2> ,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R)
{
    // eigensystem in the center - at flux i
    T h=(T).5*(U(i)(1)+U(i+1)(1));
    T sqrt_h_over_gravity=sqrt(h/gravity);
    T sqrt_gravity_over_h=0;
    if(h>min_height) sqrt_gravity_over_h=1/sqrt_h_over_gravity;
                    
    L(1,1)=(T)-.5;L(1,2)=(T).5*sqrt_h_over_gravity;
    L(2,1)=(T).5;L(2,2)=(T).5*sqrt_h_over_gravity;
    
    R(1,1)=-1;R(1,2)=sqrt_gravity_over_h;
    R(2,1)=1;R(2,2)=sqrt_gravity_over_h;
}  
//#####################################################################
template class SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_F<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_F<double>;
#endif

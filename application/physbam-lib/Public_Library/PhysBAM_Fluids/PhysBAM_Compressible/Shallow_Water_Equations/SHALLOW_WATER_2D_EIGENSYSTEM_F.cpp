//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_2D_EIGENSYSTEM_F  
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Shallow_Water_Equations/SHALLOW_WATER_2D_EIGENSYSTEM_F.h>
using namespace PhysBAM;
//#####################################################################
// Function Flux
//#####################################################################
// F(U) for i in (-2,m+3) - 3 ghost cells
template<class T> void SHALLOW_WATER_2D_EIGENSYSTEM_F<T>::
Flux(const int m,const ARRAY<VECTOR<T,3> ,VECTOR<int,1> >& U,ARRAY<VECTOR<T,3> ,VECTOR<int,1> >& F,ARRAY<VECTOR<T,3> ,VECTOR<int,1> >* U_clamped)       
{
    for(int i=-2;i<=m+3;i++){
        T one_over_h=1/U(i)(1),u=U(i)(2)*one_over_h;
        F(i)(1)=U(i)(2);                                            // h*u
        F(i)(2)=U(i)(2)*u+(T).5*gravity*sqr(U(i)(1));  // h*u^2+.5*g*h^2
        F(i)(3)=U(i)(2)*U(i)(3)*one_over_h;}              // h*u*v
}
//#####################################################################
// Function Maximum_Magnitude_Eigenvalue
//#####################################################################
// maximum magnitude eigenvalue for F(U) at point cell
template<class T> T SHALLOW_WATER_2D_EIGENSYSTEM_F<T>::
Maximum_Magnitude_Eigenvalue(const VECTOR<T,3>& U_cell)
{
    T u=U_cell(2)/U_cell(1);
    T celerity=0;if(U_cell(1) >= 0) celerity=sqrt(gravity*U_cell(1));
    return maxabs(u-celerity,u+celerity);
}
//#####################################################################
// Function Eigenvalues
//#####################################################################
// eigenvalues for F(U) at flux i and and at points i and i+1
template<class T> bool SHALLOW_WATER_2D_EIGENSYSTEM_F<T>::
Eigenvalues(const ARRAY<VECTOR<T,3> ,VECTOR<int,1> >& U,const int i,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right)
{
    bool weakly_hyperbolic=false;

    // eigenvalues on the left - at point i
    T u=U(i)(2)/U(i)(1);
    T celerity=0;if(U(i)(1) >= 0) celerity=sqrt(gravity*U(i)(1));else weakly_hyperbolic=true;
    lambda_left(1)=u-celerity;
    lambda_left(2)=u;
    lambda_left(3)=u+celerity;
        
    // eigenvalues on the right - at point i+1
    u=U(i+1)(2)/U(i+1)(1);
    celerity=0;if(U(i+1)(1) >= 0) celerity=sqrt(gravity*U(i+1)(1));else weakly_hyperbolic=true;
    lambda_right(1)=u-celerity;
    lambda_right(2)=u;
    lambda_right(3)=u+celerity;
        
    // eigenvalues in the center - at flux i
    T h=(T).5*(U(i)(1)+U(i+1)(1)),h_u=(T).5*(U(i)(2)+U(i+1)(2));
    u=h_u/h;
    celerity=0;if(h >= 0) celerity=sqrt(gravity*h);else weakly_hyperbolic=true;
    lambda(1)=u-celerity;
    lambda(2)=u;
    lambda(3)=u+celerity;

    if(weakly_hyperbolic) return false; // loss of hyperbolicity
    else return true; // eigensystem is well defined
}  
//#####################################################################
// Function Eigenvectors
//#####################################################################
// eigenvectors for F(U) at flux i between points i and i+1
template<class T> void SHALLOW_WATER_2D_EIGENSYSTEM_F<T>::
Eigenvectors(const ARRAY<VECTOR<T,3> ,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R)
{
    // eigensystem in the center - at flux i
    T h=(T).5*(U(i)(1)+U(i+1)(1)),h_u=(T).5*(U(i)(2)+U(i+1)(2)),h_v=(T).5*(U(i)(3)+U(i+1)(3));
    T one_over_h=1/h,u=h_u*one_over_h,v=h_v*one_over_h;
    T celerity=sqrt(gravity*h),one_over_2_celerity=1/(2*celerity);
                    
    L(1,1)=one_over_2_celerity*(u+celerity);L(1,2)=-one_over_2_celerity;L(1,3)=0;
    L(2,1)=-v;L(2,2)=0;L(2,3)=1;
    L(3,1)=-one_over_2_celerity*(u-celerity);L(3,2)=one_over_2_celerity;L(3,3)=0;
    
    R(1,1)=1;R(1,2)=u-celerity;R(1,3)=v;
    R(2,1)=0;R(2,2)=0;R(2,3)=1;
    R(3,1)=1;R(3,2)=u+celerity;R(3,3)=v;
}  
//#####################################################################
template class SHALLOW_WATER_2D_EIGENSYSTEM_F<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SHALLOW_WATER_2D_EIGENSYSTEM_F<double>;
#endif

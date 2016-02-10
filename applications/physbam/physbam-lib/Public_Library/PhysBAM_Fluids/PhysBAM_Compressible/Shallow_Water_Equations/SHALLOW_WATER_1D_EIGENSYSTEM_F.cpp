//#####################################################################
// Copyright 2003-2007, Eran Guendelman, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_1D_EIGENSYSTEM_F  
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Shallow_Water_Equations/SHALLOW_WATER_1D_EIGENSYSTEM_F.h>
using namespace PhysBAM;
//#####################################################################
// Function Flux
//#####################################################################
// F(U) for i in (-2,m+3) - 3 ghost cells
template<class T> void SHALLOW_WATER_1D_EIGENSYSTEM_F<T>::
Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)       
{
    for(int i=-2;i<=m+3;i++){
        T u=U(i)(2)/U(i)(1);
        F(i)(1)=U(i)(2);                                        // h*u
        F(i)(2)=U(i)(2)*u+(T).5*gravity*sqr(U(i)(1));} // h*u^2+.5*g*h^2
}
//#####################################################################
// Function Maximum_Magnitude_Eigenvalue
//#####################################################################
// maximum magnitude eigenvalue for F(U) at point cell
template<class T> T SHALLOW_WATER_1D_EIGENSYSTEM_F<T>::
Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell)
{
    T u=U_cell(2)/U_cell(1);
    T celerity=0;if(U_cell(1) >= 0) celerity=sqrt(gravity*U_cell(1));
    return maxabs(u-celerity,u+celerity);
}
//#####################################################################
// Function Eigenvalues
//#####################################################################
// eigenvalues for F(U) at flux i and and at points i and i+1
template<class T> bool SHALLOW_WATER_1D_EIGENSYSTEM_F<T>::
Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right)
{
    int weakly_hyperbolic=0;

    // eigenvalues on the left - at point i
    T u=U(i)(2)/U(i)(1);
    T celerity=0;if(U(i)(1) >= 0) celerity=sqrt(gravity*U(i)(1));else weakly_hyperbolic=1;
    lambda_left(1)=u-celerity;
    lambda_left(2)=u+celerity;

    // eigenvalues on the right - at point i+1
    u=U(i+1)(2)/U(i+1)(1);
    celerity=0;if(U(i+1)(1) >= 0) celerity=sqrt(gravity*U(i+1)(1));else weakly_hyperbolic=1;
    lambda_right(1)=u-celerity;
    lambda_right(2)=u+celerity;

    // eigenvalues in the center - at flux i
    T h=(U(i)(1)+U(i+1)(1))/2,h_u=(U(i)(2)+U(i+1)(2))/2;
    u=h_u/h;
    celerity=0;if(h >= 0) celerity=sqrt(gravity*h);else weakly_hyperbolic=1;
    lambda(1)=u-celerity;
    lambda(2)=u+celerity;

    if(weakly_hyperbolic) return false; // loss of hyperbolicity
    else return true; // eigensystem is well defined
}  
//#####################################################################
// Function Eigenvectors
//#####################################################################
// eigenvectors for F(U) at flux i between points i and i+1
template<class T> void SHALLOW_WATER_1D_EIGENSYSTEM_F<T>::
Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R)
{
    // eigensystem in the center - at flux i
    T h=(U(i)(1)+U(i+1)(1))/2,h_u=(U(i)(2)+U(i+1)(2))/2,u=h_u/h;
    T celerity=sqrt(gravity*h),one_over_2_celerity=1/(2*celerity);

    L(1,1)=one_over_2_celerity*(u+celerity);L(1,2)=-one_over_2_celerity;
    L(2,1)=-one_over_2_celerity*(u-celerity);L(2,2)=one_over_2_celerity;

    R(1,1)=1;R(1,2)=u-celerity;
    R(2,1)=1;R(2,2)=u+celerity;
}  
//#####################################################################
template class SHALLOW_WATER_1D_EIGENSYSTEM_F<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SHALLOW_WATER_1D_EIGENSYSTEM_F<double>;
#endif

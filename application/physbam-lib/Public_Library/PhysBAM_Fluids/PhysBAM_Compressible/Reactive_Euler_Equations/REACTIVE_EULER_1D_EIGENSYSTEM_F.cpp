//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Reactive_Euler_Equations/REACTIVE_EULER_1D_EIGENSYSTEM_F.h>
using namespace PhysBAM;
//#####################################################################
// Function Flux
//#####################################################################
// F(U) for i in (-2,m+3) - 3 ghost cells
template<class T> void REACTIVE_EULER_1D_EIGENSYSTEM_F<T>::
Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)       
{
    for(int i=-2;i<=m+3;i++){
        T u=U(i)(2)/U(i)(1);
        T Y=U(i)(4)/U(i)(1);
        T p=eos.p(U(i)(1),e(U(i)(1),U(i)(2),U(i)(3)),Y);
        F(i)(1)=U(i)(2);        // rho*u
        F(i)(2)=U(i)(2)*u+p;    // rho*u^2+p
        F(i)(3)=(U(i)(3)+p)*u; // (E+p)*u
        F(i)(4)=U(i)(4)*u;} // rho*Y*u
}
//#####################################################################
// Function Maximum_Magnitude_Eigenvalue
//#####################################################################
// maximum magnitude eigenvalue for F(U) at point cell
template<class T> T REACTIVE_EULER_1D_EIGENSYSTEM_F<T>::
Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell)
{
    T u=U_cell(2)/U_cell(1);
    T Y=U_cell(4)/U_cell(1);
    T sound_speed=eos.c(U_cell(1),e(U_cell(1),U_cell(2),U_cell(3)),Y);
    return maxabs(u-sound_speed,u+sound_speed);
}
//#####################################################################
// Function Eigenvalues
//#####################################################################
// eigenvalues for F(U) at flux i and and at points i and i+1
template<class T> bool REACTIVE_EULER_1D_EIGENSYSTEM_F<T>::
Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right)
{
    int cavitation=0;

    // eigenvalues on the left - at point i
    T u=U(i)(2)/U(i)(1);
    T Y=U(i)(4)/U(i)(1);
    T sound_speed=eos.c(U(i)(1),e(U(i)(1),U(i)(2),U(i)(3)),Y);
    if(sound_speed == 0) cavitation=1;
    lambda_left(1)=u-sound_speed;
    lambda_left(2)=u;
    lambda_left(3)=u;
    lambda_left(4)=u+sound_speed;
        
    // eigenvalues on the right - at point i+1
    u=U(i+1)(2)/U(i+1)(1);
    Y=U(i+1)(4)/U(i+1)(1);
    sound_speed=eos.c(U(i+1)(1),e(U(i+1)(1),U(i+1)(2),U(i+1)(3)),Y);
    if(sound_speed == 0) cavitation=1;
    lambda_right(1)=u-sound_speed;
    lambda_right(2)=u;
    lambda_right(3)=u;
    lambda_right(4)=u+sound_speed;
        
    // eigenvalues in the center - at flux i
    T rho=(U(i)(1)+U(i+1)(1))/2;
    T rho_u=(U(i)(2)+U(i+1)(2))/2;
    T E=(U(i)(3)+U(i+1)(3))/2;
    T rho_Y=(U(i)(4)+U(i+1)(4))/2;
    u=rho_u/rho;
    Y=rho_Y/rho;
    T internal_energy=e(rho,rho_u,E);
    sound_speed=eos.c(rho,internal_energy,Y);
    if(sound_speed == 0) cavitation=1;
    lambda(1)=u-sound_speed;
    lambda(2)=u;
    lambda(3)=u;
    lambda(4)=u+sound_speed;

    if(cavitation) return 0; // loss of hyperbolicity
    else return 1; // eigensystem is well defined
}  
//#####################################################################
// Function Eigenvectors
//#####################################################################
// eigenvectors for F(U) at flux i between points i and i+1
template<class T> void REACTIVE_EULER_1D_EIGENSYSTEM_F<T>::
Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R)
{
    // eigensystem in the center - at flux i
    T rho=(U(i)(1)+U(i+1)(1))/2;
    T rho_u=(U(i)(2)+U(i+1)(2))/2;
    T E=(U(i)(3)+U(i+1)(3))/2;
    T rho_Y=(U(i)(4)+U(i+1)(4))/2;
    T u=rho_u/rho;
    T Y=rho_Y/rho;
    T internal_energy=e(rho,rho_u,E);
    T sound_speed=eos.c(rho,internal_energy,Y);
        
    T q2=sqr(u);
    T h=(E+eos.p(rho,internal_energy,Y))/rho; // capital H=(E+p)/rho
    T b1=eos.p_e(rho,internal_energy,Y)/(rho*sqr(sound_speed));
    T b2=1+b1*(q2-h);
    T z=eos.e_not1-eos.e_not2;
    T b3=b1*Y*z;
        
    // some definitions to make the code faster
    T one_over_2c=1/(2*sound_speed);
    T u_over_2c=u*one_over_2c;
    T b1_over_2=b1/2;
    T b2_over_2=b2/2;
    T b3_over_2=b3/2;
    T b1_over_2_times_u=b1_over_2*u;
    T u_times_c=u*sound_speed;
                    
    L(1,1)=b2_over_2+u_over_2c+b3_over_2;
    L(1,2)=-b1_over_2_times_u-one_over_2c;
    L(1,3)=b1_over_2;
    L(1,4)=-b1_over_2*z;
    L(2,1)=1-b2-b3;
    L(2,2)=b1*u;
    L(2,3)=-b1;
    L(2,4)=b1*z; 
    L(3,1)=-Y;
    L(3,2)=0;
    L(3,3)=0;
    L(3,4)=1;
    L(4,1)=b2_over_2-u_over_2c+b3_over_2;
    L(4,2)=-b1_over_2_times_u+one_over_2c;
    L(4,3)=b1_over_2;
    L(4,4)=-b1_over_2*z;
    
    R(1,1)=1;
    R(1,2)=u-sound_speed;
    R(1,3)=h-u_times_c;
    R(1,4)=Y;
    R(2,1)=1;
    R(2,2)=u;
    R(2,3)=h-1/b1;
    R(2,4)=Y;
    R(3,1)=0;
    R(3,2)=0;
    R(3,3)=z;
    R(3,4)=1;
    R(4,1)=1;
    R(4,2)=u+sound_speed;
    R(4,3)=h+u_times_c;
    R(4,4)=Y;
}  
//#####################################################################
template class REACTIVE_EULER_1D_EIGENSYSTEM_F<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class REACTIVE_EULER_1D_EIGENSYSTEM_F<double>;
#endif

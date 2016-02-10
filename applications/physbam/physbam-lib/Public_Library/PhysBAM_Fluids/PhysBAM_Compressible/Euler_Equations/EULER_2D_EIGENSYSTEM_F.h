//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_2D_EIGENSYSTEM_F  
//##################################################################### 
#ifndef __EULER_2D_EIGENSYSTEM_F__
#define __EULER_2D_EIGENSYSTEM_F__   

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_EIGENSYSTEM.h>
namespace PhysBAM{

template<class T_input>
class EULER_2D_EIGENSYSTEM_F:public EULER_EIGENSYSTEM<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;typedef VECTOR<T,4> TV_DIMENSION;
    enum WORKAROUND1 {d=TV_DIMENSION::m};
    typedef EULER_EIGENSYSTEM<GRID<VECTOR<T_input,2> > > BASE;
public:
    using BASE::eos;

    bool only_pressure_flux;

    EULER_2D_EIGENSYSTEM_F(bool only_pressure_flux_input=false):only_pressure_flux(only_pressure_flux_input)
    {}

//#####################################################################
    void Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) PHYSBAM_OVERRIDE;
    void Flux_Using_Face_Velocity(VECTOR<int,2> range,const int face_index,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,const bool use_standard_average,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) PHYSBAM_OVERRIDE;
    T Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell) PHYSBAM_OVERRIDE;
    bool Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right) PHYSBAM_OVERRIDE; 
    void Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R) PHYSBAM_OVERRIDE; 
//#####################################################################
}; 
//#####################################################################
// Function Flux
//#####################################################################
// F(U) for i in (-2,m+3) - 3 ghost cells
template<class T> void EULER_2D_EIGENSYSTEM_F<T>::
Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)       
{
    if(only_pressure_flux){
        for(int i=-2;i<=m+3;i++){
            T p=eos->p(U(i)(1),e(U(i)(1),U(i)(2),U(i)(3),U(i)(4)));
            F(i)(1)=0;
            F(i)(2)=p;
            F(i)(3)=0;
            F(i)(4)=0;}
        return;}

    if(U_clamped){
        for(int i=-2;i<=m+3;i++){
            T u=U(i)(2)/U(i)(1),v=U(i)(3)/U(i)(1);
            T p=eos->p(U(i)(1),e(U(i)(1),U(i)(2),U(i)(3),U(i)(4)));
            F(i)(1)=(*U_clamped)(i)(2);         // rho_clamped*u
            F(i)(2)=(*U_clamped)(i)(2)*u+p;     // rho_clamped*u^2+p
            F(i)(3)=(*U_clamped)(i)(2)*v;       // rho_clamped*u*v
            F(i)(4)=((*U_clamped)(i)(4)+p)*u;}} // (E_from_rho_clamped+p)*u
    else{
        for(int i=-2;i<=m+3;i++){
            T u=U(i)(2)/U(i)(1),v=U(i)(3)/U(i)(1);
            T p=eos->p(U(i)(1),e(U(i)(1),U(i)(2),U(i)(3),U(i)(4)));
            F(i)(1)=U(i)(2);         // rho*u
            F(i)(2)=U(i)(2)*u+p;     // rho*u^2+p
            F(i)(3)=U(i)(2)*v;       // rho*u*v
            F(i)(4)=(U(i)(4)+p)*u;}} // (E+p)*u
}
//#####################################################################
// Function Flux_Using_Face_Velocity
//#####################################################################
template<class T> void EULER_2D_EIGENSYSTEM_F<T>::
Flux_Using_Face_Velocity(VECTOR<int,2> range,const int face_index,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,const bool use_standard_average,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)
{
    if(only_pressure_flux){
        for(int i=range.x;i<=range.y;i++){
            T p=eos->p(U(i)(1),e(U(i)(1),U(i)(2),U(i)(3),U(i)(4)));
            F(i)(1)=0;
            F(i)(2)=p;
            F(i)(3)=0;
            F(i)(4)=0;}
        return;}

    T average_u;
    if(use_standard_average) average_u=(U(face_index)(2)/U(face_index)(1)+U(face_index+1)(2)/U(face_index+1)(1))*(T).5;
    else average_u=(U(face_index)(2)+U(face_index+1)(2))/(U(face_index)(1)+U(face_index+1)(1));

    if(U_clamped){
        for(int i=range.x;i<=range.y;i++){
            T p=eos->p(U(i)(1),e(U(i)(1),U(i)(2),U(i)(3),U(i)(4)));
            F(i)(1)=(*U_clamped)(i)(1)*average_u;
            F(i)(2)=(*U_clamped)(i)(2)*average_u+p;
            F(i)(3)=(*U_clamped)(i)(3)*average_u;
            F(i)(4)=((*U_clamped)(i)(4)+p)*average_u;}}
    else{
        for(int i=range.x;i<=range.y;i++){
            T p=eos->p(U(i)(1),e(U(i)(1),U(i)(2),U(i)(3),U(i)(4)));
            F(i)(1)=U(i)(1)*average_u;
            F(i)(2)=U(i)(2)*average_u+p;
            F(i)(3)=U(i)(3)*average_u;
            F(i)(4)=(U(i)(4)+p)*average_u;}}
}
//#####################################################################
// Function Maximum_Magnitude_Eigenvalue
//#####################################################################
// maximum magnitude eigenvalue for F(U) at point cell
template<class T> T EULER_2D_EIGENSYSTEM_F<T>::
Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell)
{
    T u=U_cell(2)/U_cell(1);
    T sound_speed=eos->c(U_cell(1),e(U_cell(1),U_cell(2),U_cell(3),U_cell(4)));
    return maxabs(u-sound_speed,u+sound_speed);
}
//#####################################################################
// Function Eigenvalues
//#####################################################################
// eigenvalues for F(U) at flux i and and at points i and i+1
template<class T> bool EULER_2D_EIGENSYSTEM_F<T>::
Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right)
{
    bool cavitation=false;

    // eigenvalues on the left - at point i
    T u=U(i)(2)/U(i)(1);
    T sound_speed=eos->c(U(i)(1),e(U(i)(1),U(i)(2),U(i)(3),U(i)(4)));
    if(sound_speed == 0) cavitation=true;
    lambda_left(1)=u-sound_speed;
    lambda_left(2)=lambda_left(3)=u;
    lambda_left(4)=u+sound_speed;
        
    // eigenvalues on the right - at point i+1
    u=U(i+1)(2)/U(i+1)(1);
    sound_speed=eos->c(U(i+1)(1),e(U(i+1)(1),U(i+1)(2),U(i+1)(3),U(i+1)(4)));
    if(sound_speed == 0) cavitation=true;
    lambda_right(1)=u-sound_speed;
    lambda_right(2)=lambda_right(3)=u;
    lambda_right(4)=u+sound_speed;
        
    // eigenvalues in the center - at flux i
    T rho=(U(i)(1)+U(i+1)(1))/2;
    T rho_u=(U(i)(2)+U(i+1)(2))/2;
    T rho_v=(U(i)(3)+U(i+1)(3))/2;
    T E=(U(i)(4)+U(i+1)(4))/2;
    u=rho_u/rho;
    T internal_energy=e(rho,rho_u,rho_v,E);
    sound_speed=eos->c(rho,internal_energy);
    if(sound_speed == 0) cavitation=true;
    lambda(1)=u-sound_speed;
    lambda(2)=lambda(3)=u;
    lambda(4)=u+sound_speed;

    return (!cavitation); //cavitation --> loss of hyperbolicity, else well defined eigensystem
}  
//#####################################################################
// Function Eigenvectors
//#####################################################################
// eigenvectors for F(U) at flux i between points i and i+1
template<class T> void EULER_2D_EIGENSYSTEM_F<T>::
Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R)
{
    // eigensystem in the center - at flux i
    T rho=(U(i)(1)+U(i+1)(1))/2;
    T rho_u=(U(i)(2)+U(i+1)(2))/2;
    T rho_v=(U(i)(3)+U(i+1)(3))/2;
    T E=(U(i)(4)+U(i+1)(4))/2;
    T u=rho_u/rho;
    T internal_energy=e(rho,rho_u,rho_v,E);
    T sound_speed=eos->c(rho,internal_energy);
    T p=eos->p(rho,internal_energy);
        
    T v=rho_v/rho;
    T q2=sqr(u)+sqr(v);
    T h=(E+p)/rho;
    T b1=eos->p_e(rho,internal_energy)/(rho*sqr(sound_speed));
    T b2=1+b1*(q2-h);
        
    // some definitions to make the code faster
    T one_over_2c=1/(2*sound_speed);
    T u_over_2c=u*one_over_2c;
    T b1_over_2=b1/2;
    T b2_over_2=b2/2;
    T b1_over_2_times_u=b1_over_2*u;
    T b1_over_2_times_v=b1_over_2*v;
    T u_times_c=u*sound_speed;
                    
    L(1,1)=b2_over_2+u_over_2c;
    L(1,2)=-b1_over_2_times_u-one_over_2c;
    L(1,3)=-b1_over_2_times_v;
    L(1,4)=b1_over_2;
    L(2,1)=h-q2;
    L(2,2)=u;
    L(2,3)=v;
    L(2,4)=-1;
    L(3,1)=v;
    L(3,2)=0;
    L(3,3)=-1;
    L(3,4)=0;
    L(4,1)=b2_over_2-u_over_2c;
    L(4,2)=-b1_over_2_times_u+one_over_2c;
    L(4,3)=-b1_over_2_times_v;
    L(4,4)=b1_over_2;
    
    R(1,1)=1;
    R(1,2)=u-sound_speed;
    R(1,3)=v;
    R(1,4)=h-u_times_c;
    R(2,1)=b1;
    R(2,2)=b1*u;
    R(2,3)=b1*v;
    R(2,4)=b1*h-1;
    R(3,1)=0;
    R(3,2)=0;
    R(3,3)=-1;
    R(3,4)=-v;
    R(4,1)=1;
    R(4,2)=u+sound_speed;
    R(4,3)=v;
    R(4,4)=h+u_times_c;
}  
//#####################################################################
}  
#endif


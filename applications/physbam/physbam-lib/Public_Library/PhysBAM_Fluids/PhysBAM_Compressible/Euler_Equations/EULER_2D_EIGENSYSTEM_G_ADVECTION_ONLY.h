//#####################################################################
// Copyright 2007-2009, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_2D_EIGENSYSTEM_G_ADVECTION_ONLY
//#####################################################################
#ifndef __EULER_2D_EIGENSYSTEM_G_ADVECTION_ONLY__
#define __EULER_2D_EIGENSYSTEM_G_ADVECTION_ONLY__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_EIGENSYSTEM.h>
namespace PhysBAM{

template<class T_input>
class EULER_2D_EIGENSYSTEM_G_ADVECTION_ONLY:public EULER_EIGENSYSTEM<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;typedef VECTOR<T,4> TV_DIMENSION;
    enum WORKAROUND1 {d=TV_DIMENSION::m};
public:
    EULER_2D_EIGENSYSTEM_G_ADVECTION_ONLY()
    {}

    bool All_Eigenvalues_Same() PHYSBAM_OVERRIDE {return true;}

//#####################################################################
    void Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& G,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) PHYSBAM_OVERRIDE;
    void Flux_Divided_By_Velocity(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& G,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) PHYSBAM_OVERRIDE;
    T Get_Face_Velocity_Component(const int face_index,const bool use_standard_average,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U);
    void Flux_Using_Face_Velocity(VECTOR<int,2> range,const int face_index,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& G,const bool use_standard_average,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) PHYSBAM_OVERRIDE;
    T Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell) PHYSBAM_OVERRIDE;
    bool Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right) PHYSBAM_OVERRIDE;
    void Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R) PHYSBAM_OVERRIDE {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
};
//#####################################################################
// Function Flux
//#####################################################################
// G(U) for j in (-2,n+3) - 3 ghost cells
template<class T> void EULER_2D_EIGENSYSTEM_G_ADVECTION_ONLY<T>::
Flux(const int n,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& G,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)
{
    if(U_clamped){
        for(int j=-2;j<=n+3;j++){
            T u=U(j)(2)/U(j)(1),v=U(j)(3)/U(j)(1);
            G(j)(1)=(*U_clamped)(j)(3);       // rho_clamped*v
            G(j)(2)=(*U_clamped)(j)(3)*u;     // rho_clamped*u*v
            G(j)(3)=(*U_clamped)(j)(3)*v;     // rho_clamped*v^2
            G(j)(4)=(*U_clamped)(j)(4)*v;}}   // E_from_rho_clamped*v
    else{
        for(int j=-2;j<=n+3;j++){
            T u=U(j)(2)/U(j)(1),v=U(j)(3)/U(j)(1);
            G(j)(1)=U(j)(3);       // rho*v
            G(j)(2)=U(j)(3)*u;     // rho*u*v
            G(j)(3)=U(j)(3)*v;     // rho*v^2
            G(j)(4)=U(j)(4)*v;}}   // E*v
}
//#####################################################################
// Function Flux_Divided_By_Velocity
//#####################################################################
// G(U) for j in (-2,n+3) - 3 ghost cells
template<class T> void EULER_2D_EIGENSYSTEM_G_ADVECTION_ONLY<T>::
Flux_Divided_By_Velocity(const int n,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& G,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)
{
    G=(U_clamped)?*U_clamped:U;
}
//#####################################################################
// Function Get_Face_Velocity_Component
//#####################################################################
template<class T> T EULER_2D_EIGENSYSTEM_G_ADVECTION_ONLY<T>::
Get_Face_Velocity_Component(const int face_index,const bool use_standard_average,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U)
{
    T average_v;
    if(use_standard_average) average_v=(U(face_index)(3)/U(face_index)(1)+U(face_index+1)(3)/U(face_index+1)(1))*(T).5;
    else average_v=(U(face_index)(3)+U(face_index+1)(3))/(U(face_index)(1)+U(face_index+1)(1));
    return average_v;
}
//#####################################################################
// Function Flux_Using_Face_Velocity
//#####################################################################
template<class T> void EULER_2D_EIGENSYSTEM_G_ADVECTION_ONLY<T>::
Flux_Using_Face_Velocity(VECTOR<int,2> range,const int face_index,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& G,const bool use_standard_average,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)
{
    T average_v;
    if(use_standard_average) average_v=(U(face_index)(3)/U(face_index)(1)+U(face_index+1)(3)/U(face_index+1)(1))*(T).5;
    else average_v=(U(face_index)(3)+U(face_index+1)(3))/(U(face_index)(1)+U(face_index+1)(1));

    if(U_clamped){
        for(int i=range.x;i<=range.y;i++){
            G(i)(1)=(*U_clamped)(i)(1)*average_v;
            G(i)(2)=(*U_clamped)(i)(2)*average_v;
            G(i)(3)=(*U_clamped)(i)(3)*average_v;
            G(i)(4)=(*U_clamped)(i)(4)*average_v;}}
    else{
        for(int i=range.x;i<=range.y;i++){
            G(i)(1)=U(i)(1)*average_v;
            G(i)(2)=U(i)(2)*average_v;
            G(i)(3)=U(i)(3)*average_v;
            G(i)(4)=U(i)(4)*average_v;}}
}
//#####################################################################
// Function Maximum_Magnitude_Eigenvalue
//#####################################################################
// maximum magnitude eigenvalue for G(U) at point cell
template<class T> T EULER_2D_EIGENSYSTEM_G_ADVECTION_ONLY<T>::
Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell)
{
    return abs(U_cell(3)/U_cell(1));
}
//#####################################################################
// Function Eigenvalues
//#####################################################################
// eigenvalues for G(U) at flux j and and at points j and j+1
template<class T> bool EULER_2D_EIGENSYSTEM_G_ADVECTION_ONLY<T>::
Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int j,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right)
{
    bool cavitation=false;

    // eigenvalues on the left - at point j
    T v=U(j)(3)/U(j)(1);
    lambda_left(1)=lambda_left(2)=lambda_left(3)=lambda_left(4)=v;

    // eigenvalues on the right - at point j+1
    v=U(j+1)(3)/U(j+1)(1);
    lambda_right(1)=lambda_right(2)=lambda_right(3)=lambda_right(4)=v;

    // eigenvalues in the center - at flux j
    T rho=(U(j)(1)+U(j+1)(1))/2;
    T rho_v=(U(j)(3)+U(j+1)(3))/2;
    v=rho_v/rho;
    lambda(1)=lambda(2)=lambda(3)=lambda(4)=v;

    return (!cavitation); //cavitation --> loss of hyperbolicity, else well defined eigensystem
}
//#####################################################################
}
#endif

//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/POLYGON.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/GRID_LAGRANGE_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T> void GRID_LAGRANGE_2D<T>::
Euler_Step(const ARRAY<T,VECTOR<int,2> >& u,const ARRAY<T,VECTOR<int,2> >& v,const T dt)
{       
    for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){x(i,j)+=dt*u(i,j);y(i,j)+=dt*v(i,j);}
}
//#####################################################################
// Function Get_Lengths
//#####################################################################
// L1 is (1,m-1) by (1,n) while L2 is (1,m) by (1,n-1)
template<class T> void GRID_LAGRANGE_2D<T>::
Get_Lengths(ARRAY<T,VECTOR<int,2> >& L1,ARRAY<T,VECTOR<int,2> >& L2)
{
    int i,j;
    for(i=1;i<=m-1;i++) for(j=1;j<=n;j++){VECTOR<T,3> v1(x(i,j),y(i,j),0),v2(x(i+1,j),y(i+1,j),0),v3=v2-v1;L1(i,j)=v3.Magnitude();}
    for(i=1;i<=m;i++) for(j=1;j<=n-1;j++){VECTOR<T,3> v1(x(i,j),y(i,j),0),v2(x(i,j+1),y(i,j+1),0),v3=v2-v1;L2(i,j)=v3.Magnitude();}
}
//#####################################################################
// Function Get_Areas
//#####################################################################
// A is (1,m-1) by (1,n-1)
// points listed counter clockwise
template<class T> void GRID_LAGRANGE_2D<T>::
Get_Areas(ARRAY<T,VECTOR<int,2> >& A)
{
    POLYGON<VECTOR<T,2> > quadrilateral(4);

    for(int i=1;i<=m-1;i++) for(int j=1;j<=n-1;j++){
        quadrilateral.X(1)=VECTOR<T,2>(x(i,j),y(i,j));
        quadrilateral.X(2)=VECTOR<T,2>(x(i+1,j),y(i+1,j));
        quadrilateral.X(3)=VECTOR<T,2>(x(i+1,j+1),y(i+1,j+1));
        quadrilateral.X(4)=VECTOR<T,2>(x(i,j+1),y(i,j+1));
        A(i,j)=quadrilateral.Area();}  
}
//#####################################################################
// Function Get_Normals
//#####################################################################
// located at the midpoints, defined as (-y,x)/length or (y,-x)/length
// N1 points up as (1,m-1) by (1,n), N2 points right as (1,m) by (1,n-1)
template<class T> void GRID_LAGRANGE_2D<T>::
Get_Normals(ARRAY<T,VECTOR<int,2> >& N1_x,ARRAY<T,VECTOR<int,2> >& N1_y,ARRAY<T,VECTOR<int,2> >& N2_x,ARRAY<T,VECTOR<int,2> >& N2_y)
{
    int i,j;
    ARRAY<T,VECTOR<int,2> > L1(1,m-1,1,n),L2(1,m,1,n-1);Get_Lengths(L1,L2);
    for(i=1;i<=m-1;i++) for(j=1;j<=n;j++){N1_x(i,j)=-(y(i+1,j)-y(i,j))/L1(i,j);N1_y(i,j)=(x(i+1,j)-x(i,j))/L1(1,j);}
    for(i=1;i<=m;i++) for(j=1;j<=n-1;j++){N2_x(i,j)=(y(i,j+1)-y(i,j))/L2(i,j);N2_y(i,j)=-(x(i,j+1)-x(i,j))/L2(i,j);}
}
//#####################################################################
// Function Get_Centers
//#####################################################################
// C is (1,m-1) by (1,n-1)
template<class T> void GRID_LAGRANGE_2D<T>::
Get_Centers(ARRAY<T,VECTOR<int,2> >& C_x,ARRAY<T,VECTOR<int,2> >& C_y)
{
    for(int i=1;i<=m-1;i++) for(int j=1;j<=n-1;j++){
        C_x(i,j)=(x(i,j)+x(i+1,j)+x(i,j+1)+x(i+1,j+1))/4;C_y(i,j)=(y(i,j)+y(i+1,j)+y(i,j+1)+y(i+1,j+1))/4;}
}
//#####################################################################
// Function Get_Midpoints
//#####################################################################
// M1 is (1,m-1) by (1,n) and M2 is (1,m) by (1,n-1)
template<class T> void GRID_LAGRANGE_2D<T>::
Get_Midpoints(ARRAY<T,VECTOR<int,2> >& M1_x,ARRAY<T,VECTOR<int,2> >& M1_y,ARRAY<T,VECTOR<int,2> >& M2_x,ARRAY<T,VECTOR<int,2> >& M2_y)
{
    int i,j;
    for(i=1;i<=m-1;i++) for(j=1;j<=n;j++){M1_x(i,j)=(x(i,j)+x(i+1,j))/2;M1_y(i,j)=(y(i,j)+y(i+1,j))/2;}
    for(i=1;i<=m;i++) for(j=1;j<=n-1;j++){M2_x(i,j)=(x(i,j)+x(i,j+1))/2;M2_y(i,j)=(y(i,j)+y(i,j+1))/2;}
}
//#####################################################################
// Function Get_Sub_Zone_Lengths
//#####################################################################
// LL1 left, LL2 right, LL3 bottom, LL4 top: all (1,m-1) by (1,n-1)
template<class T> void GRID_LAGRANGE_2D<T>::
Get_Sub_Zone_Lengths(ARRAY<T,VECTOR<int,2> >& LL1,ARRAY<T,VECTOR<int,2> >& LL2,ARRAY<T,VECTOR<int,2> >& LL3,ARRAY<T,VECTOR<int,2> >& LL4)
{
    ARRAY<T,VECTOR<int,2> > C_x(1,m-1,1,n-1),C_y(1,m-1,1,n-1);Get_Centers(C_x,C_y);
    ARRAY<T,VECTOR<int,2> > M1_x(1,m-1,1,n),M1_y(1,m-1,1,n),M2_x(1,m,1,n-1),M2_y(1,m,1,n-1);Get_Midpoints(M1_x,M1_y,M2_x,M2_y);  
    
    for(int i=1;i<=m-1;i++) for(int j=1;j<=n-1;j++){
        VECTOR<T,3> v1(M2_x(i,j),M2_y(i,j),0),v2(C_x(i,j),C_y(i,j),0),v3=v2-v1;LL1(i,j)=v3.Magnitude();
        VECTOR<T,3> v4(C_x(i,j),C_y(i,j),0),v5(M2_x(i+1,j),M2_y(i+1,j),0),v6=v5-v4;LL2(i,j)=v6.Magnitude();
        VECTOR<T,3> v7(M1_x(i,j),M1_y(i,j),0),v8(C_x(i,j),C_y(i,j),0),v9=v8-v7;LL3(i,j)=v9.Magnitude();
        VECTOR<T,3> v10(C_x(i,j),C_y(i,j),0),v11(M1_x(i,j+1),M1_y(i,j+1),0),v12=v11-v10;LL4(i,j)=v12.Magnitude();}
}
//#####################################################################
// Function Get_Sub_Zone_Areas
//#####################################################################
// AA1 lower left, AA2 lower right, AA3 upper left, AA4 upper right: all (1,m-1) by (1,n-1)
// points listed counter clockwise
template<class T> void GRID_LAGRANGE_2D<T>::
Get_Sub_Zone_Areas(ARRAY<T,VECTOR<int,2> >& AA1,ARRAY<T,VECTOR<int,2> >& AA2,ARRAY<T,VECTOR<int,2> >& AA3,ARRAY<T,VECTOR<int,2> >& AA4)
{
    POLYGON<VECTOR<T,2> > quadrilateral(4);
    
    ARRAY<T,VECTOR<int,2> > C_x(1,m-1,1,n-1),C_y(1,m-1,1,n-1);Get_Centers(C_x,C_y);
    ARRAY<T,VECTOR<int,2> > M1_x(1,m-1,1,n),M1_y(1,m-1,1,n),M2_x(1,m,1,n-1),M2_y(1,m,1,n-1);Get_Midpoints(M1_x,M1_y,M2_x,M2_y);
    for(int i=1;i<=m-1;i++) for(int j=1;j<=n-1;j++){
        quadrilateral.X(1)=VECTOR<T,2>(x(i,j),y(i,j));
        quadrilateral.X(2)=VECTOR<T,2>(M1_x(i,j),M1_y(i,j));
        quadrilateral.X(3)=VECTOR<T,2>(C_x(i,j),C_y(i,j));
        quadrilateral.X(4)=VECTOR<T,2>(M2_x(i,j),M2_y(i,j));
        AA1(i,j)=quadrilateral.Area();
        quadrilateral.X(1)=VECTOR<T,2>(M1_x(i,j),M1_y(i,j));
        quadrilateral.X(2)=VECTOR<T,2>(x(i+1,j),y(i+1,j));
        quadrilateral.X(3)=VECTOR<T,2>(M2_x(i+1,j),M2_y(i+1,j));
        quadrilateral.X(4)=VECTOR<T,2>(C_x(i,j),C_y(i,j));
        AA2(i,j)=quadrilateral.Area();
        quadrilateral.X(1)=VECTOR<T,2>(M2_x(i,j),M2_y(i,j));
        quadrilateral.X(2)=VECTOR<T,2>(C_x(i,j),C_y(i,j));
        quadrilateral.X(3)=VECTOR<T,2>(M1_x(i,j+1),M1_y(i,j+1));
        quadrilateral.X(4)=VECTOR<T,2>(x(i,j+1),y(i,j+1));
        AA3(i,j)=quadrilateral.Area();
        quadrilateral.X(1)=VECTOR<T,2>(C_x(i,j),C_y(i,j));
        quadrilateral.X(2)=VECTOR<T,2>(M2_x(i+1,j),M2_y(i+1,j));
        quadrilateral.X(3)=VECTOR<T,2>(x(i+1,j+1),y(i+1,j+1));
        quadrilateral.X(4)=VECTOR<T,2>(M1_x(i,j+1),M1_y(i,j+1));
        AA4(i,j)=quadrilateral.Area();}  
}
//#####################################################################
// Function Get_Sub_Zone_Normals
//#####################################################################
// NN1 left, NN2 right, NN3 bottom, NN4 top: all (1,m-1) by (1,n-1)
// NN1 and NN2 point up, NN3 and NN4 point right
template<class T> void GRID_LAGRANGE_2D<T>::
Get_Sub_Zone_Normals(ARRAY<T,VECTOR<int,2> >& NN1_x,ARRAY<T,VECTOR<int,2> >& NN1_y,ARRAY<T,VECTOR<int,2> >& NN2_x,ARRAY<T,VECTOR<int,2> >& NN2_y,ARRAY<T,VECTOR<int,2> >& NN3_x,ARRAY<T,VECTOR<int,2> >& NN3_y,ARRAY<T,VECTOR<int,2> >& NN4_x,ARRAY<T,VECTOR<int,2> >& NN4_y)
{
    ARRAY<T,VECTOR<int,2> > C_x(1,m-1,1,n-1),C_y(1,m-1,1,n-1);Get_Centers(C_x,C_y);
    ARRAY<T,VECTOR<int,2> > M1_x(1,m-1,1,n),M1_y(1,m-1,1,n),M2_x(1,m,1,n-1),M2_y(1,m,1,n-1);Get_Midpoints(M1_x,M1_y,M2_x,M2_y);
    ARRAY<T,VECTOR<int,2> > L1(1,m-1,1,n-1),L2(1,m-1,1,n-1),L3(1,m-1,1,n-1),L4(1,m-1,1,n-1);Get_Sub_Zone_Lengths(L1,L2,L3,L4);
    for(int i=1;i<=m-1;i++) for(int j=1;j<=n-1;j++){
        NN1_x(i,j)=-(C_y(i,j)-M2_y(i,j))/L1(i,j);NN1_y(i,j)=(C_x(i,j)-M2_x(i,j))/L1(i,j);
        NN2_x(i,j)=-(M2_y(i+1,j)-C_y(i,j))/L2(i,j);NN2_y(i,j)=(M2_x(i+1,j)-C_x(i,j))/L2(i,j);
        NN3_x(i,j)=(C_y(i,j)-M1_y(i,j))/L3(i,j);NN3_y(i,j)=-(C_x(i,j)-M1_x(i,j))/L3(i,j);
        NN4_x(i,j)=(M1_y(i,j)-C_y(i,j))/L4(i,j);NN4_y(i,j)=-(M1_x(i,j)-C_x(i,j))/L4(i,j);}    
}
//#####################################################################
template class GRID_LAGRANGE_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_LAGRANGE_2D<double>;
#endif

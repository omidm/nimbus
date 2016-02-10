//#####################################################################
// Copyright 2005-2007, Kevin Der, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SIMPLEX_INTERACTIONS
//##################################################################### 
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/SIMPLEX_INTERACTIONS.h>
using namespace PhysBAM;
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool SIMPLEX_INTERACTIONS<T>::
Intersection(const VECTOR<VECTOR<T,2>,2>& segment_1,const VECTOR<VECTOR<T,2>,2>& segment_2)
{
    T points[4][2],m2p[16],d3p[16];
    int plus,minus;
    unsigned int mask;

    points[0][0]=segment_1[1].x;points[0][1]=segment_1[1].y;points[1][0]=segment_1[2].x;points[1][1]=segment_1[2].y;
    points[2][0]=segment_2[1].x;points[2][1]=segment_2[1].y;points[3][0]=segment_2[2].x;points[3][1]=segment_2[2].y;
    for(unsigned int i=0;i<=2;i++) for(unsigned int j=i+1;j<=3;j++) m2p[(1<<i)|(1<<j)]=points[i][0]*points[j][1]-points[i][1]*points[j][0];
    // Compute signed areas for all triangles formed by three input points
    for(unsigned int i=0;i<=1;i++) for(unsigned int j=i+1;j<=2;j++) for(unsigned int k=j+1;k<=3;k++) 
        d3p[(1<<i)|(1<<j)|(1<<k)]=m2p[(1<<j)|(1<<k)]-m2p[(1<<i)|(1<<k)]+m2p[(1<<i)|(1<<j)]; 
    // Separator line must contain a vertex from each segment
    for(unsigned int i=0;i<=1;i++) for(unsigned int j=2;j<=3;j++){
        plus=minus=0;mask=(1<<i)|(1<<j);
        for(unsigned int k=0;k<i;k++)    if(d3p[(1<<k)|mask]>(T)0) plus++;  else if(d3p[(1<<k)|mask]<(T)0) minus++;
        for(unsigned int k=i+1;k<=1;k++) if(d3p[(1<<k)|mask]>(T)0) minus++; else if(d3p[(1<<k)|mask]<(T)0) plus++;
        for(unsigned int k=2;k<j;k++)    if(d3p[(1<<k)|mask]>(T)0) plus++;  else if(d3p[(1<<k)|mask]<(T)0) minus++;
        for(unsigned int k=j+1;k<=3;k++) if(d3p[(1<<k)|mask]>(T)0) minus++; else if(d3p[(1<<k)|mask]<(T)0) plus++;
        if (plus==0||minus==0) return false;}
    // No separator line found, objects are intersecting
    return true;
}
//#####################################################################
// Function Two_Segment_Intersection_Barycentric_Coordinates
//#####################################################################
template<class T> bool SIMPLEX_INTERACTIONS<T>::
Two_Segment_Intersection_Barycentric_Coordinates(const VECTOR<VECTOR<T,2>,2>& segment1,const VECTOR<VECTOR<T,2>,2>& segment2,VECTOR<T,1>& weights1,VECTOR<T,1>& weights2)
{
    MATRIX<T,2> A(segment1[2]-segment1[1],segment2[1]-segment2[2]);if(A.Determinant()==0) return false;
    VECTOR<T,2> weights=A.Solve_Linear_System(segment2[1]-segment1[1]);
    weights1.x=weights[1];weights2.x=weights[2];
    return weights1.x>=0 && weights2.x>=0 && weights1.x<=1 && weights2.x<=1;
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool SIMPLEX_INTERACTIONS<T>::
Intersection(const VECTOR<VECTOR<T,2>,3>& triangle,const VECTOR<T,2>& point)
{
    T points[4][2],m2p[16],d3p[16];
    int plus,minus;
    unsigned int mask;

    points[0][0]=triangle[1].x;points[0][1]=triangle[1].y;points[1][0]=triangle[2].x;points[1][1]=triangle[2].y;points[2][0]=triangle[3].x;points[2][1]=triangle[3].y;
    points[3][0]=point.x;points[3][1]=point.y;
    for(unsigned int i=0;i<=2;i++) for(unsigned int j=i+1;j<=3;j++) m2p[(1<<i)|(1<<j)]=points[i][0]*points[j][1]-points[i][1]*points[j][0];
    // Compute signed areas for all triangles formed by three input points
    for(unsigned int i=0;i<=1;i++) for(unsigned int j=i+1;j<=2;j++) for(unsigned int k=j+1;k<=3;k++) 
        d3p[(1<<i)|(1<<j)|(1<<k)]=m2p[(1<<j)|(1<<k)]-m2p[(1<<i)|(1<<k)]+m2p[(1<<i)|(1<<j)];
    // Separator line must contain the lone point
    for(unsigned int i=0;i<=2;i++){
        plus=minus=0;mask=(1<<i|1<<3);
        for(unsigned int k=0;k<i;k++)    if(d3p[(1<<k)|mask]>(T)0) plus++;  else if(d3p[(1<<k)|mask]<(T)0) minus++;
        for(unsigned int k=i+1;k<=2;k++) if(d3p[(1<<k)|mask]>(T)0) minus++; else if(d3p[(1<<k)|mask]<(T)0) plus++;
        if (plus==0||minus==0) return false;}
    // No separator line found, objects are intersecting
    return true;
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool SIMPLEX_INTERACTIONS<T>::
Intersection(const VECTOR<VECTOR<T,2>,3>& triangle,const VECTOR<VECTOR<T,2>,2>& segment)
{
    T points[5][2],m2p[64],d3p[64];
    int plus,minus;
    unsigned int mask;

    points[0][0]=triangle[1].x;points[0][1]=triangle[1].y;points[1][0]=triangle[2].x;points[1][1]=triangle[2].y;points[2][0]=triangle[3].x;points[2][1]=triangle[3].y;
    points[3][0]=segment[1].x;points[3][1]=segment[1].y;points[4][0]=segment[2].x;points[4][1]=segment[2].y;
    for(unsigned int i=0;i<=3;i++) for(unsigned int j=i+1;j<=4;j++) m2p[(1<<i)|(1<<j)]=points[i][0]*points[j][1]-points[i][1]*points[j][0];
    // Compute signed areas for all triangles formed by three input points
    for(unsigned int i=0;i<=2;i++) for(unsigned int j=i+1;j<=3;j++) for(unsigned int k=j+1;k<=4;k++) 
        d3p[(1<<i)|(1<<j)|(1<<k)]=m2p[(1<<j)|(1<<k)]-m2p[(1<<i)|(1<<k)]+m2p[(1<<i)|(1<<j)]; 
    // Separator line must contain a vertex from the triangle and the segment
    for(unsigned int i=0;i<=2;i++) for(unsigned int j=3;j<=4;j++){
        plus=minus=0;mask=(1<<i)|(1<<j);
        for(unsigned int k=0;k<i;k++)    if(d3p[(1<<k)|mask]>(T)0) plus++;  else if(d3p[(1<<k)|mask]<(T)0) minus++;
        for(unsigned int k=i+1;k<=2;k++) if(d3p[(1<<k)|mask]>(T)0) minus++; else if(d3p[(1<<k)|mask]<(T)0) plus++;
        for(unsigned int k=3;k<j;k++)    if(d3p[(1<<k)|mask]>(T)0) plus++;  else if(d3p[(1<<k)|mask]<(T)0) minus++;
        for(unsigned int k=j+1;k<=4;k++) if(d3p[(1<<k)|mask]>(T)0) minus++; else if(d3p[(1<<k)|mask]<(T)0) plus++;
        if (plus==0||minus==0) return false;}
    // No separator line found, objects are intersecting
    return true;
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T> bool SIMPLEX_INTERACTIONS<T>::
Intersection(const VECTOR<VECTOR<T,2>,3>& triangle_1,const VECTOR<VECTOR<T,2>,3>& triangle_2)
{
    T points[6][2],m2p[64],d3p[64];
    int plus,minus;
    unsigned int mask;

    points[0][0]=triangle_1[1].x;points[0][1]=triangle_1[1].y;points[1][0]=triangle_1[2].x;points[1][1]=triangle_1[2].y;points[2][0]=triangle_1[3].x;points[2][1]=triangle_1[3].y;
    points[3][0]=triangle_2[1].x;points[3][1]=triangle_2[1].y;points[4][0]=triangle_2[2].x;points[4][1]=triangle_2[2].y;points[5][0]=triangle_2[3].x;points[5][1]=triangle_2[3].y;
    for(unsigned int i=0;i<=4;i++) for(unsigned int j=i+1;j<=5;j++) m2p[(1<<i)|(1<<j)]=points[i][0]*points[j][1]-points[i][1]*points[j][0];
    // Compute signed areas for all triangles formed by three input points
    for(unsigned int i=0;i<=3;i++) for(unsigned int j=i+1;j<=4;j++) for(unsigned int k=j+1;k<=5;k++) 
        d3p[(1<<i)|(1<<j)|(1<<k)]=m2p[(1<<j)|(1<<k)]-m2p[(1<<i)|(1<<k)]+m2p[(1<<i)|(1<<j)]; 
    // Separator line must contain a vertex from each triangle
    for(unsigned int i=0;i<=2;i++) for(unsigned int j=3;j<=5;j++){
        plus=minus=0;mask=(1<<i)|(1<<j);
        for(unsigned int k=0;k<i;k++)    if(d3p[(1<<k)|mask]>(T)0) plus++;  else if(d3p[(1<<k)|mask]<(T)0) minus++;
        for(unsigned int k=i+1;k<=2;k++) if(d3p[(1<<k)|mask]>(T)0) minus++; else if(d3p[(1<<k)|mask]<(T)0) plus++;
        for(unsigned int k=3;k<j;k++)    if(d3p[(1<<k)|mask]>(T)0) plus++;  else if(d3p[(1<<k)|mask]<(T)0) minus++;
        for(unsigned int k=j+1;k<=5;k++) if(d3p[(1<<k)|mask]>(T)0) minus++; else if(d3p[(1<<k)|mask]<(T)0) plus++;
        if (plus==0||minus==0) return false;}
    // No separator line found, objects are intersecting
    return true;
}
//#####################################################################
// Function Intersection
//#####################################################################
template <class T> bool SIMPLEX_INTERACTIONS<T>::
Intersection(const VECTOR<VECTOR<T,3>,3>& triangle_1,const VECTOR<VECTOR<T,3>,3>& triangle_2)
{
    T points[6][3],m3p[64],d4p[64];
    int plus,minus;
    unsigned int mask;

    points[0][0]=triangle_1[1].x;points[0][1]=triangle_1[1].y;points[0][2]=triangle_1[1].z;points[1][0]=triangle_1[2].x;points[1][1]=triangle_1[2].y;points[1][2]=triangle_1[2].z;
    points[2][0]=triangle_1[3].x;points[2][1]=triangle_1[3].y;points[2][2]=triangle_1[3].z;
    points[3][0]=triangle_2[1].x;points[3][1]=triangle_2[1].y;points[3][2]=triangle_2[1].z;points[4][0]=triangle_2[2].x;points[4][1]=triangle_2[2].y;points[4][2]=triangle_2[2].z;
    points[5][0]=triangle_2[3].x;points[5][1]=triangle_2[3].y;points[5][2]=triangle_2[3].z;
    for (unsigned int i=0;i<=3;i++) for (unsigned int j=i+1;j<=4;j++) for (unsigned int k=j+1;k<=5;k++)
        m3p[(1<<i)|(1<<j)|(1<<k)]=
            +points[i][0]*points[j][1]*points[k][2]+points[i][1]*points[j][2]*points[k][0]+points[i][2]*points[j][0]*points[k][1]
            -points[i][0]*points[j][2]*points[k][1]-points[i][1]*points[j][0]*points[k][2]-points[i][2]*points[j][1]*points[k][0];
    // Compute signed volumes for all tetrahedra formed by four input points
    for (unsigned int i=0;i<=2;i++) for (unsigned int j=i+1;j<=3;j++) for (unsigned int k=j+1;k<=4;k++) for (unsigned int l=k+1;l<=5;l++)
        d4p[(1<<i)|(1<<j)|(1<<k)|(1<<l)]=
            -m3p[(1<<j)|(1<<k)|(1<<l)]+m3p[(1<<i)|(1<<k)|(1<<l)]
            -m3p[(1<<i)|(1<<j)|(1<<l)]+m3p[(1<<i)|(1<<j)|(1<<k)];
    // Case 1 : Separator plane contains two vertices on first triangle
    for (unsigned int i=0;i<=2;i++)             // The vertex of the 1st triangle not on the plane
        for (unsigned int j=3;j<=5;j++) {       // The vertex of the 2nd triangle on the plane
            plus=minus=0;
            mask=((1<<j)|0x7)&(~(1<<i));
            switch(i) { case 0: if (d4p[0x01|mask]>(T)0) plus++;  else if (d4p[0x01|mask]<(T)0) minus++; break;
                        case 1: if (d4p[0x02|mask]>(T)0) minus++; else if (d4p[0x02|mask]<(T)0) plus++; break;
                        case 2: if (d4p[0x04|mask]>(T)0) plus++;  else if (d4p[0x04|mask]<(T)0) minus++; break; }    
            for(unsigned int k=3;k<j;k++)    if (d4p[(1<<k)|mask]>(T)0) minus++; else if (d4p[(1<<k)|mask]<(T)0) plus++;
            for(unsigned int k=j+1;k<=5;k++) if (d4p[(1<<k)|mask]>(T)0) plus++;  else if (d4p[(1<<k)|mask]<(T)0) minus++;
            if (plus==0||minus==0) return false; }
    // Case 2 : Separator plane contains two vertices on the second triangle
    for (unsigned int i=0;i<=2;i++)             // The vertex of the 1st triangle on the plane
        for (unsigned int j=3;j<=5;j++) {       // The vertex of the 2nd triangle not on the plane
            plus=minus=0;
            mask=((1<<i)|0x38)&(~(1<<j));
            for(unsigned int k=0;k<i;k++)    if (d4p[(1<<k)|mask]>(T)0) plus++;  else if (d4p[(1<<k)|mask]<(T)0) minus++;
            for(unsigned int k=i+1;k<=2;k++) if (d4p[(1<<k)|mask]>(T)0) minus++; else if (d4p[(1<<k)|mask]<(T)0) plus++;
            switch(j) { case 3: if (d4p[0x08|mask]>(T)0) plus++;  else if (d4p[0x08|mask]<(T)0) minus++; break;
                        case 4: if (d4p[0x10|mask]>(T)0) minus++; else if (d4p[0x10|mask]<(T)0) plus++; break;
                        case 5: if (d4p[0x20|mask]>(T)0) plus++;  else if (d4p[0x20|mask]<(T)0) minus++; break; }    
            if (plus==0||minus==0) return false; }
    // No separator plane found, objects are intersecting
    return true;
}
//#####################################################################
// Function Intersection
//#####################################################################
template <class T> bool SIMPLEX_INTERACTIONS<T>::
Intersection(const VECTOR<VECTOR<T,3>,3>& triangle,const VECTOR<VECTOR<T,3>,2>& segment)
{
    T points[5][3],m3p[32],d4p[32];
    int plus,minus;
    unsigned int mask;

    points[0][0]=triangle[1].x;points[0][1]=triangle[1].y;points[0][2]=triangle[1].z;points[1][0]=triangle[2].x;points[1][1]=triangle[2].y;points[1][2]=triangle[2].z;
    points[2][0]=triangle[3].x;points[2][1]=triangle[3].y;points[2][2]=triangle[3].z;
    points[3][0]=segment[1].x;points[3][1]=segment[1].y;points[3][2]=segment[1].z;points[4][0]=segment[2].x;points[4][1]=segment[2].y;points[4][2]=segment[2].z;
    for (unsigned int i=0;i<=2;i++) for (unsigned int j=i+1;j<=3;j++) for (unsigned int k=j+1;k<=4;k++)
        m3p[(1<<i)|(1<<j)|(1<<k)]=
            +points[i][0]*points[j][1]*points[k][2]+points[i][1]*points[j][2]*points[k][0]+points[i][2]*points[j][0]*points[k][1]
            -points[i][0]*points[j][2]*points[k][1]-points[i][1]*points[j][0]*points[k][2]-points[i][2]*points[j][1]*points[k][0];
    // Compute signed volumes for all tetrahedra formed by four input points
    for (unsigned int i=0;i<=1;i++) for (unsigned int j=i+1;j<=2;j++) for (unsigned int k=j+1;k<=3;k++) for (unsigned int l=k+1;l<=4;l++)
        d4p[(1<<i)|(1<<j)|(1<<k)|(1<<l)]=
            -m3p[(1<<j)|(1<<k)|(1<<l)]+m3p[(1<<i)|(1<<k)|(1<<l)]
            -m3p[(1<<i)|(1<<j)|(1<<l)]+m3p[(1<<i)|(1<<j)|(1<<k)];
    // Separator plane must contain two triangle vertices and one segment vertex 
    for (unsigned int i=0;i<=1;i++)             // 1st triangle vertex on plane
        for (unsigned int j=i+1;j<=2;j++)       // 2nd triangle vertex on plane
            for (unsigned int k=3;k<=4;k++) {   // Segment vertex on plane
                plus=minus=0;
                mask=(1<<i)|(1<<j)|(1<<k);
                for (unsigned int l=0;l<i;l++)    if (d4p[(1<<l)|mask]>(T)0) plus++;  else if (d4p[(1<<l)|mask]<(T)0) minus++;
                for (unsigned int l=i+1;l<j;l++)  if (d4p[(1<<l)|mask]>(T)0) minus++; else if (d4p[(1<<l)|mask]<(T)0) plus++;
                for (unsigned int l=j+1;l<=2;l++) if (d4p[(1<<l)|mask]>(T)0) plus++;  else if (d4p[(1<<l)|mask]<(T)0) minus++;
                for (unsigned int l=3;l<k;l++)    if (d4p[(1<<l)|mask]>(T)0) minus++; else if (d4p[(1<<l)|mask]<(T)0) plus++;
                for (unsigned int l=k+1;l<=4;l++) if (d4p[(1<<l)|mask]>(T)0) plus++;  else if (d4p[(1<<l)|mask]<(T)0) minus++;
                if (plus==0||minus==0) return false; }
    // No separator plane found, objects are intersecting
    return true;
}
//#####################################################################
// Function Intersection_Fraction
//#####################################################################
template <class T> T SIMPLEX_INTERACTIONS<T>::
Intersection_Fraction(const VECTOR<VECTOR<T,3>,3>& triangle,const VECTOR<VECTOR<T,3>,2>& segment)
{
    assert(segment[1]!=segment[2]);
    T volume1=TETRAHEDRON<T>(triangle[1],triangle[2],triangle[3],segment[1]).Signed_Volume(),volume2=TETRAHEDRON<T>(triangle[1],triangle[2],triangle[3],segment[2]).Signed_Volume();
    return volume1/(volume1-volume2);
}
//#####################################################################
// Function Intersection
//#####################################################################
template <class T> bool SIMPLEX_INTERACTIONS<T>::
Intersection(const VECTOR<VECTOR<T,3>,4>& tet,const VECTOR<VECTOR<T,3>,2>& segment)
{
    T points[6][3],m3p[64],d4p[64];
    int plus,minus;
    unsigned int mask;

    points[0][0]=tet[1].x;points[0][1]=tet[1].y;points[0][2]=tet[1].z;points[1][0]=tet[2].x;points[1][1]=tet[2].y;points[1][2]=tet[2].z;
    points[2][0]=tet[3].x;points[2][1]=tet[3].y;points[2][2]=tet[3].z;points[3][0]=tet[4].x;points[3][1]=tet[4].y;points[3][2]=tet[4].z;
    points[4][0]=segment[1].x;points[4][1]=segment[1].y;points[4][2]=segment[1].z;points[5][0]=segment[2].x;points[5][1]=segment[2].y;points[5][2]=segment[2].z;
    for (unsigned int i=0;i<=3;i++) for (unsigned int j=i+1;j<=4;j++) for (unsigned int k=j+1;k<=5;k++)
        m3p[(1<<i)|(1<<j)|(1<<k)]=
            +points[i][0]*points[j][1]*points[k][2]+points[i][1]*points[j][2]*points[k][0]+points[i][2]*points[j][0]*points[k][1]
            -points[i][0]*points[j][2]*points[k][1]-points[i][1]*points[j][0]*points[k][2]-points[i][2]*points[j][1]*points[k][0];
    // Compute signed volumes for all tetrahedra formed by four input points
    for (unsigned int i=0;i<=2;i++) for (unsigned int j=i+1;j<=3;j++) for (unsigned int k=j+1;k<=4;k++) for (unsigned int l=k+1;l<=5;l++)
        d4p[(1<<i)|(1<<j)|(1<<k)|(1<<l)]=
            -m3p[(1<<j)|(1<<k)|(1<<l)]+m3p[(1<<i)|(1<<k)|(1<<l)]
            -m3p[(1<<i)|(1<<j)|(1<<l)]+m3p[(1<<i)|(1<<j)|(1<<k)];
    // Separator plane must contain two tetrahedron vertices and one segment vertex 
    for (unsigned int i=0;i<=2;i++)             // 1st tetrahedron vertex on plane
        for (unsigned int j=i+1;j<=3;j++)       // 2nd tetrahedron vertex on plane
            for (unsigned int k=4;k<=5;k++) {   // Segment vertex on plane
                plus=minus=0;
                mask=(1<<i)|(1<<j)|(1<<k);
                for (unsigned int l=0;l<i;l++)    if (d4p[(1<<l)|mask]>(T)0) plus++;  else if (d4p[(1<<l)|mask]<(T)0) minus++;
                for (unsigned int l=i+1;l<j;l++)  if (d4p[(1<<l)|mask]>(T)0) minus++; else if (d4p[(1<<l)|mask]<(T)0) plus++;
                for (unsigned int l=j+1;l<=3;l++) if (d4p[(1<<l)|mask]>(T)0) plus++;  else if (d4p[(1<<l)|mask]<(T)0) minus++;
                for (unsigned int l=4;l<k;l++)    if (d4p[(1<<l)|mask]>(T)0) minus++; else if (d4p[(1<<l)|mask]<(T)0) plus++;
                for (unsigned int l=k+1;l<=5;l++) if (d4p[(1<<l)|mask]>(T)0) plus++;  else if (d4p[(1<<l)|mask]<(T)0) minus++;
                if (plus==0||minus==0) return false; }
    // No separator plane found, objects are intersecting
    return true;
}
//#####################################################################
// Function Intersection
//#####################################################################
template <class T> bool SIMPLEX_INTERACTIONS<T>::
Intersection(const VECTOR<VECTOR<T,3>,4>& tet,const VECTOR<VECTOR<T,3>,3>& triangle)
{
    T points[7][3],m3p[128],d4p[128];
    int plus,minus;
    unsigned int mask;

    points[0][0]=tet[1].x;points[0][1]=tet[1].y;points[0][2]=tet[1].z;points[1][0]=tet[2].x;points[1][1]=tet[2].y;points[1][2]=tet[2].z;
    points[2][0]=tet[3].x;points[2][1]=tet[3].y;points[2][2]=tet[3].z;points[3][0]=tet[4].x;points[3][1]=tet[4].y;points[3][2]=tet[4].z;
    points[4][0]=triangle[1].x;points[4][1]=triangle[1].y;points[4][2]=triangle[1].z;points[5][0]=triangle[2].x;points[5][1]=triangle[2].y;points[5][2]=triangle[2].z;
    points[6][0]=triangle[3].x;points[6][1]=triangle[3].y;points[6][2]=triangle[3].z;
    for (unsigned int i=0;i<=4;i++) for (unsigned int j=i+1;j<=5;j++) for (unsigned int k=j+1;k<=6;k++)
        m3p[(1<<i)|(1<<j)|(1<<k)]=
            +points[i][0]*points[j][1]*points[k][2]+points[i][1]*points[j][2]*points[k][0]+points[i][2]*points[j][0]*points[k][1]
            -points[i][0]*points[j][2]*points[k][1]-points[i][1]*points[j][0]*points[k][2]-points[i][2]*points[j][1]*points[k][0];
    // Compute signed volumes for all tetrahedra formed by four input points
    for (unsigned int i=0;i<=3;i++) for (unsigned int j=i+1;j<=4;j++) for (unsigned int k=j+1;k<=5;k++) for (unsigned int l=k+1;l<=6;l++)
        d4p[(1<<i)|(1<<j)|(1<<k)|(1<<l)]=
            -m3p[(1<<j)|(1<<k)|(1<<l)]+m3p[(1<<i)|(1<<k)|(1<<l)]
            -m3p[(1<<i)|(1<<j)|(1<<l)]+m3p[(1<<i)|(1<<j)|(1<<k)];
    // Case 1 : Entire tetrahedron on same half-space of triangle plane
    plus=minus=0;mask=0x70;
    for (unsigned int i=0;i<=3;i++) if (d4p[(1<<i)|mask]>(T)0) plus++; else if (d4p[(1<<i)|mask]<(T)0) minus++;
    if (plus==0||minus==0) return false;
    // Case 2 : Separator plane contains two triangle vertices
    for (unsigned int i=4;i<=6;i++)             // The 3rd triangle vertex
        for (unsigned int j=0;j<=3;j++) {       // The tetrahedron vertex on the plane
            plus=minus=0;
            mask=((1<<j)|0x70)&(~(1<<i));
            for (unsigned int k=0;k<j;k++)    if (d4p[(1<<k)|mask]>(T)0) plus++;  else if (d4p[(1<<k)|mask]<(T)0) minus++;
            for (unsigned int k=j+1;k<=3;k++) if (d4p[(1<<k)|mask]>(T)0) minus++; else if (d4p[(1<<k)|mask]<(T)0) plus++;
            switch(i) { case 4: if (d4p[0x10|mask]>(T)0) plus++;  else if (d4p[0x10|mask]<(T)0) minus++; break;
                        case 5: if (d4p[0x20|mask]>(T)0) minus++; else if (d4p[0x20|mask]<(T)0) plus++; break;
                        case 6: if (d4p[0x40|mask]>(T)0) plus++;  else if (d4p[0x40|mask]<(T)0) minus++; break; } 
            if (plus==0||minus==0) return false; }
    // Case 3 : Separator plane contains two tetraherdon vertices
    for (unsigned int i=0;i<=2;i++)             // 1st tetrahedron vertex on plane
        for (unsigned int j=i+1;j<=3;j++)       // 2nd tetrahedron vertex on plane
            for (unsigned int k=4;k<=6;k++) {   // Triangle vertex on plane
                plus=minus=0;
                mask=(1<<i)|(1<<j)|(1<<k);
                for (unsigned int l=0;l<i;l++)    if (d4p[(1<<l)|mask]>(T)0) plus++;  else if (d4p[(1<<l)|mask]<(T)0) minus++;
                for (unsigned int l=i+1;l<j;l++)  if (d4p[(1<<l)|mask]>(T)0) minus++; else if (d4p[(1<<l)|mask]<(T)0) plus++;
                for (unsigned int l=j+1;l<=3;l++) if (d4p[(1<<l)|mask]>(T)0) plus++;  else if (d4p[(1<<l)|mask]<(T)0) minus++;
                for (unsigned int l=4;l<k;l++)    if (d4p[(1<<l)|mask]>(T)0) minus++; else if (d4p[(1<<l)|mask]<(T)0) plus++;
                for (unsigned int l=k+1;l<=6;l++) if (d4p[(1<<l)|mask]>(T)0) plus++;  else if (d4p[(1<<l)|mask]<(T)0) minus++;
                if (plus==0||minus==0) return false; }
    // No separator plane found, objects are intersecting
    return true;
}
//#####################################################################
// Function Intersection
//#####################################################################
template <class T> bool SIMPLEX_INTERACTIONS<T>::
Intersection(const VECTOR<VECTOR<T,3>,4>& tetrahedron_1,const VECTOR<VECTOR<T,3>,4>& tetrahedron_2)
{
    TETRAHEDRON<T> tet(tetrahedron_1);
    if (tet.Inside(tetrahedron_2[1])||tet.Inside(tetrahedron_2[2])||tet.Inside(tetrahedron_2[3])||tet.Inside(tetrahedron_2[4])) return true;
    if (Intersection(tetrahedron_2,tetrahedron_1.Remove_Index(4)) || Intersection(tetrahedron_2,tetrahedron_1.Remove_Index(1)) 
        || Intersection(tetrahedron_2,VECTOR<VECTOR<T,3>,3>(tetrahedron_1[3],tetrahedron_1[4],tetrahedron_1[1])) 
        || Intersection(tetrahedron_2,VECTOR<VECTOR<T,3>,3>(tetrahedron_1[4],tetrahedron_1[1],tetrahedron_1[2]))) return true;
    return false;
}
//#####################################################################
// Function Three_Triangle_Intersection_Barycentric_Coordinates
//#####################################################################
template <class T> static VECTOR<T,2>
Three_Triangle_Intersection_Barycentric_Coordinates_Helper(const VECTOR<VECTOR<T,3>,3>& triangle1,const VECTOR<VECTOR<T,3>,3>& triangle2,const VECTOR<VECTOR<T,3>,3>& triangle3)
{
    T vol_x2=VECTOR<T,3>::Triple_Product(triangle2.x-triangle1.x,triangle2.y-triangle1.x,triangle2.z-triangle1.x);
    T vol_y2=VECTOR<T,3>::Triple_Product(triangle2.x-triangle1.y,triangle2.y-triangle1.y,triangle2.z-triangle1.y);
    T vol_z2=VECTOR<T,3>::Triple_Product(triangle2.x-triangle1.z,triangle2.y-triangle1.z,triangle2.z-triangle1.z);
    T vol_x3=VECTOR<T,3>::Triple_Product(triangle3.x-triangle1.x,triangle3.y-triangle1.x,triangle3.z-triangle1.x);
    T vol_y3=VECTOR<T,3>::Triple_Product(triangle3.x-triangle1.y,triangle3.y-triangle1.y,triangle3.z-triangle1.y);
    T vol_z3=VECTOR<T,3>::Triple_Product(triangle3.x-triangle1.z,triangle3.y-triangle1.z,triangle3.z-triangle1.z);
    return MATRIX<T,2>(vol_z2-vol_x2,vol_z3-vol_x3,vol_z2-vol_y2,vol_z3-vol_y3).Solve_Linear_System(VECTOR<T,2>(vol_z2,vol_z3));
}
template <class T> bool SIMPLEX_INTERACTIONS<T>::
Three_Triangle_Intersection_Barycentric_Coordinates(const VECTOR<VECTOR<T,3>,3>& triangle1,const VECTOR<VECTOR<T,3>,3>& triangle2,const VECTOR<VECTOR<T,3>,3>& triangle3,
        VECTOR<T,2>& weights1,VECTOR<T,2>& weights2,VECTOR<T,2>& weights3)
{
    weights1=Three_Triangle_Intersection_Barycentric_Coordinates_Helper(triangle1,triangle2,triangle3);
    weights2=Three_Triangle_Intersection_Barycentric_Coordinates_Helper(triangle2,triangle3,triangle1);
    weights3=Three_Triangle_Intersection_Barycentric_Coordinates_Helper(triangle3,triangle1,triangle2);
    if(weights1.Max()>(T)1||weights1.Min()<(T)0||weights1.Sum()>(T)1||
        weights2.Max()>(T)1||weights2.Min()<(T)0||weights2.Sum()>(T)1||
        weights3.Max()>(T)1||weights3.Min()<(T)0||weights3.Sum()>(T)1) return false;
    return true;
}
//#####################################################################
// Function Triangle_Segment_Intersection_Fraction
//#####################################################################
template <class T> void SIMPLEX_INTERACTIONS<T>::
Triangle_Segment_Intersection_Barycentric_Coordinates(const VECTOR<VECTOR<T,3>,3>& triangle,const VECTOR<VECTOR<T,3>,2>& segment,VECTOR<T,2>& triangle_weights,T& segment_weight)
{
    VECTOR<T,3> weights=MATRIX<T,3>(triangle.x-triangle.z,triangle.y-triangle.z,segment.y-segment.x).Robust_Solve_Linear_System(VECTOR<T,3>(segment.y-triangle.z));
    triangle_weights=VECTOR<T,2>(weights.x,weights.y);segment_weight=weights.z;
}
//####################################################################
template class SIMPLEX_INTERACTIONS<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SIMPLEX_INTERACTIONS<double>;
#endif

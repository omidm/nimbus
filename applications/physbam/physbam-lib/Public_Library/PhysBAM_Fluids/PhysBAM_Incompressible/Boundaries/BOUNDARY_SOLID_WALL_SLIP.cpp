//#####################################################################
// Copyright 2004, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_QUADTREE_MESH.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Boundaries/BOUNDARY_SOLID_WALL_SLIP.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BOUNDARY_SOLID_WALL_SLIP<TV>::
BOUNDARY_SOLID_WALL_SLIP(const bool left_constant_extrapolation_input,const bool right_constant_extrapolation_input,const bool bottom_constant_extrapolation_input,
    const bool top_constant_extrapolation_input,const bool front_constant_extrapolation_input,const bool back_constant_extrapolation_input)
{
    Set_Constant_Extrapolation(left_constant_extrapolation_input,right_constant_extrapolation_input,bottom_constant_extrapolation_input,top_constant_extrapolation_input,
        front_constant_extrapolation_input,back_constant_extrapolation_input);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BOUNDARY_SOLID_WALL_SLIP<TV>::
~BOUNDARY_SOLID_WALL_SLIP()
{
}
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV> void BOUNDARY_SOLID_WALL_SLIP<TV>::
Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAY<TV,VECTOR<int,2> >& V,ARRAY<TV,VECTOR<int,2> >& V_ghost,const T dt,const T time,const int number_of_ghost_cells)
{
    if(V.length != 1) PHYSBAM_FUNCTION_IS_NOT_DEFINED();
    int i,j,m=grid.m,n=grid.n;

    for(i=1;i<=m;i++) for(j=1;j<=n;j++) V_ghost(i,j)=V(i,j); // interior

    if(left_constant_extrapolation) Fill_Left_Ghost_Cells(grid,V_ghost,time);
    else for(j=1;j<=n;j++){
        V_ghost(0,j).x=-V_ghost(2,j).x;V_ghost(-1,j).x=-V_ghost(3,j).x;V_ghost(-2,j).x=-V_ghost(4,j).x;
        V_ghost(0,j).y=V_ghost(2,j).y;V_ghost(-1,j).y=V_ghost(3,j).y;V_ghost(-2,j).y=V_ghost(4,j).y;}  
    if(right_constant_extrapolation) Fill_Right_Ghost_Cells(grid,V_ghost,time);
    else for(j=1;j<=n;j++){
        V_ghost(m+1,j).x=-V_ghost(m-1,j).x;V_ghost(m+2,j).x=-V_ghost(m-2,j).x;V_ghost(m+3,j).x=-V_ghost(m-3,j).x;
        V_ghost(m+1,j).y=V_ghost(m-1,j).y;V_ghost(m+2,j).y=V_ghost(m-2,j).y;V_ghost(m+3,j).y=V_ghost(m-3,j).y;}  
    if(bottom_constant_extrapolation) Fill_Bottom_Ghost_Cells(grid,V_ghost,time);
    else for(i=-2;i<=m+3;i++){
        V_ghost(i,0).x=V_ghost(i,2).x;V_ghost(i,-1).x=V_ghost(i,3).x;V_ghost(i,-2).x=V_ghost(i,4).x;
        V_ghost(i,0).y=-V_ghost(i,2).y;V_ghost(i,-1).y=-V_ghost(i,3).y;V_ghost(i,-2).y=-V_ghost(i,4).y;}  
    if(top_constant_extrapolation) Fill_Top_Ghost_Cells(grid,V_ghost,time);
    else for(i=-2;i<=m+3;i++){
        V_ghost(i,n+1).x=V_ghost(i,n-1).x;V_ghost(i,n+2).x=V_ghost(i,n-2).x;V_ghost(i,n+3).x=V_ghost(i,n-3).x;
        V_ghost(i,n+1).y=-V_ghost(i,n-1).y;V_ghost(i,n+2).y=-V_ghost(i,n-2).y;V_ghost(i,n+3).y=-V_ghost(i,n-3).y;}  
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class TV> void BOUNDARY_SOLID_WALL_SLIP<TV>::
Apply_Boundary_Condition(const GRID<TV>& grid,ARRAY<TV,VECTOR<int,2> >& V,const T time) 
{
    if(V.length != 1) PHYSBAM_FUNCTION_IS_NOT_DEFINED();
    int i,j,m=grid.m,n=grid.n;
    if(!left_constant_extrapolation)for(j=1;j<=n;j++)V(1,j).x=0;
    if(!right_constant_extrapolation)for(j=1;j<=n;j++)V(m,j).x=0;  
    if(!bottom_constant_extrapolation)for(i=1;i<=m;i++)V(i,1).y=0;
    if(!top_constant_extrapolation)for(i=1;i<=m;i++)V(i,n).y=0;
}
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV> void BOUNDARY_SOLID_WALL_SLIP<TV>::
Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAY<TV,VECTOR<int,3> >& V,ARRAY<TV,VECTOR<int,3> >& V_ghost,const T dt,const T time,const int number_of_ghost_cells)
{
    if(V.length != 1) PHYSBAM_FUNCTION_IS_NOT_DEFINED();
    int i,j,ij,m=grid.m,n=grid.n,mn=grid.mn;

    for(i=1;i<=m;i++) for(j=1;j<=n;j++) for(ij=1;ij<=mn;ij++) V_ghost(i,j,ij)=V(i,j,ij); // interior

    if(left_constant_extrapolation) Fill_Left_Ghost_Cells(grid,V_ghost,time);
    else for(j=1;j<=n;j++) for(ij=1;ij<=mn;ij++){
        V_ghost(0,j,ij)[1]=-V_ghost(2,j,ij)[1];V_ghost(-1,j,ij)[1]=-V_ghost(3,j,ij)[1];V_ghost(-2,j,ij)[1]=-V_ghost(4,j,ij)[1];
        V_ghost(0,j,ij)[2]=V_ghost(2,j,ij)[2];V_ghost(-1,j,ij)[2]=V_ghost(3,j,ij)[2];V_ghost(-2,j,ij)[2]=V_ghost(4,j,ij)[2];
        V_ghost(0,j,ij)[3]=V_ghost(2,j,ij)[3];V_ghost(-1,j,ij)[3]=V_ghost(3,j,ij)[3];V_ghost(-2,j,ij)[3]=V_ghost(4,j,ij)[3];}  
    if(right_constant_extrapolation) Fill_Right_Ghost_Cells(grid,V_ghost,time);
    else for(j=1;j<=n;j++) for(ij=1;ij<=mn;ij++){
        V_ghost(m+1,j,ij)[1]=-V_ghost(m-1,j,ij)[1];V_ghost(m+2,j,ij)[1]=-V_ghost(m-2,j,ij)[1];V_ghost(m+3,j,ij)[1]=-V_ghost(m-3,j,ij)[1];
        V_ghost(m+1,j,ij)[2]=V_ghost(m-1,j,ij)[2];V_ghost(m+2,j,ij)[2]=V_ghost(m-2,j,ij)[2];V_ghost(m+3,j,ij)[2]=V_ghost(m-3,j,ij)[2];
        V_ghost(m+1,j,ij)[3]=V_ghost(m-1,j,ij)[3];V_ghost(m+2,j,ij)[3]=V_ghost(m-2,j,ij)[3];V_ghost(m+3,j,ij)[3]=V_ghost(m-3,j,ij)[3];}  
    if(bottom_constant_extrapolation) Fill_Bottom_Ghost_Cells(grid,V_ghost,time);
    else for(i=-2;i<=m+3;i++) for(ij=1;ij<=mn;ij++){
        V_ghost(i,0,ij)[1]=V_ghost(i,2,ij)[1];V_ghost(i,-1,ij)[1]=V_ghost(i,3,ij)[1];V_ghost(i,-2,ij)[1]=V_ghost(i,4,ij)[1];
        V_ghost(i,0,ij)[2]=-V_ghost(i,2,ij)[2];V_ghost(i,-1,ij)[2]=-V_ghost(i,3,ij)[2];V_ghost(i,-2,ij)[2]=-V_ghost(i,4,ij)[2];
        V_ghost(i,0,ij)[3]=V_ghost(i,2,ij)[3];V_ghost(i,-1,ij)[3]=V_ghost(i,3,ij)[3];V_ghost(i,-2,ij)[3]=V_ghost(i,4,ij)[3];}
    if(top_constant_extrapolation) Fill_Top_Ghost_Cells(grid,V_ghost,time);
    else for(i=-2;i<=m+3;i++) for(ij=1;ij<=mn;ij++){
        V_ghost(i,n+1,ij)[1]=V_ghost(i,n-1,ij)[1];V_ghost(i,n+2,ij)[1]=V_ghost(i,n-2,ij)[1];V_ghost(i,n+3,ij)[1]=V_ghost(i,n-3,ij)[1];
        V_ghost(i,n+1,ij)[2]=-V_ghost(i,n-1,ij)[2];V_ghost(i,n+2,ij)[2]=-V_ghost(i,n-2,ij)[2];V_ghost(i,n+3,ij)[2]=-V_ghost(i,n-3,ij)[2];
        V_ghost(i,n+1,ij)[3]=V_ghost(i,n-1,ij)[3];V_ghost(i,n+2,ij)[3]=V_ghost(i,n-2,ij)[3];V_ghost(i,n+3,ij)[3]=V_ghost(i,n-3,ij)[3];}
    if(front_constant_extrapolation) Fill_Front_Ghost_Cells(grid,V_ghost,time);
    else for(i=-2;i<=m+3;i++) for(j=-2;j<=n+3;j++){
        V_ghost(i,j,0)[1]=V_ghost(i,j,2)[1];V_ghost(i,j,-1)[1]=V_ghost(i,j,3)[1];V_ghost(i,j,-2)[1]=V_ghost(i,j,4)[1];
        V_ghost(i,j,0)[2]=V_ghost(i,j,2)[2];V_ghost(i,j,-1)[2]=V_ghost(i,j,3)[2];V_ghost(i,j,-2)[2]=V_ghost(i,j,4)[2];
        V_ghost(i,j,0)[3]=-V_ghost(i,j,2)[3];V_ghost(i,j,-1)[3]=-V_ghost(i,j,3)[3];V_ghost(i,j,-2)[3]=-V_ghost(i,j,4)[3];}
    if(back_constant_extrapolation) Fill_Back_Ghost_Cells(grid,V_ghost,time);
    else for(i=-2;i<=m+3;i++) for(j=-2;j<=n+3;j++){
        V_ghost(i,j,mn+1)[1]=V_ghost(i,j,mn-1)[1];V_ghost(i,j,mn+2)[1]=V_ghost(i,j,mn-2)[1];V_ghost(i,j,mn+3)[1]=V_ghost(i,j,mn-3)[1];
        V_ghost(i,j,mn+1)[2]=V_ghost(i,j,mn-1)[2];V_ghost(i,j,mn+2)[2]=V_ghost(i,j,mn-2)[2];V_ghost(i,j,mn+3)[2]=V_ghost(i,j,mn-3)[2];
        V_ghost(i,j,mn+1)[3]=-V_ghost(i,j,mn-1)[3];V_ghost(i,j,mn+2)[3]=-V_ghost(i,j,mn-2)[3];V_ghost(i,j,mn+3)[3]=-V_ghost(i,j,mn-3)[3];}
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class TV> void BOUNDARY_SOLID_WALL_SLIP<TV>::
Apply_Boundary_Condition(const GRID<TV>& grid,ARRAY<TV,VECTOR<int,3> >& V,const T time) 
{
    if(V.length != 1) PHYSBAM_FUNCTION_IS_NOT_DEFINED();
    int i,j,ij,m=grid.m,n=grid.n,mn=grid.mn;
    if(!left_constant_extrapolation)for(j=1;j<=n;j++)for(ij=1;ij<=mn;ij++)V(1,j,ij)[1]=0;
    if(!right_constant_extrapolation)for(j=1;j<=n;j++)for(ij=1;ij<=mn;ij++)V(m,j,ij)[1]=0;
    if(!bottom_constant_extrapolation)for(i=1;i<=m;i++)for(ij=1;ij<=mn;ij++)V(i,1,ij)[2]=0;
    if(!top_constant_extrapolation)for(i=1;i<=m;i++)for(ij=1;ij<=mn;ij++)V(i,n,ij)[2]=0;
    if(!front_constant_extrapolation)for(i=1;i<=m;i++)for(j=1;j<=n;j++)V(i,j,1)[3]=0;
    if(!back_constant_extrapolation)for(i=1;i<=m;i++)for(j=1;j<=n;j++)V(i,j,mn)[3]=0;
}
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class TV> void BOUNDARY_SOLID_WALL_SLIP<TV>::
Apply_Boundary_Condition(const QUADTREE_GRID<T>& grid,ARRAY<TV>& V,const T time)
{
    if(!left_constant_extrapolation) for(DYADIC_GRID_ITERATOR_NODE<QUADTREE_GRID<T> > iterator(grid,grid.Map_Individual_Side_Boundary_Nodes(1));iterator.Valid();iterator.Next()) V(iterator.Node_Index())[1]=0;
    if(!right_constant_extrapolation) for(DYADIC_GRID_ITERATOR_NODE<QUADTREE_GRID<T> > iterator(grid,grid.Map_Individual_Side_Boundary_Nodes(2));iterator.Valid();iterator.Next()) V(iterator.Node_Index())[1]=0;
    if(!bottom_constant_extrapolation) for(DYADIC_GRID_ITERATOR_NODE<QUADTREE_GRID<T> > iterator(grid,grid.Map_Individual_Side_Boundary_Nodes(3));iterator.Valid();iterator.Next()) V(iterator.Node_Index())[2]=0;
    if(!top_constant_extrapolation) for(DYADIC_GRID_ITERATOR_NODE<QUADTREE_GRID<T> > iterator(grid,grid.Map_Individual_Side_Boundary_Nodes(4));iterator.Valid();iterator.Next()) V(iterator.Node_Index())[2]=0;
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class TV> void BOUNDARY_SOLID_WALL_SLIP<TV>::
Apply_Boundary_Condition(const OCTREE_GRID<T>& grid,ARRAY<TV>& V,const T time)
{
    if(!left_constant_extrapolation) for(DYADIC_GRID_ITERATOR_NODE<OCTREE_GRID<T> > iterator(grid,grid.Map_Individual_Side_Boundary_Nodes(1));iterator.Valid();iterator.Next()) V(iterator.Node_Index())[1]=0;
    if(!right_constant_extrapolation) for(DYADIC_GRID_ITERATOR_NODE<OCTREE_GRID<T> > iterator(grid,grid.Map_Individual_Side_Boundary_Nodes(2));iterator.Valid();iterator.Next()) V(iterator.Node_Index())[1]=0;
    if(!bottom_constant_extrapolation) for(DYADIC_GRID_ITERATOR_NODE<OCTREE_GRID<T> > iterator(grid,grid.Map_Individual_Side_Boundary_Nodes(3));iterator.Valid();iterator.Next()) V(iterator.Node_Index())[2]=0;
    if(!top_constant_extrapolation) for(DYADIC_GRID_ITERATOR_NODE<OCTREE_GRID<T> > iterator(grid,grid.Map_Individual_Side_Boundary_Nodes(4));iterator.Valid();iterator.Next()) V(iterator.Node_Index())[2]=0;
    if(!front_constant_extrapolation) for(DYADIC_GRID_ITERATOR_NODE<OCTREE_GRID<T> > iterator(grid,grid.Map_Individual_Side_Boundary_Nodes(5));iterator.Valid();iterator.Next()) V(iterator.Node_Index())[3]=0;
    if(!back_constant_extrapolation) for(DYADIC_GRID_ITERATOR_NODE<OCTREE_GRID<T> > iterator(grid,grid.Map_Individual_Side_Boundary_Nodes(6));iterator.Valid();iterator.Next()) V(iterator.Node_Index())[3]=0;
}
#endif

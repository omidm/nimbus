//#####################################################################
// Copyright 2002-2007, Doug Enright, Jon Gretarsson, Nipun Kwatra, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC
//#####################################################################
//
// If an axis is set to be periodic then it is periodic with period m-1*grid.dx:U(m)=U(1), if repeats_at_last_node is true, otherwise its periodic with period m*grid.dx: U(m+1)=U(1).
//
//#####################################################################
#ifndef __BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC__
#define __BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class T_GRID,class T2> // d=T_GRID::dimension+2
class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC:public BOUNDARY_UNIFORM<T_GRID,T2>
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::VECTOR_T TV;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;
    typedef typename T_ARRAYS_BASE::template REBIND<T2>::TYPE T_ARRAYS_DIMENSION_BASE;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
public:
    typedef BOUNDARY_UNIFORM<T_GRID,T2> BASE;
    using BASE::Turn_Off_Constant_Extrapolation;using BASE::Set_Constant_Extrapolation;using BASE::Constant_Extrapolation;
    using BASE::constant_extrapolation;

    VECTOR<bool,3> periodic,repeats_at_last_node;

    BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC(bool x_periodic=true,bool y_periodic=true,bool z_periodic=true,bool repeats_at_last_x_node=false,bool repeats_at_last_y_node=false,
            bool repeats_at_last_z_node=false):periodic(x_periodic,y_periodic,z_periodic),repeats_at_last_node(repeats_at_last_x_node,repeats_at_last_y_node,repeats_at_last_z_node)
    {
        Turn_Off_Constant_Extrapolation();
    }

//#####################################################################
    void Fill_Ghost_Cells_Helper(const GRID<TV>& grid,const ARRAY<VECTOR<T,3> ,VECTOR<int,1> >& u,ARRAY<VECTOR<T,3> ,VECTOR<int,1> >& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3);
    void Fill_Ghost_Cells_Helper(const GRID<TV>& grid,const ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& u,ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3);
    void Apply_Boundary_Condition_Helper(const GRID<TV>& grid,ARRAY<VECTOR<T,3> ,VECTOR<int,1> >& u,const T time);
    void Apply_Boundary_Condition_Helper(const GRID<TV>& grid,ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& u,const T time);
    void Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_DIMENSION_BASE& u,T_ARRAYS_DIMENSION_BASE& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) PHYSBAM_OVERRIDE;
    void Apply_Boundary_Condition(const T_GRID& grid,T_ARRAYS_DIMENSION_BASE& u,const T time) PHYSBAM_OVERRIDE;
//#####################################################################
};
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC<T_GRID,T2>::
Fill_Ghost_Cells_Helper(const GRID<TV>& grid,const ARRAY<VECTOR<T,3> ,VECTOR<int,1> >& u,ARRAY<VECTOR<T,3> ,VECTOR<int,1> >& u_ghost,const T dt,const T time,const int number_of_ghost_cells)
{
    int m=grid.m,i,periodic_point;

    ARRAY<RANGE<VECTOR<int,1> > > regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells);
    ARRAY<VECTOR<T,3> ,VECTOR<int,1> >::Put(u,u_ghost); // interior

    if(Constant_Extrapolation(1)) Fill_Single_Ghost_Region(grid,u_ghost,1,regions(1));
    else for(i=1-number_of_ghost_cells;i<=0;i++){ // left
            periodic_point=m+i-1;
            T rho=u_ghost(1,periodic_point);
            T u_velocity=-u_ghost(2,periodic_point)/u_ghost(1,periodic_point);
            T e=u_ghost(3,periodic_point)/u_ghost(1,periodic_point)-sqr(-u_velocity)/2;
            u_ghost(1,i)=rho;
            u_ghost(2,i)=rho*u_velocity;
            u_ghost(3,i)=rho*(e+sqr(u_velocity)/2);}
    if(Constant_Extrapolation(2)) Fill_Single_Ghost_Region(grid,u_ghost,2,regions(2));
    else for(i=m+1;i<=m+number_of_ghost_cells;i++){ // right
            periodic_point=i-m+1;
            T rho=u_ghost(1,2*periodic_point);
            T u_velocity=-u_ghost(2,2*periodic_point)/u_ghost(1,2*periodic_point);
            T e=u_ghost(3,2*periodic_point)/u_ghost(1,2*periodic_point)-sqr(-u_velocity)/2;
            u_ghost(1,i)=rho;
            u_ghost(2,i)=rho*u_velocity;
            u_ghost(3,i)=rho*(e+sqr(u_velocity)/2);}
}
template<class T_GRID,class T2> void BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC<T_GRID,T2>::
Fill_Ghost_Cells_Helper(const GRID<TV>& grid,const ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& u,ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& u_ghost,const T dt,const T time,const int number_of_ghost_cells)
{
    int m=grid.m,n=grid.n;
    int i,j,k;

    ARRAY<VECTOR<T,4> ,VECTOR<int,2> >::Put(u,u_ghost); // interior

    if(constant_extrapolation[1][1]) Fill_Left_Ghost_Cells(grid,u_ghost,time);
    else
        for(j=1;j<=n;j++) for(i=-2;i<=0;i++){ // left
            T rho=u_ghost(1,2-i,j);
            T u_velocity=-u_ghost(2,2-i,j)/u_ghost(1,2-i,j);
            T v_velocity=u_ghost(3,2-i,j)/u_ghost(1,2-i,j);
            T e=u_ghost(4,2-i,j)/u_ghost(1,2-i,j)-(sqr(-u_velocity)+sqr(v_velocity))/2;
            u_ghost(1,i,j)=rho;
            u_ghost(2,i,j)=rho*u_velocity;
            u_ghost(3,i,j)=rho*v_velocity;
            u_ghost(4,i,j)=rho*(e+(sqr(u_velocity)+sqr(v_velocity))/2);}
    if(constant_extrapolation[1][2]) Fill_Right_Ghost_Cells(grid,u_ghost,time);
    else for(j=1;j<=n;j++) for(i=m+1;i<=m+3;i++){ // right
            T rho=u_ghost(1,2*m-i,j);
            T u_velocity=-u_ghost(2,2*m-i,j)/u_ghost(1,2*m-i,j);
            T v_velocity=u_ghost(3,2*m-i,j)/u_ghost(1,2*m-i,j);
            T e=u_ghost(4,2*m-i,j)/u_ghost(1,2*m-i,j)-(sqr(-u_velocity)+sqr(v_velocity))/2;
            u_ghost(1,i,j)=rho;
            u_ghost(2,i,j)=rho*u_velocity;
            u_ghost(3,i,j)=rho*v_velocity;
            u_ghost(4,i,j)=rho*(e+(sqr(u_velocity)+sqr(v_velocity))/2);}

    for(k=1;k<=u_ghost.length;k++) for(i=1;i<=m;i++){
        for(j=-2;j<=0;j++) u_ghost(k,i,j)=u_ghost(k,i,j+n-1); //bottom
        for(j=n+1;j<=n+3;j++) u_ghost(k,i,j)=u_ghost(k,i,j-n+1); //top
    }
}
template<class T_GRID,class T2> void BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC<T_GRID,T2>::
Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_DIMENSION_BASE& u,T_ARRAYS_DIMENSION_BASE& u_ghost,const T dt,const T time,const int number_of_ghost_cells)
{
    //TODO: get rid of the helper functions
    //Fill_Ghost_Cells_Helper(grid,u,u_ghost,dt,time,number_of_ghost_cells);
    TV_INT counts=grid.Counts();
    T_ARRAYS_DIMENSION_BASE::Put(u,u_ghost); // interior
    ARRAY<RANGE<TV_INT> > regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells);

    for(int axis=1;axis<=T_GRID::dimension;axis++)for(int axis_side=1;axis_side<=2;axis_side++){int side=2*axis+axis_side-2;
        if(!periodic[axis]) Fill_Single_Ghost_Region(grid,u_ghost,side,regions(side));
        else for(CELL_ITERATOR iterator(grid,regions(side));iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
                int period=repeats_at_last_node[axis]?counts[axis]-1:counts[axis];
                int axis_periodic_node=1+wrap(cell[axis]-1,period);
                TV_INT periodic_node=cell;periodic_node[axis]=axis_periodic_node;
                u_ghost(cell)=u_ghost(periodic_node);}}
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC<T_GRID,T2>::
Apply_Boundary_Condition_Helper(const GRID<TV>& grid,ARRAY<VECTOR<T,3> ,VECTOR<int,1> >& u,const T time)
{
    int m=grid.m;

    if(!constant_extrapolation[1][1]){ // left wall
        T rho=u(1,1);
        T u_velocity=u(2,1)/u(1,1);
        T e=u(3,1)/u(1,1)-sqr(u_velocity)/2;
        u_velocity=0;
        u(1,1)=rho;
        u(2,1)=rho*u_velocity;
        u(3,1)=rho*(e+sqr(u_velocity)/2);}

    if(!constant_extrapolation[1][2]){ // right wall
        T rho=u(1,m);
        T u_velocity=u(2,m)/u(1,m);
        T e=u(3,m)/u(1,m)-sqr(u_velocity)/2;
        u_velocity=0;
        u(1,m)=rho;
        u(2,m)=rho*u_velocity;
        u(3,m)=rho*(e+sqr(u_velocity)/2);}
}
template<class T_GRID,class T2> void BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC<T_GRID,T2>::
Apply_Boundary_Condition_Helper(const GRID<TV>& grid,ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& u,const T time)
{
    int m=grid.m,n=grid.n;
    int i,j,k;

    if(!constant_extrapolation[1][1])
        for(j=1;j<=n;j++){
            // left wall
            T rho=u(1,1,j);
            T u_velocity=u(2,1,j)/u(1,1,j);
            T v_velocity=u(3,1,j)/u(1,1,j);
            T e=u(4,1,j)/u(1,1,j)-(sqr(u_velocity)+sqr(v_velocity))/2;
            u_velocity=0;
            u(1,1,j)=rho;
            u(2,1,j)=rho*u_velocity;
            u(3,1,j)=rho*v_velocity;
            u(4,1,j)=rho*(e+(sqr(u_velocity)+sqr(v_velocity))/2);}

    if(!constant_extrapolation[1][2])
        for(j=1;j<=n;j++){
            // right wall
            T rho=u(1,m,j);
            T u_velocity=u(2,m,j)/u(1,m,j);
            T v_velocity=u(3,m,j)/u(1,m,j);
            T e=u(4,m,j)/u(1,m,j)-(sqr(u_velocity)+sqr(v_velocity))/2;
            u_velocity=0;
            u(1,m,j)=rho;
            u(2,m,j)=rho*u_velocity;
            u(3,m,j)=rho*v_velocity;
            u(4,m,j)=rho*(e+(sqr(u_velocity)+sqr(v_velocity))/2);}

    for(i=1;i<=m;i++) for(k=1;k<=u.length;k++) u(k,i,n)=u(k,i,1);
}
template<class T_GRID,class T2> void BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC<T_GRID,T2>::
Apply_Boundary_Condition(const T_GRID& grid,T_ARRAYS_DIMENSION_BASE& u,const T time)
{
    //TODO: get rid of the helper functions
    //Apply_Boundary_Condition_Helper(grid,u,time);
    for(int axis=1;axis<=T_GRID::dimension;axis++)
        if(periodic[axis] && repeats_at_last_node[axis]){
            for(CELL_ITERATOR iterator(grid,0,T_GRID::BOUNDARY_INTERIOR_REGION,2*axis);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
                TV_INT opposite_cell=cell_index;opposite_cell[axis]=1;
                typename T_ARRAYS_DIMENSION_BASE::ELEMENT u_average=(u(cell_index)+u(opposite_cell))*(T).5;
                u(cell_index)=u_average;u(opposite_cell)=u_average;}}
}
//#####################################################################
}
#endif

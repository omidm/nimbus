//#####################################################################
// Copyright 2002-2005, Doug Enright, Ronald Fedkiw, Geoffrey Irving, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T> void LEVELSET_ADVECTION_2D<T>::
Euler_Step(const ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocity,const T dt,const T time,const int number_of_ghost_cells)
{
    DEBUG_UTILITIES::Debug_Breakpoint();
    
    GRID<TV>& grid=levelset->grid;
    T_BOUNDARY_SCALAR* boundary=levelset->boundary;
    T_ARRAYS_SCALAR& phi=levelset->phi;
    
    assert(grid.Is_MAC_Grid());
    int m=grid.counts.x,n=grid.counts.y;
    ARRAY<T,VECTOR<int,2> > phi_ghost(grid.Domain_Indices(number_of_ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,number_of_ghost_cells); 
    
    if(levelset->curvature_motion){ // do curvature first - based on phi^n
        T one_over_two_dx=1/(2*grid.dX.x),one_over_two_dy=1/(2*grid.dX.y);
        bool curvature_defined=(levelset->curvature!=0);levelset->Compute_Curvature(time);
        for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){
            T phix=(phi_ghost(i+1,j)-phi_ghost(i-1,j))*one_over_two_dx,phiy=(phi_ghost(i,j+1)-phi_ghost(i,j-1))*one_over_two_dy;
            phi(i,j)-=dt*levelset->sigma*(*levelset->curvature)(i,j)*sqrt(sqr(phix)+sqr(phiy));}
        boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,number_of_ghost_cells); 
        if(!curvature_defined){delete levelset->curvature;levelset->curvature=0;}}

    advection->Update_Advection_Equation_Cell(grid,phi,phi_ghost,face_velocity,*boundary,dt,time);
    boundary->Apply_Boundary_Condition(grid,phi,time+dt); 
}
//#####################################################################
// Functions Reinitialize
//#####################################################################
template<class T> void LEVELSET_ADVECTION_2D<T>::
Reinitialize(const int time_steps,const T time)
{
    GRID<TV>& grid=levelset->grid;
    T_ARRAYS_SCALAR& phi=levelset->phi;

    int m=grid.counts.x,n=grid.counts.y;
    
    ARRAY<T,VECTOR<int,2> > sign_phi(1,m,1,n); // smeared out sign function
    T epsilon=sqr(grid.dX.Max());
    for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) sign_phi(i,j)=phi(i,j)/sqrt(sqr(phi(i,j))+epsilon);

    T dt=reinitialization_cfl*grid.min_dX;
    RUNGEKUTTA<ARRAY<T,VECTOR<int,2> > > rungekutta(phi); 
    rungekutta.Set_Grid_And_Boundary_Condition(grid,*levelset->boundary);
    rungekutta.Set_Order(reinitialization_runge_kutta_order);
    rungekutta.Set_Time(time);
    rungekutta.Pseudo_Time();
    for(int k=1;k<=time_steps;k++){
        rungekutta.Start(dt);
        for(int kk=1;kk<=rungekutta.order;kk++){Euler_Step_Of_Reinitialization(sign_phi,dt,time);rungekutta.Main();}
    } 
}
//#####################################################################
// Functions Euler_Step_Of_Reinitialization
//#####################################################################
template<class T> void LEVELSET_ADVECTION_2D<T>::
Euler_Step_Of_Reinitialization(const ARRAY<T,VECTOR<int,2> >& sign_phi,const T dt,const T time)
{
    GRID<TV>& grid=levelset->grid;
    T_BOUNDARY_SCALAR* boundary=levelset->boundary;
    T_ARRAYS_SCALAR& phi=levelset->phi;
    
    int i,j;int m=grid.counts.x,n=grid.counts.y;T dx=grid.dX.x,dy=grid.dX.y;
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,2> > phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,ghost_cells);
    ARRAY<T,VECTOR<int,2> > rhs(1,m,1,n);
        
    ARRAY<T,VECTOR<int,1> > phi_1d_x(1-ghost_cells,m+ghost_cells),phix_minus(1,m),phix_plus(1,m);
    for(j=1;j<=n;j++){
        for(i=1-ghost_cells;i<=m+ghost_cells;i++) phi_1d_x(i)=phi_ghost(i,j);
        if(reinitialization_spatial_order == 5) HJ_WENO(m,dx,phi_1d_x,phix_minus,phix_plus);
        else HJ_ENO(reinitialization_spatial_order,m,dx,phi_1d_x,phix_minus,phix_plus);
        for(i=1;i<=m;i++)
            if(LEVELSET_UTILITIES<T>::Sign(phi(i,j)) < 0) rhs(i,j)=sqr(max(-phix_minus(i),phix_plus(i),(T)0));
            else rhs(i,j)=sqr(max(phix_minus(i),-phix_plus(i),(T)0));}

    ARRAY<T,VECTOR<int,1> > phi_1d_y(1-ghost_cells,n+ghost_cells),phiy_minus(1,n),phiy_plus(1,n);
    for(i=1;i<=m;i++){
        for(j=1-ghost_cells;j<=n+ghost_cells;j++) phi_1d_y(j)=phi_ghost(i,j);
        if(reinitialization_spatial_order == 5) HJ_WENO(n,dy,phi_1d_y,phiy_minus,phiy_plus);
        else HJ_ENO(reinitialization_spatial_order,n,dy,phi_1d_y,phiy_minus,phiy_plus);
        for(j=1;j<=n;j++)
            if(LEVELSET_UTILITIES<T>::Sign(phi(i,j)) < 0) rhs(i,j)+=sqr(max(-phiy_minus(j),phiy_plus(j),(T)0));
            else rhs(i,j)+=sqr(max(phiy_minus(j),-phiy_plus(j),(T)0));}

    for(i=1;i<=m;i++) for(j=1;j<=n;j++){
        phi(i,j)-=dt*sign_phi(i,j)*(sqrt(rhs(i,j))-1);
        if(LEVELSET_UTILITIES<T>::Interface(phi_ghost(i,j),phi(i,j))) phi(i,j)=LEVELSET_UTILITIES<T>::Sign(phi_ghost(i,j))*levelset->small_number*grid.min_dX;} 
       
    boundary->Apply_Boundary_Condition(grid,phi,time); // time not incremented, pseudo-time
}
//#####################################################################
template class LEVELSET_ADVECTION_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_ADVECTION_2D<double>;
#endif

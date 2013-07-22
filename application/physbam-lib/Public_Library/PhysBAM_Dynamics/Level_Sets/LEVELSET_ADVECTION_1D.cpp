//#####################################################################
// Copyright 2002-2005, Doug Enright, Ronald Fedkiw, Geoffrey Irving, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_1D.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T> void LEVELSET_ADVECTION_1D<T>::
Euler_Step(const ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocity,const T dt,const T time,const int number_of_ghost_cells)
{
    DEBUG_UTILITIES::Debug_Breakpoint();
    
    GRID<TV>& grid=levelset->grid;
    T_BOUNDARY_SCALAR* boundary=levelset->boundary;
    T_ARRAYS_SCALAR& phi=levelset->phi;
    
    assert(grid.Is_MAC_Grid());
    ARRAY<T,VECTOR<int,1> > phi_ghost(grid.Domain_Indices(number_of_ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,number_of_ghost_cells);

    advection->Update_Advection_Equation_Cell(grid,phi,phi_ghost,face_velocity,*boundary,dt,time);
    boundary->Apply_Boundary_Condition(grid,phi,time+dt);
}
//#####################################################################
// Functions Reinitialize
//#####################################################################
template<class T> void LEVELSET_ADVECTION_1D<T>::
Reinitialize(const int time_steps,const T time)
{
    GRID<TV>& grid=levelset->grid;
    T_ARRAYS_SCALAR& phi=levelset->phi;

    int m=grid.counts.x;T dx=grid.dX.x;
    
    ARRAY<T,VECTOR<int,1> > sign_phi(1,m); // smeared out sign function
    T epsilon=sqr(dx);
    for(int i=1;i<=m;i++) sign_phi(i)=phi(i)/sqrt(sqr(phi(i))+epsilon);

    T dt=reinitialization_cfl*dx;
    RUNGEKUTTA<ARRAY<T,VECTOR<int,1> > > rungekutta(phi);
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
template<class T> void LEVELSET_ADVECTION_1D<T>::
Euler_Step_Of_Reinitialization(const ARRAY<T,VECTOR<int,1> >& sign_phi,const T dt,const T time)
{
    GRID<TV>& grid=levelset->grid;
    T_BOUNDARY_SCALAR* boundary=levelset->boundary;
    T_ARRAYS_SCALAR& phi=levelset->phi;
    
    int i;int m=grid.counts.x;T dx=grid.dX.x; 
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,1> > phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,ghost_cells);
    ARRAY<T,VECTOR<int,1> > rhs(1,m);
    
    ARRAY<T,VECTOR<int,1> > phi_1d_x(1-ghost_cells,m+ghost_cells),phix_minus(1,m),phix_plus(1,m);
    for(i=1-ghost_cells;i<=m+ghost_cells;i++) phi_1d_x(i)=phi_ghost(i);
    if(reinitialization_spatial_order == 5) HJ_WENO(m,dx,phi_1d_x,phix_minus,phix_plus);
    else HJ_ENO(reinitialization_spatial_order,m,dx,phi_1d_x,phix_minus,phix_plus);
    for(i=1;i<=m;i++)
        if(LEVELSET_UTILITIES<T>::Sign(phi(i)) < 0) rhs(i)=sqr(max(-phix_minus(i),phix_plus(i),(T)0));
        else rhs(i)=sqr(max(phix_minus(i),-phix_plus(i),(T)0));

    for(i=1;i<=m;i++){
        phi(i)-=dt*sign_phi(i)*(sqrt(rhs(i))-1);
        if(LEVELSET_UTILITIES<T>::Interface(phi_ghost(i),phi(i))) phi(i)=LEVELSET_UTILITIES<T>::Sign(phi_ghost(i))*levelset->small_number*dx;} 
       
    boundary->Apply_Boundary_Condition(grid,phi,time); // time not incremented - pseudo-time
}
//#####################################################################
template class LEVELSET_ADVECTION_1D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_ADVECTION_1D<double>;
#endif

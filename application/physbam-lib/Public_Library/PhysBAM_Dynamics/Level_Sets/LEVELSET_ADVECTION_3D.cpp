//#####################################################################
// Copyright 2002-2005, Doug Enright, Ronald Fedkiw, Geoffrey Irving, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_3D.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T> void LEVELSET_ADVECTION_3D<T>::
Euler_Step(const ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocity,const T dt,const T time,const int number_of_ghost_cells)
{
    GRID<TV>& grid=levelset->grid;
    Euler_Step_Subset(face_velocity,1,grid.counts.x,1,grid.counts.y,1,grid.counts.z,dt,time,number_of_ghost_cells);
}
//#####################################################################
// Function Euler_Step_Cell
//#####################################################################
template<class T> void LEVELSET_ADVECTION_3D<T>::
Euler_Step_Cell(const ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocity,int i,int j,int ij,const T dt,const T time,const int number_of_ghost_cells)
{       
    Euler_Step_Subset(face_velocity,i,i,j,j,ij,ij,dt,time,number_of_ghost_cells);
}
//#####################################################################
// Function Euler_Step_Subset
//#####################################################################
template<class T> void LEVELSET_ADVECTION_3D<T>::
Euler_Step_Subset(const ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocity,int m_start,int m_end,int n_start,int n_end,int mn_start,int mn_end,const T dt,const T time,const int number_of_ghost_cells)
{       
    DEBUG_UTILITIES::Debug_Breakpoint();
    
    GRID<TV>& grid=levelset->grid;
    T_BOUNDARY_SCALAR* boundary=levelset->boundary;
    T_ARRAYS_SCALAR& phi=levelset->phi;
    
    assert(grid.Is_MAC_Grid() && advection);
    assert(grid.Is_MAC_Grid());
    ARRAY<T,VECTOR<int,3> > phi_ghost(grid.Domain_Indices(number_of_ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,number_of_ghost_cells);
    int z=phi_ghost.counts.z,yz=phi_ghost.counts.y*z;

    if(levelset->curvature_motion){ // do curvature first - based on phi^n
        T one_over_two_dx=1/(2*grid.dX.x),one_over_two_dy=1/(2*grid.dX.y),one_over_two_dz=1/(2*grid.dX.z);
        bool curvature_defined=levelset->curvature!=0;levelset->Compute_Curvature(time);
        TV_INT i;
        for(i.x=m_start;i.x<=m_end;i.x++) for(i.y=n_start;i.y<=n_end;i.y++) for(i.z=mn_start;i.z<=mn_end;i.z++){
            int index=phi_ghost.Standard_Index(i);
            T phix=(phi_ghost.array(index+yz)-phi_ghost.array(index-yz))*one_over_two_dx,
              phiy=(phi_ghost.array(index+z)-phi_ghost.array(index-z))*one_over_two_dy,
              phiz=(phi_ghost.array(index+1)-phi_ghost.array(index-1))*one_over_two_dz;
            phi(i)-=dt*levelset->sigma*(*levelset->curvature)(i)*sqrt(sqr(phix)+sqr(phiy)+sqr(phiz));}
        boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,number_of_ghost_cells);
        if(!curvature_defined){delete levelset->curvature;levelset->curvature=0;}}

    advection->Update_Advection_Equation_Cell(grid,phi,phi_ghost,face_velocity,*boundary,dt,time);
    boundary->Apply_Boundary_Condition(grid,phi,time+dt);
}
//#####################################################################
// Functions Reinitialize
//#####################################################################
template<class T> void LEVELSET_ADVECTION_3D<T>::
Reinitialize(const int time_steps,const T time)
{
    GRID<TV>& grid=levelset->grid;
    T_ARRAYS_SCALAR& phi=levelset->phi;

    int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z;

    ARRAY<T,VECTOR<int,3> > sign_phi(1,m,1,n,1,mn); // smeared out sign function
    T epsilon=sqr(grid.dX.Max());
    TV_INT i;
    for(i.x=1;i.x<=m;i.x++) for(i.y=1;i.y<=n;i.y++) for(i.z=1;i.z<=mn;i.z++) sign_phi(i)=phi(i)/sqrt(sqr(phi(i))+epsilon);

    T dt=reinitialization_cfl*grid.min_dX;
    RUNGEKUTTA<ARRAY<T,VECTOR<int,3> > > rungekutta(phi);
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
template<class T> void LEVELSET_ADVECTION_3D<T>::
Euler_Step_Of_Reinitialization(const ARRAY<T,VECTOR<int,3> >& sign_phi,const T dt,const T time)
{
    GRID<TV>& grid=levelset->grid;
    T_BOUNDARY_SCALAR* boundary=levelset->boundary;
    T_ARRAYS_SCALAR& phi=levelset->phi;
    
    TV_INT i;int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z;T dx=grid.dX.x,dy=grid.dX.y,dz=grid.dX.z;
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,3> > phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,ghost_cells);
    ARRAY<T,VECTOR<int,3> > rhs(1,m,1,n,1,mn);

    ARRAY<T,VECTOR<int,1> > phi_1d_x(1-ghost_cells,m+ghost_cells),phix_minus(1,m),phix_plus(1,m);
    for(i.y=1;i.y<=n;i.y++) for(i.z=1;i.z<=mn;i.z++){
        for(i.x=1-ghost_cells;i.x<=m+ghost_cells;i.x++) phi_1d_x(i.x)=phi_ghost(i);
        if(reinitialization_spatial_order == 5) HJ_WENO(m,dx,phi_1d_x,phix_minus,phix_plus);
        else HJ_ENO(reinitialization_spatial_order,m,dx,phi_1d_x,phix_minus,phix_plus);
        for(i.x=1;i.x<=m;i.x++)
            if(LEVELSET_UTILITIES<T>::Sign(phi(i)) < 0) rhs(i)=sqr(max(-phix_minus(i.x),phix_plus(i.x),(T)0));
            else rhs(i)=sqr(max(phix_minus(i.x),-phix_plus(i.x),(T)0));}

    ARRAY<T,VECTOR<int,1> > phi_1d_y(1-ghost_cells,n+ghost_cells),phiy_minus(1,n),phiy_plus(1,n);
    for(i.x=1;i.x<=m;i.x++) for(i.z=1;i.z<=mn;i.z++){
        for(i.y=1-ghost_cells;i.y<=n+ghost_cells;i.y++) phi_1d_y(i.y)=phi_ghost(i);
        if(reinitialization_spatial_order == 5) HJ_WENO(n,dy,phi_1d_y,phiy_minus,phiy_plus);
        else HJ_ENO(reinitialization_spatial_order,n,dy,phi_1d_y,phiy_minus,phiy_plus);
        for(i.y=1;i.y<=n;i.y++)
            if(LEVELSET_UTILITIES<T>::Sign(phi(i)) < 0) rhs(i)+=sqr(max(-phiy_minus(i.y),phiy_plus(i.y),(T)0));
            else rhs(i)+=sqr(max(phiy_minus(i.y),-phiy_plus(i.y),(T)0));}

    ARRAY<T,VECTOR<int,1> > phi_1d_z(1-ghost_cells,mn+ghost_cells),phiz_minus(1,mn),phiz_plus(1,mn);
    for(i.x=1;i.x<=m;i.x++) for(i.y=1;i.y<=n;i.y++){
        for(i.z=1-ghost_cells;i.z<=mn+ghost_cells;i.z++) phi_1d_z(i.z)=phi_ghost(i);
        if(reinitialization_spatial_order == 5) HJ_WENO(mn,dz,phi_1d_z,phiz_minus,phiz_plus);
        else HJ_ENO(reinitialization_spatial_order,mn,dz,phi_1d_z,phiz_minus,phiz_plus);
        for(i.z=1;i.z<=mn;i.z++)
            if(LEVELSET_UTILITIES<T>::Sign(phi(i)) < 0) rhs(i)+=sqr(max(-phiz_minus(i.z),phiz_plus(i.z),(T)0));
            else rhs(i)+=sqr(max(phiz_minus(i.z),-phiz_plus(i.z),(T)0));}

    for(i.x=1;i.x<=m;i.x++) for(i.y=1;i.y<=n;i.y++) for(i.z=1;i.z<=mn;i.z++){
        phi(i)-=dt*sign_phi(i)*(sqrt(rhs(i))-1);
        if(LEVELSET_UTILITIES<T>::Interface(phi_ghost(i),phi(i))) phi(i)=LEVELSET_UTILITIES<T>::Sign(phi_ghost(i))*levelset->small_number*grid.min_dX;}

    boundary->Apply_Boundary_Condition(grid,phi,time); // time not incremented, pseudo-time
}
//#####################################################################
template class LEVELSET_ADVECTION_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_ADVECTION_3D<double>;
#endif

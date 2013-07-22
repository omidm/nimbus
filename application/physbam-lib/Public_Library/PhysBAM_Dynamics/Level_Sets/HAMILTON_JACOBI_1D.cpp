//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/HAMILTON_JACOBI_1D.h>
#include <PhysBAM_Dynamics/Level_Sets/HAMILTONIAN_1D.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
// providing time is optional
template<class T> void HAMILTON_JACOBI_1D<T>::
Euler_Step(const T dt,const T time)
{       
    int i;int m=grid.counts.x;
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,1> > phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,ghost_cells);
        
    // find phx_plus and phix_minus
    ARRAY<T,VECTOR<int,1> > phix_minus(1,m),phix_plus(1,m);
    Calculate_Derivatives(phi_ghost,phix_minus,phix_plus);
        
    if(LF_viscosity){
        T phix_min=min(phix_minus(1),phix_plus(1)),phix_max=max(phix_minus(1),phix_plus(1));
        for(i=1;i<=m;i++){phix_min=min(phix_min,phix_minus(i),phix_plus(i));phix_max=max(phix_max,phix_minus(i),phix_plus(i));}
        for(i=1;i<=m;i++){
            T phix_ave=(phix_minus(i)+phix_plus(i))/2,phix_difference=(phix_plus(i)-phix_minus(i))/2;
            phi(i)-=dt*(hamiltonian.H(phix_ave,i,time)-hamiltonian.Maxabs_H1(phix_min,phix_max,i,time)*phix_difference);}}
    else // LLF and LLLF are the same in 1D
        for(i=1;i<=m;i++){
            T phix_min=min(phix_minus(i),phix_plus(i)),phix_max=max(phix_minus(i),phix_plus(i));
            T phix_ave=(phix_minus(i)+phix_plus(i))/2,phix_difference=(phix_plus(i)-phix_minus(i))/2;
            phi(i)-=dt*(hamiltonian.H(phix_ave,i,time)-hamiltonian.Maxabs_H1(phix_min,phix_max,i,time)*phix_difference);}
                        
    boundary->Apply_Boundary_Condition(grid,phi,time+dt); 
}
//#####################################################################
// Function Calculate_Derivatives
//#####################################################################
template<class T> void HAMILTON_JACOBI_1D<T>::
Calculate_Derivatives(ARRAY<T,VECTOR<int,1> >& phi_ghost,ARRAY<T,VECTOR<int,1> >& phix_minus,ARRAY<T,VECTOR<int,1> >& phix_plus)
{
    int m=grid.counts.x;T dx=grid.dX.x;
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,1> > phi_1d_x(1-ghost_cells,m+ghost_cells);
    for(int i=1-ghost_cells;i<=m+ghost_cells;i++) phi_1d_x(i)=phi_ghost(i);
    if(spatial_order == 5) HJ_WENO(m,dx,phi_1d_x,phix_minus,phix_plus);
    else HJ_ENO(spatial_order,m,dx,phi_1d_x,phix_minus,phix_plus);
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T> T HAMILTON_JACOBI_1D<T>::
CFL(const T time)
{
    int i; int m=grid.counts.x;T dx=grid.dX.x;
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,1> > phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
        
    // find phx_plus and phix_minus
    ARRAY<T,VECTOR<int,1> > phix_minus(1,m),phix_plus(1,m);
    Calculate_Derivatives(phi_ghost,phix_minus,phix_plus);
    
    T maxabs_H1;
    if(LF_viscosity){
        T phix_min=min(phix_minus(1),phix_plus(1)),phix_max=max(phix_minus(1),phix_plus(1));
        for(i=1;i<=m;i++){phix_min=min(phix_min,phix_minus(i),phix_plus(i));phix_max=max(phix_max,phix_minus(i),phix_plus(i));}
        maxabs_H1=hamiltonian.Maxabs_H1(phix_min,phix_max,1,time);
        for(i=1;i<=m;i++) maxabs_H1=max(maxabs_H1,hamiltonian.Maxabs_H1(phix_min,phix_max,i,time));}
    else{ // LLF and LLLF are the same in 1D
        T phix_min=min(phix_minus(1),phix_plus(1)),phix_max=max(phix_minus(1),phix_plus(1));
        maxabs_H1=hamiltonian.Maxabs_H1(phix_min,phix_max,1,time);
        for(i=1;i<=m;i++){
            phix_min=min(phix_minus(i),phix_plus(i)),phix_max=max(phix_minus(i),phix_plus(i));
            maxabs_H1=max(maxabs_H1,hamiltonian.Maxabs_H1(phix_min,phix_max,i,time));}}

    T dt_convect=maxabs_H1/dx;
    dt_convect=max(dt_convect,1/max_time_step); // avoids division by zero
    return 1/dt_convect;
}
//#####################################################################
template class HAMILTON_JACOBI_1D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class HAMILTON_JACOBI_1D<double>;
#endif

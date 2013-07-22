//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Dynamics/Level_Sets/HAMILTON_JACOBI_2D.h>
#include <PhysBAM_Dynamics/Level_Sets/HAMILTONIAN_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
// providing time is optional
template<class T> void HAMILTON_JACOBI_2D<T>::
Euler_Step(const T dt,const T time)
{       
    int i,j;int m=grid.counts.x,n=grid.counts.y;T dx=grid.dX.x,dy=grid.dX.y;
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,2> > phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,ghost_cells);
        
    // do curvature first - based on phi^n
    if(curvature_motion){ 
        bool curvature_defined=(curvature!=0);Compute_Curvature(time);
        for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){
            T phix=(phi_ghost(i+1,j)-phi_ghost(i-1,j))/(2*dx);
            T phiy=(phi_ghost(i,j+1)-phi_ghost(i,j-1))/(2*dy);
            phi(i,j)-=dt*sigma*(*curvature)(i,j)*sqrt(sqr(phix)+sqr(phiy));}
        if(!curvature_defined){delete curvature;curvature=0;}}
    
    // find phx_plus and phix_minus
    ARRAY<T,VECTOR<int,2> > phix_minus(1,m,1,n),phix_plus(1,m,1,n),phiy_minus(1,m,1,n),phiy_plus(1,m,1,n);
    Calculate_Derivatives(phi_ghost,phix_minus,phix_plus,phiy_minus,phiy_plus);

    if(LF_viscosity){
        T phix_min=min(phix_minus(1,1),phix_plus(1,1)),phix_max=max(phix_minus(1,1),phix_plus(1,1)),
                    phiy_min=min(phiy_minus(1,1),phiy_plus(1,1)),phiy_max=max(phiy_minus(1,1),phiy_plus(1,1));
        for(i=1;i<=m;i++) for(j=1;j<=n;j++){
            phix_min=min(phix_min,phix_minus(i,j),phix_plus(i,j));phix_max=max(phix_max,phix_minus(i,j),phix_plus(i,j));
            phiy_min=min(phiy_min,phiy_minus(i,j),phiy_plus(i,j));phiy_max=max(phiy_max,phiy_minus(i,j),phiy_plus(i,j));}
        for(i=1;i<=m;i++) for(j=1;j<=n;j++){
            T phix_ave=(phix_minus(i,j)+phix_plus(i,j))/2,phix_difference=(phix_plus(i,j)-phix_minus(i,j))/2,
                        phiy_ave=(phiy_minus(i,j)+phiy_plus(i,j))/2,phiy_difference=(phiy_plus(i,j)-phiy_minus(i,j))/2;
            phi(i,j)-=dt*(hamiltonian.H(phix_ave,phiy_ave,i,j,time)
                               -hamiltonian.Maxabs_H1(phix_min,phix_max,phiy_min,phiy_max,i,j,time)*phix_difference
                               -hamiltonian.Maxabs_H2(phix_min,phix_max,phiy_min,phiy_max,i,j,time)*phiy_difference);}}
    else if(LLF_viscosity){
        T phix_min=min(phix_minus(1,1),phix_plus(1,1)),phix_max=max(phix_minus(1,1),phix_plus(1,1)),
                    phiy_min=min(phiy_minus(1,1),phiy_plus(1,1)),phiy_max=max(phiy_minus(1,1),phiy_plus(1,1));
        for(i=1;i<=m;i++) for(j=1;j<=n;j++){
            phix_min=min(phix_min,phix_minus(i,j),phix_plus(i,j));phix_max=max(phix_max,phix_minus(i,j),phix_plus(i,j));
            phiy_min=min(phiy_min,phiy_minus(i,j),phiy_plus(i,j));phiy_max=max(phiy_max,phiy_minus(i,j),phiy_plus(i,j));}
        for(i=1;i<=m;i++) for(j=1;j<=n;j++){
            T phix_min_local=min(phix_minus(i,j),phix_plus(i,j)),phix_max_local=max(phix_minus(i,j),phix_plus(i,j)),
                        phiy_min_local=min(phiy_minus(i,j),phiy_plus(i,j)),phiy_max_local=max(phiy_minus(i,j),phiy_plus(i,j));
            T phix_ave=(phix_minus(i,j)+phix_plus(i,j))/2,phix_difference=(phix_plus(i,j)-phix_minus(i,j))/2,
                        phiy_ave=(phiy_minus(i,j)+phiy_plus(i,j))/2,phiy_difference=(phiy_plus(i,j)-phiy_minus(i,j))/2;
            phi(i,j)-=dt*(hamiltonian.H(phix_ave,phiy_ave,i,j,time)
                               -hamiltonian.Maxabs_H1(phix_min_local,phix_max_local,phiy_min,phiy_max,i,j,time)*phix_difference
                               -hamiltonian.Maxabs_H2(phix_min,phix_max,phiy_min_local,phiy_max_local,i,j,time)*phiy_difference);}}
    else if(LLLF_viscosity)
        for(i=1;i<=m;i++) for(j=1;j<=n;j++){
            T phix_min_local=min(phix_minus(i,j),phix_plus(i,j)),phix_max_local=max(phix_minus(i,j),phix_plus(i,j)),
                        phiy_min_local=min(phiy_minus(i,j),phiy_plus(i,j)),phiy_max_local=max(phiy_minus(i,j),phiy_plus(i,j));
            T phix_ave=(phix_minus(i,j)+phix_plus(i,j))/2,phix_difference=(phix_plus(i,j)-phix_minus(i,j))/2,
                        phiy_ave=(phiy_minus(i,j)+phiy_plus(i,j))/2,phiy_difference=(phiy_plus(i,j)-phiy_minus(i,j))/2;
            phi(i,j)-=dt*(hamiltonian.H(phix_ave,phiy_ave,i,j,time)
                               -hamiltonian.Maxabs_H1(phix_min_local,phix_max_local,phiy_min_local,phiy_max_local,i,j,time)*phix_difference
                               -hamiltonian.Maxabs_H2(phix_min_local,phix_max_local,phiy_min_local,phiy_max_local,i,j,time)*phiy_difference);}
                        
    boundary->Apply_Boundary_Condition(grid,phi,time+dt); 
}
//#####################################################################
// Function Calculate_Derivatives
//#####################################################################
template<class T> void HAMILTON_JACOBI_2D<T>::
Calculate_Derivatives(ARRAY<T,VECTOR<int,2> >& phi_ghost,ARRAY<T,VECTOR<int,2> >& phix_minus,ARRAY<T,VECTOR<int,2> >& phix_plus,ARRAY<T,VECTOR<int,2> >& phiy_minus,ARRAY<T,VECTOR<int,2> >& phiy_plus)
{
    int i,j;int m=grid.counts.x,n=grid.counts.y;T dx=grid.dX.x,dy=grid.dX.y;
    int ghost_cells=3;
    
    // x-direction
    ARRAY<T,VECTOR<int,1> > phi_1d_x(1-ghost_cells,m+ghost_cells),phix_minus_1d(1,m),phix_plus_1d(1,m); 
    for(j=1;j<=n;j++){
        for(i=1-ghost_cells;i<=m+ghost_cells;i++) phi_1d_x(i)=phi_ghost(i,j);
        if(spatial_order == 5) HJ_WENO(m,dx,phi_1d_x,phix_minus_1d,phix_plus_1d);
        else HJ_ENO(spatial_order,m,dx,phi_1d_x,phix_minus_1d,phix_plus_1d);
        for(i=1;i<=m;i++){phix_minus(i,j)=phix_minus_1d(i);phix_plus(i,j)=phix_plus_1d(i);}}

    // y-direction
    ARRAY<T,VECTOR<int,1> > phi_1d_y(1-ghost_cells,n+ghost_cells),phiy_minus_1d(1,n),phiy_plus_1d(1,n); 
    for(i=1;i<=m;i++){
        for(j=1-ghost_cells;j<=n+ghost_cells;j++) phi_1d_y(j)=phi_ghost(i,j);
        if(spatial_order == 5) HJ_WENO(n,dy,phi_1d_y,phiy_minus_1d,phiy_plus_1d);
        else HJ_ENO(spatial_order,n,dy,phi_1d_y,phiy_minus_1d,phiy_plus_1d);
        for(j=1;j<=n;j++){phiy_minus(i,j)=phiy_minus_1d(j);phiy_plus(i,j)=phiy_plus_1d(j);}}
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T> T HAMILTON_JACOBI_2D<T>::
CFL(const T time)
{
    int i,j;int m=grid.counts.x,n=grid.counts.y;T dx=grid.dX.x,dy=grid.dX.y;
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,2> > phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
        
    // find phx_plus and phix_minus
    ARRAY<T,VECTOR<int,2> > phix_minus(1,m,1,n),phix_plus(1,m,1,n),phiy_minus(1,m,1,n),phiy_plus(1,m,1,n);
    Calculate_Derivatives(phi_ghost,phix_minus,phix_plus,phiy_minus,phiy_plus);
    
    T maxabs_H1=0,maxabs_H2=0;
    if(LF_viscosity){
        T phix_min=min(phix_minus(1,1),phix_plus(1,1)),phix_max=max(phix_minus(1,1),phix_plus(1,1)),
                    phiy_min=min(phiy_minus(1,1),phiy_plus(1,1)),phiy_max=max(phiy_minus(1,1),phiy_plus(1,1));
        for(i=1;i<=m;i++) for(j=1;j<=n;j++){
            phix_min=min(phix_min,phix_minus(i,j),phix_plus(i,j));phix_max=max(phix_max,phix_minus(i,j),phix_plus(i,j));
            phiy_min=min(phiy_min,phiy_minus(i,j),phiy_plus(i,j));phiy_max=max(phiy_max,phiy_minus(i,j),phiy_plus(i,j));}
        maxabs_H1=hamiltonian.Maxabs_H1(phix_min,phix_max,phiy_min,phiy_max,1,1,time);
        maxabs_H2=hamiltonian.Maxabs_H2(phix_min,phix_max,phiy_min,phiy_max,1,1,time);
        for(i=1;i<=m;i++) for(j=1;j<=n;j++){
            maxabs_H1=max(maxabs_H1,hamiltonian.Maxabs_H1(phix_min,phix_max,phiy_min,phiy_max,i,j,time));
            maxabs_H2=max(maxabs_H2,hamiltonian.Maxabs_H2(phix_min,phix_max,phiy_min,phiy_max,i,j,time));}}
    else if(LLF_viscosity){
        T phix_min=min(phix_minus(1,1),phix_plus(1,1)),phix_max=max(phix_minus(1,1),phix_plus(1,1)),
                    phiy_min=min(phiy_minus(1,1),phiy_plus(1,1)),phiy_max=max(phiy_minus(1,1),phiy_plus(1,1));
        for(i=1;i<=m;i++) for(j=1;j<=n;j++){
            phix_min=min(phix_min,phix_minus(i,j),phix_plus(i,j));phix_max=max(phix_max,phix_minus(i,j),phix_plus(i,j));
            phiy_min=min(phiy_min,phiy_minus(i,j),phiy_plus(i,j));phiy_max=max(phiy_max,phiy_minus(i,j),phiy_plus(i,j));}
        maxabs_H1=hamiltonian.Maxabs_H1(phix_minus(1,1),phix_plus(1,1),phiy_min,phiy_max,1,1,time);
        maxabs_H2=hamiltonian.Maxabs_H2(phix_min,phix_max,phiy_minus(1,1),phiy_plus(1,1),1,1,time);
        for(i=1;i<=m;i++) for(j=1;j<=n;j++){
            T phix_min_local=min(phix_minus(i,j),phix_plus(i,j)),phix_max_local=max(phix_minus(i,j),phix_plus(i,j)),
                        phiy_min_local=min(phiy_minus(i,j),phiy_plus(i,j)),phiy_max_local=max(phiy_minus(i,j),phiy_plus(i,j));
            maxabs_H1=max(maxabs_H1,hamiltonian.Maxabs_H1(phix_min_local,phix_max_local,phiy_min,phiy_max,i,j,time));
            maxabs_H2=max(maxabs_H2,hamiltonian.Maxabs_H2(phix_min,phix_max,phiy_min_local,phiy_max_local,i,j,time));}}
    else if(LLLF_viscosity){
        maxabs_H1=hamiltonian.Maxabs_H1(phix_minus(1,1),phix_plus(1,1),phiy_minus(1,1),phiy_plus(1,1),1,1,time);
        maxabs_H2=hamiltonian.Maxabs_H2(phix_minus(1,1),phix_plus(1,1),phiy_minus(1,1),phiy_plus(1,1),1,1,time);
        for(i=1;i<=m;i++) for(j=1;j<=n;j++){
            T phix_min_local=min(phix_minus(i,j),phix_plus(i,j)),phix_max_local=max(phix_minus(i,j),phix_plus(i,j)),
                        phiy_min_local=min(phiy_minus(i,j),phiy_plus(i,j)),phiy_max_local=max(phiy_minus(i,j),phiy_plus(i,j));
            maxabs_H1=max(maxabs_H1,hamiltonian.Maxabs_H1(phix_min_local,phix_max_local,phiy_min_local,phiy_max_local,i,j,time));
            maxabs_H2=max(maxabs_H2,hamiltonian.Maxabs_H2(phix_min_local,phix_max_local,phiy_min_local,phiy_max_local,i,j,time));}}

    T dt_convection=maxabs_H1/dx+maxabs_H2/dy;
    T dt_curvature=0;if(curvature_motion) dt_curvature=sigma*(2/sqr(dx)+2/sqr(dy)); // note that sigma is negative
    T dt_overall=dt_convection+dt_curvature;
    return 1/max(dt_overall,1/max_time_step);
}
//#####################################################################
template class HAMILTON_JACOBI_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class HAMILTON_JACOBI_2D<double>;
#endif

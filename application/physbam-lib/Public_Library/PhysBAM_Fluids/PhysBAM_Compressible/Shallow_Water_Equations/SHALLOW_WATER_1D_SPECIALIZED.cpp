//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Eran Guendelman, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_1D_SPECIALIZED  
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Shallow_Water_Equations/SHALLOW_WATER_1D_SPECIALIZED.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T> void SHALLOW_WATER_1D_SPECIALIZED<T>::
Euler_Step(const T dt,const T time)
{   
    int m=grid.counts.x;
    int ghost_cells=3;

    ARRAY<TV_DIMENSION,VECTOR<int,1> > U_ghost(1-ghost_cells,m+ghost_cells);boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time,ghost_cells);
    ARRAY<bool,VECTOR<int,1> > psi(1,m,false);psi.Fill(true); // no cut out grids

    static ARRAY<bool,VECTOR<int,1> > zero_height(1-ghost_cells,m+ghost_cells);zero_height.Fill(false);
    for(int i=1-ghost_cells;i<=grid.counts.x+ghost_cells;i++) if(U_ghost(i)(1)<=min_height) zero_height(i)=true;

    ARRAY<T,VECTOR<int,1> > ground_ghost(1-ghost_cells,grid.counts.x+ghost_cells);
    if(ground) BOUNDARY_UNIFORM<GRID<TV>,T>().Fill_Ghost_Cells(grid,*ground,ground_ghost,dt,time,ghost_cells);
    for(int i=1-ghost_cells;i<=grid.counts.x+ghost_cells;i++) eta_ghost(i)=U_ghost(i)(1)+ground_ghost(i);

    T_FACE_ARRAYS_BOOL psi_N(grid.Get_MAC_Grid_At_Regular_Positions());
    T_FACE_ARRAYS_SCALAR face_velocities(grid.Get_MAC_Grid_At_Regular_Positions());
    conservation->Save_Fluxes();
    VECTOR<EIGENSYSTEM<T,VECTOR<T,2> >*,1> eigensystem(&eigensystem_F);
    conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystem,eigensystem,psi_N,face_velocities);

    ARRAY_VIEW<TV_DIMENSION,VECTOR<int,1> >& old_flux=conservation->fluxes.Component(1);
    for(int i=0;i<=grid.counts.x;i++){
        bool zero_left=zero_height(i),zero_right=zero_height(i+1),update_flux=false;
        ARRAY<T,VECTOR<int,1> > new_flux_left(1,2),new_flux_right(1,2);new_flux_left(1)=new_flux_right(1)=old_flux(i+1)(1);new_flux_left(2)=new_flux_right(2)=old_flux(i+1)(2);

        if(zero_left && zero_right){ // no mass flux between empty cells
            new_flux_left(1)=new_flux_right(1)=0;
            update_flux=true;}
        else if(eta_ghost(i) < ground_ghost(i+1) || (zero_left && !zero_right)){ // cliff to right OR only flow to left via wetting
            T extra_momentum_flux=0;
            if(U_ghost(i+1)(2)<0 && !zero_right){ // water flowing off of cliff
                new_flux_left(1)=new_flux_right(1)=U_ghost(i+1)(1)*U_ghost(i+1)(2);
                extra_momentum_flux=(T).5*sqr(U_ghost(i+1)(2));}
            else new_flux_left(1)=new_flux_right(1)=0;
            new_flux_left(2)=gravity*eta_ghost(i)+extra_momentum_flux;
            new_flux_right(2)=gravity*(T).5*(eta_ghost(i)+eta_ghost(i+1))+extra_momentum_flux;
            update_flux=true;}
        else if(ground_ghost(i) > eta_ghost(i+1) || (!zero_left && zero_right)){ // cliff to left OR only flow right via wetting
            T extra_momentum_flux=0;
            if(U_ghost(i)(2)>0 && !zero_left){ // water flowing off of cliff
                new_flux_left(1)=new_flux_right(1)=U_ghost(i)(1)*U_ghost(i)(2);
                extra_momentum_flux=(T).5*sqr(U_ghost(i)(2));}
            else new_flux_left(1)=new_flux_right(1)=0;
            new_flux_right(2)=gravity*eta_ghost(i+1)+extra_momentum_flux;
            new_flux_left(2)=gravity*(T).5*(eta_ghost(i)+eta_ghost(i+1))+extra_momentum_flux;
            update_flux=true;}

        if(update_flux){ 
            if(i>0) for(int k=1;k<=d;k++) U(i)(k)+=dt*(old_flux(i+1)(k)-new_flux_left(k))*grid.one_over_dX.x; // update cell on left
            if(i<grid.counts.x) for(int k=1;k<=d;k++) U(i+1)(k)-=dt*(old_flux(i+1)(k)-new_flux_right(k))*grid.one_over_dX.x;} // update cell on right

        // for debugging
        postprocessed_flux(i)(1)=new_flux_left(1);postprocessed_flux(i)(2)=new_flux_left(2);postprocessed_flux(i)(3)=new_flux_right(1);postprocessed_flux(i)(4)=new_flux_right(2);}

    T two_min_height=(T)2.01*min_height;
    for(int i=1;i<=grid.counts.x;i++){
        if(U_ghost(i)(1)<=two_min_height) U(i)(2)=0; // correct for where we have zero fluxes due to small average h
        if(U(i)(1)<0) U(i)(1)=0;}

    boundary->Apply_Boundary_Condition(grid,U,time+dt); 
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T> T SHALLOW_WATER_1D_SPECIALIZED<T>::
CFL()
{
    T max_speed=0;
    for(int i=1;i<=grid.counts.x;i++){
        T u=U(i)(2),celerity=sqrt(gravity*U(i)(1));
        max_speed=max(max_speed,abs(u)+celerity);}
    T dt_convect=max_speed*grid.one_over_dX.x;
    return 1/dt_convect;
}
//#####################################################################
template class SHALLOW_WATER_1D_SPECIALIZED<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SHALLOW_WATER_1D_SPECIALIZED<double>;
#endif

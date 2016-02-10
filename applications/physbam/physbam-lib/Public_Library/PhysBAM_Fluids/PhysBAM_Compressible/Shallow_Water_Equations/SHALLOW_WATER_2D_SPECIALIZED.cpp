//#####################################################################
// Copyright 2002-2004, Ron Fedkiw, Eran Guendelman, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_2D_SPECIALIZED  
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Shallow_Water_Equations/SHALLOW_WATER_2D_SPECIALIZED.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T> void SHALLOW_WATER_2D_SPECIALIZED<T>::
Euler_Step(const T dt,const T time)
{  
    int m=grid.counts.x,n=grid.counts.y;
    int ghost_cells=3;

    ARRAY<TV_DIMENSION,VECTOR<int,2> > U_ghost(1-ghost_cells,m+ghost_cells,1-ghost_cells,n+ghost_cells);boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time,ghost_cells);
    ARRAY<bool,VECTOR<int,2> > psi(1,m,1,n);psi.Fill(true); // no cut out grids

    static ARRAY<bool,VECTOR<int,2> > zero_height(1-ghost_cells,m+ghost_cells,1-ghost_cells,n+ghost_cells);zero_height.Fill(false);
    for(int i=1-ghost_cells;i<=grid.counts.x+ghost_cells;i++) for(int j=1-ghost_cells;j<=grid.counts.y+ghost_cells;j++) if(U_ghost(i,j)(1)<=min_height) zero_height(i,j)=true;

    for(int i=1-ghost_cells;i<=grid.counts.x+ghost_cells;i++) for(int j=1-ghost_cells;j<=grid.counts.y+ghost_cells;j++) eta_ghost(i,j)=U_ghost(i,j)(1)+(*ground_ghost)(i,j);

    internal_conservation.Save_Fluxes();conservation->Save_Fluxes();
    internal_conservation.Update_Conservation_Law_For_Specialized_Shallow_Water_Equations(grid,U,U_ghost,psi,dt,eigensystem_F,eigensystem_G,*conservation);

    if(epsilon_u || epsilon_v) for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++){
        T u_yy=sqr(grid.one_over_dX.y)*(U_ghost(i,j+1)(2)-2*U_ghost(i,j)(2)+U_ghost(i,j-1)(2));
        T v_xx=sqr(grid.one_over_dX.x)*(U_ghost(i+1,j)(3)-2*U_ghost(i,j)(3)+U_ghost(i-1,j)(3));
        U(i,j)(2)+=dt*epsilon_u*abs(U(i,j)(3))*u_yy;
        U(i,j)(3)+=dt*epsilon_v*abs(U(i,j)(2))*v_xx;}

    // x fluxes
    ARRAY_VIEW<TV_DIMENSION,VECTOR<int,2> >& old_flux_x=internal_conservation.fluxes.Component(1);
    for(int i=0;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++){
        bool zero_left=zero_height(i,j),zero_right=zero_height(i+1,j),update_flux=false;
        ARRAY<T,VECTOR<int,1> > new_flux_left(1,3),new_flux_right(1,3);for(int k=1;k<=3;k++) new_flux_left(k)=new_flux_right(k)=old_flux_x(i+1,j)(k);

        if(zero_left && zero_right){ // no mass flux between empty cells
            new_flux_left(1)=new_flux_right(1)=0;
            update_flux=true;}        
        else if(eta_ghost(i,j) < (*ground_ghost)(i+1,j) || (zero_left && !zero_right)){ // cliff to right OR only flow to left via wetting
            T extra_momentum_flux=0;
            if(U_ghost(i+1,j)(2)<0 && !zero_right){ // water flowing off of cliff
                new_flux_left(1)=new_flux_right(1)=U_ghost(i+1,j)(1)*U_ghost(i+1,j)(2);
                extra_momentum_flux=(T).5*sqr(U_ghost(i+1,j)(2));}
            else new_flux_left(1)=new_flux_right(1)=0;
            new_flux_left(2)=gravity*eta_ghost(i,j)+extra_momentum_flux;
            new_flux_right(2)=gravity*(T).5*(eta_ghost(i,j)+eta_ghost(i+1,j))+extra_momentum_flux;
            update_flux=true;}
        else if((*ground_ghost)(i,j) > eta_ghost(i+1,j) || (!zero_left && zero_right)){ // cliff to left OR only flow right via wetting
            T extra_momentum_flux=0;
            if(U_ghost(i,j)(2)>0 && !zero_left){ // water flowing off of cliff
                new_flux_left(1)=new_flux_right(1)=U_ghost(i,j)(1)*U_ghost(i,j)(2);
                extra_momentum_flux=(T).5*sqr(U_ghost(i,j)(2));}
            else new_flux_left(1)=new_flux_right(1)=0;
            new_flux_right(2)=gravity*eta_ghost(i+1,j)+extra_momentum_flux;
            new_flux_left(2)=gravity*(T).5*(eta_ghost(i,j)+eta_ghost(i+1,j))+extra_momentum_flux;
            update_flux=true;}

        if(update_flux){
            if(i>0) for(int k=1;k<=d;k++) U(i,j)(k)+=dt*(old_flux_x(i+1,j)(k)-new_flux_left(k))*grid.one_over_dX.x; // update cell on left
            if(i<grid.counts.x) for(int k=1;k<=d;k++) U(i+1,j)(k)-=dt*(old_flux_x(i+1,j)(k)-new_flux_right(k))*grid.one_over_dX.x;} // update cell on right

        // for debugging
        postprocessed_flux_x(i,j)(1)=new_flux_left(1);postprocessed_flux_x(i,j)(2)=new_flux_left(2);postprocessed_flux_x(i,j)(3)=new_flux_left(3);
        postprocessed_flux_x(i,j)(4)=new_flux_right(1);postprocessed_flux_x(i,j)(5)=new_flux_right(2);postprocessed_flux_x(i,j)(6)=new_flux_right(3);}

    // y fluxes
    ARRAY_VIEW<TV_DIMENSION,VECTOR<int,2> >& old_flux_y=internal_conservation.fluxes.Component(2);
    for(int i=1;i<=grid.counts.x;i++) for(int j=0;j<=grid.counts.y;j++){
        bool zero_left=zero_height(i,j),zero_right=zero_height(i,j+1),update_flux=false;
        ARRAY<T,VECTOR<int,1> > new_flux_left(1,3),new_flux_right(1,3);for(int k=1;k<=3;k++) new_flux_left(k)=new_flux_right(k)=old_flux_y(i,j+1)(k);
        
        if(zero_left && zero_right){ // no mass flux between empty cells
            new_flux_left(1)=new_flux_right(1)=0;
            update_flux=true;}        
        else if(eta_ghost(i,j) < (*ground_ghost)(i,j+1) || (zero_left && !zero_right)){ // cliff to right OR only flow to left via wetting
            T extra_momentum_flux=0;
            if(U_ghost(i,j+1)(3)<0 && !zero_right){ // water flowing off of cliff
                new_flux_left(1)=new_flux_right(1)=U_ghost(i,j+1)(1)*U_ghost(i,j+1)(3);
                extra_momentum_flux=(T).5*sqr(U_ghost(i,j+1)(3));}
            else new_flux_left(1)=new_flux_right(1)=0;
            new_flux_left(3)=gravity*eta_ghost(i,j)+extra_momentum_flux;
            new_flux_right(3)=gravity*(T).5*(eta_ghost(i,j)+eta_ghost(i,j+1))+extra_momentum_flux;
            update_flux=true;}
        else if((*ground_ghost)(i,j) > eta_ghost(i,j+1) || (!zero_left && zero_right)){ // cliff to left OR only flow right via wetting
            T extra_momentum_flux=0;
            if(U_ghost(i,j)(3)>0 && !zero_left){ // water flowing off of cliff
                new_flux_left(1)=new_flux_right(1)=U_ghost(i,j)(1)*U_ghost(i,j)(3);
                extra_momentum_flux=(T).5*sqr(U_ghost(i,j)(3));}
            else new_flux_left(1)=new_flux_right(1)=0;
            new_flux_right(3)=gravity*eta_ghost(i,j+1)+extra_momentum_flux;
            new_flux_left(3)=gravity*(T).5*(eta_ghost(i,j)+eta_ghost(i,j+1))+extra_momentum_flux;
            update_flux=true;}
        
        if(update_flux){
            if(j>0) for(int k=1;k<=d;k++) U(i,j)(k)+=dt*(old_flux_y(i,j+1)(k)-new_flux_left(k))*grid.one_over_dX.y; // update cell on left
            if(j<grid.counts.y) for(int k=1;k<=d;k++) U(i,j+1)(k)-=dt*(old_flux_y(i,j+1)(k)-new_flux_right(k))*grid.one_over_dX.y;} // update cell on right
        
        // for debugging
        postprocessed_flux_y(i,j)(1)=new_flux_left(1);postprocessed_flux_y(i,j)(2)=new_flux_left(2);postprocessed_flux_y(i,j)(3)=new_flux_left(3);
        postprocessed_flux_y(i,j)(4)=new_flux_right(1);postprocessed_flux_y(i,j)(5)=new_flux_right(2);postprocessed_flux_y(i,j)(6)=new_flux_right(3);}
    
    T two_min_height=(T)2.01*min_height;
    for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++){
        if(U_ghost(i,j)(1)<=two_min_height) U(i,j)(2)=U(i,j)(3)=0; // correct for where we have zero fluxes due to small average h
        if(U(i,j)(1)<0) U(i,j)(1)=0;}

    boundary->Apply_Boundary_Condition(grid,U,time+dt); 
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T> T SHALLOW_WATER_2D_SPECIALIZED<T>::
CFL()
{
    T max_x_speed=0,max_y_speed=0;
    for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++){
        T u=U(i,j)(2),v=U(i,j)(3),celerity=sqrt(gravity*U(i,j)(1));
        max_x_speed=max(max_x_speed,abs(u)+celerity);
        max_y_speed=max(max_y_speed,abs(v)+celerity);}
    T dt_convect=max_x_speed*grid.one_over_dX.x+max_y_speed*grid.one_over_dX.y;
    return 1/dt_convect;
}
//#####################################################################
template class SHALLOW_WATER_2D_SPECIALIZED<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SHALLOW_WATER_2D_SPECIALIZED<double>;
#endif

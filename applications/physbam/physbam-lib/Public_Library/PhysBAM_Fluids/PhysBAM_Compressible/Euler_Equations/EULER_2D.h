//#####################################################################
// Copyright 2002, 2003, Doug Enright, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_2D  
//##################################################################### 
//
// Input an EOS class, that is a class that inherits the virtual base class EOS.
// Input a GRID_2D class.
// Input U as 4 by (1,m) by (1,n) for mass, momentum, and energy.
//
// Use Set_Up_Cut_Out_Grid(psi) to define a cut out grid with psi as (1,m) by (1,n).
// When psi=true, solve the equaitions. 
// When psi=false, do NOT solve the equations.
//
//#####################################################################
#ifndef __EULER_2D__
#define __EULER_2D__

#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_2D_EIGENSYSTEM_F.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_2D_EIGENSYSTEM_G.h>
namespace PhysBAM{

template<class T_input>
class EULER_2D:public EULER<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<T,4> TV_DIMENSION;
protected:
    using EULER<GRID<TV> >::cut_out_grid;using EULER<GRID<TV> >::max_time_step;
public:
    using EULER<GRID<TV> >::boundary;using EULER<GRID<TV> >::conservation;using EULER<GRID<TV> >::eos;

    GRID<TV> grid;
    ARRAY<TV_DIMENSION,VECTOR<int,2> >& U; // mass, momentum, and energy
    ARRAY<bool,VECTOR<int,2> >* psi_pointer; // defines cut out grid
    EULER_2D_EIGENSYSTEM_F<T> eigensystem_F;
    EULER_2D_EIGENSYSTEM_G<T> eigensystem_G;

    EULER_2D(ARRAY<TV_DIMENSION,VECTOR<int,2> >& U_input)
        :U(U_input)
    {}
    
    void Set_Up_Cut_Out_Grid(ARRAY<bool,VECTOR<int,2> >& psi_input)
    {psi_pointer=&psi_input;cut_out_grid=true;}

    void Set_Custom_Equation_Of_State(EOS<T>& eos_input)
    {eigensystem_F.Set_Custom_Equation_Of_State(eos_input);eigensystem_G.Set_Custom_Equation_Of_State(eos_input);EULER<GRID<TV> >::Set_Custom_Equation_Of_State(eos_input);}
    
    void Initialize_Domain(const int m, const int n, const T xmin, const T xmax, const T ymin, const T ymax)
    {grid.Initialize(m,n,xmin,xmax,ymin,ymax);}
    
//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    T CFL();
//#####################################################################
};
//#####################################################################
// Function Euler_Step
//#####################################################################
// providing time is optional
template<class T> void EULER_2D<T>::
Euler_Step(const T dt,const T time)
{  
    int m=grid.m,n=grid.n;
    
    ARRAY<TV_DIMENSION,VECTOR<int,2> > U_ghost(-2,m+3,-2,n+3);
    boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time);
    
    if(cut_out_grid) conservation->Update_Conservation_Law(grid,U,U_ghost,*psi_pointer,dt,eigensystem_F,eigensystem_G);  
    else{ // not a cut out grid
        ARRAY<bool,VECTOR<int,2> > psi(1,m,1,n);psi.Fill(true);
        conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystem_F,eigensystem_G);}

    boundary->Apply_Boundary_Condition(grid,U,time+dt); 
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T> T EULER_2D<T>::
CFL()
{
    int m=grid.m,n=grid.n;T dx=grid.dx,dy=grid.dy;
    
    ARRAY<T,VECTOR<int,2> > u_minus_c(1,m,1,n),u_plus_c(1,m,1,n),v_minus_c(1,m,1,n),v_plus_c(1,m,1,n);
    for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){
        if(!cut_out_grid || (cut_out_grid && (*psi_pointer)(i,j))){
            T u=U(2,i,j)/U(1,i,j),v=U(3,i,j)/U(1,i,j);
            T sound_speed=eos->c(U(1,i,j),e(U(1,i,j),U(2,i,j),U(3,i,j),U(4,i,j)));
            u_minus_c(i,j)=u-sound_speed;u_plus_c(i,j)=u+sound_speed;
            v_minus_c(i,j)=v-sound_speed;v_plus_c(i,j)=v+sound_speed;}}
    T dt_convect=max(u_minus_c.Maxabs(),u_plus_c.Maxabs())/dx+max(v_minus_c.Maxabs(),v_plus_c.Maxabs())/dy;
    dt_convect=max(dt_convect,1/max_time_step);
    return 1/dt_convect;
}           
//#####################################################################     
}
#endif

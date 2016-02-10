//#####################################################################
// Copyright 2002, 2003, Doug Enright, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_1D  
//##################################################################### 
//
// Input an EOS class, that is a class that inherits the virtual base class EOS.
// Input a GRID_1D class.
// Input U as 3 by (1,m) for mass, momentum, and energy.
//
// Use Set_Up_Cut_Out_Grid(psi) to define a cut out grid with psi as (1,m).
// When psi=true, solve the equaitions. 
// When psi=false, do NOT solve the equations.
//
//#####################################################################
#ifndef __EULER_1D__
#define __EULER_1D__

#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_1D_EIGENSYSTEM_F.h>
namespace PhysBAM{

template<class T_input>
class EULER_1D:public EULER<GRID<VECTOR<T_input,1> > >
{
    typedef T_input T;typedef VECTOR<T,1> TV;typedef VECTOR<T,3> TV_DIMENSION;
protected:
    using EULER<GRID<TV> >::cut_out_grid;using EULER<GRID<TV> >::max_time_step;
public:
    using EULER<GRID<TV> >::boundary;using EULER<GRID<TV> >::conservation;using EULER<GRID<TV> >::eos;
    
    GRID<TV> grid;
    ARRAY<TV_DIMENSION,VECTOR<int,1> >& U;         // mass, momentum, and energy
    ARRAY<bool,VECTOR<int,1> >* psi_pointer; // defines cut out grid
    EULER_1D_EIGENSYSTEM_F<T> eigensystem_F;

    EULER_1D(ARRAY<TV_DIMENSION,VECTOR<int,1> >& U_input)
        :U(U_input)
    {}
    
    void Set_Up_Cut_Out_Grid(ARRAY<bool,VECTOR<int,1> >& psi_input)
    {psi_pointer=&psi_input;cut_out_grid=true;}

    void Set_Custom_Equation_Of_State(EOS<T>& eos_input)
    {eigensystem_F.Set_Custom_Equation_Of_State(eos_input);EULER<GRID<TV> >::Set_Custom_Equation_Of_State(eos_input);}
    
    void Initialize_Domain(const int m, const T xmin, const T xmax)
    {grid.Initialize(m,xmin,xmax);}
    
//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    T CFL();
//#####################################################################
};
//#####################################################################
// Function Euler_Step
//#####################################################################
// providing time is optional
template<class T> void EULER_1D<T>::
Euler_Step(const T dt,const T time)
{  
    int m=grid.m;
    
    ARRAY<TV_DIMENSION,VECTOR<int,1> > U_ghost(-2,m+3);
    boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time);
    
    if(cut_out_grid) conservation->Update_Conservation_Law(grid,U,U_ghost,*psi_pointer,dt,eigensystem_F);   
    else{ // not a cut out grid
        ARRAY<bool,VECTOR<int,1> > psi(1,m);psi.Fill(true);
        conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystem_F);}

    boundary->Apply_Boundary_Condition(grid,U,time+dt); 
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T> T EULER_1D<T>::
CFL()
{
    int m=grid.m;T dx=grid.dx;
    
    ARRAY<T,VECTOR<int,1> > u_minus_c(1,m),u_plus_c(1,m);
    for(int i=1;i<=m;i++){
        if(!cut_out_grid || (cut_out_grid && (*psi_pointer)(i))){
            T u=U(i)(2)/U(i)(1);
            T sound_speed=eos->c(U(i)(1),e(U(i)(1),U(i)(2),U(i)(3)));
            u_minus_c(i)=u-sound_speed;u_plus_c(i)=u+sound_speed;}}
    T dt_convect=max(u_minus_c.Maxabs(),u_plus_c.Maxabs())/dx;
    dt_convect=max(dt_convect,1/max_time_step);
    return 1/dt_convect;
}
//#####################################################################
}
#endif

//#####################################################################
// Copyright 2002-2003, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Math_Tools/maxabs.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Reactive_Euler_Equations/REACTIVE_EULER_1D.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
// providing time is optional
template<class T> void REACTIVE_EULER_1D<T>::
Euler_Step(const T dt,const T time)
{  
    int m=grid.counts.x;
    int ghost_cells=3;
    
    ARRAY<TV_DIMENSION,VECTOR<int,1> > U_ghost(1-ghost_cells,m+ghost_cells);
    boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time,ghost_cells);
    
    T_FACE_ARRAYS_BOOL psi_N(grid.Get_MAC_Grid_At_Regular_Positions());
    T_FACE_ARRAYS_SCALAR face_velocities(grid.Get_MAC_Grid_At_Regular_Positions());
    VECTOR<EIGENSYSTEM<T,VECTOR<T,4> >*,1> eigensystem(&eigensystem_F);
    if(cut_out_grid) conservation->Update_Conservation_Law(grid,U,U_ghost,*psi_pointer,dt,eigensystem,eigensystem,psi_N,face_velocities);
    else{ // not a cut out grid
        ARRAY<bool,VECTOR<int,1> > psi(1,m);psi.Fill(1);
        conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystem,eigensystem,psi_N,face_velocities);}

    boundary->Apply_Boundary_Condition(grid,U,time+dt); 
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T> T REACTIVE_EULER_1D<T>::
CFL()
{
    int m=grid.counts.x;T dx=grid.dX.x;
    
    ARRAY<T,VECTOR<int,1> > u_minus_c(1,m),u_plus_c(1,m);
    for(int i=1;i<=m;i++){
        if(!cut_out_grid || (cut_out_grid && (*psi_pointer)(i)==1)){
            T u=U(i)(2)/U(i)(1);
            T Y=U(i)(4)/U(i)(1);
            T sound_speed=eos.c(U(i)(1),e(U(i)(1),U(i)(2),U(i)(3)),Y);
            u_minus_c(i)=u-sound_speed;u_plus_c(i)=u+sound_speed;}}
    T dt_convect=max(u_minus_c.Maxabs(),u_plus_c.Maxabs())/dx;
    return 1/dt_convect;
}
//#####################################################################
template class REACTIVE_EULER_1D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class REACTIVE_EULER_1D<double>;
#endif

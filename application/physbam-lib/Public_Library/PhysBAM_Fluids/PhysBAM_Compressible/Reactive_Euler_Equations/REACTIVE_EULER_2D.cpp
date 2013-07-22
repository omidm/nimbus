//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Duc Nguyen.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Reactive_Euler_Equations/REACTIVE_EULER_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
// providing time is optional
template<class T> void REACTIVE_EULER_2D<T>::
Euler_Step(const T dt,const T time)
{  
    int m=grid.counts.x,n=grid.counts.y;
    int ghost_cells=3;
    
    ARRAY<TV_DIMENSION,VECTOR<int,2> > U_ghost(1-ghost_cells,m+ghost_cells,1-ghost_cells,n+ghost_cells);
    boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time,ghost_cells);

    T_FACE_ARRAYS_BOOL psi_N(grid.Get_MAC_Grid_At_Regular_Positions());
    T_FACE_ARRAYS_SCALAR face_velocities(grid.Get_MAC_Grid_At_Regular_Positions());
    VECTOR<EIGENSYSTEM<T,VECTOR<T,5> >*,2> eigensystem(&eigensystem_F,&eigensystem_G);
    if(cut_out_grid) 
        conservation->Update_Conservation_Law(grid,U,U_ghost,*psi_pointer,dt,eigensystem,eigensystem,psi_N,face_velocities);
    else{ // not a cut out grid
        ARRAY<bool,VECTOR<int,2> > psi(1,m,1,n);psi.Fill(1);
        conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystem,eigensystem,psi_N,face_velocities);}

    boundary->Apply_Boundary_Condition(grid,U,time+dt); 
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T> T REACTIVE_EULER_2D<T>::
CFL()
{
    int m=grid.counts.x,n=grid.counts.y;T dx=grid.dX.x,dy=grid.dX.y;
    
    ARRAY<T,VECTOR<int,2> > u_minus_c(1,m,1,n),u_plus_c(1,m,1,n),v_minus_c(1,m,1,n),v_plus_c(1,m,1,n);
    for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){
        if(!cut_out_grid || (cut_out_grid && (*psi_pointer)(i,j)==1)){
            T u=U(i,j)(2)/U(i,j)(1),v=U(i,j)(3)/U(i,j)(1);
            T Y=U(i,j)(5)/U(i,j)(1);
            T sound_speed=eos.c(U(i,j)(1),e(U(i,j)(1),U(i,j)(2),U(i,j)(3),U(i,j)(4)),Y);
            u_minus_c(i,j)=u-sound_speed;u_plus_c(i,j)=u+sound_speed;
            v_minus_c(i,j)=v-sound_speed;v_plus_c(i,j)=v+sound_speed;}}
    T dt_convect=max(u_minus_c.Maxabs(),u_plus_c.Maxabs())/dx+max(v_minus_c.Maxabs(),v_plus_c.Maxabs())/dy;
    return 1/dt_convect;
}
//#####################################################################
template class REACTIVE_EULER_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class REACTIVE_EULER_2D<double>;
#endif

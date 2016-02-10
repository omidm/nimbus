//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Math_Tools/maxabs.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_3D.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T> void EULER_3D<T>::
Euler_Step(const T dt,const T time)
{  
    int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z;
    int ghost_cells=3;
    
    ARRAY<TV_DIMENSION,VECTOR<int,3> > U_ghost(1-ghost_cells,m+ghost_cells,1-ghost_cells,n+ghost_cells,1-ghost_cells,mn+ghost_cells);
    boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time,ghost_cells);
    
    T_FACE_ARRAYS_BOOL psi_N(grid.Get_MAC_Grid_At_Regular_Positions());
    T_FACE_ARRAYS_SCALAR face_velocities(grid.Get_MAC_Grid_At_Regular_Positions());
    VECTOR<EIGENSYSTEM<T,VECTOR<T,5> >*,3> eigensystem(&eigensystem_F,&eigensystem_G,&eigensystem_H);
    if(psi_pointer) 
        conservation->Update_Conservation_Law(grid,U,U_ghost,*psi_pointer,dt,eigensystem,eigensystem,psi_N,face_velocities);
    else{ // not a cut out grid
        ARRAY<bool,VECTOR<int,3> > psi(1,m,1,n,1,mn);psi.Fill(1);
        conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystem,eigensystem,psi_N,face_velocities);}

    boundary->Apply_Boundary_Condition(grid,U,time+dt); 
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T> T EULER_3D<T>::
CFL()
{
    int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z;T dx=grid.dX.x,dy=grid.dX.y,dz=grid.dX.z;
    
    ARRAY<T,VECTOR<int,3> > u_minus_c(1,m,1,n,1,mn),u_plus_c(1,m,1,n,1,mn),v_minus_c(1,m,1,n,1,mn),
             v_plus_c(1,m,1,n,1,mn),w_minus_c(1,m,1,n,1,mn),w_plus_c(1,m,1,n,1,mn);
    for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) for(int ij=1;ij<=mn;ij++){
        if(!psi_pointer || (*psi_pointer)(i,j,ij)==1){
            T u=U(i,j,ij)(2)/U(i,j,ij)(1),v=U(i,j,ij)(3)/U(i,j,ij)(1),
                   w=U(i,j,ij)(4)/U(i,j,ij)(1);
            T sound_speed=eos->c(U(i,j,ij)(1),e(U(i,j,ij)(1),U(i,j,ij)(2),U(i,j,ij)(3),
                                     U(i,j,ij)(4),U(i,j,ij)(5)));
            u_minus_c(i,j,ij)=u-sound_speed;u_plus_c(i,j,ij)=u+sound_speed;
            v_minus_c(i,j,ij)=v-sound_speed;v_plus_c(i,j,ij)=v+sound_speed;
            w_minus_c(i,j,ij)=w-sound_speed;w_plus_c(i,j,ij)=w+sound_speed;}}
    T dt_convect=max(u_minus_c.Maxabs(),u_plus_c.Maxabs())/dx+max(v_minus_c.Maxabs(),v_plus_c.Maxabs())/dy+max(w_minus_c.Maxabs(),w_plus_c.Maxabs())/dz;
    return 1/dt_convect;
}              
//#####################################################################
template class EULER_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EULER_3D<double>;
#endif

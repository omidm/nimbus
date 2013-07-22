//#####################################################################
// Copyright 2002-2004, Eran Guendelman, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_2D  
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Shallow_Water_Equations/SHALLOW_WATER_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T> void SHALLOW_WATER_2D<T>::
Euler_Step(const T dt,const T time)
{  
    int m=grid.counts.x,n=grid.counts.y;
    int ghost_cells=3;

    // make sure things'll work in conservation law solver
    for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++) if (U(i,j)(1) < min_height){U(i,j)(1)=min_height;U(i,j)(2)=U(i,j)(3)=0;}

    ARRAY<TV_DIMENSION,VECTOR<int,2> > U_ghost(1-ghost_cells,m+ghost_cells,1-ghost_cells,n+ghost_cells);boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time,ghost_cells);

    ARRAY<bool,VECTOR<int,2> > psi(1,m,1,n);psi.Fill(true); // no cut out grids
    T_FACE_ARRAYS_BOOL psi_N(grid.Get_MAC_Grid_At_Regular_Positions());
    T_FACE_ARRAYS_SCALAR face_velocities(grid.Get_MAC_Grid_At_Regular_Positions());
    VECTOR<EIGENSYSTEM<T,VECTOR<T,3> >*,2> eigensystem(&eigensystem_F,&eigensystem_G);
    conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystem,eigensystem,psi_N,face_velocities);

    boundary->Apply_Boundary_Condition(grid,U,time+dt); 
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T> T SHALLOW_WATER_2D<T>::
CFL()
{
    T max_x_speed=0,max_y_speed=0;
    for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++){
        T one_over_h=1/U(i,j)(1);
        T u=U(i,j)(2)*one_over_h,v=U(i,j)(3)*one_over_h,celerity=sqrt(gravity*U(i,j)(1));
        max_x_speed=max(max_x_speed,abs(u)+celerity);
        max_y_speed=max(max_y_speed,abs(v)+celerity);}
    T dt_convect=max_x_speed*grid.one_over_dX.x+max_y_speed*grid.one_over_dX.y;
    return 1/dt_convect;
}
//#####################################################################
template class SHALLOW_WATER_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SHALLOW_WATER_2D<double>;
#endif

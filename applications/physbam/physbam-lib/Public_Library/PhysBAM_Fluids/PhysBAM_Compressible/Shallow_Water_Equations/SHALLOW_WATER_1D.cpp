//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Eran Guendelman, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_1D  
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Shallow_Water_Equations/SHALLOW_WATER_1D.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T> void SHALLOW_WATER_1D<T>::
Euler_Step(const T dt,const T time)
{  
    int m=grid.counts.x;
    int ghost_cells=3;
    ARRAY<TV_DIMENSION,VECTOR<int,1> > U_ghost(1-ghost_cells,m+ghost_cells);boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time,ghost_cells);

    // make sure things'll work in conservation law solver
    for(int i=1;i<=grid.counts.x;i++) if (U(i)(1) < min_height){U(i)(1)=min_height;U(i)(2)=0;}
    for(int i=1-ghost_cells;i<=m+ghost_cells;i++) if(U_ghost(i)(1) < min_height){U_ghost(i)(1)=min_height;U_ghost(i)(2)=0;}

    ARRAY<bool,VECTOR<int,1> > psi(1,m);psi.Fill(true); // no cut out grids
    T_FACE_ARRAYS_BOOL psi_N(grid.Get_MAC_Grid_At_Regular_Positions());
    T_FACE_ARRAYS_SCALAR face_velocities(grid.Get_MAC_Grid_At_Regular_Positions());
    VECTOR<EIGENSYSTEM<T,VECTOR<T,2> >*,1> eigensystem(&eigensystem_F);
    conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystem,eigensystem,psi_N,face_velocities);
    boundary->Apply_Boundary_Condition(grid,U,time+dt); 
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T> T SHALLOW_WATER_1D<T>::
CFL()
{
    T max_speed=0;
    for(int i=1;i<=grid.counts.x;i++){
        T u=(U(i)(1)!=0)?(U(i)(2)/U(i)(1)):0,celerity=sqrt(gravity*U(i)(1));
        max_speed=max(max_speed,abs(u)+celerity);}
    T dt_convect=max_speed*grid.one_over_dX.x;
    return 1/dt_convect;
}
//#####################################################################
template class SHALLOW_WATER_1D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SHALLOW_WATER_1D<double>;
#endif

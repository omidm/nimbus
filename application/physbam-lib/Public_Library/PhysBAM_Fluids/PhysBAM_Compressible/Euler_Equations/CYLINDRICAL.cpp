//#####################################################################
// Copyright 2002-2003, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/CYLINDRICAL.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
// providing time is optional
template<class T> void CYLINDRICAL<T>::
Euler_Step(const T dt,const T time)
{  
    int m=grid.counts.x,n=grid.counts.y;
    int i,j;
    int ghost_cells=3;

    ARRAY<TV_DIMENSION,VECTOR<int,2> > U_ghost(1-ghost_cells,m+ghost_cells,1-ghost_cells,n+ghost_cells);
    boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time,ghost_cells);
    
    // evaluate source terms
    ARRAY<TV_DIMENSION,VECTOR<int,2> > S(1,m,1,n);
    for(i=1;i<=m;i++) for(j=1;j<=n;j++){
        T rho=U(i,j)(1);
        T u=U(i,j)(2)/U(i,j)(1);
        T v=U(i,j)(3)/U(i,j)(1);
        T e=U(i,j)(4)/U(i,j)(1)-(sqr(u)+sqr(v))/2;
        T rho_u=U(i,j)(2);
        T coefficient=-1/grid.Axis_X(i,1);
        S(i,j)(1)=coefficient*rho_u;
        S(i,j)(2)=coefficient*rho_u*u;
        S(i,j)(3)=coefficient*rho_u*v;
        S(i,j)(4)=coefficient*(U(i,j)(4)+eos->p(rho,e))*u;}
    
    T_FACE_ARRAYS_BOOL psi_N(grid.Get_MAC_Grid_At_Regular_Positions());
    T_FACE_ARRAYS_SCALAR face_velocities(grid.Get_MAC_Grid_At_Regular_Positions());
    VECTOR<EIGENSYSTEM<T,VECTOR<T,4> >*,2> eigensystem(&eigensystem_F,&eigensystem_G);
    if(cut_out_grid) conservation->Update_Conservation_Law(grid,U,U_ghost,*psi_pointer,dt,eigensystem,eigensystem,psi_N,face_velocities);
    else{ // not a cut out grid
        ARRAY<bool,VECTOR<int,2> > psi(1,m,1,n);psi.Fill(1);
        conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystem,eigensystem,psi_N,face_velocities);}

    // add source terms
    for(i=1;i<=m;i++) for(j=1;j<=n;j++) U(i,j)+=dt*S(i,j);
    
    boundary->Apply_Boundary_Condition(grid,U,time+dt); 
}
//#####################################################################
template class CYLINDRICAL<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class CYLINDRICAL<double>;
#endif

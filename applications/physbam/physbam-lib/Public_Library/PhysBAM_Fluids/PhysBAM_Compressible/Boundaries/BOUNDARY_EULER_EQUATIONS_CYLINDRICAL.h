//#####################################################################
// Copyright 2002-2003, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_EULER_EQUATIONS_CYLINDRICAL
//#####################################################################
//
// Assumes a MAC grid straddling the centerline on the left, applying reflection.
// Extrapolation is applied to the right, top and bottom.
//
//#####################################################################
#ifndef __BOUNDARY_EULER_EQUATIONS_CYLINDRICAL__
#define __BOUNDARY_EULER_EQUATIONS_CYLINDRICAL__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
namespace PhysBAM{

template<class T_input>
class BOUNDARY_EULER_EQUATIONS_CYLINDRICAL:public BOUNDARY_UNIFORM<GRID<VECTOR<T_input,2> >,VECTOR<T_input,4> >
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<T,4> TV_DIMENSION;typedef VECTOR<int,2> TV_INT;typedef ARRAYS_ND_BASE<VECTOR<TV_DIMENSION,2> > T_ARRAYS_T2;
    enum {d=4};
public:
    BOUNDARY_EULER_EQUATIONS_CYLINDRICAL()
    {}

//#####################################################################
    void Fill_Ghost_Cells(const GRID<TV>& grid,const T_ARRAYS_T2& u,T_ARRAYS_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) PHYSBAM_OVERRIDE;
    void Apply_Boundary_Condition(const GRID<TV>& grid,T_ARRAYS_T2& u,const T time) PHYSBAM_OVERRIDE {} // do nothing
//#####################################################################
};
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class T> void BOUNDARY_EULER_EQUATIONS_CYLINDRICAL<T>::
Fill_Ghost_Cells(const GRID<TV>& grid,const T_ARRAYS_T2& u,T_ARRAYS_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells)
{
    int m=grid.counts.x,n=grid.counts.y;
    int i,j;

    T_ARRAYS_T2::Put(u,u_ghost); // interior

    // left
    for(i=-2;i<=0;i++) for(j=1;j<=n;j++){
        T rho=u_ghost(1-i,j)(1);
        T u_velocity=-u_ghost(1-i,j)(2)/u_ghost(1-i,j)(1);
        T v_velocity=u_ghost(1-i,j)(3)/u_ghost(1-i,j)(1);
        T e=u_ghost(1-i,j)(4)/u_ghost(1-i,j)(1)-(sqr(u_velocity)+sqr(v_velocity))/2;
        u_ghost(i,j)(1)=rho;
        u_ghost(i,j)(2)=rho*u_velocity;
        u_ghost(i,j)(3)=rho*v_velocity;
        u_ghost(i,j)(4)=rho*(e+(sqr(u_velocity)+sqr(v_velocity))/2);}

    // constant extrapolation
    for(j=1;j<=n;j++) u_ghost(m+3,j)=u_ghost(m+2,j)=u_ghost(m+1,j)=u_ghost(m,j); // right
    for(i=-2;i<=m+3;i++){
        u_ghost(i,-2)=u_ghost(i,-1)=u_ghost(i,0)=u_ghost(i,1);           // bottom
        u_ghost(i,n+3)=u_ghost(i,n+2)=u_ghost(i,n+1)=u_ghost(i,n);} // top
}
//#####################################################################
}
#endif

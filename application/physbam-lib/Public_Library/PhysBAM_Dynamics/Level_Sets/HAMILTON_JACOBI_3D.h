//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HAMILTON_JACOBI_3D  
//##################################################################### 
//
// Input a GRID_3D class.
// Input phi as (1,m) by (1,n) by (1,mn).
//
//#####################################################################
#ifndef __HAMILTON_JACOBI_3D__
#define __HAMILTON_JACOBI_3D__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Dynamics/Level_Sets/HAMILTON_JACOBI.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_3D.h>
namespace PhysBAM{

template<class T> class HAMILTONIAN_3D;

template<class T_input>
class HAMILTON_JACOBI_3D:public HAMILTON_JACOBI,public LEVELSET_3D<GRID<VECTOR<T_input,3> > >,LEVELSET_ADVECTION_3D<T_input>
{
    typedef T_input T;typedef VECTOR<T,3> TV;
public:
    typedef LEVELSET_3D<GRID<TV> > BASE;
    using BASE::grid;using BASE::phi;using BASE::boundary;using BASE::max_time_step;
    using BASE::curvature;using BASE::curvature_motion;using BASE::sigma;

    HAMILTONIAN_3D<T>& hamiltonian;

    HAMILTON_JACOBI_3D(HAMILTONIAN_3D<T>& hamiltonian_input,GRID<TV>& grid_input,ARRAY<T,VECTOR<int,3> >& phi_input) 
        :LEVELSET_3D<GRID<TV> >(grid_input,phi_input),LEVELSET_ADVECTION_3D<T>((LEVELSET_3D<GRID<TV> >*)this),hamiltonian(hamiltonian_input)
    {}
    
//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    void Calculate_Derivatives(ARRAY<T,VECTOR<int,3> >& phi_ghost,ARRAY<T,VECTOR<int,3> >& phix_minus,ARRAY<T,VECTOR<int,3> >& phix_plus,ARRAY<T,VECTOR<int,3> >& phiy_minus,ARRAY<T,VECTOR<int,3> >& phiy_plus,ARRAY<T,VECTOR<int,3> >& phiz_minus,ARRAY<T,VECTOR<int,3> >& phiz_plus);
    T CFL(const T time=0);
//#####################################################################
};   
}
#endif

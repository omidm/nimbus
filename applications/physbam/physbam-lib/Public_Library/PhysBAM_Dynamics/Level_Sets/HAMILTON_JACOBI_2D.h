//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HAMILTON_JACOBI_2D  
//##################################################################### 
//
// Input a GRID_2D class.
// Input phi as (1,m) by (1,n).
//
//#####################################################################
#ifndef __HAMILTON_JACOBI_2D__
#define __HAMILTON_JACOBI_2D__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Dynamics/Level_Sets/HAMILTON_JACOBI.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_2D.h>
namespace PhysBAM{

template<class T> class HAMILTONIAN_2D;

template<class T_input>
class HAMILTON_JACOBI_2D:public HAMILTON_JACOBI,public LEVELSET_2D<GRID<VECTOR<T_input,2> > >,public LEVELSET_ADVECTION_2D<T_input>
{
    typedef T_input T;typedef VECTOR<T,2> TV;
public:
    typedef LEVELSET_2D<GRID<TV> > BASE;
    using BASE::grid;using BASE::phi;using BASE::boundary;using BASE::max_time_step;
    using BASE::curvature;using BASE::curvature_motion;using BASE::sigma;

    HAMILTONIAN_2D<T>& hamiltonian;

    HAMILTON_JACOBI_2D(HAMILTONIAN_2D<T>& hamiltonian_input,GRID<VECTOR<T,2> >& grid_input,ARRAY<T,VECTOR<int,2> >& phi_input) 
        :LEVELSET_2D<GRID<TV> >(grid_input,phi_input),LEVELSET_ADVECTION_2D<T>((LEVELSET_2D<GRID<TV> >*)this),hamiltonian(hamiltonian_input)
    {}
    
//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    void Calculate_Derivatives(ARRAY<T,VECTOR<int,2> >& phi_ghost,ARRAY<T,VECTOR<int,2> >& phix_minus,ARRAY<T,VECTOR<int,2> >& phix_plus,ARRAY<T,VECTOR<int,2> >& phiy_minus,ARRAY<T,VECTOR<int,2> >& phiy_plus);
    T CFL(const T time=0);
//#####################################################################
};   
}
#endif

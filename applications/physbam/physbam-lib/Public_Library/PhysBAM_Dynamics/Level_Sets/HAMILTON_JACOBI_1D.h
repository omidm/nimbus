//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HAMILTON_JACOBI_1D  
//##################################################################### 
//
// Input a GRID_1D class.
// Input phi as (1,m).
//
//#####################################################################
#ifndef __HAMILTON_JACOBI_1D__
#define __HAMILTON_JACOBI_1D__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Dynamics/Level_Sets/HAMILTON_JACOBI.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_1D.h>
namespace PhysBAM{

template<class T> class HAMILTONIAN_1D;

template<class T_input>
class HAMILTON_JACOBI_1D:public HAMILTON_JACOBI,public LEVELSET_1D<T_input>,public LEVELSET_ADVECTION_1D<T_input>
{
    typedef T_input T;
public:
    typedef LEVELSET_1D<T> BASE;
    using BASE::grid;using BASE::phi;using BASE::boundary;using BASE::max_time_step;

    HAMILTONIAN_1D<T>& hamiltonian;
    
    HAMILTON_JACOBI_1D(HAMILTONIAN_1D<T>& hamiltonian_input,GRID<VECTOR<T,1> >& grid_input,ARRAY<T,VECTOR<int,1> >& phi_input) 
        :LEVELSET_1D<T>(grid_input,phi_input),LEVELSET_ADVECTION_1D<T>((LEVELSET_1D<T>*)this),hamiltonian(hamiltonian_input)
    {}
    
//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    void Calculate_Derivatives(ARRAY<T,VECTOR<int,1> >& phi_ghost,ARRAY<T,VECTOR<int,1> >& phix_minus,ARRAY<T,VECTOR<int,1> >& phix_plus);
    T CFL(const T time=0);
//#####################################################################
};
}    
#endif

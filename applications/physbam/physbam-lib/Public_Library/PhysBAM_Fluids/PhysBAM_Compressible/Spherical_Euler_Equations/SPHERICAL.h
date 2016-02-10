//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPHERICAL  
//##################################################################### 
//
// Input an EOS class, that is a class that inherits the virtual base class EOS.
// Input a GRID_1D class which MUST BE a MAC grid to avoid division by zero at the origin!
// Input U as 3 by (1,m) for mass, momentum, and energy.
//
// Inherits EULER_1D.
//
// Use Set_Up_Cut_Out_Grid(psi) to define a cut out grid with psi as (1,m).
// When psi=1, solve the equaitions. 
// When psi=0, do NOT solve the equations.
//
//#####################################################################
#ifndef __SPHERICAL__
#define __SPHERICAL__

#include <PhysBAM_Fluids/PhysBAM_Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_SPHERICAL.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_1D.h>
namespace PhysBAM{

template<class T_input>
class SPHERICAL:public EULER_1D<T_input>
{
    typedef T_input T;typedef VECTOR<T,1> TV;typedef VECTOR<T,3> TV_DIMENSION;
    typedef typename  GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
public:
    using EULER_1D<T>::boundary;using EULER_1D<T>::conservation;using EULER_1D<T>::eos;using EULER_1D<T>::grid;using EULER_1D<T>::U;
    using EULER_1D<T>::cut_out_grid;using EULER_1D<T>::psi_pointer;using EULER_1D<T>::eigensystem_F;

private:
    BOUNDARY_EULER_EQUATIONS_SPHERICAL<T> boundary_default;
public:
    SPHERICAL(EOS<T>& eos_input,GRID<VECTOR<T,1> >& grid_input,ARRAY<TV_DIMENSION,VECTOR<int,1> >& U_input) 
        :EULER_1D<T>(U_input)
    {
        Set_Custom_Equation_Of_State(eos_input);
        boundary=&boundary_default;
    }
    
//#####################################################################
    void Euler_Step(const T dt,const T time=0);
//#####################################################################
};
}    
#endif

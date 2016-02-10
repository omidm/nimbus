//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CYLINDRICAL  
//##################################################################### 
//
// Input an EOS class, that is a class that inherits the virtual base class EOS.
// Input a GRID_2D class which MUST BE a MAC grid to avoid division by zero on the left wall!
// Input U as 4 by (1,m) by (1,n) for mass, momentum, and energy.
//
// Inherits EULER_2D.
//
// Use Set_Up_Cut_Out_Grid(psi) to define a cut out grid with psi as (1,m) by (1,n).
// When psi=1, solve the equaitions. 
// When psi=0, do NOT solve the equations.
//
//#####################################################################
#ifndef __CYLINDRICAL__
#define __CYLINDRICAL__

#include <PhysBAM_Fluids/PhysBAM_Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_CYLINDRICAL.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_2D.h>
namespace PhysBAM{

template<class T_input>
class CYLINDRICAL:public EULER_2D<T_input>
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<T,4> TV_DIMENSION;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
public:
    using EULER_2D<T>::boundary;using EULER_2D<T>::conservation;using EULER_2D<T>::eos;using EULER_2D<T>::grid;using EULER_2D<T>::U;
    using EULER_2D<T>::cut_out_grid;using EULER_2D<T>::psi_pointer;using EULER_2D<T>::eigensystem_F;using EULER_2D<T>::eigensystem_G;

private:
    BOUNDARY_EULER_EQUATIONS_CYLINDRICAL<T> boundary_default;
public:
    CYLINDRICAL(EOS<T>& eos_input,GRID<VECTOR<T,2> >& grid_input,ARRAY<TV_DIMENSION,VECTOR<int,2> >& U_input) 
        :EULER_2D<T>(U_input)
    {
        boundary=&boundary_default;
    }
    
//#####################################################################
    void Euler_Step(const T dt,const T time=0);
//#####################################################################
};   
}
#endif

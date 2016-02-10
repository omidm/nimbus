//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REACTIVE_EULER_1D  
//##################################################################### 
//
// Input an EOS_REACTIVE class, that is a class that inherits the virtual base class EOS_REACTIVE.
// Input a GRID_1D class.
// Input U as 4 by (1,m) for mass, momentum, energy, and mass fraction density
//
// Use Set_Up_Cut_Out_Grid(psi) to define a cut out grid with psi as (1,m).
// When psi=1, solve the equaitions. 
// When psi=0, do NOT solve the equations.
//
//#####################################################################
#ifndef __REACTIVE_EULER_1D__
#define __REACTIVE_EULER_1D__

#include <PhysBAM_Fluids/PhysBAM_Compressible/Reactive_Euler_Equations/REACTIVE_EULER.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Reactive_Euler_Equations/REACTIVE_EULER_1D_EIGENSYSTEM_F.h>
namespace PhysBAM{

template<class T_input>
class REACTIVE_EULER_1D:public REACTIVE_EULER<GRID<VECTOR<T_input,1> > >
{   
    typedef T_input T;typedef VECTOR<T,1> TV;typedef VECTOR<T,4> TV_DIMENSION;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
protected:
    typedef REACTIVE_EULER<GRID<TV> > BASE;
    using BASE::cut_out_grid;using BASE::boundary;using BASE::eos;using BASE::conservation;

    GRID<TV>& grid;
    ARRAY<TV_DIMENSION,VECTOR<int,1> >& U;         // mass, momentum, and energy
    ARRAY<bool,VECTOR<int,1> >* psi_pointer; // defines cut out grid
    REACTIVE_EULER_1D_EIGENSYSTEM_F<T> eigensystem_F;

public:
    REACTIVE_EULER_1D(REACTIVE_EOS<T>& eos_input,GRID<VECTOR<T,1> >& grid_input,ARRAY<TV_DIMENSION,VECTOR<int,1> >& U_input)  
        :REACTIVE_EULER<GRID<TV> >(eos_input),grid(grid_input),U(U_input),eigensystem_F(eos_input)
    {}
    
    void Set_Up_Cut_Out_Grid(ARRAY<bool,VECTOR<int,1> >& psi_input)
    {psi_pointer=&psi_input;cut_out_grid=1;}
    
//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    T CFL();
//#####################################################################
};   
}
#endif

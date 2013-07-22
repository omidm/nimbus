//#####################################################################
// Copyright 2002-2003, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REACTIVE_EULER_2D  
//##################################################################### 
//
// Input an EOS_REACTIVE class, that is a class that inherits the virtual base class EOS_REACTIVE.
// Input a GRID_1D class.
// Input U as 4 by (1,m) by (1,n) for mass, momentum, energy, and mass fraction density
//
// Use Set_Up_Cut_Out_Grid(psi) to define a cut out grid with psi as (1,m) by (1,n).
// When psi=1, solve the equaitions. 
// When psi=0, do NOT solve the equations.
//
//#####################################################################
#ifndef __REACTIVE_EULER_2D__
#define __REACTIVE_EULER_2D__

#include <PhysBAM_Fluids/PhysBAM_Compressible/Reactive_Euler_Equations/REACTIVE_EULER.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Reactive_Euler_Equations/REACTIVE_EULER_2D_EIGENSYSTEM_F.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Reactive_Euler_Equations/REACTIVE_EULER_2D_EIGENSYSTEM_G.h>
namespace PhysBAM{

template<class T_input>
class REACTIVE_EULER_2D:public REACTIVE_EULER<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<T,5> TV_DIMENSION;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename REBIND<T_FACE_ARRAYS_SCALAR,bool>::TYPE T_FACE_ARRAYS_BOOL;
protected:
    typedef REACTIVE_EULER<GRID<TV> > BASE;
    using BASE::cut_out_grid;using BASE::boundary;using BASE::eos;using BASE::conservation;

    GRID<TV>& grid;
    ARRAY<TV_DIMENSION,VECTOR<int,2> >& U;         // mass, momentum, and energy
    ARRAY<bool,VECTOR<int,2> >* psi_pointer; // defines cut out grid
    REACTIVE_EULER_2D_EIGENSYSTEM_F<T> eigensystem_F;
    REACTIVE_EULER_2D_EIGENSYSTEM_G<T> eigensystem_G;

public:
    REACTIVE_EULER_2D(REACTIVE_EOS<T>& eos_input,GRID<VECTOR<T,2> >& grid_input,ARRAY<TV_DIMENSION,VECTOR<int,2> >& U_input)  
        :REACTIVE_EULER<GRID<TV> >(eos_input),grid(grid_input),U(U_input),eigensystem_F(eos_input),eigensystem_G(eos_input)
    {}
    
    void Set_Up_Cut_Out_Grid(ARRAY<bool,VECTOR<int,2> >& psi_input)
    {psi_pointer=&psi_input;cut_out_grid=1;}
    
//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    T CFL();
//#####################################################################
};   
}
#endif

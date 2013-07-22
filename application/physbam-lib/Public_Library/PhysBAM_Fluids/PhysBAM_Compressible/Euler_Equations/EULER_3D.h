//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_3D  
//##################################################################### 
//
// Input an EOS class, that is a class that inherits the virtual base class EOS.
// Input a GRID_3D class.
// Input U as 5 by (1,m) by (1,n) by (1,mn) for mass, momentum, and energy.
//
// Use Set_Up_Cut_Out_Grid(psi) to define a cut out grid with psi as (1,m) by (1,n) by (1,mn).
// When psi=1, solve the equaitions. 
// When psi=0, do NOT solve the equations.
//
//#####################################################################
#ifndef __EULER_3D__
#define __EULER_3D__

#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_3D_EIGENSYSTEM_F.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_3D_EIGENSYSTEM_G.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_3D_EIGENSYSTEM_H.h>
namespace PhysBAM{

template<class T_input>
class EULER_3D:public EULER<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;typedef VECTOR<T,3> TV;typedef VECTOR<T,5> TV_DIMENSION;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename REBIND<T_FACE_ARRAYS_SCALAR,bool>::TYPE T_FACE_ARRAYS_BOOL;
protected:
    using EULER<GRID<TV> >::boundary;using EULER<GRID<TV> >::conservation;using EULER<GRID<TV> >::eos;

    GRID<TV>& grid;
    ARRAY<TV_DIMENSION,VECTOR<int,3> >& U;             // mass, momentum, and energy
    ARRAY<bool,VECTOR<int,3> >* psi_pointer; // defines cut out grid
    EULER_3D_EIGENSYSTEM_F<T> eigensystem_F;
    EULER_3D_EIGENSYSTEM_G<T> eigensystem_G;
    EULER_3D_EIGENSYSTEM_H<T> eigensystem_H;
    
public:
    EULER_3D(EOS<T>& eos_input,GRID<TV>& grid_input,ARRAY<TV_DIMENSION,VECTOR<int,3> >& U_input)  
        :grid(grid_input),U(U_input),psi_pointer(0),eigensystem_F(),eigensystem_G(),eigensystem_H()
    {
        Set_Custom_Equation_Of_State(eos_input);
    }
    
    void Set_Up_Cut_Out_Grid(ARRAY<bool,VECTOR<int,3> >& psi_input)
    {psi_pointer=&psi_input;}
    
//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    T CFL();
//#####################################################################
};
}    
#endif

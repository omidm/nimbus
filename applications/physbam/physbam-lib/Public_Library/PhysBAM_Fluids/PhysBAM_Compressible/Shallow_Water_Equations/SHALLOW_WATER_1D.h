//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Eran Guendelman, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_1D  
//##################################################################### 
//
// Input U as 2 by (1,m) for h and h*u
//
//#####################################################################
#ifndef __SHALLOW_WATER_1D__
#define __SHALLOW_WATER_1D__

#include <PhysBAM_Fluids/PhysBAM_Compressible/Shallow_Water_Equations/SHALLOW_WATER.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Shallow_Water_Equations/SHALLOW_WATER_1D_EIGENSYSTEM_F.h>
namespace PhysBAM{

template<class T_input>
class SHALLOW_WATER_1D:public SHALLOW_WATER<GRID<VECTOR<T_input,1> > >
{
    typedef T_input T;typedef VECTOR<T,1> TV;typedef VECTOR<T,2> TV_DIMENSION;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
public:
    using SHALLOW_WATER<GRID<TV> >::boundary;using SHALLOW_WATER<GRID<TV> >::conservation;

    T gravity;
    T min_height;
protected:
    GRID<TV>& grid;
    ARRAY<TV_DIMENSION,VECTOR<int,1> >& U; // h and h*u
    SHALLOW_WATER_1D_EIGENSYSTEM_F<T> eigensystem_F;

public:
    SHALLOW_WATER_1D(GRID<TV>& grid_input,ARRAY<TV_DIMENSION,VECTOR<int,1> >& U_input,T gravity_input=9.8,T min_height_input=1e-3)
        :gravity(gravity_input),min_height(min_height_input),grid(grid_input),U(U_input),eigensystem_F(gravity_input)
    {}

//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    T CFL();
//#####################################################################
};   
}
#endif

//#####################################################################
// Copyright 2002-2004, Eran Guendelman, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_2D  
//##################################################################### 
//
// Input U as 3 by (1,m) by (1,n) for h, h*u, and h*v.
//
//#####################################################################
#ifndef __SHALLOW_WATER_2D__
#define __SHALLOW_WATER_2D__

#include <PhysBAM_Fluids/PhysBAM_Compressible/Shallow_Water_Equations/SHALLOW_WATER.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Shallow_Water_Equations/SHALLOW_WATER_2D_EIGENSYSTEM_F.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Shallow_Water_Equations/SHALLOW_WATER_2D_EIGENSYSTEM_G.h>
namespace PhysBAM{

template<class T_input>
class SHALLOW_WATER_2D:public SHALLOW_WATER<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<T,3> TV_DIMENSION;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename REBIND<T_FACE_ARRAYS_SCALAR,bool>::TYPE T_FACE_ARRAYS_BOOL;
    T_FACE_ARRAYS_SCALAR face_velocities;
public:
    using SHALLOW_WATER<GRID<TV> >::boundary;using SHALLOW_WATER<GRID<TV> >::conservation;

    T gravity;
    T min_height;
protected:
    GRID<TV>& grid;
    ARRAY<TV_DIMENSION,VECTOR<int,2> >& U; // h, h*u, and h*v
    SHALLOW_WATER_2D_EIGENSYSTEM_F<T> eigensystem_F;
    SHALLOW_WATER_2D_EIGENSYSTEM_G<T> eigensystem_G;

public:
    SHALLOW_WATER_2D(GRID<TV>& grid_input,ARRAY<TV_DIMENSION,VECTOR<int,2> >& U_input,T gravity_input=9.8,T min_height_input=1e-3)
        :gravity(gravity_input),min_height(min_height_input),grid(grid_input),U(U_input),eigensystem_F(gravity_input),eigensystem_G(gravity_input)
    {}
    
//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    T CFL();
//#####################################################################
};   
}
#endif

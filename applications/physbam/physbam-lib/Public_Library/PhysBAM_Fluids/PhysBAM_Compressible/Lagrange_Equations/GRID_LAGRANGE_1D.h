//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_LAGRANGE_1D  
//#####################################################################
//
// One dimensional Lagrangian grid. 
//
//#####################################################################
#ifndef __GRID_LAGRANGE_1D__
#define __GRID_LAGRANGE_1D__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/GRID_LAGRANGE.h>
namespace PhysBAM{

template<class T>
class GRID_LAGRANGE_1D:public GRID_LAGRANGE<T>
{
public:
    int m; // # of points: dimension 1
    ARRAY<T,VECTOR<int,1> >& x; // x spatial location

    GRID_LAGRANGE_1D(const int m_input,ARRAY<T,VECTOR<int,1> >& x_input)
        :m(m_input),x(x_input)
    {}

//#####################################################################
    void Euler_Step(const ARRAY<T,VECTOR<int,1> >& u,const T dt);
    void Get_Lengths(ARRAY<T,VECTOR<int,1> >& L);
    void Get_Midpoints(ARRAY<T,VECTOR<int,1> >& M);
//#####################################################################
};
}
#endif

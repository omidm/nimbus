//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_LAGRANGE_2D
//#####################################################################
#ifndef __GRID_LAGRANGE_2D__
#define __GRID_LAGRANGE_2D__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/GRID_LAGRANGE.h>
namespace PhysBAM{

template<class T>
class GRID_LAGRANGE_2D:public GRID_LAGRANGE<T>
{
public:
    int m;       // # of points: dimension 1
    int n;       // # of points: dimension 2
    ARRAY<T,VECTOR<int,2> >& x; // x spatial location
    ARRAY<T,VECTOR<int,2> >& y; // y spatial location

    GRID_LAGRANGE_2D(const int m_input,const int n_input,ARRAY<T,VECTOR<int,2> >& x_input,ARRAY<T,VECTOR<int,2> >& y_input)
        :m(m_input),n(n_input),x(x_input),y(y_input)
    {}

//#####################################################################
    void Euler_Step(const ARRAY<T,VECTOR<int,2> >& u,const ARRAY<T,VECTOR<int,2> >& v,const T dt);
    void Get_Lengths(ARRAY<T,VECTOR<int,2> >& L1,ARRAY<T,VECTOR<int,2> >& L2);
    void Get_Areas(ARRAY<T,VECTOR<int,2> >& A);
    void Get_Normals(ARRAY<T,VECTOR<int,2> >& N1_x,ARRAY<T,VECTOR<int,2> >& N1_y,ARRAY<T,VECTOR<int,2> >& N2_x,ARRAY<T,VECTOR<int,2> >& N2_y);
    void Get_Centers(ARRAY<T,VECTOR<int,2> >& C_x,ARRAY<T,VECTOR<int,2> >& C_y);
    void Get_Midpoints(ARRAY<T,VECTOR<int,2> >& M1_x,ARRAY<T,VECTOR<int,2> >& M1_y,ARRAY<T,VECTOR<int,2> >& M2_x,ARRAY<T,VECTOR<int,2> >& M2_y);
    void Get_Sub_Zone_Lengths(ARRAY<T,VECTOR<int,2> >& LL1,ARRAY<T,VECTOR<int,2> >& LL2,ARRAY<T,VECTOR<int,2> >& LL3,ARRAY<T,VECTOR<int,2> >& LL4);
    void Get_Sub_Zone_Areas(ARRAY<T,VECTOR<int,2> >& AA1,ARRAY<T,VECTOR<int,2> >& AA2,ARRAY<T,VECTOR<int,2> >& AA3,ARRAY<T,VECTOR<int,2> >& AA4);
    void Get_Sub_Zone_Normals(ARRAY<T,VECTOR<int,2> >& NN1_x,ARRAY<T,VECTOR<int,2> >& NN1_y,ARRAY<T,VECTOR<int,2> >& NN2_x,ARRAY<T,VECTOR<int,2> >& NN2_y,ARRAY<T,VECTOR<int,2> >& NN3_x,ARRAY<T,VECTOR<int,2> >& NN3_y,ARRAY<T,VECTOR<int,2> >& NN4_x,ARRAY<T,VECTOR<int,2> >& NN4_y);
//#####################################################################
};
}
#endif

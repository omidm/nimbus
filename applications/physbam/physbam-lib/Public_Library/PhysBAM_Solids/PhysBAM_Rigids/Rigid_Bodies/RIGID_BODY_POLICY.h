//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_POLICY 
//#####################################################################
#ifndef __RIGID_BODY_POLICY__
#define __RIGID_BODY_POLICY__

namespace PhysBAM{

template<class TV> struct RIGID_BODY_POLICY;
template<class T,int m_input,int n_input> class MATRIX;
template<class T,int d> class DIAGONAL_MATRIX;
template<class T,int d> class SYMMETRIC_MATRIX;
template<class T,int d> class VECTOR;

//#####################################################################
// 1D
//#####################################################################
template<class T>
struct RIGID_BODY_POLICY<VECTOR<T,1> >
{
    typedef MATRIX<T,0,0> INERTIA_TENSOR;
    typedef MATRIX<T,0,0> WORLD_SPACE_INERTIA_TENSOR;
};
//#####################################################################
// 2D
//#####################################################################
template<class T>
struct RIGID_BODY_POLICY<VECTOR<T,2> >
{
    typedef MATRIX<T,1,1> INERTIA_TENSOR;
    typedef MATRIX<T,1,1> WORLD_SPACE_INERTIA_TENSOR;
};
//#####################################################################
// 3D
//#####################################################################
template<class T>
struct RIGID_BODY_POLICY<VECTOR<T,3> >
{
    typedef DIAGONAL_MATRIX<T,3> INERTIA_TENSOR;
    typedef SYMMETRIC_MATRIX<T,3> WORLD_SPACE_INERTIA_TENSOR;
};
//#####################################################################
}
#endif

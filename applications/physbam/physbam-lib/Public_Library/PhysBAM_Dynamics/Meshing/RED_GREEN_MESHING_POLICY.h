//#####################################################################
// Copyright 2006, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RED_GREEN_MESHING_POLICY 
//#####################################################################
#ifndef __RED_GREEN_MESHING_POLICY__
#define __RED_GREEN_MESHING_POLICY__

#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
namespace PhysBAM{

template<class TV> class RED_GREEN_TRIANGLES;
template<class T> class RED_GREEN_TETRAHEDRA;

template<class TV> struct RED_GREEN_MESHING_POLICY;

//#####################################################################
// 2D
//#####################################################################
template<class T>
struct RED_GREEN_MESHING_POLICY<VECTOR<T,2> >
{
    typedef RED_GREEN_TRIANGLES<VECTOR<T,2> > RED_GREEN_SIMPLICES;
};
//#####################################################################
// 3D
//#####################################################################
template<class T>
struct RED_GREEN_MESHING_POLICY<VECTOR<T,3> >
{
    typedef RED_GREEN_TETRAHEDRA<T> RED_GREEN_SIMPLICES;
};
//#####################################################################
}
#endif

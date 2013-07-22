//#####################################################################
// Copyright 2006-2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FORCES_POLICY 
//#####################################################################
#ifndef __SOLIDS_FORCES_POLICY__
#define __SOLIDS_FORCES_POLICY__

#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class TV> class SOLIDS_FORCES;
template<class T> class LINEAR_ALTITUDE_SPRINGS_2D;
template<class T> class LINEAR_ALTITUDE_SPRINGS_S3D;
template<class T> class LINEAR_ALTITUDE_SPRINGS_3D;

template<class TV,int d> struct SOLIDS_FORCES_POLICY;

//#####################################################################
// 2D
//#####################################################################
template<class T>
struct SOLIDS_FORCES_POLICY<VECTOR<T,2>,2>
{
    typedef LINEAR_ALTITUDE_SPRINGS_2D<T> LINEAR_ALTITUDE_SPRINGS;
};
//#####################################################################
// S3D
//#####################################################################
template<class T>
struct SOLIDS_FORCES_POLICY<VECTOR<T,3>,2>
{
    typedef LINEAR_ALTITUDE_SPRINGS_S3D<T> LINEAR_ALTITUDE_SPRINGS;
};
//#####################################################################
// 3D
//#####################################################################
template<class T>
struct SOLIDS_FORCES_POLICY<VECTOR<T,3>,3>
{
    typedef LINEAR_ALTITUDE_SPRINGS_3D<T> LINEAR_ALTITUDE_SPRINGS;
};
//#####################################################################
}
#endif

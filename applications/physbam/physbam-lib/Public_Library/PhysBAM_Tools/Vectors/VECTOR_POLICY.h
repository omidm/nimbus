//#####################################################################
// Copyright 2011, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VECTOR_POLICY
//#####################################################################
#ifndef __VECTOR_POLICY__
#define __VECTOR_POLICY__

#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template <class TV> struct VECTOR_POLICY;

template <class T>
struct VECTOR_POLICY<VECTOR<T,1> >{
    static const int DIMMINUSONE = 0;
};

template <class T>
struct VECTOR_POLICY<VECTOR<T,2> >{
    static const int DIMMINUSONE = 1;
};

template <class T>
struct VECTOR_POLICY<VECTOR<T,3> >{
    static const int DIMMINUSONE = 2;
};

}
#endif

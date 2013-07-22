//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIAGONALIZED_SEMI_IMPLICIT_ELEMENT
//##################################################################### 
#ifndef __DIAGONALIZED_SEMI_IMPLICIT_ELEMENT__
#define __DIAGONALIZED_SEMI_IMPLICIT_ELEMENT__

#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class TV,int d>
class DIAGONALIZED_SEMI_IMPLICIT_ELEMENT
{
    typedef typename TV::SCALAR T;
public:
    VECTOR<int,d+1> nodes;
    T Bm_scale;
    DIAGONAL_MATRIX<T,d> W;
    MATRIX<T,d> Dm_inverse;
    T dt_cfl,time;

//#####################################################################
};
}
#endif

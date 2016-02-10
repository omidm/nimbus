//#####################################################################
// Copyright 2011, Valeria Nikolaenko
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SCALAR_COMPONENT
//#####################################################################
#ifndef __SCALAR_COMPONENT__
#define __SCALAR_COMPONENT__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <climits>

namespace PhysBAM {

template<class TV>
class SCALAR_COMPONENT {
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
public:
    virtual ARRAY<T,TV_INT>* Get_Scalar_Values_Simulated() = 0;
    virtual GRID<TV>* Get_Grid() = 0;
    virtual int Get_Upsample_Scale() = 0;
};

}

#endif



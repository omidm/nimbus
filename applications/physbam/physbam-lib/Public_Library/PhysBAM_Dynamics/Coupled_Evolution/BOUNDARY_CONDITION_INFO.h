//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_CONDITION_INFO
//#####################################################################
#ifndef __BOUNDARY_CONDITION_INFO__
#define __BOUNDARY_CONDITION_INFO__
#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{
template<class TV>
struct BOUNDARY_CONDITION_INFO
{
    typedef VECTOR<int,TV::m> TV_INT;typedef typename TV::SCALAR T;
    int type;
    T theta,value;
    FACE_INDEX<TV::m> f;
    int in_side;
};
}
#endif

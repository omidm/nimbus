#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_WRAPPER_FIRE_RLE
//#####################################################################
#ifndef __ADVECTION_WRAPPER_FIRE_RLE__
#define __ADVECTION_WRAPPER_FIRE_RLE__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_NESTED_ADVECTION>
class ADVECTION_WRAPPER_FIRE_RLE:public ADVECTION<T_GRID,T2>
{
    typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::PROJECTION T_PROJECTION;typedef typename T_GRID::LEVELSET_MULTIPLE T_LEVELSET_MULTIPLE;
public:
    T_NESTED_ADVECTION& nested_advection;
    const T_PROJECTION& projection;
    const ARRAY<T>& phi_n_plus_one;

    ADVECTION_WRAPPER_FIRE_RLE(T_NESTED_ADVECTION& nested_advection_input,const T_PROJECTION& projection_input,const ARRAY<T>& phi_n_plus_one_input)
        :nested_advection(nested_advection_input),projection(projection_input),phi_n_plus_one(phi_n_plus_one_input)
    {
        PHYSBAM_NOT_IMPLEMENTED();
    }

    ADVECTION_WRAPPER_FIRE_RLE(T_NESTED_ADVECTION& nested_advection_input,const T_PROJECTION& projection_input,const T_LEVELSET_MULTIPLE& levelset_multiple_input)
        :nested_advection(nested_advection_input),projection(projection_input)
    {
        PHYSBAM_NOT_IMPLEMENTED();
    }

//#####################################################################
};
}
#endif
#endif 

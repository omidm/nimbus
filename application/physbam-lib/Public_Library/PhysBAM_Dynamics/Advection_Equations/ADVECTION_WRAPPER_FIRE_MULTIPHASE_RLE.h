#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_WRAPPER_FIRE_MULTIPHASE_RLE
//#####################################################################
#ifndef __ADVECTION_WRAPPER_FIRE_MULTIPHASE_RLE__
#define __ADVECTION_WRAPPER_FIRE_MULTIPHASE_RLE__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/LEVELSET_POLICY_RLE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_POLICY.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_NESTED_ADVECTION>
class ADVECTION_WRAPPER_FIRE_MULTIPHASE_RLE:public ADVECTION<T_GRID,T2>
{
    typedef typename INCOMPRESSIBLE_POLICY<T_GRID>::PROJECTION T_PROJECTION;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET_MULTIPLE T_LEVELSET_MULTIPLE;
public:

    ADVECTION_WRAPPER_FIRE_MULTIPHASE_RLE(T_NESTED_ADVECTION& nested_advection_input,const T_PROJECTION& projection_input,const T_LEVELSET_MULTIPLE& levelset_multiple_input)
    {
        PHYSBAM_NOT_IMPLEMENTED();
    }

//#####################################################################
};
}
#endif
#endif

//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace SIGNED_DISTANCE
//##################################################################### 
#ifndef __BOX_SIGNED_DISTANCE__
#define __BOX_SIGNED_DISTANCE__
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>

namespace PhysBAM{
namespace SIGNED_DISTANCE{
template<class TV,class T_GRID,class T_ARRAY> void Calculate(BOX<TV>& box,const T_GRID& grid,T_ARRAY& array);
};
};
#endif

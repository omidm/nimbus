//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header FIRE_INTERPOLATION_FORWARD
//#####################################################################
#ifndef __FIRE_INTERPOLATION_FORWARD__
#define __FIRE_INTERPOLATION_FORWARD__

#include <PhysBAM_Geometry/Grids_Dyadic_Interpolation_Collidable/INTERPOLATION_COLLIDABLE_DYADIC_FORWARD.h>
#include <PhysBAM_Geometry/Grids_RLE_Interpolation_Collidable/INTERPOLATION_COLLIDABLE_RLE_FORWARD.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/INTERPOLATION_COLLIDABLE_UNIFORM_FORWARD.h>
namespace PhysBAM{

template<class T_GRID> class FACE_LOOKUP_FIRE_DYADIC;
template<class T_GRID> class FACE_LOOKUP_FIRE_MULTIPHASE_DYADIC;
template<class T_GRID> class FACE_LOOKUP_FIRE_MULTIPHASE_RLE;
template<class T_GRID> class FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM;
template<class T_GRID> class FACE_LOOKUP_FIRE_RLE;

}
#endif

//#####################################################################
// Copyright 2006, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header INCOMPRESSIBLE_FORWARD
//#####################################################################
#ifndef __INCOMPRESSIBLE_FORWARD__
#define __INCOMPRESSIBLE_FORWARD__

#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD_FORWARD.h>
namespace PhysBAM{

template<class T_GRID> class INCOMPRESSIBLE_DYADIC;
template<class T> class INCOMPRESSIBLE_QUADTREE;
template<class T> class INCOMPRESSIBLE_OCTREE;
template<class T_GRID> class INCOMPRESSIBLE_UNIFORM;
template<class T_GRID> class INCOMPRESSIBLE_MULTIPHASE_UNIFORM;
template<class T_GRID> class INCOMPRESSIBLE_RLE;
template<class T_GRID> class SPH_EVOLUTION_UNIFORM;
template<class T_GRID> class SPH_CALLBACKS;
template<class T_GRID> class PROJECTION_DYNAMICS_UNIFORM;
template<class T_GRID> class PROJECTION_DYADIC;
template<class T_GRID> class PROJECTION_RLE;

}
#endif

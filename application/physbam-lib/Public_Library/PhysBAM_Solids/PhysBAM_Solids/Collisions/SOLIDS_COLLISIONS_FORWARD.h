//#####################################################################
// Copyright 2006-2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header SOLIDS_COLLISIONS_FORWARD
//#####################################################################
#ifndef __SOLIDS_COLLISIONS_FORWARD__
#define __SOLIDS_COLLISIONS_FORWARD__

#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLES_COLLISIONS_FORWARD.h>
namespace PhysBAM{

template<class T> class TETRAHEDRON_COLLISION_BODY;
template<class TV> class COLLISION_PENALTY_FORCES;
template<class TV> class RIGID_DEFORMABLE_COLLISIONS;

template<class TV> struct PRECOMPUTE_PROJECT_POINT_FACE;
template<class TV> struct PRECOMPUTE_PROJECT_EDGE_EDGE;

}
#endif

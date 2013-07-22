//#####################################################################
// Copyright 2006-2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header DEFORMABLES_COLLISIONS_FORWARD
//#####################################################################
#ifndef __DEFORMABLES_COLLISIONS_FORWARD__
#define __DEFORMABLES_COLLISIONS_FORWARD__

#include <PhysBAM_Geometry/Collisions/COLLISIONS_GEOMETRY_FORWARD.h>
namespace PhysBAM{

template<class TV> struct POINT_FACE_REPULSION_PAIR;
template<class TV> struct EDGE_EDGE_REPULSION_PAIR;
template<class TV> class SEGMENT_REPULSIONS_AND_COLLISIONS_GEOMETRY;
template<class TV> class TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY;
template<class TV> class SEGMENT_REPULSIONS;
template<class TV> class TRIANGLE_REPULSIONS;
template<class TV> class SEGMENT_COLLISIONS;
template<class TV> class TRIANGLE_COLLISIONS;

}
#endif

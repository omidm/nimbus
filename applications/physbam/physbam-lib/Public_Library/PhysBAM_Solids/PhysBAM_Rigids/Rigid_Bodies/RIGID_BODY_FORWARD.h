//#####################################################################
// Copyright 2007-2008, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header RIGID_BODY_FORWARD
//#####################################################################
#ifndef __RIGID_BODY_FORWARD__
#define __RIGID_BODY_FORWARD__

namespace PhysBAM{

template<class TV> class RIGID_BODY;
template<class TV> class RIGID_BODY_STATE;
template<class TV> class RIGID_BODY_COLLISIONS;
template<class TV> class RIGID_BODY_COLLECTION;
template<class TV,bool world_space=false> class RIGID_BODY_MASS;
template<class TV> class RIGID_BODY_CONTACT_GRAPH;
template<class TV> class RIGID_BODY_INTERSECTIONS;
template<class TV> struct RIGID_BODY_PARTICLE_INTERSECTION;
template<class TV> class RIGID_BODY_EVOLUTION;
template<class TV> class RIGID_DEFORMABLE_COLLISIONS;
class RIGID_BODY_SKIP_COLLISION_CHECK;
class RIGID_BODY_COLLISION_MANAGER;

}
#endif

//#####################################################################
// Copyright 2009, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __GEOMETRY_PARTICLES_FORWARD__
#define __GEOMETRY_PARTICLES_FORWARD__

#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD_FORWARD.h>
namespace PhysBAM{
template<class TV> class GEOMETRY_PARTICLES;
const ATTRIBUTE_ID ATTRIBUTE_ID_V(2);
const ATTRIBUTE_ID ATTRIBUTE_ID_ROTATION(3);
const ATTRIBUTE_ID ATTRIBUTE_ID_ANGULAR_VELOCITY(4);
const ATTRIBUTE_ID ATTRIBUTE_ID_RIGID_GEOMETRY(5);
const ATTRIBUTE_ID ATTRIBUTE_ID_STRUCTURE_IDS(6);
const ATTRIBUTE_ID ATTRIBUTE_ID_COLLIDABLE(30);
const ATTRIBUTE_ID ATTRIBUTE_ID_COLOR(31);
const ATTRIBUTE_ID ATTRIBUTE_ID_DISPLAY_SIZE(32);
const ATTRIBUTE_ID ATTRIBUTE_ID_RADIUS(15);
const ATTRIBUTE_ID ATTRIBUTE_ID_BRADIUS(50);
const ATTRIBUTE_ID ATTRIBUTE_ID_BMASS(51);
const ATTRIBUTE_ID ATTRIBUTE_ID_PRESSURE(52);
const ATTRIBUTE_ID ATTRIBUTE_ID_RADIAL_VELOCITY(53);
const ATTRIBUTE_ID ATTRIBUTE_ID_WEIGHTS(54);
const ATTRIBUTE_ID ATTRIBUTE_ID_ANGULAR(34);
}
#endif

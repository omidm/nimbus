//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RIGIDS_PARTICLES_FORWARD__
#define __RIGIDS_PARTICLES_FORWARD__

#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY_PARTICLES;

const ATTRIBUTE_ID ATTRIBUTE_ID_RIGID_MASS(8);
const ATTRIBUTE_ID ATTRIBUTE_ID_ANGULAR_MOMENTUM(9);
const ATTRIBUTE_ID ATTRIBUTE_ID_RIGID_INERTIA_TENSOR(23);
const ATTRIBUTE_ID ATTRIBUTE_ID_KINEMATIC(24);

}
#endif

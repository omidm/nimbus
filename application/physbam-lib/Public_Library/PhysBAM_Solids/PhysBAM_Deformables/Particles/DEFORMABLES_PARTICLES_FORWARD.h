//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __DEFORMABLES_PARTICLES_FORWARD__
#define __DEFORMABLES_PARTICLES_FORWARD__

#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
namespace PhysBAM{

template<class TV> class PARTICLES;

const ATTRIBUTE_ID ATTRIBUTE_ID_MASS(7);
const ATTRIBUTE_ID ATTRIBUTE_ID_ONE_OVER_MASS(19);
const ATTRIBUTE_ID ATTRIBUTE_ID_EFFECTIVE_MASS(21);
const ATTRIBUTE_ID ATTRIBUTE_ID_ONE_OVER_EFFECTIVE_MASS(22);

}
#endif

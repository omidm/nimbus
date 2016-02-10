//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PARTICLES_FORWARD__
#define __PARTICLES_FORWARD__

#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLES_PARTICLES_FORWARD.h>
namespace PhysBAM{
void Initialize_Particles();

template<class TV> class SPH_PARTICLES;
template<class TV> class PARTICLE_LEVELSET_PARTICLES;
template<class TV> class PARTICLE_LEVELSET_REMOVED_PARTICLES;
template<class TV> class VORTICITY_PARTICLES;
template<class TV> class FLUID_PARTICLES;
template<class TV> class CODIMENSIONAL_FLUID_PARTICLES;

const ATTRIBUTE_ID ATTRIBUTE_ID_E(10);
const ATTRIBUTE_ID ATTRIBUTE_ID_RHO(11);
const ATTRIBUTE_ID ATTRIBUTE_ID_AGE(12);
const ATTRIBUTE_ID ATTRIBUTE_ID_PHI(13);
const ATTRIBUTE_ID ATTRIBUTE_ID_GRAD_PHI(14);
const ATTRIBUTE_ID ATTRIBUTE_ID_VORTICITY(16);
const ATTRIBUTE_ID ATTRIBUTE_ID_MATERIAL_VOLUME(17);
const ATTRIBUTE_ID ATTRIBUTE_ID_QUANTIZED_COLLISION_DISTANCE(18);
const ATTRIBUTE_ID ATTRIBUTE_ID_VOLUME(35);
const ATTRIBUTE_ID ATTRIBUTE_ID_FORCE(36);
const ATTRIBUTE_ID ATTRIBUTE_ID_VOLUME_RESIDUAL(37);
const ATTRIBUTE_ID ATTRIBUTE_ID_MOMENTUM_RESIDUAL(38);
const ATTRIBUTE_ID ATTRIBUTE_ID_ON_BOUNDARY(39);
const ATTRIBUTE_ID ATTRIBUTE_ID_NORMAL(40);
const ATTRIBUTE_ID ATTRIBUTE_ID_TYPE(41);
const ATTRIBUTE_ID ATTRIBUTE_ID_NODE(42);
}
#endif

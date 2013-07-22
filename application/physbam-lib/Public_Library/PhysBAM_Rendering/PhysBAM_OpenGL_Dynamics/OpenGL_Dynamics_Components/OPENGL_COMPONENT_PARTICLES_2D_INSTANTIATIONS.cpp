//#####################################################################
// Copyright 2009, Nipun Kwatra
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Particles/VORTICITY_PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_2D_DEFINITIONS.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/SPH_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/FLUID_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/BUBBLE_PARTICLES.h>

template class OPENGL_COMPONENT_PARTICLES_2D<float,SPH_PARTICLES<VECTOR<float,2> >,float>;
template class OPENGL_COMPONENT_PARTICLES_2D<float,BUBBLE_PARTICLES<VECTOR<float,2> >,float>;
template class OPENGL_COMPONENT_PARTICLES_2D<float,VORTICITY_PARTICLES<VECTOR<float,2> >,float>;
template class OPENGL_COMPONENT_PARTICLES_2D<float,PARTICLE_LEVELSET_PARTICLES<VECTOR<float,2> >,float>;
template class OPENGL_COMPONENT_PARTICLES_2D<float,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,2> >,float>;
template class OPENGL_COMPONENT_PARTICLES_2D<float,GEOMETRY_PARTICLES<VECTOR<float,2> >,float>;
template class OPENGL_COMPONENT_PARTICLES_2D<float,FLUID_PARTICLES<VECTOR<float,2> >,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_PARTICLES_2D<double,SPH_PARTICLES<VECTOR<double,2> >,double>;
template class OPENGL_COMPONENT_PARTICLES_2D<double,BUBBLE_PARTICLES<VECTOR<double,2> >,double>;
template class OPENGL_COMPONENT_PARTICLES_2D<double,VORTICITY_PARTICLES<VECTOR<double,2> >,double>;
template class OPENGL_COMPONENT_PARTICLES_2D<double,PARTICLE_LEVELSET_PARTICLES<VECTOR<double,2> >,double>;
template class OPENGL_COMPONENT_PARTICLES_2D<double,PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,2> >,double>;
template class OPENGL_COMPONENT_PARTICLES_2D<double,GEOMETRY_PARTICLES<VECTOR<double,2> >,double>;
template class OPENGL_COMPONENT_PARTICLES_2D<double,FLUID_PARTICLES<VECTOR<double,2> >,double>;
#endif

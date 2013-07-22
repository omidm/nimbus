//#####################################################################
// Copyright 2007-2008, Geoffrey Irving, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_CONNECTIVITY
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fragments/PARTICLE_CONNECTIVITY.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLE_CONNECTIVITY<TV>::
PARTICLE_CONNECTIVITY(const PARTICLES<TV>& particles,const RIGID_BODY_COLLECTION<TV>& rigid_body_collection)
    :particles_number(particles.array_collection->Size()),rigid_body_collection(rigid_body_collection),union_find(particles_number+rigid_body_collection.rigid_body_particle.array_collection->Size())
{}
//#####################################################################
// Function Exclude_Particle
//#####################################################################
template<class TV> bool PARTICLE_CONNECTIVITY<TV>::
Exclude_Particle(const int i) const
{   // TODO(jontg): This is a terrible reason to add dependencies to the rats-nest of dependencies that is RIGID_BODIES
    if(i<=particles_number) return false;
    if(!rigid_body_collection.Is_Active(i-particles_number)) return true;
    return !rigid_body_collection.Rigid_Body(i-particles_number).Is_Simulated();
}
//#####################################################################
// Function Union_All_Rigid_Body_Particles
//#####################################################################
template<class TV> void PARTICLE_CONNECTIVITY<TV>::
Union_All_Rigid_Body_Particles()
{
    Union(IDENTITY_ARRAY<>(rigid_body_collection.rigid_body_particle.array_collection->Size())+particles_number);
}
//#####################################################################
template class PARTICLE_CONNECTIVITY<VECTOR<float,1> >;
template class PARTICLE_CONNECTIVITY<VECTOR<float,2> >;
template class PARTICLE_CONNECTIVITY<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PARTICLE_CONNECTIVITY<VECTOR<double,1> >;
template class PARTICLE_CONNECTIVITY<VECTOR<double,2> >;
template class PARTICLE_CONNECTIVITY<VECTOR<double,3> >;
#endif
}

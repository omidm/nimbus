//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
using namespace PhysBAM;
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void GRAVITY<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    if(!gravity) return;
    TV acceleration=gravity*downward_direction;
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();F(p)+=particles.mass(p)*acceleration;}
    for(ELEMENT_ITERATOR iterator(force_rigid_body_particles);iterator.Valid();iterator.Next()){
        int p=iterator.Data();rigid_F(p).linear+=rigid_body_collection.rigid_body_particle.mass(p)*acceleration;}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR GRAVITY<TV>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    TV acceleration=gravity*downward_direction;
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
        potential_energy-=particles.mass(p)*TV::Dot_Product(particles.X(p),acceleration);}
    for(ELEMENT_ITERATOR iterator(force_rigid_body_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
        potential_energy-=rigid_body_collection.rigid_body_particle.mass(p)*TV::Dot_Product(rigid_body_collection.rigid_body_particle.X(p),acceleration);}
    return potential_energy;
}
//#####################################################################
template class GRAVITY<VECTOR<float,1> >;
template class GRAVITY<VECTOR<float,2> >;
template class GRAVITY<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRAVITY<VECTOR<double,1> >;
template class GRAVITY<VECTOR<double,2> >;
template class GRAVITY<VECTOR<double,3> >;
#endif

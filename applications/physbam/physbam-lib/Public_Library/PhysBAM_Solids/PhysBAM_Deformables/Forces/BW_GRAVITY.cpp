//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BW_GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void BW_GRAVITY<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    if(!gravity) return;
    TV acceleration=gravity*downward_direction;
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){
        int p=iterator.Data();F(p)+=particles.mass(p)*acceleration;}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR BW_GRAVITY<TV>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    TV acceleration=gravity*downward_direction;
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
        potential_energy-=particles.mass(p)*TV::Dot_Product(particles.X(p),acceleration);}
    return potential_energy;
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
//TODO: This needs to be fixed if particles are deleted
template<class TV> void BW_GRAVITY<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    force_particles.Update(Get_Particle_List(IDENTITY_ARRAY<>(particles.array_collection->Size())),particle_is_simulated);
}
//#####################################################################
template class BW_GRAVITY<VECTOR<float,1> >;
template class BW_GRAVITY<VECTOR<float,2> >;
template class BW_GRAVITY<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BW_GRAVITY<VECTOR<double,1> >;
template class BW_GRAVITY<VECTOR<double,2> >;
template class BW_GRAVITY<VECTOR<double,3> >;
#endif

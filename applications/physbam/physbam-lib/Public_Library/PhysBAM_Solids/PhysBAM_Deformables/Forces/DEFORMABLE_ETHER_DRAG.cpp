//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLE_ETHER_DRAG.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void DEFORMABLE_ETHER_DRAG<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int k=iterator.Data();
        if(use_spatially_varying_wind){
            if(spatially_varying_wind_domain.Lazy_Inside(particles.X(k))) F(k)+=spatially_varying_wind_viscosity*particles.mass(k)*Spatially_Varying_Wind_Velocity(particles.X(k));
            else if(use_constant_wind) F(k)+=constant_wind_viscosity*particles.mass(k)*constant_wind;}
        else if(use_constant_wind) F(k)+=constant_wind_viscosity*particles.mass(k)*constant_wind;}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void DEFORMABLE_ETHER_DRAG<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int k=iterator.Data();
        if(use_spatially_varying_wind){
            if(spatially_varying_wind_domain.Lazy_Inside(particles.X(k))) F(k)-=spatially_varying_wind_viscosity*particles.mass(k)*V(k);
            else if(use_constant_wind) F(k)-=constant_wind_viscosity*particles.mass(k)*V(k);}
        else if(use_constant_wind) F(k)-=constant_wind_viscosity*particles.mass(k)*V(k);}
}
//#####################################################################
template class DEFORMABLE_ETHER_DRAG<VECTOR<float,2> >;
template class DEFORMABLE_ETHER_DRAG<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DEFORMABLE_ETHER_DRAG<VECTOR<double,2> >;
template class DEFORMABLE_ETHER_DRAG<VECTOR<double,3> >;
#endif

//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/ETHER_DRAG.h>
using namespace PhysBAM;
template<class T_GRID> ETHER_DRAG<T_GRID>::
ETHER_DRAG(PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,ARRAY<int>* influenced_particles_input,
    ARRAY<int>* influenced_rigid_body_particles_input,T dynamic_ether_viscosity,T angular_viscosity)
    :POINTWISE_FORCE<TV>(particles_input,rigid_body_collection_input,influenced_particles_input,influenced_rigid_body_particles_input),use_constant_wind(dynamic_ether_viscosity!=0),
    constant_wind_viscosity(dynamic_ether_viscosity),constant_wind_angular_viscosity(angular_viscosity),use_spatially_varying_wind(false),spatially_varying_wind_viscosity(0)
{
}
template<class T_GRID> ETHER_DRAG<T_GRID>::
ETHER_DRAG(PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,const bool influence_all_particles_input,
    const bool influence_all_rigid_body_particles_input,T dynamic_ether_viscosity,T angular_viscosity)
    :POINTWISE_FORCE<TV>(particles_input,rigid_body_collection_input,influence_all_particles_input,influence_all_rigid_body_particles_input),use_constant_wind(dynamic_ether_viscosity!=0),
    constant_wind_viscosity(dynamic_ether_viscosity),constant_wind_angular_viscosity(angular_viscosity),use_spatially_varying_wind(false),spatially_varying_wind_viscosity(0)
{
}
template<class T_GRID> ETHER_DRAG<T_GRID>::
~ETHER_DRAG()
{
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T_GRID> void ETHER_DRAG<T_GRID>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int k=iterator.Data();
        if(use_spatially_varying_wind){
            if(spatially_varying_wind_domain.Lazy_Inside(particles.X(k))) F(k)+=spatially_varying_wind_viscosity*particles.mass(k)*Spatially_Varying_Wind_Velocity(particles.X(k));
            else if(use_constant_wind) F(k)+=constant_wind_viscosity*particles.mass(k)*constant_wind;}
        else if(use_constant_wind) F(k)+=constant_wind_viscosity*particles.mass(k)*constant_wind;}
    for(ELEMENT_ITERATOR iterator(force_rigid_body_particles);iterator.Valid();iterator.Next()){int k=iterator.Data();
        if(use_spatially_varying_wind){
            if(spatially_varying_wind_domain.Lazy_Inside(rigid_body_collection.rigid_body_particle.X(k)))
                rigid_F(k).linear+=spatially_varying_wind_viscosity*rigid_body_collection.rigid_body_particle.mass(k)*
                    Spatially_Varying_Wind_Velocity(rigid_body_collection.rigid_body_particle.X(k));
            else if(use_constant_wind) rigid_F(k).linear+=constant_wind_viscosity*rigid_body_collection.rigid_body_particle.mass(k)*constant_wind;}
        else if(use_constant_wind) rigid_F(k).linear+=constant_wind_viscosity*rigid_body_collection.rigid_body_particle.mass(k)*constant_wind;}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class T_GRID> void ETHER_DRAG<T_GRID>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int k=iterator.Data();
        if(use_spatially_varying_wind){
            if(spatially_varying_wind_domain.Lazy_Inside(particles.X(k))) F(k)-=spatially_varying_wind_viscosity*particles.mass(k)*V(k);
            else if(use_constant_wind) F(k)-=constant_wind_viscosity*particles.mass(k)*V(k);}
        else if(use_constant_wind) F(k)-=constant_wind_viscosity*particles.mass(k)*V(k);}
    for(ELEMENT_ITERATOR iterator(force_rigid_body_particles);iterator.Valid();iterator.Next()){int k=iterator.Data();
        if(use_spatially_varying_wind){
            if(spatially_varying_wind_domain.Lazy_Inside(rigid_body_collection.rigid_body_particle.X(k)))
                rigid_F(k).linear-=spatially_varying_wind_viscosity*rigid_body_collection.rigid_body_particle.mass(k)*rigid_V(k).linear;
            else if(use_constant_wind) rigid_F(k).linear-=constant_wind_viscosity*rigid_body_collection.rigid_body_particle.mass(k)*rigid_V(k).linear;}
        else if(use_constant_wind) rigid_F(k).linear-=constant_wind_viscosity*rigid_body_collection.rigid_body_particle.mass(k)*rigid_V(k).linear;
        if(constant_wind_angular_viscosity) rigid_F(k).angular-=constant_wind_angular_viscosity*rigid_V(k).angular;}
}
template<class T_GRID> void ETHER_DRAG<T_GRID>::
Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const
{
}
template<class T_GRID> void ETHER_DRAG<T_GRID>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
}
template<class T_GRID> void ETHER_DRAG<T_GRID>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
}
//#####################################################################
template class ETHER_DRAG<GRID<VECTOR<float,2> > >;
template class ETHER_DRAG<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ETHER_DRAG<GRID<VECTOR<double,2> > >;
template class ETHER_DRAG<GRID<VECTOR<double,3> > >;
#endif

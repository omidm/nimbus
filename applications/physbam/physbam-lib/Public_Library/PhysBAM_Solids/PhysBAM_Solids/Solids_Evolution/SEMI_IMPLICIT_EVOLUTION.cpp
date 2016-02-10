//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEMI_IMPLICIT_EVOLUTION
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SEMI_IMPLICIT_EVOLUTION.h>
using namespace PhysBAM;
#if 0
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class T> void SEMI_IMPLICIT_EVOLUTION<T>::
Advance_One_Time_Step(const T dt,const T time)
{
    PARTICLES<TV>& particles=deformable_object.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.binding_list;
    EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities=*deformable_object.example_forces_and_velocities;

    // do precomputation if necessary
    for(int k=1;k<=forces.m;k++){
        FINITE_VOLUME<TV,3>* fvm=dynamic_cast<FINITE_VOLUME<TV,3>*>(forces(k));if(!fvm) continue;
        if(!fvm->semi_implicit_data) fvm->Semi_Implicit_Impulse_Precomputation(time,FLT_MAX,0,true);}
    // update position
    Save_Position();
    Euler_Step_Position(dt,time);
    // update embedded positions
    binding_list.Clamp_Particles_To_Embedded_Positions();
    // update velocity with impulses
    for(int k=1;k<=forces.m;k++){
        FINITE_VOLUME<TV,3>* fvm=dynamic_cast<FINITE_VOLUME<TV,3>*>(forces(k));if(!fvm) continue;
        fvm->Add_Semi_Implicit_Impulses(dt,0);}
    example_forces_and_velocities.Add_External_Impulses(particles.V,time,dt);
    // collision body collisions
    if(solids_parameters.perform_collision_body_collisions){
        int interactions=deformable_object.collisions.Adjust_Nodes_For_Collision_Body_Collisions(binding_list,solid_body_collection.deformable_body_collection.soft_bindings,X_save,dt,0);
        if(interactions) {std::stringstream ss;ss<<"collision body collisions = "<<interactions<<std::endl;LOG::filecout(ss.str());}}
    // positions just changed, so update again
    example_forces_and_velocities.Update_Time_Varying_Material_Properties(time);
    // updated embedded velocities
    binding_list.Clamp_Particles_To_Embedded_Velocities();
    // output
    {std::stringstream ss;ss<<"maximum velocity = "<<ARRAY<TV>::Maximum_Magnitude(particles.V.Subset(dynamic_particles))<<std::endl;LOG::filecout(ss.str());}
}
//#####################################################################
template class SEMI_IMPLICIT_EVOLUTION<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SEMI_IMPLICIT_EVOLUTION<double>;
#endif
#endif

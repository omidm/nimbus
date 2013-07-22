//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/IDENTITY_ARRAY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGID_POINTWISE_FORCE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
using namespace PhysBAM;
//#####################################################################
// Function Get_Rigid_Body_Particle_List
//#####################################################################
template<class TV> template<class T_ARRAY> ARRAY<int> RIGID_POINTWISE_FORCE<TV>::
Get_Rigid_Body_Particle_List(const T_ARRAY& array)
{
    {std::stringstream ss;ss<<"attempting to add "<<array.Size()<<" rigid body particles"<<std::endl;LOG::filecout(ss.str());}
    ARRAY<int> list;
    for(int i=1;i<=array.Size();i++){
        int p=array(i);if(rigid_body_collection.Is_Active(p) && !rigid_body_collection.Rigid_Body(p).Has_Infinite_Inertia()) list.Append(p);}
    return list;
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void RIGID_POINTWISE_FORCE<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated)
{
    if(influence_all_rigid_body_particles){
        ARRAY<int> all_rigid=Get_Rigid_Body_Particle_List(IDENTITY_ARRAY<>(rigid_body_collection.rigid_body_particle.array_collection->Size()));
        force_rigid_body_particles.Update(all_rigid,particle_is_simulated);}
    else if(influenced_rigid_body_particles) force_rigid_body_particles.Update(*influenced_rigid_body_particles,particle_is_simulated);
}
//#####################################################################
template class RIGID_POINTWISE_FORCE<VECTOR<float,1> >;
template class RIGID_POINTWISE_FORCE<VECTOR<float,2> >;
template class RIGID_POINTWISE_FORCE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_POINTWISE_FORCE<VECTOR<double,1> >;
template class RIGID_POINTWISE_FORCE<VECTOR<double,2> >;
template class RIGID_POINTWISE_FORCE<VECTOR<double,3> >;
#endif

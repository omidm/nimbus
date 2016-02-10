//#####################################################################
// Copyright 2007-2008, Michael Lentine, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/IDENTITY_ARRAY.h>
#include <PhysBAM_Geometry/Topology/TETRAHEDRON_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/POINTWISE_FORCE.h>
using namespace PhysBAM;
template<class TV> POINTWISE_FORCE<TV>::
POINTWISE_FORCE(PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,ARRAY<int>* influenced_particles_input,
    ARRAY<int>* influenced_rigid_body_particles_input)
    :SOLIDS_FORCES<TV>(particles_input,rigid_body_collection_input),influenced_particles(influenced_particles_input),
    influenced_rigid_body_particles(influenced_rigid_body_particles_input),
    need_destroy_influenced_particles(false),need_destroy_influenced_rigid_body_particles(false),influence_all_particles(false),influence_all_rigid_body_particles(false),mpi_solids(0)
{
}
template<class TV> POINTWISE_FORCE<TV>::
POINTWISE_FORCE(PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,const bool influence_all_particles_input,
    const bool influence_all_rigid_body_particles_input)
    :SOLIDS_FORCES<TV>(particles_input,rigid_body_collection_input),influenced_particles(0),influenced_rigid_body_particles(0),need_destroy_influenced_particles(true),
    need_destroy_influenced_rigid_body_particles(true),influence_all_particles(influence_all_particles_input),influence_all_rigid_body_particles(influence_all_rigid_body_particles_input),
    mpi_solids(0)
{
}
template<class TV> template<class T_MESH> POINTWISE_FORCE<TV>::
POINTWISE_FORCE(PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,const T_MESH& mesh,ARRAY<int>* influenced_rigid_body_particles_input)
    :SOLIDS_FORCES<TV>(particles_input,rigid_body_collection_input),influenced_particles(new ARRAY<int>),
    influenced_rigid_body_particles(influenced_rigid_body_particles_input),
    need_destroy_influenced_particles(true),need_destroy_influenced_rigid_body_particles(false),influence_all_particles(false),influence_all_rigid_body_particles(false)
{
    mesh.elements.Flattened().Get_Unique(*influenced_particles);
}
template<class TV> POINTWISE_FORCE<TV>::
~POINTWISE_FORCE()
{
    if(need_destroy_influenced_particles) delete influenced_particles;
    if(need_destroy_influenced_rigid_body_particles) delete influenced_rigid_body_particles;
}
//#####################################################################
// Function Get_Rigid_Body_Particle_List
//#####################################################################
template<class TV> template<class T_ARRAY> ARRAY<int> POINTWISE_FORCE<TV>::
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
template<class TV> void POINTWISE_FORCE<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,const ARRAY<bool>& rigid_particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    if(influence_all_particles) force_particles.Update(Get_Particle_List(IDENTITY_ARRAY<>(particles.array_collection->Size())),particle_is_simulated);
    else if(influenced_particles) force_particles.Update(*influenced_particles,particle_is_simulated);

    if(influence_all_rigid_body_particles){
        ARRAY<int> all_rigid=Get_Rigid_Body_Particle_List(IDENTITY_ARRAY<>(rigid_body_collection.rigid_body_particle.array_collection->Size()));
        force_rigid_body_particles.Update(all_rigid,rigid_particle_is_simulated);}
    else if(influenced_rigid_body_particles) force_rigid_body_particles.Update(*influenced_rigid_body_particles,rigid_particle_is_simulated);
}
//#####################################################################
template class POINTWISE_FORCE<VECTOR<float,1> >;
template class POINTWISE_FORCE<VECTOR<float,2> >;
template class POINTWISE_FORCE<VECTOR<float,3> >;
template POINTWISE_FORCE<VECTOR<float,3> >::POINTWISE_FORCE(PARTICLES<VECTOR<float,3> >&,RIGID_BODY_COLLECTION<VECTOR<float,3> >&,TETRAHEDRON_MESH const&,ARRAY<int,int>*);
template POINTWISE_FORCE<VECTOR<float,3> >::POINTWISE_FORCE(PARTICLES<VECTOR<float,3> >&,RIGID_BODY_COLLECTION<VECTOR<float,3> >&,TRIANGLE_MESH const&,ARRAY<int,int>*);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class POINTWISE_FORCE<VECTOR<double,1> >;
template class POINTWISE_FORCE<VECTOR<double,2> >;
template class POINTWISE_FORCE<VECTOR<double,3> >;
template POINTWISE_FORCE<VECTOR<double,3> >::POINTWISE_FORCE(PARTICLES<VECTOR<double,3> >&,RIGID_BODY_COLLECTION<VECTOR<double,3> >&,TETRAHEDRON_MESH const&,ARRAY<int,int>*);
template POINTWISE_FORCE<VECTOR<double,3> >::POINTWISE_FORCE(PARTICLES<VECTOR<double,3> >&,RIGID_BODY_COLLECTION<VECTOR<double,3> >&,TRIANGLE_MESH const&,ARRAY<int,int>*);
#endif

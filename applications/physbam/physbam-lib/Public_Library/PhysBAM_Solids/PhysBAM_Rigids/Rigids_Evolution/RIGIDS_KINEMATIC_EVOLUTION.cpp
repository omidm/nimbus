//#####################################################################
// Copyright 2006-2009, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigids_Evolution/RIGIDS_KINEMATIC_EVOLUTION.h>
using namespace PhysBAM;
template<class TV> RIGIDS_KINEMATIC_EVOLUTION<TV>::
RIGIDS_KINEMATIC_EVOLUTION(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,bool use_kinematic_keyframes_input)
    :BASE(rigid_body_collection_input.rigid_geometry_collection,use_kinematic_keyframes_input),rigid_body_collection(rigid_body_collection_input)
{
}
template<class TV> RIGIDS_KINEMATIC_EVOLUTION<TV>::
~RIGIDS_KINEMATIC_EVOLUTION()
{
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class TV> void RIGIDS_KINEMATIC_EVOLUTION<TV>::
Set_External_Velocities(TV& V,T_SPIN& angular_velocity,const T time,const int id)
{
    assert(rigid_body_collection.rigid_body_particle.kinematic(id));
    BASE::Set_External_Velocities(V,angular_velocity,time,id);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
template<class TV> void RIGIDS_KINEMATIC_EVOLUTION<TV>::
Set_Kinematic_Velocities(TV& V,T_SPIN& angular_velocity,const T frame_dt,const T time,const int id)
{
    RIGID_BODY<TV>* rigid_body=&rigid_body_collection.Rigid_Body(id);
    assert(rigid_body_collection.rigid_body_particle.kinematic(id));
    int index=rigid_body->particle_index;
    int new_id=id;
    if(rigid_body_collection.rigid_body_cluster_bindings.Is_Parent(index)){
        typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER* cluster=rigid_body_collection.rigid_body_cluster_bindings.reverse_bindings.Get(index);
        rigid_body=&rigid_body_collection.Rigid_Body(cluster->infinite_body);
        new_id=rigid_body->particle_index;}
    BASE::Set_Kinematic_Velocities(V,angular_velocity,frame_dt,time,new_id);
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
template<class TV> void RIGIDS_KINEMATIC_EVOLUTION<TV>::
Set_External_Positions(TV& X,ROTATION<TV>& rotation,const T time,const int id)
{
    RIGID_BODY<TV>* rigid_body=&rigid_body_collection.Rigid_Body(id);
    assert(rigid_body_collection.rigid_body_particle.kinematic(id));
    int new_id=id;
    int index=rigid_body->particle_index;
    if(rigid_body_collection.rigid_body_cluster_bindings.Is_Parent(index)){
        typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER* cluster=rigid_body_collection.rigid_body_cluster_bindings.reverse_bindings.Get(index);
        rigid_body=&rigid_body_collection.Rigid_Body(cluster->infinite_body);
        new_id=rigid_body->particle_index;}
    BASE::Set_External_Positions(X,rotation,time,new_id);
}
//#####################################################################
template class RIGIDS_KINEMATIC_EVOLUTION<VECTOR<float,1> >;
template class RIGIDS_KINEMATIC_EVOLUTION<VECTOR<float,2> >;
template class RIGIDS_KINEMATIC_EVOLUTION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGIDS_KINEMATIC_EVOLUTION<VECTOR<double,1> >;
template class RIGIDS_KINEMATIC_EVOLUTION<VECTOR<double,2> >;
template class RIGIDS_KINEMATIC_EVOLUTION<VECTOR<double,3> >;
#endif

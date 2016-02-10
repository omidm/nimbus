//#####################################################################
// Copyright 2007, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_FRACTURE_QUASISTATICS_FORCES
//#####################################################################
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_IMPULSE_ACCUMULATOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/RIGID_BODY_FRACTURE_OBJECT_3D.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/RIGID_FRACTURE_QUASISTATICS_FORCES.h>

using namespace PhysBAM;
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void RIGID_FRACTURE_QUASISTATICS_FORCES<T>::
Initialize(const RIGID_BODY_FRACTURE_OBJECT_3D<T>& fracture_object_input,const VECTOR<int,3>& null_space_nodes_input)
{
    initialized=true;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=fracture_object_input.rigid_body_collection;
    COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR<TV>* collision_body_impulse_accumulator=
        rigid_body_collection.rigid_geometry_collection.collision_body_list->Get_Collision_Geometry(fracture_object_input.particle_index)->impulse_accumulator;

    fracture_object=&fracture_object_input;impulses=dynamic_cast<RIGID_BODY_IMPULSE_ACCUMULATOR<TV,3>&>(*collision_body_impulse_accumulator).accumulated_node_impulses;

    const PARTICLES<TV>& particles=fracture_object->particles;
    null_space_nodes=null_space_nodes_input;

    null_space_position=fracture_object->Frame()*particles.X(null_space_nodes(1));
    TV x1_x2_direction=fracture_object->Frame()*particles.X(null_space_nodes(2))-null_space_position,x1_x3_direction=fracture_object->Frame()*particles.X(null_space_nodes(3))-null_space_position;
    direction1=x1_x2_direction.Normalized();

    //make sure x1_x2 and x1_x3 are not collinear
    while(direction1==x1_x3_direction.Normalized()&& null_space_nodes(3)<particles.array_collection->Size()){null_space_nodes(3)++;x1_x3_direction=fracture_object->Frame()*particles.X(null_space_nodes(3))-null_space_position;}
    direction2=TV::Cross_Product(x1_x2_direction,x1_x3_direction).Normalized();

    // New way of removing null space
    TV q1=TV::Cross_Product(x1_x3_direction,x1_x2_direction).Normalized();
    TV q2=TV::Cross_Product(q1,x1_x2_direction).Normalized();
    TV q3=TV::Cross_Product(q1,q2);
    dx2_transform=MATRIX<T,3>::Outer_Product(q3,q3);
    dx3_transform=MATRIX<T,3>::Identity_Matrix()-MATRIX<T,3>::Outer_Product(q1,q1);
}
//#####################################################################
// Function Zero_Out_Enslaved_Position_Nodes
//#####################################################################
template<class T> void RIGID_FRACTURE_QUASISTATICS_FORCES<T>::
Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> dX,const T time)
{
    if(!initialized) return;
    dX(null_space_nodes(1))=TV();
    dX(null_space_nodes(2)).Projected_On_Unit_Direction(direction2);
    dX(null_space_nodes(3)).Projected_Orthogonal_To_Unit_Direction(direction1);
}
//#####################################################################
// Function Add_External_Forces
//#####################################################################
template<class T> void RIGID_FRACTURE_QUASISTATICS_FORCES<T>::
Add_External_Forces(ARRAY_VIEW<TV> F,const T time)
{
    if(!initialized) return;
    for(int i=1;i<=fracture_object->rigid_to_deformable_particles.m;i++)F(fracture_object->rigid_to_deformable_particles(i))+=fracture_object->World_Space_Vector((*impulses)(i));
    //F+=*impulses;
}
//#####################################################################
template class RIGID_FRACTURE_QUASISTATICS_FORCES<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_FRACTURE_QUASISTATICS_FORCES<double>;
#endif

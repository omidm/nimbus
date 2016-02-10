//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Joseph Teran, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COLLISION_PENALTY_FORCES
//#####################################################################
#ifndef __COLLISION_PENALTY_FORCES__
#define __COLLISION_PENALTY_FORCES__

#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TETRAHEDRON_COLLISION_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

using ::std::exp;

template<class TV_input>
class COLLISION_PENALTY_FORCES:public DEFORMABLES_FORCES<TV_input>
{
    typedef TV_input TV;typedef typename TV::SCALAR T;
    typedef DEFORMABLES_FORCES<TV_input> BASE;
    typedef typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX T_SYMMETRIC_MATRIX;
public:
    using BASE::particles;

    COLLISION_GEOMETRY_COLLECTION<TV>* collision_body_list;
    ARRAY<bool,COLLISION_GEOMETRY_ID> skip_collision_body;
    ARRAY<int> check_collision; // TODO: If Append is being called on this, why isn't this a ARRAY?
    ARRAY<TV> collision_force;
    ARRAY<T_SYMMETRIC_MATRIX> collision_force_derivative;
    T stiffness;
    T separation_parameter;
    T self_collision_reciprocity_factor;
    COLLISION_GEOMETRY_ID collision_body_list_id;
    bool perform_self_collision;

    COLLISION_PENALTY_FORCES(PARTICLES<TV>& particles)
        :DEFORMABLES_FORCES<TV>(particles),collision_body_list_id(0)
    {
        Set_Stiffness();Set_Separation_Parameter();Set_Self_Collision_Reciprocity_Factor();Set_Perform_Self_Collision();
    }

    ~COLLISION_PENALTY_FORCES()
    {}

    void Set_Perform_Self_Collision(const bool perform_self_collision_input=true)
    {perform_self_collision=perform_self_collision_input;}

    void Set_Stiffness(const T stiffness_input=(T)1e4)
    {stiffness=stiffness_input;}

    void Set_Separation_Parameter(const T separation_parameter_input=(T)1e-4)
    {separation_parameter=separation_parameter_input;}

    void Set_Collision_Body_List(COLLISION_GEOMETRY_COLLECTION<TV>& collision_body_list_input)
    {collision_body_list=&collision_body_list_input;skip_collision_body.Resize(collision_body_list->collision_bodies.m);ARRAYS_COMPUTATIONS::Fill(skip_collision_body,false);}

    template<class T_MESH>
    void Set_Boundary_Only_Collisions(T_MESH& mesh)
    {ARRAY<bool>* old_node_on_boundary=mesh.node_on_boundary;mesh.node_on_boundary=0;
    mesh.Initialize_Node_On_Boundary();for(int p=1;p<=particles.array_collection->Size();p++) if((*mesh.node_on_boundary)(p)) check_collision.Append(p);
    collision_force.Resize(check_collision.m);collision_force_derivative.Resize(check_collision.m);
    delete mesh.node_on_boundary;mesh.node_on_boundary=old_node_on_boundary;}

    void Set_Collision_Body_List_ID(const COLLISION_GEOMETRY_ID id)
    {collision_body_list_id=id;}

    void Set_Self_Collision_Reciprocity_Factor(const T self_collision_reciprocity_factor_input=(T)2)
    {self_collision_reciprocity_factor=self_collision_reciprocity_factor_input;}

    void Check_Collision(const int particle)
    {check_collision.Append_Unique(particle);} // TODO: not very efficient.

    void Omit_Collision(const int particle)
    {int index;if(check_collision.Find(particle,index)) check_collision.Remove_Index(index);} // TODO: not very efficient.

    void Resize_Collision_Arrays_From_Check_Collision()
    {collision_force.Resize(check_collision.m);collision_force_derivative.Resize(check_collision.m);}

    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE
    {}

    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE
    {}

    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE // Currently works only with a single fragment
    {if(collision_body_list_id && typeid((*collision_body_list)(collision_body_list_id))==typeid(TETRAHEDRON_COLLISION_BODY<T>)){
        TETRAHEDRON_COLLISION_BODY<T>& collision_body=(TETRAHEDRON_COLLISION_BODY<T>&)((*collision_body_list)(collision_body_list_id));
        collision_body.tetrahedralized_volume.hierarchy->Update_Boxes(collision_body.collision_thickness);
        collision_body.tetrahedralized_volume.triangulated_surface->hierarchy->Update_Boxes();
        collision_body.tetrahedralized_volume.triangulated_surface->Update_Vertex_Normals();}}

    void Update_Forces_And_Derivatives() // Currently works only with a single fragment
    {for(int p=1;p<=check_collision.m;p++){
        int index=check_collision(p);collision_force(p)=TV();collision_force_derivative(p)=T_SYMMETRIC_MATRIX();
        for(COLLISION_GEOMETRY_ID r(1);r<=collision_body_list->bodies.m;r++) if(collision_body_list->Is_Active(r)){
            COLLISION_GEOMETRY<TV>& collision_body=*collision_body_list->bodies(r);
            if(!skip_collision_body(r) && (perform_self_collision || collision_body.collision_geometry_id!=collision_body_list_id)){
                int collision_body_particle_index=0;if(collision_body_list_id==collision_body.collision_geometry_id) collision_body_particle_index=index;
                T phi_value;int aggregate=-1;TV normal=collision_body.Implicit_Geometry_Extended_Normal(particles.X(index),phi_value,aggregate,collision_body_particle_index);
                T scaled_stiffness=stiffness;if(collision_body_list_id==collision_body.collision_geometry_id) scaled_stiffness*=self_collision_reciprocity_factor;
                if(phi_value<=0){
                    collision_force(p)+=stiffness*(-phi_value+separation_parameter)*normal;
                    collision_force_derivative(p)-=scaled_stiffness*T_SYMMETRIC_MATRIX::Outer_Product(normal);}
                else if(phi_value<collision_body.collision_thickness){
                    collision_force(p)+=stiffness*separation_parameter*exp(-phi_value/separation_parameter)*normal;
                    collision_force_derivative(p)-=scaled_stiffness*exp(-phi_value/separation_parameter)*T_SYMMETRIC_MATRIX::Outer_Product(normal);}}}}}

    // Currently works only with a single fragment
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE
    {for(int p=1;p<=check_collision.m;p++) F(check_collision(p))+=collision_force(p);}

    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}

    void Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const PHYSBAM_OVERRIDE // Currently works only with a single fragment
    {for(int p=1;p<=check_collision.m;p++) dF(check_collision(p))+=collision_force_derivative(p)*dX(check_collision(p));}

    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}

//#####################################################################
};
}
#endif

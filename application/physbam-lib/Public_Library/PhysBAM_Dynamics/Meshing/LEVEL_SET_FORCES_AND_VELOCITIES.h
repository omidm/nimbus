//#####################################################################
// Copyright 2002-2007, Christopher Allocco, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVEL_SET_FORCES_AND_VELOCITIES
//#####################################################################
#ifndef __LEVEL_SET_FORCES_AND_VELOCITIES__
#define __LEVEL_SET_FORCES_AND_VELOCITIES__

#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
namespace PhysBAM{

template<class TV>
class LEVEL_SET_FORCES_AND_VELOCITIES:public EXAMPLE_FORCES_AND_VELOCITIES<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m>::OBJECT T_MESH_OBJECT;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_BOUNDARY_OBJECT;
    typedef EXAMPLE_FORCES_AND_VELOCITIES<TV> BASE;
    using BASE::Add_External_Forces;using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes; // silence -Woverloaded-virtual
public:
    T_MESH_OBJECT& mesh_object;
    IMPLICIT_OBJECT<TV>& implicit_object;
    T force_attraction_coefficient,velocity_attraction_coefficient;
    bool use_external_forces,use_external_velocities;
    bool use_external_velocities_normal_to_boundary;
    bool allow_tangential_velocity_slip;

    LEVEL_SET_FORCES_AND_VELOCITIES(T_MESH_OBJECT& mesh_object,IMPLICIT_OBJECT<TV>& implicit_object)
        :mesh_object(mesh_object),implicit_object(implicit_object)
    {
        Set_Force_Attraction_Coefficient();Set_Velocity_Attraction_Coefficient();
        Use_External_Forces();
        Use_External_Velocities_Normal_To_Boundary();
        Allow_Tangential_Velocity_Slip();
    }

    ~LEVEL_SET_FORCES_AND_VELOCITIES()
    {}

    void Set_Force_Attraction_Coefficient(const T coefficient=.1)
    {force_attraction_coefficient=coefficient;}

    void Set_Velocity_Attraction_Coefficient(const T coefficient=.1)
    {velocity_attraction_coefficient=coefficient;}

    void Use_External_Forces()
    {use_external_forces=true;use_external_velocities=false;}

    void Use_External_Velocities()
    {use_external_forces=false;use_external_velocities=true;}

    void Use_External_Velocities_Normal_To_Boundary()
    {use_external_velocities_normal_to_boundary=true;}

    void Use_External_Velocities_Towards_Zero_Isocontour()
    {use_external_velocities_normal_to_boundary=false;}

    void Allow_Tangential_Velocity_Slip(const bool allow_tangential_velocity_slip_input=true)
    {allow_tangential_velocity_slip=allow_tangential_velocity_slip_input;}

//#####################################################################
// Function Add_External_Forces
//#####################################################################
void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE
{
    if(!use_external_forces) return;

    PARTICLES<TV>& particles=dynamic_cast<PARTICLES<TV>&>(mesh_object.particles);
    const ARRAY<int>& boundary_nodes=*mesh_object.mesh.boundary_nodes;
    T_BOUNDARY_OBJECT& boundary_object=mesh_object.Get_Boundary_Object();

    boundary_object.Update_Vertex_Normals();
    for(int i=1;i<=boundary_nodes.m;i++){int p=boundary_nodes(i);
        F(p)-=force_attraction_coefficient*sqrt(particles.mass(p))*implicit_object(particles.X(p))*(*boundary_object.vertex_normals)(p);}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(!use_external_velocities) return;

    GEOMETRY_PARTICLES<TV>& particles=mesh_object.particles;
    const ARRAY<int>& boundary_nodes=*mesh_object.mesh.boundary_nodes;

    if(allow_tangential_velocity_slip){
        if(use_external_velocities_normal_to_boundary){
            T_BOUNDARY_OBJECT& boundary_object=mesh_object.Get_Boundary_Object();
            boundary_object.Update_Vertex_Normals();
            for(int i=1;i<=boundary_nodes.m;i++){int p=boundary_nodes(i);
                TV N=(*boundary_object.vertex_normals)(p);
                V(p)-=(TV::Dot_Product(V(p),N)+velocity_attraction_coefficient*implicit_object(particles.X(p)))*N;}}
        else // use external velocities normal to level set
            for(int i=1;i<=boundary_nodes.m;i++){int p=boundary_nodes(i);
                TV N=implicit_object.Normal(particles.X(p));
                V(p)-=(TV::Dot_Product(V(p),N)+velocity_attraction_coefficient*implicit_object(particles.X(p)))*N;}}
    else{
        if(use_external_velocities_normal_to_boundary){
            T_BOUNDARY_OBJECT& boundary_object=mesh_object.Get_Boundary_Object();
            boundary_object.Update_Vertex_Normals();
            for(int i=1;i<=boundary_nodes.m;i++){int p=boundary_nodes(i);
                TV N=(*boundary_object.vertex_normals)(p);
                V(p)=-velocity_attraction_coefficient*implicit_object(particles.X(p))*N;}}
        else // use external velocities normal to level set
            for(int i=1;i<=boundary_nodes.m;i++){int p=boundary_nodes(i);
                TV N=implicit_object.Normal(particles.X(p));
                V(p)=-velocity_attraction_coefficient*implicit_object(particles.X(p))*N;}}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(!use_external_velocities) return;

    GEOMETRY_PARTICLES<TV>& particles=mesh_object.particles;
    const ARRAY<int>& boundary_nodes=*mesh_object.mesh.boundary_nodes;

    if(allow_tangential_velocity_slip){
        if(use_external_velocities_normal_to_boundary){
            T_BOUNDARY_OBJECT& boundary_object=mesh_object.Get_Boundary_Object();
            boundary_object.Update_Vertex_Normals();
            for(int i=1;i<=boundary_nodes.m;i++){int p=boundary_nodes(i);
                TV N=(*boundary_object.vertex_normals)(p);
                V(p)-=TV::Dot_Product(V(p),N)*N;}}
        else // use external velocities normal to level set
            for(int i=1;i<=boundary_nodes.m;i++){int p=boundary_nodes(i);
                TV N=implicit_object.Normal(particles.X(p));
                V(p)-=TV::Dot_Product(V(p),N)*N;}}
    else for(int i=1;i<=boundary_nodes.m;i++){int p=boundary_nodes(i);
        V(p)=TV();}
}
//#####################################################################
};
}
#endif



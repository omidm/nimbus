//#####################################################################
// Copyright 2006-2009, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_STANDARD_TESTS
//#####################################################################
#include <PhysBAM_Tools/Ordinary_Differential_Equations/EXAMPLE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/SEGMENTED_CURVE_2D_SIGNED_DISTANCE.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/MASS_PROPERTIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Bindings/RIGID_BODY_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLIDS_STANDARD_TESTS<TV>::
SOLIDS_STANDARD_TESTS(EXAMPLE<TV>& example_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input)
    :BASE(example_input,solid_body_collection_input.deformable_body_collection),RIGIDS_STANDARD_TESTS<TV>(example_input.data_directory,example_input.stream_type,solid_body_collection_input.rigid_body_collection),
    solid_body_collection(solid_body_collection_input)
{
}
//#####################################################################
// Function Add_Gravity
//#####################################################################
template<class TV> void SOLIDS_STANDARD_TESTS<TV>::
Add_Gravity()
{
    // add gravity on all deformable particles and on all rigid body particles
    solid_body_collection.Add_Force(new GRAVITY<TV>(solid_body_collection.deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
}
//#####################################################################
// Function Bind_Particles_In_Rigid_Body
//#####################################################################
template<class TV> void SOLIDS_STANDARD_TESTS<TV>::
Bind_Particles_In_Rigid_Body(RIGID_BODY<TV>& rigid_body)
{
    Bind_Particles_In_Rigid_Body(rigid_body,IDENTITY_ARRAY<>(solid_body_collection.deformable_body_collection.particles.array_collection->Size()));
}
//#####################################################################
// Function Bind_Unbound_Particles_In_Rigid_Body
//#####################################################################
template<class TV> template<class T_ARRAY> void SOLIDS_STANDARD_TESTS<TV>::
Bind_Unbound_Particles_In_Rigid_Body(RIGID_BODY<TV>& rigid_body,const T_ARRAY& particle_array)
{
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    HASHTABLE<int> exempt_particles;
    ARRAY<int> binding_candidates;
    for(int i=1;i<=binding_list.bindings.m;i++) if(!dynamic_cast<RIGID_BODY_BINDING<TV>*>(binding_list.bindings(i))) exempt_particles.Set_All(binding_list.bindings(i)->Parents());
    for(int i=1;i<=particle_array.Size();i++) if(exempt_particles.Set(particle_array(i))) binding_candidates.Append(particle_array(i));  // Ignore duplicates
    Bind_Particles_In_Rigid_Body(rigid_body,binding_candidates);
}
//#####################################################################
// Function Bind_Particles_In_Rigid_Body
//#####################################################################
template<class TV> template<class T_ARRAY> void SOLIDS_STANDARD_TESTS<TV>::
Bind_Particles_In_Rigid_Body(RIGID_BODY<TV>& rigid_body,const T_ARRAY& particle_array)
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    for(typename T_ARRAY::ELEMENT i(1);i<=particle_array.Size();i++){int p=particle_array(i);
        if(rigid_body.implicit_object->Inside(particles.X(p)))
            binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(particles,p,rigid_body_collection,rigid_body.particle_index,rigid_body.Object_Space_Point(particles.X(p))));}
}
//#####################################################################
// Function PD_Curl
//#####################################################################
template<class TV> void SOLIDS_STANDARD_TESTS<TV>::
PD_Curl(const T scale,const TV shift,const ROTATION<TV> orient,const T k_p,const int number_of_joints,const bool parent_static,const T friction)
{
    PHYSBAM_ASSERT(scale>0);
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    RIGID_BODY<TV> *parent_body=0,*child_body=0;
    T cheight=(T)0;

    // Create first body
    parent_body=&Add_Rigid_Body("miniplank25wide2",scale,friction);
    parent_body->X()=shift;
    parent_body->Rotation()=orient;
    parent_body->Set_Name("parent");
    parent_body->Set_Mass(50);
    parent_body->is_static=parent_static;

    // Add children and joints
    T desired_x=(T)two_pi/(T)(number_of_joints+1);
    for(int i=0;i<number_of_joints;i++){
        cheight+=scale*(T)1.25;
        child_body=&Add_Rigid_Body("miniplank25wide2",scale,friction);
        child_body->X()=orient.Rotate(TV(cheight,0,0))+shift;
        child_body->Rotation()=orient;
        child_body->Set_Coefficient_Of_Restitution((T)0.5);
        child_body->Set_Name(STRING_UTILITIES::string_sprintf("child_%d",i));
        child_body->Set_Mass(50);

        ROTATION<TV> desired_rotation(desired_x,TV(1,0,0));

        JOINT<TV>* joint=new POINT_JOINT<TV>();arb.joint_mesh.Add_Articulation(child_body->particle_index-1,child_body->particle_index,joint);
        JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
        joint_function->Set_k_p(k_p);joint_function->Set_Target_Angle(desired_rotation);
        joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(scale*(T).625,0,0),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
        joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-scale*(T).625,0,0),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));

        parent_body=child_body;}
}
//#####################################################################
// Function Create_From_Tetrahedralized_Volume
//#####################################################################
template<class TV> RIGID_BODY<TV>* SOLIDS_STANDARD_TESTS<TV>::
Create_Rigid_Body_From_Tetrahedralized_Volume(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const T density,
    const T cell_size,const int subdivision_loops,bool perform_manifold_check,const bool (*create_levelset_test)(TETRAHEDRALIZED_VOLUME<T>&),
    const bool use_levelset_maker,const int levels_of_octree)
{
    PHYSBAM_ASSERT(density>0);
    RIGID_BODY<TV>* rigid_body=new RIGID_BODY<TV>(rigid_body_collection,true);
    if(!tetrahedralized_volume.triangulated_surface) tetrahedralized_volume.Initialize_Triangulated_Surface();
    MASS_PROPERTIES<TV> mass_properties(*tetrahedralized_volume.triangulated_surface,true);
    mass_properties.Set_Density(density);rigid_body->Mass().mass=mass_properties.Mass();
    FRAME<TV> frame_local;
    mass_properties.Transform_To_Object_Frame(frame_local,rigid_body->Mass().inertia_tensor,dynamic_cast<PARTICLES<TV>&>(tetrahedralized_volume.particles));
    rigid_body->Set_Frame(frame_local);
    rigid_body->Initialize_From_Tetrahedralized_Volume_And_Triangulated_Surface(tetrahedralized_volume,*tetrahedralized_volume.triangulated_surface,cell_size,subdivision_loops,
        create_levelset_test,use_levelset_maker,levels_of_octree);
    return rigid_body;
}
//#####################################################################
// Function Create_From_Fracture_Tetrahedralized_Volume
//#####################################################################
template<class TV> RIGID_BODY<TV>* SOLIDS_STANDARD_TESTS<TV>::
Create_Rigid_Body_From_Fracture_Tetrahedralized_Volume(EMBEDDED_MATERIAL_SURFACE<TV,3>& embedded_material_surface,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,
    const T density,const T cell_size,const int subdivision_loops,const bool perform_manifold_check,const bool (*create_levelset_test)(TETRAHEDRALIZED_VOLUME<T>&),
    const bool use_levelset_maker,const int levels_of_octree)
{
    PHYSBAM_ASSERT(density>0);
    RIGID_BODY<TV>* rigid_body=new RIGID_BODY<TV>(rigid_body_collection,true);
    if(!embedded_material_surface.material_surface_mesh.elements.m) embedded_material_surface.Create_Material_Surface();
    TRIANGULATED_SURFACE<T>& material_surface=embedded_material_surface.material_surface;
    if(perform_manifold_check){
        ARRAY<int> non_manifold_nodes;material_surface.mesh.Non_Manifold_Nodes(non_manifold_nodes);
        if(non_manifold_nodes.m) return 0;}
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=embedded_material_surface.embedded_object.simplicial_object;
    MASS_PROPERTIES<TV> mass_properties(material_surface,false);
    mass_properties.Set_Density(density);rigid_body->Mass().mass=mass_properties.Mass();
    FRAME<TV> frame_local;
    mass_properties.Transform_To_Object_Frame(frame_local,rigid_body->Mass().inertia_tensor,dynamic_cast<PARTICLES<TV>&>(tetrahedralized_volume.particles));
    rigid_body->Set_Frame(frame_local);
    embedded_material_surface.embedded_object.Update_Embedded_Particle_Positions();
    rigid_body->Initialize_From_Tetrahedralized_Volume_And_Triangulated_Surface(tetrahedralized_volume,material_surface,cell_size,subdivision_loops,create_levelset_test,
        use_levelset_maker,levels_of_octree);
    return rigid_body;
}
//#####################################################################
// Function Create_From_Triangulated_Surface
//#####################################################################
template<class TV> RIGID_BODY<TV>* SOLIDS_STANDARD_TESTS<TV>::
Create_Rigid_Body_From_Triangulated_Surface(TRIANGULATED_SURFACE<T>& triangulated_surface,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const T density)
{
    PHYSBAM_ASSERT(density>0);
    RIGID_BODY<TV>* rigid_body=new RIGID_BODY<TV>(rigid_body_collection,true);
    MASS_PROPERTIES<TV> mass_properties(triangulated_surface,true);
    mass_properties.Set_Density(density);rigid_body->Mass().mass=mass_properties.Mass();
    FRAME<TV> frame_local;
    mass_properties.Transform_To_Object_Frame(frame_local,rigid_body->Mass().inertia_tensor,dynamic_cast<PARTICLES<TV>&>(triangulated_surface.particles));
    rigid_body->Set_Frame(frame_local);
    rigid_body->Add_Structure(triangulated_surface);
    return rigid_body;
}
//#####################################################################
// Function Create_From_Triangulated_Area
//#####################################################################
template<class TV> RIGID_BODY<TV>* SOLIDS_STANDARD_TESTS<TV>::
Create_Rigid_Body_From_Triangulated_Area(TRIANGULATED_AREA<T>& triangulated_area,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const T density,const T cell_size,
    const bool move_body_to_center_of_mass,const bool move_only_mesh_particles)
{
    PHYSBAM_ASSERT(density>0);
    RIGID_BODY<TV>* rigid_body=new RIGID_BODY<TV>(rigid_body_collection,true);

    if(!triangulated_area.segmented_curve) triangulated_area.Initialize_Segmented_Curve();
    MASS_PROPERTIES<TV> mass_properties(*triangulated_area.segmented_curve,true);
    mass_properties.Set_Density(density);
    VECTOR<T,2> center_of_mass=mass_properties.Center();
    rigid_body->Mass().mass=mass_properties.Mass();
    rigid_body->Mass().inertia_tensor=mass_properties.Inertia_Tensor();
    
    // Translate triangulated area to be in object space (centered at center of mass)
    if(move_only_mesh_particles){
        // mark which nodes are used
        ARRAY<bool> node_is_used(triangulated_area.mesh.number_nodes);
        for(int i=1;i<=triangulated_area.mesh.elements.m;i++){
            int node1,node2,node3;triangulated_area.mesh.elements(i).Get(node1,node2,node3);
            node_is_used(node1)=node_is_used(node2)=node_is_used(node3)=true;}        
        // shift only used particles
        for(int i=1;i<=triangulated_area.mesh.number_nodes;i++) if(node_is_used(i)) triangulated_area.particles.X(i)-=center_of_mass;}
    else triangulated_area.particles.X-=center_of_mass;

    if(move_body_to_center_of_mass) rigid_body->Frame().t=center_of_mass;

    rigid_body->Add_Structure(triangulated_area);

    // Make a copy of the boundary segmented curve for the rigid body
    SEGMENTED_CURVE_2D<T>* segmented_curve=SEGMENTED_CURVE_2D<T>::Create();
    segmented_curve->mesh.Initialize_Mesh(triangulated_area.segmented_curve->mesh);
    segmented_curve->particles.Initialize(triangulated_area.segmented_curve->particles);
    segmented_curve->Discard_Valence_Zero_Particles_And_Renumber();
    rigid_body->Add_Structure(*segmented_curve);

    LEVELSET_IMPLICIT_OBJECT<TV>* implicit_curve=LEVELSET_IMPLICIT_OBJECT<TV>::Create();
    implicit_curve->levelset.grid=GRID<TV>::Create_Grid_Given_Cell_Size(*segmented_curve->bounding_box,cell_size,false,5);
    SIGNED_DISTANCE::Calculate(*segmented_curve,implicit_curve->levelset.grid,implicit_curve->levelset.phi);

    rigid_body->Add_Structure(*implicit_curve);

    return rigid_body;
}
//#####################################################################
#define INSTANTIATION_HELPER2_1D(T) \
    template SOLIDS_STANDARD_TESTS<VECTOR<T,1> >::SOLIDS_STANDARD_TESTS(EXAMPLE<VECTOR<T,1> >&,SOLID_BODY_COLLECTION<VECTOR<T,1> >&); \
    template void SOLIDS_STANDARD_TESTS<VECTOR<T,1> >::Add_Gravity();

#define INSTANTIATION_HELPER2_ALL(T,d) \
    template SOLIDS_STANDARD_TESTS<VECTOR<T,d> >::SOLIDS_STANDARD_TESTS(EXAMPLE<VECTOR<T,d> >&,SOLID_BODY_COLLECTION<VECTOR<T,d> >&); \
    template void SOLIDS_STANDARD_TESTS<VECTOR<T,d> >::Bind_Particles_In_Rigid_Body(RIGID_BODY<VECTOR<T,d> >& rigid_body); \
    template void SOLIDS_STANDARD_TESTS<VECTOR<T,d> >::Bind_Unbound_Particles_In_Rigid_Body(RIGID_BODY<VECTOR<T,d> >& rigid_body,const ARRAY_VIEW<int>& particles_array); \
    template void SOLIDS_STANDARD_TESTS<VECTOR<T,d> >::Bind_Unbound_Particles_In_Rigid_Body(RIGID_BODY<VECTOR<T,d> >& rigid_body,const ARRAY<int>& particles_array); \
    template void SOLIDS_STANDARD_TESTS<VECTOR<T,d> >::Add_Gravity();

#define INSTANTIATION_HELPER2(T) \
    INSTANTIATION_HELPER2_1D(T) \
    INSTANTIATION_HELPER2_ALL(T,3) \
    template void SOLIDS_STANDARD_TESTS<VECTOR<T,3> >::PD_Curl(const T,const VECTOR<T,3>,const ROTATION<VECTOR<T,3> >,const T,const int,const bool,const T);

INSTANTIATION_HELPER2(float);
template SOLIDS_STANDARD_TESTS<VECTOR<float,2> >::SOLIDS_STANDARD_TESTS(EXAMPLE<VECTOR<float,2> >&,SOLID_BODY_COLLECTION<VECTOR<float,2> >&);
template void SOLIDS_STANDARD_TESTS<VECTOR<float,2> >::Add_Gravity();
template void SOLIDS_STANDARD_TESTS<VECTOR<float,2> >::Bind_Unbound_Particles_In_Rigid_Body<ARRAY<int,int> >(RIGID_BODY<VECTOR<float,2> >&,ARRAY<int,int> const&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER2(double);
template SOLIDS_STANDARD_TESTS<VECTOR<double,2> >::SOLIDS_STANDARD_TESTS(EXAMPLE<VECTOR<double,2> >&,SOLID_BODY_COLLECTION<VECTOR<double,2> >&);
template void SOLIDS_STANDARD_TESTS<VECTOR<double,2> >::Add_Gravity();
template void SOLIDS_STANDARD_TESTS<VECTOR<double,2> >::Bind_Unbound_Particles_In_Rigid_Body<ARRAY<int,int> >(RIGID_BODY<VECTOR<double,2> >&,ARRAY<int,int> const&);
#endif

//#####################################################################
// Copyright 2004-2006, Zhaosheng Bao, Geoffrey Irving, Neil Molino, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLES_OF_MATERIAL
//##################################################################### 
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/TRIANGLES_OF_MATERIAL.h>
using namespace PhysBAM;
//#####################################################################
// Register this class as read-write
//#####################################################################
namespace {
bool Register_Triangles_Of_Material(){
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<TRIANGLES_OF_MATERIAL<VECTOR<float,2> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<TRIANGLES_OF_MATERIAL<VECTOR<float,3> > >();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<TRIANGLES_OF_MATERIAL<VECTOR<double,2> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<TRIANGLES_OF_MATERIAL<VECTOR<double,3> > >();
#endif
    return true;
}
static bool registered=Register_Triangles_Of_Material();
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> TRIANGLES_OF_MATERIAL<TV>::
TRIANGLES_OF_MATERIAL(T_EMBEDDED_OBJECT& embedded_object_input)
    :EMBEDDED_MATERIAL_SURFACE<TV,2>(embedded_object_input)
{
    PHYSBAM_ASSERT(registered);
}
//#####################################################################
// Function Perturb_Node_For_Collision_Freeness
//#####################################################################
template<class TV> void TRIANGLES_OF_MATERIAL<TV>::
Perturb_Nodes_For_Collision_Freeness(const T perturb_amount)
{
    // perturb those that have not been previously perturbed
    for(int node=1;node<=embedded_object.embedded_particles.active_indices.m;node++) if(!previously_perturbed(node)){
        int parent_1,parent_2;embedded_object.parent_particles(node).Get(parent_1,parent_2);
        ARRAY<int> triangles_on_edge;embedded_object.simplicial_object.mesh.Triangles_On_Edge(parent_1,parent_2,&triangles_on_edge);
        bool parent1_in_material=false,parent2_in_material=false,more_than_one_other_node_in_material=false;int other_node_in_material=0;
        for(int t=1;t<=triangles_on_edge.m;t++){
            int triangle=triangles_on_edge(t);
            if(embedded_object.Node_In_Simplex_Is_Material(parent_1,triangle)) parent1_in_material=true;
            if(embedded_object.Node_In_Simplex_Is_Material(parent_2,triangle)) parent2_in_material=true;
            int other_node=embedded_object.simplicial_object.mesh.Other_Node(parent_1,parent_2,triangle);
            if(embedded_object.Node_In_Simplex_Is_Material(other_node,triangle) && other_node == embedded_object.Diamond_Node(t)){
                if(!other_node_in_material) other_node_in_material=other_node;else more_than_one_other_node_in_material=true;}}
        if(parent1_in_material){
            if(!parent2_in_material){
                TV perturb_direction=(particles.X(parent_1)-embedded_particles.X(node)).Normalized();
                embedded_particles.X(node)+=perturb_amount*perturb_direction;previously_perturbed(node)=true;}}
        else if(parent2_in_material){
            TV perturb_direction=(particles.X(parent_2)-embedded_particles.X(node)).Normalized();
            embedded_particles.X(node)+=perturb_amount*perturb_direction;previously_perturbed(node)=true;}
        else if(other_node_in_material && !more_than_one_other_node_in_material){
            TV perturb_direction=(particles.X(other_node_in_material)-embedded_particles.X(node)).Normalized();           
            embedded_particles.X(node)+=perturb_amount*perturb_direction;previously_perturbed(node)=true;}}
}
//#####################################################################
// Function Construct_Material_Surface_Mesh
//#####################################################################
template<class TV> void TRIANGLES_OF_MATERIAL<TV>::
Construct_Material_Surface_Mesh() 
{
    bool embedded_segments_in_triangle_defined=embedded_object.embedded_subelements_in_parent_element!=0;
    if(!embedded_segments_in_triangle_defined) embedded_object.Initialize_Embedded_Subelements_In_Parent_Element(); 
    for(int t=1;t<=embedded_object.simplicial_object.mesh.elements.m;t++){
        bool i_is_material,j_is_material,k_is_material;embedded_object.node_in_simplex_is_material(t).Get(i_is_material,j_is_material,k_is_material);
        if(i_is_material && j_is_material && k_is_material){
            Add_To_Material_Surface_Mesh_Face_Triangle(t);continue;}
        VECTOR<int,2> emb_segs=embedded_object.Embedded_Subelements_In_Element(t);
        if(emb_segs[2]){
            Add_To_Material_Surface_Mesh_Subquadrilateral_Containg_Diamond_Node(t);
            Add_To_Material_Surface_Mesh_Corner_Triangles(t);} 
        else if(emb_segs[1]){
            Add_To_Material_Surface_Mesh_Isolated_Node_Subtriangle(t);
            Add_To_Material_Surface_Mesh_Subquadrilateral_Opposite_Isolated_Node(t);}}
    if(!embedded_segments_in_triangle_defined){
        delete embedded_object.embedded_subelements_in_parent_element_index;embedded_object.embedded_subelements_in_parent_element_index=0;
        delete embedded_object.number_of_embedded_subelements_in_parent_element;embedded_object.number_of_embedded_subelements_in_parent_element=0;
        delete embedded_object.embedded_subelements_in_parent_element;embedded_object.embedded_subelements_in_parent_element=0;}
}
//#####################################################################
// Function Add_To_Material_Surface_Mesh_Triangle
//#####################################################################
template<class TV> void TRIANGLES_OF_MATERIAL<TV>::
Add_To_Material_Surface_Mesh_Triangle(const int x1,const int x2,const int x3)
{
    material_surface_mesh.elements.Append(VECTOR<int,3>(x1,x2,x3));
}
//#####################################################################
// Function Add_To_Material_Surface_Mesh_Quad
//#####################################################################
template<class TV> void TRIANGLES_OF_MATERIAL<TV>::
Add_To_Material_Surface_Mesh_Quad(const int x1,const int x2,const int x3,const int x4)
{
    material_surface_mesh.elements.Append(VECTOR<int,3>(x1,x2,x3));
    material_surface_mesh.elements.Append(VECTOR<int,3>(x1,x3,x4));
}
//#####################################################################
// Function Add_To_Material_Surface_Mesh_Subquadrilateral_Containing_Diamond_Node             
//#####################################################################
template<class TV> void TRIANGLES_OF_MATERIAL<TV>::
Add_To_Material_Surface_Mesh_Face_Triangle(const int triangle)
{
    int i,j,k;embedded_object.simplicial_object.mesh.elements(triangle).Get(i,j,k);
    int ij=embedded_object.Particle_Embedded_On_Segment(i,j),ik=embedded_object.Particle_Embedded_On_Segment(i,k),jk=embedded_object.Particle_Embedded_On_Segment(j,k);
    if(ij && ik && jk){
        Add_To_Material_Surface_Mesh_Triangle(i,ij,ik);Add_To_Material_Surface_Mesh_Triangle(ij,j,jk);
        Add_To_Material_Surface_Mesh_Triangle(ik,jk,k);Add_To_Material_Surface_Mesh_Triangle(ik,ij,jk);}
    else if(!ij && ik && jk){Add_To_Material_Surface_Mesh_Triangle(k,ik,jk);Add_To_Material_Surface_Mesh_Quad(i,j,jk,ik);}
    else if(ij && !ik && jk){Add_To_Material_Surface_Mesh_Triangle(j,jk,ij);Add_To_Material_Surface_Mesh_Quad(i,ij,jk,k);}
    else if(ij && ik && !jk){Add_To_Material_Surface_Mesh_Triangle(i,ij,ik);Add_To_Material_Surface_Mesh_Quad(ij,j,k,ik);}
    else if(ij && !ik && !jk){Add_To_Material_Surface_Mesh_Triangle(j,k,ij);Add_To_Material_Surface_Mesh_Triangle(ij,k,i);}
    else if(!ij && ik && !jk){Add_To_Material_Surface_Mesh_Triangle(i,j,ik);Add_To_Material_Surface_Mesh_Triangle(ik,j,k);}
    else if(!ij && !ik && jk){Add_To_Material_Surface_Mesh_Triangle(i,j,jk);Add_To_Material_Surface_Mesh_Triangle(jk,k,i);}
    else Add_To_Material_Surface_Mesh_Triangle(i,j,k);
}
//#####################################################################
// Function Add_To_Material_Surface_Mesh_Subquadrilateral_Containing_Diamond_Node             
//#####################################################################
template<class TV> void TRIANGLES_OF_MATERIAL<TV>::
Add_To_Material_Surface_Mesh_Subquadrilateral_Containg_Diamond_Node(const int triangle)
{
    int i,j,k,diamond_node=embedded_object.Diamond_Node(triangle);
    if(!embedded_object.Node_In_Simplex_Is_Material(diamond_node,triangle))return;
    embedded_object.simplicial_object.mesh.elements(triangle).Get(i,j,k);
    assert(embedded_object.Number_Of_Embedded_Subelements_In_Element(triangle)==2);
    int ij=embedded_object.Particle_Embedded_On_Segment(i,j),ik=embedded_object.Particle_Embedded_On_Segment(i,k),jk=embedded_object.Particle_Embedded_On_Segment(j,k);
    if(diamond_node == i) Add_To_Material_Surface_Mesh_Diamond_Quad(i,ij,jk,ik);
    else if (diamond_node == j) Add_To_Material_Surface_Mesh_Diamond_Quad(j,jk,ik,ij);
    else if (diamond_node == k) Add_To_Material_Surface_Mesh_Diamond_Quad(k,ik,ij,jk);
}
//#####################################################################
// Function Add_To_Material_Surface_Mesh_Diamond_Quad
//#####################################################################
template<class TV> void TRIANGLES_OF_MATERIAL<TV>::
Add_To_Material_Surface_Mesh_Diamond_Quad(const int tri_node,const int emb_node1,const int emb_node2,const int emb_node3)
{
    Add_To_Material_Surface_Mesh_Triangle(tri_node,emb_node1,emb_node2); 
    Add_To_Material_Surface_Mesh_Triangle(tri_node,emb_node2,emb_node3); 
}
//#####################################################################
// Function Add_To_Material_Surface_Mesh_Corner_Triangles
//#####################################################################
template<class TV> void TRIANGLES_OF_MATERIAL<TV>::
Add_To_Material_Surface_Mesh_Corner_Triangles(const int triangle) 
{
    int i,j,k,diamond_node=embedded_object.Diamond_Node(triangle);embedded_object.simplicial_object.mesh.elements(triangle).Get(i,j,k);
    assert(embedded_object.Number_Of_Embedded_Subelements_In_Element(triangle)==2);
    int ij=embedded_object.Particle_Embedded_On_Segment(i,j),ik=embedded_object.Particle_Embedded_On_Segment(i,k),jk=embedded_object.Particle_Embedded_On_Segment(j,k);
    if(diamond_node == i){
        if(embedded_object.Node_In_Simplex_Is_Material(j,triangle))Add_To_Material_Surface_Mesh_Corner_Triangle(ij,j,jk);
        if(embedded_object.Node_In_Simplex_Is_Material(k,triangle))Add_To_Material_Surface_Mesh_Corner_Triangle(jk,k,ik);}
    else if (diamond_node == j){
        if(embedded_object.Node_In_Simplex_Is_Material(i,triangle))Add_To_Material_Surface_Mesh_Corner_Triangle(ik,i,ij);
        if(embedded_object.Node_In_Simplex_Is_Material(k,triangle))Add_To_Material_Surface_Mesh_Corner_Triangle(jk,k,ik);} 
    else if (diamond_node == k){
        if(embedded_object.Node_In_Simplex_Is_Material(i,triangle))Add_To_Material_Surface_Mesh_Corner_Triangle(ik,i,ij);
        if(embedded_object.Node_In_Simplex_Is_Material(j,triangle))Add_To_Material_Surface_Mesh_Corner_Triangle(ij,j,jk);}
}
//#####################################################################
// Function Add_To_Material_Surface_Mesh_Corner_Triangles
//#####################################################################
template<class TV> void TRIANGLES_OF_MATERIAL<TV>::
Add_To_Material_Surface_Mesh_Corner_Triangle(const int emb_node1,const int tri_node,const int emb_node2) 
{
    Add_To_Material_Surface_Mesh_Triangle(emb_node1,tri_node,emb_node2); 
}
//#####################################################################
// Function Add_To_Material_Surface_Mesh_Isolated_Node_Subtriangle
//#####################################################################
template<class TV> void TRIANGLES_OF_MATERIAL<TV>::
Add_To_Material_Surface_Mesh_Isolated_Node_Subtriangle(const int triangle)
{
    int isolated_node=embedded_object.Isolated_Node(triangle);
    if(!embedded_object.Node_In_Simplex_Is_Material(isolated_node,triangle)) return;
    int i,j,k;embedded_object.simplicial_object.mesh.elements(triangle).Get(i,j,k);
    VECTOR<int,2> emb_segments=embedded_object.Embedded_Subelements_In_Element(triangle);
    assert(emb_segments[1] && !emb_segments[2]);
    int curve_particle1,curve_particle2;embedded_object.embedded_mesh.elements(emb_segments[1]).Get(curve_particle1,curve_particle2);
    int embedded_curve_particle1=embedded_object.embedded_particles.subset_index_from_point_cloud_index(curve_particle1);
    if(isolated_node == i){
        if(embedded_object.Is_Parent(j,embedded_curve_particle1)) Add_To_Material_Surface_Mesh_Subtriangle(curve_particle1,curve_particle2,i);
        else Add_To_Material_Surface_Mesh_Subtriangle(curve_particle2,curve_particle1,i);}
    else if(isolated_node == j){
        if(embedded_object.Is_Parent(k,embedded_curve_particle1)) Add_To_Material_Surface_Mesh_Subtriangle(curve_particle1,curve_particle2,j);
        else Add_To_Material_Surface_Mesh_Subtriangle(curve_particle2,curve_particle1,j);}
    else{assert(isolated_node == k);
        if(embedded_object.Is_Parent(i,embedded_curve_particle1)) Add_To_Material_Surface_Mesh_Subtriangle(curve_particle1,curve_particle2,k);
        else Add_To_Material_Surface_Mesh_Subtriangle(curve_particle2,curve_particle1,k);}
}
//#####################################################################
// Function Add_To_Material_Surface_Mesh_Subtriangle
//#####################################################################
template<class TV> void TRIANGLES_OF_MATERIAL<TV>::
Add_To_Material_Surface_Mesh_Subtriangle(const int curve_particle1,const int curve_particle2,const int triangle_particle)
{
    Add_To_Material_Surface_Mesh_Triangle(triangle_particle,curve_particle1,curve_particle2); 
}
//#####################################################################
// Function Add_To_Material_Surface_Mesh_Subquadrilateral_Opposite_Isolated_Node
//#####################################################################
template<class TV> void TRIANGLES_OF_MATERIAL<TV>::
Add_To_Material_Surface_Mesh_Subquadrilateral_Opposite_Isolated_Node(const int triangle)
{
    int i,j,k,isolated_node=embedded_object.Isolated_Node(triangle);embedded_object.simplicial_object.mesh.elements(triangle).Get(i,j,k);
    VECTOR<int,2> emb_segments=embedded_object.Embedded_Subelements_In_Element(triangle);
    assert(emb_segments[1] && !emb_segments[2]);
    int curve_particle1,curve_particle2;embedded_object.embedded_mesh.elements(emb_segments[1]).Get(curve_particle1,curve_particle2);
    int embedded_curve_particle1=embedded_object.embedded_particles.subset_index_from_point_cloud_index(curve_particle1);
    if(isolated_node == i){
        if(!embedded_object.Node_In_Simplex_Is_Material(j,triangle)) return;
        if(embedded_object.Is_Parent(j,embedded_curve_particle1)) Add_To_Material_Surface_Mesh_Subquadrilateral(curve_particle1,curve_particle2,j,k);
        else Add_To_Material_Surface_Mesh_Subquadrilateral(curve_particle2,curve_particle1,j,k);}
    else if(isolated_node == j){
        if(!embedded_object.Node_In_Simplex_Is_Material(i,triangle)) return;
        if(embedded_object.Is_Parent(k,embedded_curve_particle1)) Add_To_Material_Surface_Mesh_Subquadrilateral(curve_particle1,curve_particle2,k,i);
        else Add_To_Material_Surface_Mesh_Subquadrilateral(curve_particle2,curve_particle1,k,i);}
    else{assert(isolated_node == k);
        if(!embedded_object.Node_In_Simplex_Is_Material(i,triangle)) return;
        if(embedded_object.Is_Parent(i,embedded_curve_particle1)) Add_To_Material_Surface_Mesh_Subquadrilateral(curve_particle1,curve_particle2,i,j);
        else Add_To_Material_Surface_Mesh_Subquadrilateral(curve_particle2,curve_particle1,i,j);}
}
//#####################################################################
// Function Add_To_Material_Surface_Mesh_Subquadrilateral
//#####################################################################
template<class TV> void TRIANGLES_OF_MATERIAL<TV>::
Add_To_Material_Surface_Mesh_Subquadrilateral(const int curve_particle1,const int curve_particle2,const int triangle_particle1,const int triangle_particle2)
{
    int jk=embedded_object.Particle_Embedded_On_Segment(triangle_particle1,triangle_particle2);
    if(jk){
        Add_To_Material_Surface_Mesh_Triangle(curve_particle1,jk,curve_particle2);
        Add_To_Material_Surface_Mesh_Triangle(curve_particle1,triangle_particle1,jk);
        Add_To_Material_Surface_Mesh_Triangle(curve_particle2,jk,triangle_particle2);} 
    else{
        Add_To_Material_Surface_Mesh_Triangle(curve_particle1,triangle_particle1,triangle_particle2); 
        Add_To_Material_Surface_Mesh_Triangle(curve_particle1,triangle_particle2,curve_particle2);}
}
//#####################################################################
template class TRIANGLES_OF_MATERIAL<VECTOR<float,2> >;
template class TRIANGLES_OF_MATERIAL<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TRIANGLES_OF_MATERIAL<VECTOR<double,2> >;
template class TRIANGLES_OF_MATERIAL<VECTOR<double,3> >;
#endif

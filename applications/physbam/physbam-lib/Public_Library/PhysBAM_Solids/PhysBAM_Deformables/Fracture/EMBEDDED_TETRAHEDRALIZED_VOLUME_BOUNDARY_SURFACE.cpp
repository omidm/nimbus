//#####################################################################
// Copyright 2003-2006, Bao, Ron Fedkiw, Geoffrey Irving, Josh, Neil Molino, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license  contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/exchange_sort.h>
#include <PhysBAM_Tools/Math_Tools/permutation.h>
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Register this class as read-write
//#####################################################################
namespace {
bool Register_Embedded_Tetrahedralized_Volume_Boundary_Surface(){
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<float> >();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<double> >();
#endif
    return true;
}
static bool registered=Register_Embedded_Tetrahedralized_Volume_Boundary_Surface();
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& embedded_object)
    :EMBEDDED_MATERIAL_SURFACE<TV,3>(embedded_object)
{
    PHYSBAM_ASSERT(registered);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
~EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE()
{}
//#####################################################################
// Function Create_Material_Surface_From_Manifold_Embedded_Surface
//#####################################################################
template<class T> void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Create_Material_Surface_From_Manifold_Embedded_Surface(const bool verbose)
{
    material_surface_mesh.Clean_Memory();material_surface.Clean_Memory();
    TETRAHEDRON_MESH& mesh=embedded_object.simplicial_object.mesh;
    embedded_object.Update_Embedded_Particle_Positions(); // TODO: remove these?
    material_surface_mesh.number_nodes=particles.array_collection->Size();
    material_surface_mesh.elements.Append_Elements(embedded_object.embedded_object.mesh.elements); // add embedded elements to boundary mesh
    ARRAY<bool> node_is_material(mesh.number_nodes);
    for(int t=1;t<=mesh.elements.m;t++) node_is_material.Subset(mesh.elements(t))=embedded_object.node_in_simplex_is_material(t);
    bool boundary_mesh_defined=mesh.boundary_mesh!=0;if(!boundary_mesh_defined) mesh.Initialize_Boundary_Mesh();
    // TODO: the following loop only does something when part of the simulation boundary is real boundary, but this function doesn't handle that case correctly
    for(int t=1;t<=mesh.boundary_mesh->elements.m;t++){ // add pure material tetrahedron faces to boundary mesh
        VECTOR<int,3>& triangle=mesh.boundary_mesh->elements(t);
        if(node_is_material(triangle.x) && node_is_material(triangle.y) && node_is_material(triangle.z)) material_surface_mesh.elements.Append(triangle);}
    if(!boundary_mesh_defined){delete mesh.boundary_mesh;mesh.boundary_mesh=0;}
    if(verbose) {std::stringstream ss;ss<<"Triangles On Material Surface Mesh: "<<material_surface_mesh.elements.m<<std::endl;LOG::filecout(ss.str());}
}
//#####################################################################
// Function Conservative_Perturb_Nodes_For_Collision_Freeness
//#####################################################################
template<class T> void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Conservative_Perturb_Nodes_For_Collision_Freeness(const T perturb_amount,const ARRAY<bool>& particle_on_surface) // TODO: not currently called, so should be dealt with later
{
    embedded_object.Initialize_Orientation_Index_If_Necessary();

    // perturb those that have not been previously perturbed
    for(int node=1;node<=embedded_object.embedded_particles.active_indices.m;node++) if(!previously_perturbed(node) && particle_on_surface(node)){
        const VECTOR<int,2>& parents=embedded_object.parent_particles(node);
        ARRAY<int> tetrahedrons_on_edge;embedded_object.simplicial_object.mesh.Tetrahedrons_On_Edge(parents,tetrahedrons_on_edge);
        bool parent1_ever_material_by_itself=false,parent2_ever_material_by_itself=false,more_than_one_other_node_on_boundary=false;int other_node_on_boundary1=0,other_node_on_boundary2=0;
        for(int t=1;t<=tetrahedrons_on_edge.m;t++){
            int tetrahedron=tetrahedrons_on_edge(t);
            bool parent1_material_in_current_tetrahedron=embedded_object.Node_In_Simplex_Is_Material(parents[1],tetrahedron);
            bool parent2_material_in_current_tetrahedron=embedded_object.Node_In_Simplex_Is_Material(parents[2],tetrahedron);
            if(!parent1_material_in_current_tetrahedron || !parent2_material_in_current_tetrahedron){
                parent1_ever_material_by_itself=parent1_ever_material_by_itself || parent1_material_in_current_tetrahedron;
                parent2_ever_material_by_itself=parent2_ever_material_by_itself || parent2_material_in_current_tetrahedron;
                int other_node1,other_node2;embedded_object.simplicial_object.mesh.Other_Two_Nodes(parents[1],parents[2],tetrahedron,other_node1,other_node2);
                if((embedded_object.Node_In_Simplex_Is_Material(other_node1,tetrahedron) || embedded_object.Node_In_Simplex_Is_Material(other_node2,tetrahedron))
                        && Center_Octahedron_In_Material(t)){
                    if(!other_node_on_boundary1 && !other_node_on_boundary2){other_node_on_boundary1=other_node1;other_node_on_boundary2=other_node2;}
                    else more_than_one_other_node_on_boundary=true;}}}
        if(parent1_ever_material_by_itself){
            if(!parent2_ever_material_by_itself){
                VECTOR<T,3> perturb_direction=(particles.X(parents[1])-embedded_particles.X(node)).Normalized();
                embedded_particles.X(node)+=perturb_amount*perturb_direction;previously_perturbed(node)=true;}}
        else if(parent2_ever_material_by_itself){
            VECTOR<T,3> perturb_direction=(particles.X(parents[2])-embedded_particles.X(node)).Normalized();
            embedded_particles.X(node)+=perturb_amount*perturb_direction;previously_perturbed(node)=true;}
        else if(!more_than_one_other_node_on_boundary){
            if(other_node_on_boundary1 && other_node_on_boundary2){
                VECTOR<T,3> opposite_midpoint((T).5*(particles.X(other_node_on_boundary1)+particles.X(other_node_on_boundary2)));
                VECTOR<T,3> perturb_direction=(opposite_midpoint-embedded_particles.X(node)).Normalized();
                embedded_particles.X(node)+=perturb_amount*perturb_direction;previously_perturbed(node)=true;}
            else if(other_node_on_boundary1 && !other_node_on_boundary2){
                 VECTOR<T,3> perturb_direction=(particles.X(other_node_on_boundary1)-embedded_particles.X(node)).Normalized();
                embedded_particles.X(node)+=perturb_amount*perturb_direction;previously_perturbed(node)=true;}
            else if(!other_node_on_boundary1 && other_node_on_boundary2){
                VECTOR<T,3> perturb_direction=(particles.X(other_node_on_boundary2)-embedded_particles.X(node)).Normalized();
                embedded_particles.X(node)+=perturb_amount*perturb_direction;previously_perturbed(node)=true;}}}
}
//#####################################################################
// Function Perturb_Nodes_For_Collision_Freeness
//#####################################################################
template<class T> void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Perturb_Nodes_For_Collision_Freeness(const T perturb_amount) // used to be called Aggressive_
{
    embedded_object.Initialize_Orientation_Index_If_Necessary();

    for(int node=1;node<=embedded_particles.active_indices.m;node++){
        const VECTOR<int,2>& parents=embedded_object.parent_particles(node);
        ARRAY<int> tetrahedrons_on_edge;embedded_object.simplicial_object.mesh.Tetrahedrons_On_Edge(parents,tetrahedrons_on_edge);
        for(int t=1;t<=tetrahedrons_on_edge.m;t++){
            int tetrahedron=tetrahedrons_on_edge(t);
            bool parent1_material_in_current_tetrahedron=embedded_object.Node_In_Simplex_Is_Material(parents[1],tetrahedron);
            bool parent2_material_in_current_tetrahedron=embedded_object.Node_In_Simplex_Is_Material(parents[2],tetrahedron);
            int other_node1,other_node2;embedded_object.simplicial_object.mesh.Other_Two_Nodes(parents[1],parents[2],tetrahedron,other_node1,other_node2);
            bool other_node1_material_in_current_tetrahedron=embedded_object.Node_In_Simplex_Is_Material(other_node1,tetrahedron);
            bool other_node2_material_in_current_tetrahedron=embedded_object.Node_In_Simplex_Is_Material(other_node2,tetrahedron);
            if(parent1_material_in_current_tetrahedron){
                    if(!parent2_material_in_current_tetrahedron){
                        VECTOR<T,3> perturb_direction=(particles.X(parents[1])-embedded_particles.X(node)).Normalized();
                        embedded_particles.X(node)+=perturb_amount*perturb_direction;previously_perturbed(node)=true;}}
            else if(parent2_material_in_current_tetrahedron){
                    VECTOR<T,3> perturb_direction=(particles.X(parents[2])-embedded_particles.X(node)).Normalized();
                    embedded_particles.X(node)+=perturb_amount*perturb_direction;previously_perturbed(node)=true;}
            else if(other_node1_material_in_current_tetrahedron && other_node2_material_in_current_tetrahedron){
                        VECTOR<T,3> opposite_midpoint((T).5*(particles.X(other_node1)+particles.X(other_node2)));
                        VECTOR<T,3> perturb_direction=(opposite_midpoint-embedded_particles.X(node)).Normalized();
                        embedded_particles.X(node)+=perturb_amount*perturb_direction;previously_perturbed(node)=true;}
            else if(other_node1_material_in_current_tetrahedron && !other_node2_material_in_current_tetrahedron){
                VECTOR<T,3> perturb_direction=(particles.X(other_node1)-embedded_particles.X(node)).Normalized();
                embedded_particles.X(node)+=perturb_amount*perturb_direction;previously_perturbed(node)=true;}
            else if(!other_node1_material_in_current_tetrahedron && other_node2_material_in_current_tetrahedron){
                VECTOR<T,3> perturb_direction=(particles.X(other_node2)-embedded_particles.X(node)).Normalized();
                embedded_particles.X(node)+=perturb_amount*perturb_direction;previously_perturbed(node)=true;}}}
}
//#####################################################################
// Function Center_Octahedron_In_Material
//#####################################################################
template<class T> bool EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Center_Octahedron_In_Material(const int tetrahedron)
{
    VECTOR<int,4> emb_tris=embedded_object.Embedded_Subelements_In_Element(tetrahedron);
    VECTOR<bool,4> is_material=permute_four(embedded_object.node_in_simplex_is_material(tetrahedron),embedded_object.orientation_index(tetrahedron));
    int i,j,k,l;permute_four(embedded_object.simplicial_object.mesh.elements(tetrahedron),embedded_object.orientation_index(tetrahedron)).Get(i,j,k,l);
    int ij=embedded_object.Particle_Embedded_On_Segment(i,j),ik=embedded_object.Particle_Embedded_On_Segment(i,k),il=embedded_object.Particle_Embedded_On_Segment(i,l),
        jk=embedded_object.Particle_Embedded_On_Segment(j,k),jl=embedded_object.Particle_Embedded_On_Segment(j,l),kl=embedded_object.Particle_Embedded_On_Segment(k,l);
    if(emb_tris[1] && emb_tris[2] && !emb_tris[3] && embedded_object.embedded_object.mesh.Node_In_Triangle(ij,emb_tris[1])) return true;
    else if(emb_tris[1] && emb_tris[2] && emb_tris[3] && !emb_tris[4] && embedded_object.Node_Separated_By_Embedded_Subelement(emb_tris[1])
        && embedded_object.Node_Separated_By_Embedded_Subelement(emb_tris[2]) && embedded_object.Node_Separated_By_Embedded_Subelement(emb_tris[3])
        && ((is_material[3] && is_material[4])
            || (is_material[4] && embedded_object.embedded_object.mesh.Triangle(kl,ik,jk))
            || (is_material[3] && embedded_object.embedded_object.mesh.Triangle(kl,jl,il)))) return true;
    return false;
}
//#####################################################################
// Function Construct_Material_Surface_Mesh
//#####################################################################
template<class T> void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Construct_Material_Surface_Mesh()
{
    embedded_object.Initialize_Orientation_Index_If_Necessary();
    if(!embedded_object.simplicial_object.mesh.incident_elements) embedded_object.simplicial_object.mesh.Initialize_Incident_Elements();
    if(!embedded_object.embedded_object.mesh.incident_elements) embedded_object.embedded_object.mesh.Initialize_Incident_Elements();
    if(!embedded_object.embedded_subelements_in_parent_element) embedded_object.Initialize_Embedded_Subelements_In_Parent_Element();
    for(int t=1;t<=embedded_object.simplicial_object.mesh.elements.m;t++){
        bool i_is_material,j_is_material,k_is_material,l_is_material;
        embedded_object.node_in_simplex_is_material(t).Get(i_is_material,j_is_material,k_is_material,l_is_material);
        if(i_is_material && j_is_material && k_is_material && l_is_material){
            int i,j,k,l;embedded_object.simplicial_object.mesh.elements(t).Get(i,j,k,l);
            Add_To_Material_Surface_Tetrahedron(i,j,k,l);continue;}
        VECTOR<int,4> emb_tris=embedded_object.Embedded_Subelements_In_Element(t);
        int i,j,k,l;permute_four(embedded_object.simplicial_object.mesh.elements(t),embedded_object.orientation_index(t)).Get(i,j,k,l);
        if(emb_tris[1] && !emb_tris[2]) Add_To_Material_Surface_Subtetrahedron_And_Subprism(t,emb_tris[1],i,j,k,l);
        else if(emb_tris[1] && emb_tris[2] && !emb_tris[3]){
            if(embedded_object.embedded_object.mesh.Node_In_Triangle(embedded_object.Particle_Embedded_On_Segment(i,j),emb_tris[1]))
                Add_To_Material_Surface_Subtetrahedrons_And_Wrick(t,emb_tris[1],emb_tris[2],i,j,k,l);
            else Add_To_Material_Surface_Wedge_On_Either_Side(t,emb_tris[1],emb_tris[2],i,j,k,l);}
        else if(emb_tris[1] && emb_tris[2] && emb_tris[3] && !emb_tris[4]){
            int isolated_node1=embedded_object.Node_Separated_By_Embedded_Subelement(emb_tris[1]);
            int isolated_node2=embedded_object.Node_Separated_By_Embedded_Subelement(emb_tris[2]);
            int isolated_node3=embedded_object.Node_Separated_By_Embedded_Subelement(emb_tris[3]);
            if(isolated_node1 && isolated_node2 && isolated_node3)Add_To_Material_Surface_Subtetrahedron_Tets_And_Oct_Plus_Tet(t,i,j,k,l);
            else Add_To_Material_Surface_Subtetrahedron_And_Wedge_And_Half_Oct_Plus_Tet(t,i,j,k,l);}
        else if(emb_tris[1] && emb_tris[2] && emb_tris[3] && emb_tris[4])Add_To_Material_Surface_Subtetrahedron_And_Wedge_And_Half_Oct_Plus_Tet(t,i,j,k,l);}
    Merge_Or_Cancel_Duplicate_Triangles();
}  
//#####################################################################
// Function Merge_Or_Cancel_Duplicate_Triangles
//#####################################################################
template<class T> void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Merge_Or_Cancel_Duplicate_Triangles()
{
    HASHTABLE<VECTOR<int,3>,int> surface_hash_table(60);
    for(int t=1;t<=material_surface_mesh.elements.m;t++){
        VECTOR<int,3>& triangle=material_surface_mesh.elements(t);
        VECTOR<int,3> sorted_triangle=triangle.Sorted();
        int t2;
        if(!surface_hash_table.Get(sorted_triangle,t2))
            surface_hash_table.Insert(sorted_triangle,t);
        else{
            VECTOR<int,3>& triangle2=material_surface_mesh.elements(t2);
            if(triangle2.x){
                if(!TRIANGLE_MESH::Equivalent_Oriented_Triangles(triangle.x,triangle.y,triangle.z,triangle2.x,triangle2.y,triangle2.z)) triangle2.x=0; // the two triangles cancel out
                else PHYSBAM_FATAL_ERROR();} // TODO: when we transformed this code, we didn't know whether this assert would trigger (if it does, feel free to remove it)
            triangle.x=0;}} // never need two copies of a triangle, so mark this one for deletion
    material_surface_mesh.Delete_Elements_With_Missing_Nodes();
}
//#####################################################################
// Function Add_To_Candidate_Material_Surface_Tetrahedron
//#####################################################################
template<class T> void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Add_To_Material_Surface_Tetrahedron(const int i,const int j,const int k,const int l) 
{
    Add_To_Material_Surface_Tetrahedron_Face(i,k,j,true);Add_To_Material_Surface_Tetrahedron_Face(l,i,j,true);
    Add_To_Material_Surface_Tetrahedron_Face(l,j,k,true);Add_To_Material_Surface_Tetrahedron_Face(l,k,i,true);
}
//#####################################################################
// Function Add_To_Material_Surface_Tetrahedron_Face
//#####################################################################
template<class T> void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Add_To_Material_Surface_Tetrahedron_Face(const int i,const int j,const int l,const bool is_clockwise)
{
    int ij=embedded_object.Particle_Embedded_On_Segment(i,j),il=embedded_object.Particle_Embedded_On_Segment(i,l),jl=embedded_object.Particle_Embedded_On_Segment(j,l);
    if(ij && il && jl){
        Add_To_Material_Surface_Triangle(i,ij,il,is_clockwise);Add_To_Material_Surface_Triangle(ij,j,jl,is_clockwise);
        Add_To_Material_Surface_Triangle(il,jl,l,is_clockwise);Add_To_Material_Surface_Triangle(il,ij,jl,is_clockwise);}   
    else if(!ij && il && jl){Add_To_Material_Surface_Triangle(l,il,jl,is_clockwise);Add_To_Material_Surface_Planar_Quad(i,j,jl,il,is_clockwise);}
    else if(ij && !il && jl){Add_To_Material_Surface_Triangle(j,jl,ij,is_clockwise);Add_To_Material_Surface_Planar_Quad(i,ij,jl,l,is_clockwise);}
    else if(ij && il && !jl){Add_To_Material_Surface_Triangle(i,ij,il,is_clockwise);Add_To_Material_Surface_Planar_Quad(ij,j,l,il,is_clockwise);}
    else if(ij && !il && !jl){Add_To_Material_Surface_Triangle(j,l,ij,is_clockwise);Add_To_Material_Surface_Triangle(ij,l,i,is_clockwise);}
    else if(!ij && il && !jl){Add_To_Material_Surface_Triangle(i,j,il,is_clockwise);Add_To_Material_Surface_Triangle(il,j,l,is_clockwise);}
    else if(!ij && !il && jl){Add_To_Material_Surface_Triangle(i,j,jl,is_clockwise);Add_To_Material_Surface_Triangle(jl,l,i,is_clockwise);}
    else Add_To_Material_Surface_Triangle(i,j,l,is_clockwise);
}
//#####################################################################
// Function Add_To_Material_Surface_Triangle
//#####################################################################
template<class T> void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Add_To_Material_Surface_Triangle(int x1,int x2,int x3,const bool is_clockwise)
{
    if(!is_clockwise) exchange(x2,x3);
    material_surface_mesh.elements.Append(VECTOR<int,3>(x1,x2,x3));
}
//#####################################################################
// Function Add_To_Material_Surface_Quad
//#####################################################################
template<class T> void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Add_To_Material_Surface_Quad_Cut(const int il,const int jl,const int jk,const int ik,const bool is_clockwise)
{
    Add_To_Material_Surface_Triangle(il,jl,ik,is_clockwise);Add_To_Material_Surface_Triangle(jl,jk,ik,is_clockwise);
}
//#####################################################################
// Function Add_To_Material_Surface_Quad
//#####################################################################
template<class T> void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Add_To_Material_Surface_Planar_Quad(const int x1,const int x2,const int x3,const int x4,const bool is_clockwise) 
{           
    if(min(x1,x3)<min(x2,x4)){Add_To_Material_Surface_Triangle(x1,x2,x3,is_clockwise);Add_To_Material_Surface_Triangle(x3,x4,x1,is_clockwise);}
    else{Add_To_Material_Surface_Triangle(x1,x2,x4,is_clockwise);Add_To_Material_Surface_Triangle(x2,x3,x4,is_clockwise);}
}
//#####################################################################
// Function Add_To_Material_Surface_Subtetrahedron_And_Subprism
//#####################################################################
template<class T> void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Add_To_Material_Surface_Subtetrahedron_And_Subprism(const int tetrahedron,const int embedded_triangle1,const int i,const int j,const int k,const int l)
{    
    assert((embedded_object.embedded_object.mesh.elements(embedded_triangle1)==
        VECTOR<int,3>(embedded_object.Particle_Embedded_On_Segment(i,j),embedded_object.Particle_Embedded_On_Segment(i,k),embedded_object.Particle_Embedded_On_Segment(i,l))));
    bool is_clockwise=permutation_of_four_is_even(embedded_object.orientation_index(tetrahedron));
    if(embedded_object.Node_In_Simplex_Is_Material(i,tetrahedron)) Add_To_Material_Surface_Subtetrahedron(i,j,k,l,is_clockwise);
    if(embedded_object.Node_In_Simplex_Is_Material(j,tetrahedron) && embedded_object.Node_In_Simplex_Is_Material(k,tetrahedron) && 
        embedded_object.Node_In_Simplex_Is_Material(l,tetrahedron)) Add_To_Material_Surface_Subprism(i,j,k,l,is_clockwise);
}
//#####################################################################
// Function Add_To_Material_Surface_Subtetrahedron
//#####################################################################
template<class T>
void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>:: 
Add_To_Material_Surface_Subtetrahedron(const int i,const int j,const int k,const int l,const bool is_clockwise)
{ 
    int ij=embedded_object.Particle_Embedded_On_Segment(i,j),ik=embedded_object.Particle_Embedded_On_Segment(i,k),il=embedded_object.Particle_Embedded_On_Segment(i,l);
    Add_To_Material_Surface_Triangle(ij,ik,il,is_clockwise);Add_To_Material_Surface_Triangle(i,ik,ij,is_clockwise);
    Add_To_Material_Surface_Triangle(i,il,ik,is_clockwise);Add_To_Material_Surface_Triangle(i,ij,il,is_clockwise);
}
//#####################################################################
// Function Add_To_Material_Surface_Subprism
//#####################################################################
template<class T> void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Add_To_Material_Surface_Subprism(const int i,const int j,const int k,const int l,const bool is_clockwise) 
{
    int ij=embedded_object.Particle_Embedded_On_Segment(i,j),ik=embedded_object.Particle_Embedded_On_Segment(i,k),il=embedded_object.Particle_Embedded_On_Segment(i,l),
        jk=embedded_object.Particle_Embedded_On_Segment(j,k),kl=embedded_object.Particle_Embedded_On_Segment(k,l),jl=embedded_object.Particle_Embedded_On_Segment(j,l);
    Add_To_Material_Surface_Triangle(il,ik,ij,is_clockwise);Add_To_Material_Surface_Tetrahedron_Face(j,k,l,is_clockwise);    
    if(!jk)Add_To_Material_Surface_Planar_Quad(k,j,ij,ik,is_clockwise);
    else{Add_To_Material_Surface_Triangle(j,ij,jk,is_clockwise);Add_To_Material_Surface_Triangle(ij,ik,jk,is_clockwise);Add_To_Material_Surface_Triangle(ik,k,jk,is_clockwise);}
    if(!kl)Add_To_Material_Surface_Planar_Quad(l,k,ik,il,is_clockwise);
    else{Add_To_Material_Surface_Triangle(k,ik,kl,is_clockwise);Add_To_Material_Surface_Triangle(ik,il,kl,is_clockwise);Add_To_Material_Surface_Triangle(il,l,kl,is_clockwise);}
    if(!jl)Add_To_Material_Surface_Planar_Quad(j,l,il,ij,is_clockwise);
    else{Add_To_Material_Surface_Triangle(l,il,jl,is_clockwise);Add_To_Material_Surface_Triangle(il,ij,jl,is_clockwise);Add_To_Material_Surface_Triangle(ij,j,jl,is_clockwise);}
}
//#####################################################################
// Function Add_To_Material_Surface_Subtetrahedrons_And_Wrick
//#####################################################################
template<class T> void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Add_To_Material_Surface_Subtetrahedrons_And_Wrick(const int tetrahedron,const int embedded_triangle1,const int embedded_triangle2,const int i,const int j,const int k,const int l) 
{
    bool is_clockwise=permutation_of_four_is_even(embedded_object.orientation_index(tetrahedron));
    if(embedded_object.Node_In_Simplex_Is_Material(i,tetrahedron))Add_To_Material_Surface_Subtetrahedron(i,j,k,l,is_clockwise);
    if(embedded_object.Node_In_Simplex_Is_Material(j,tetrahedron))Add_To_Material_Surface_Subtetrahedron(j,i,l,k,is_clockwise);
    if(embedded_object.Node_In_Simplex_Is_Material(k,tetrahedron) && embedded_object.Node_In_Simplex_Is_Material(l,tetrahedron))
       Add_To_Material_Surface_Wrick(i,j,k,l,is_clockwise);
}
//#####################################################################
// Function Add_To_Material_Surface_Wrick
//#####################################################################
template<class T> void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Add_To_Material_Surface_Wrick(const int i,const int j,const int k,const int l,const bool is_clockwise) 
{
    int ij=embedded_object.Particle_Embedded_On_Segment(i,j),ik=embedded_object.Particle_Embedded_On_Segment(i,k),jk=embedded_object.Particle_Embedded_On_Segment(j,k),
        jl=embedded_object.Particle_Embedded_On_Segment(j,l),il=embedded_object.Particle_Embedded_On_Segment(i,l),kl=embedded_object.Particle_Embedded_On_Segment(k,l);
    Add_To_Material_Surface_Triangle(l,il,jl,is_clockwise);Add_To_Material_Surface_Triangle(jl,il,ij,is_clockwise);
    Add_To_Material_Surface_Triangle(ik,k,jk,is_clockwise);Add_To_Material_Surface_Triangle(ik,jk,ij,is_clockwise);
    if(!kl) Add_To_Material_Surface_Planar_Quad(il,l,k,ik,is_clockwise);
    else{Add_To_Material_Surface_Triangle(l,kl,il,is_clockwise);Add_To_Material_Surface_Triangle(kl,ik,il,is_clockwise);Add_To_Material_Surface_Triangle(kl,k,ik,is_clockwise);}
    if(!kl) Add_To_Material_Surface_Planar_Quad(k,l,jl,jk,is_clockwise);
    else{Add_To_Material_Surface_Triangle(l,jl,kl,is_clockwise);Add_To_Material_Surface_Triangle(jl,jk,kl,is_clockwise);Add_To_Material_Surface_Triangle(jk,k,kl,is_clockwise);}
    Add_To_Material_Surface_Triangle(ij,jk,jl,is_clockwise);Add_To_Material_Surface_Triangle(ij,il,ik,is_clockwise);
}
//#####################################################################
// Function Add_To_Material_Surface_Wedge_On_Either_Side
//#####################################################################
template<class T> void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Add_To_Material_Surface_Wedge_On_Either_Side(const int tetrahedron,const int embedded_triangle1,const int embedded_triangle2,const int i,const int j,const int k,const int l) 
{
    VECTOR<int,4> surface_particles;
    embedded_object.embedded_object.mesh.elements(embedded_triangle1).Get(surface_particles[1],surface_particles[2],surface_particles[3]);
    embedded_object.embedded_object.mesh.elements(embedded_triangle2).Get(surface_particles[3],surface_particles[4],surface_particles[1]);
    assert((surface_particles==VECTOR<int,4>(embedded_object.Particle_Embedded_On_Segment(i,k),embedded_object.Particle_Embedded_On_Segment(k,j),
        embedded_object.Particle_Embedded_On_Segment(j,l),embedded_object.Particle_Embedded_On_Segment(l,i))));
    bool is_clockwise=permutation_of_four_is_even(embedded_object.orientation_index(tetrahedron));
    if(embedded_object.Node_In_Simplex_Is_Material(i,tetrahedron) && embedded_object.Node_In_Simplex_Is_Material(j,tetrahedron))
        Add_To_Material_Surface_Subwedge(i,j,k,l,is_clockwise);
    if(embedded_object.Node_In_Simplex_Is_Material(k,tetrahedron) && embedded_object.Node_In_Simplex_Is_Material(l,tetrahedron))        
        Add_To_Material_Surface_Subwedge(k,l,i,j,is_clockwise);
}
//#####################################################################
// Function Add_To_Material_Surface_Subwedge
//#####################################################################
template<class T> void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Add_To_Material_Surface_Subwedge(const int i,const int j,const int k,const int l,const bool is_clockwise)
{
    int ij=embedded_object.Particle_Embedded_On_Segment(i,j),ik=embedded_object.Particle_Embedded_On_Segment(i,k),il=embedded_object.Particle_Embedded_On_Segment(i,l),
        jk=embedded_object.Particle_Embedded_On_Segment(j,k),jl=embedded_object.Particle_Embedded_On_Segment(j,l);
    Add_To_Material_Surface_Quad_Cut(il,jl,jk,ik,is_clockwise);
    if(!ij) Add_To_Material_Surface_Planar_Quad(ik,jk,j,i,is_clockwise);
    else{Add_To_Material_Surface_Triangle(jk,j,ij,is_clockwise);Add_To_Material_Surface_Triangle(ij,ik,jk,is_clockwise);Add_To_Material_Surface_Triangle(ij,i,ik,is_clockwise);}    
    if(!ij) Add_To_Material_Surface_Planar_Quad(jl,il,i,j,is_clockwise);
    else{Add_To_Material_Surface_Triangle(j,jl,ij,is_clockwise);Add_To_Material_Surface_Triangle(jl,il,ij,is_clockwise);Add_To_Material_Surface_Triangle(il,i,ij,is_clockwise);}    
    Add_To_Material_Surface_Triangle(i,il,ik,is_clockwise);
    Add_To_Material_Surface_Triangle(jk,jl,j,is_clockwise);
}
//#####################################################################
// Function Add_To_Material_Surface_Subtetrahedron_Tets_And_Oct_Plus_Tet
//#####################################################################
template<class T> void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Add_To_Material_Surface_Subtetrahedron_Tets_And_Oct_Plus_Tet(const int tetrahedron,const int i,const int j,const int k,const int l) 
{
    int ik=embedded_object.Particle_Embedded_On_Segment(i,k),il=embedded_object.Particle_Embedded_On_Segment(i,l),jk=embedded_object.Particle_Embedded_On_Segment(j,k),
        jl=embedded_object.Particle_Embedded_On_Segment(j,l),kl=embedded_object.Particle_Embedded_On_Segment(k,l);
    bool is_clockwise=permutation_of_four_is_even(embedded_object.orientation_index(tetrahedron));
    if(embedded_object.Node_In_Simplex_Is_Material(i,tetrahedron)) Add_To_Material_Surface_Subtetrahedron(i,j,k,l,is_clockwise);
    if(embedded_object.Node_In_Simplex_Is_Material(j,tetrahedron)) Add_To_Material_Surface_Subtetrahedron(j,i,l,k,is_clockwise);
    if(embedded_object.Node_In_Simplex_Is_Material(k,tetrahedron) && !embedded_object.Node_In_Simplex_Is_Material(l,tetrahedron)){        
        if(embedded_object.embedded_object.mesh.Triangle(kl,ik,jk))
            Add_To_Material_Surface_Subtetrahedron(k,l,i,j,is_clockwise);
        else{
            assert(embedded_object.embedded_object.mesh.Triangle(kl,jl,il));
            Add_To_Material_Surface_Oct_Plus_Tet(l,k,j,i,is_clockwise);}}
    if(!embedded_object.Node_In_Simplex_Is_Material(k,tetrahedron) && embedded_object.Node_In_Simplex_Is_Material(l,tetrahedron)){
        if(embedded_object.embedded_object.mesh.Triangle(kl,jl,il))
            Add_To_Material_Surface_Subtetrahedron(l,k,j,i,is_clockwise);
        else{
            assert(embedded_object.embedded_object.mesh.Triangle(kl,ik,jk));
            Add_To_Material_Surface_Oct_Plus_Tet(k,l,i,j,is_clockwise);}}
    if(embedded_object.Node_In_Simplex_Is_Material(k,tetrahedron) && embedded_object.Node_In_Simplex_Is_Material(l,tetrahedron))
        Add_To_Material_Surface_Wrick(i,j,k,l,is_clockwise);
}
//#####################################################################
// Function Add_To_Material_Surface_Oct_Plus_Tet
//#####################################################################
template<class T> void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Add_To_Material_Surface_Oct_Plus_Tet(const int i,const int j,const int k,const int l,const bool is_clockwise)
{
    int ij=embedded_object.Particle_Embedded_On_Segment(i,j),ik=embedded_object.Particle_Embedded_On_Segment(i,k),il=embedded_object.Particle_Embedded_On_Segment(i,l),
        jk=embedded_object.Particle_Embedded_On_Segment(j,k),jl=embedded_object.Particle_Embedded_On_Segment(j,l),kl=embedded_object.Particle_Embedded_On_Segment(k,l);
    Add_To_Material_Surface_Triangle(il,jl,kl,is_clockwise);Add_To_Material_Surface_Triangle(il,kl,ik,is_clockwise);Add_To_Material_Surface_Triangle(kl,jk,ik,is_clockwise);
    Add_To_Material_Surface_Triangle(ij,il,ik,is_clockwise);Add_To_Material_Surface_Triangle(j,ij,jk,is_clockwise);Add_To_Material_Surface_Triangle(ij,ik,jk,is_clockwise);
    Add_To_Material_Surface_Triangle(jl,il,ij,is_clockwise);Add_To_Material_Surface_Triangle(ij,j,jl,is_clockwise);
    Add_To_Material_Surface_Triangle(j,jk,jl,is_clockwise);Add_To_Material_Surface_Triangle(jk,kl,jl,is_clockwise);
}
//#####################################################################
// Function Add_To_Material_Surface_Half_Oct_Plus_Tet
//#####################################################################
template<class T> void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Add_To_Material_Surface_Half_Oct_Plus_Tet(const int i,const int j,const int k,const int l,const bool is_clockwise)
{
    int ij=embedded_object.Particle_Embedded_On_Segment(i,j),ik=embedded_object.Particle_Embedded_On_Segment(i,k),il=embedded_object.Particle_Embedded_On_Segment(i,l),
        jk=embedded_object.Particle_Embedded_On_Segment(j,k),jl=embedded_object.Particle_Embedded_On_Segment(j,l);
    Add_To_Material_Surface_Quad_Cut(il,jl,jk,ik,is_clockwise);Add_To_Material_Surface_Triangle(j,ij,jk,is_clockwise);Add_To_Material_Surface_Triangle(ij,ik,jk,is_clockwise);
    Add_To_Material_Surface_Triangle(jl,il,ij,is_clockwise);Add_To_Material_Surface_Triangle(ij,j,jl,is_clockwise);Add_To_Material_Surface_Triangle(ij,il,ik,is_clockwise);
    Add_To_Material_Surface_Triangle(jk,jl,j,is_clockwise);
}
//#####################################################################
// Function Add_To_Material_Surface_Subtetrahedron_And_Wedge_And_Half_Oct_Plus_Tet
//#####################################################################
template<class T> void EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>::
Add_To_Material_Surface_Subtetrahedron_And_Wedge_And_Half_Oct_Plus_Tet(const int tetrahedron,const int i,const int j,const int k,const int l) 
{
    int ik=embedded_object.Particle_Embedded_On_Segment(i,k),il=embedded_object.Particle_Embedded_On_Segment(i,l),jk=embedded_object.Particle_Embedded_On_Segment(j,k),
        jl=embedded_object.Particle_Embedded_On_Segment(j,l),kl=embedded_object.Particle_Embedded_On_Segment(k,l);
    bool is_clockwise=permutation_of_four_is_even(embedded_object.orientation_index(tetrahedron));
    if(embedded_object.Node_In_Simplex_Is_Material(i,tetrahedron))
        Add_To_Material_Surface_Subtetrahedron(i,j,k,l,is_clockwise);
    if(embedded_object.Node_In_Simplex_Is_Material(k,tetrahedron) && embedded_object.Node_In_Simplex_Is_Material(l,tetrahedron))        
        Add_To_Material_Surface_Subwedge(k,l,i,j,is_clockwise);
    if(embedded_object.Node_In_Simplex_Is_Material(k,tetrahedron) && !embedded_object.Node_In_Simplex_Is_Material(l,tetrahedron)){        
        if(embedded_object.embedded_object.mesh.Triangle(kl,ik,jk))
            Add_To_Material_Surface_Subtetrahedron(k,l,i,j,is_clockwise);
        else{
            assert(embedded_object.embedded_object.mesh.Triangle(kl,jl,il));
            Add_To_Material_Surface_Half_Oct_Plus_Tet(l,k,j,i,is_clockwise);}}
    if(!embedded_object.Node_In_Simplex_Is_Material(k,tetrahedron) && embedded_object.Node_In_Simplex_Is_Material(l,tetrahedron)){
        if(embedded_object.embedded_object.mesh.Triangle(kl,jl,il))
            Add_To_Material_Surface_Subtetrahedron(l,k,j,i,is_clockwise);
        else{
            assert(embedded_object.embedded_object.mesh.Triangle(kl,ik,jk));
            Add_To_Material_Surface_Half_Oct_Plus_Tet(k,l,i,j,is_clockwise);}}
    if(embedded_object.Node_In_Simplex_Is_Material(j,tetrahedron)) Add_To_Material_Surface_Half_Oct_Plus_Tet(i,j,k,l,is_clockwise);
}
//#####################################################################
template class EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<double>;
#endif

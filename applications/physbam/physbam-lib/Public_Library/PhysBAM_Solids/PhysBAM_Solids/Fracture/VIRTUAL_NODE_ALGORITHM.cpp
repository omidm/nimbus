//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Neil Molino, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VIRTUAL_NODE_ALGORITHM
//#####################################################################
#include <PhysBAM_Tools/Arrays/IDENTITY_ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TRIANGULATED_OBJECT.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/VIRTUAL_NODE_ALGORITHM.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/VIRTUAL_NODES.h>
namespace PhysBAM{
namespace VIRTUAL_NODE_ALGORITHM{
//#####################################################################
// Function Add_Element_If_Necessary
//#####################################################################
template<class TV,int d> void
Add_Element_If_Necessary(EMBEDDED_OBJECT<TV,d>& embedded_object,const VECTOR<int,d+1>& nodes,const int center_node)
{
    SIMPLEX_MESH<d>& mesh=embedded_object.simplicial_object.mesh;
    if(int element=mesh.Simplex(nodes)){embedded_object.node_in_simplex_is_material(element)(center_node)=true;return;}
    mesh.elements.Append(nodes);
    for(int i=1;i<=nodes.m;i++) (*mesh.incident_elements)(nodes[i]).Append(mesh.elements.m);
    embedded_object.embedded_subelements_in_parent_element_index->Append(0);
    embedded_object.node_in_simplex_is_material.Append(VECTOR<bool,d+1>());
    assert(embedded_object.node_in_simplex_is_material.m==mesh.elements.m);
    embedded_object.node_in_simplex_is_material(mesh.elements.m)(center_node)=true;
}
//#####################################################################
// Function Rebuild_Embedded_Object
//#####################################################################
template<class TV,int d> void
Rebuild_Embedded_Object(EMBEDDED_OBJECT<TV,d>& embedded_object,ARRAY<int>& map_to_old_particles,ARRAY<int>& map_to_old_embedded_particles,
    ARRAY<int>& map_to_old_simplices,const bool verbose)
{
    typedef typename MESH_POLICY<d>::MESH T_MESH;
    T_MESH& mesh=embedded_object.simplicial_object.mesh;

    // initialize some auxiliary structures
    mesh.Initialize_Neighbor_Nodes();mesh.Initialize_Incident_Elements();mesh.Initialize_Segment_Mesh();mesh.segment_mesh->Initialize_Incident_Elements();
    embedded_object.embedded_mesh.Initialize_Incident_Elements();

    // build connectivity
    VIRTUAL_NODES virtual_nodes;
    Construct_Virtual_Nodes(embedded_object,map_to_old_particles,virtual_nodes);

    // save old topology (TODO: invert this a bit for speed)
    T_MESH old_mesh(mesh);
    typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d>::OBJECT old_simplicial_object(old_mesh,embedded_object.particles);
    typename EMBEDDING_POLICY<TV,d>::EMBEDDED_OBJECT old_embedded_object(old_simplicial_object);
    old_embedded_object.embedded_particles.Initialize_Subset(embedded_object.embedded_particles);
    old_embedded_object.embedded_mesh.Initialize_Mesh(embedded_object.embedded_mesh);
    old_embedded_object.node_in_simplex_is_material=embedded_object.node_in_simplex_is_material;
    old_embedded_object.parent_particles=embedded_object.parent_particles;
    old_embedded_object.interpolation_fraction=embedded_object.interpolation_fraction;
    old_embedded_object.Initialize_Parents_To_Embedded_Particles_Hash_Table();
    old_mesh.Initialize_Incident_Elements();
    old_embedded_object.Initialize_Embedded_Subelements_In_Parent_Element();

    // clean out structures from current configuration
    mesh.Clean_Memory();embedded_object.simplicial_object.Clean_Memory();
    embedded_object.embedded_mesh.Clean_Memory();embedded_object.node_in_simplex_is_material.Clean_Memory();embedded_object.Delete_Auxiliary_Structures();
    embedded_object.parents_to_embedded_particles_hash_table=new HASHTABLE<VECTOR<int,2>,int>;
    embedded_object.embedded_particles.Update_Subset_Index_From_Element_Index();
    embedded_object.Initialize_Embedded_Subelements_In_Parent_Element();

    mesh.incident_elements=new ARRAY<ARRAY<int> >(mesh.number_nodes);
    for(int node=1;node<=old_mesh.number_nodes;node++) for(int i=1;i<=(*old_mesh.incident_elements)(node).m;i++){
        int element=(*old_mesh.incident_elements)(node)(i);
        VECTOR<int,d+1> nodes=old_mesh.elements(element);int center_node=nodes.Find(node);
        if(old_embedded_object.node_in_simplex_is_material(element)(center_node)){
            VECTOR<int,d+1> donor_nodes;
            for(int i=1;i<=nodes.m;i++) donor_nodes[i]=virtual_nodes.Donor_Node(node,nodes[i]);
            Add_Element_If_Necessary<TV,d>(embedded_object,donor_nodes,center_node);}}
    embedded_object.simplicial_object.mesh.Initialize_Adjacent_Elements();// we want this for biasing and stress smoothing

    // make a map from the new elements to the old ones
    map_to_old_simplices.Resize(mesh.elements.m,false,false);
    for(int t=1;t<=mesh.elements.m;t++) map_to_old_simplices(t)=old_mesh.Simplex(VECTOR<int,d+1>::Map(map_to_old_particles,mesh.elements(t)));

    // add surface particles to current embedded_object, and copy state from old embedded particles
    // make a map from the new embedded particles to the old ones - map to 0, if it doesn't exist
    // reuse each old embedded particle once
    map_to_old_embedded_particles=IDENTITY_ARRAY<>(embedded_object.embedded_particles.active_indices.m);
    ARRAY<bool> embedded_node_already_used(embedded_object.embedded_particles.active_indices.m);
    mesh.Initialize_Segment_Mesh();mesh.segment_mesh->Initialize_Incident_Elements();
    for(int s=1;s<=mesh.segment_mesh->elements.m;s++){
        VECTOR<int,2> nodes=mesh.segment_mesh->elements(s);
        VECTOR<int,2> old_nodes(map_to_old_particles.Subset(nodes));
        int old_embedded_particle=old_embedded_object.Embedded_Particle_On_Segment(old_nodes);
        if(!old_embedded_particle) continue;
        if(old_embedded_object.parent_particles(old_embedded_particle)!=old_nodes) exchange(nodes[1],nodes[2]); // swap if necessary to preserve node order
        if(embedded_node_already_used(old_embedded_particle)){ // add a new particle and copy state from old particle
            int p=embedded_object.Add_Embedded_Particle(nodes,old_embedded_object.interpolation_fraction(old_embedded_particle),false);
            embedded_object.embedded_particles.point_cloud.array_collection->Copy_Element(*embedded_object.embedded_particles.point_cloud.array_collection,old_embedded_particle,p);
            map_to_old_embedded_particles.Append(old_embedded_particle);}
        else{ // reuse existing embedded particle
            embedded_node_already_used(old_embedded_particle)=true;
            embedded_object.parent_particles(old_embedded_particle)=nodes;
            embedded_object.parents_to_embedded_particles_hash_table->Insert(nodes.Sorted(),old_embedded_particle);}}
    delete mesh.segment_mesh;mesh.segment_mesh=0;
    embedded_object.Initialize_Parents_To_Embedded_Particles_Hash_Table();
    embedded_object.Initialize_Embedded_Children();
    assert(!embedded_node_already_used.Number_False());

    // add surface elements to current embedded_object
    embedded_object.embedded_mesh.Initialize_Incident_Elements();
    for(int t=1;t<=mesh.elements.m;t++){
        int old_element=map_to_old_simplices(t);
        VECTOR<int,2*d-2> old_emb_subelements=old_embedded_object.Embedded_Subelements_In_Element(old_element);
        for(int i=1;i<=old_emb_subelements.m && old_emb_subelements[i];i++)
            Add_Embedded_Subelement(embedded_object,old_embedded_object,old_emb_subelements[i],t,old_element);}

    // remap orientation_index if necessary
    if(embedded_object.orientation_index.m) embedded_object.orientation_index=ARRAY<int>(embedded_object.orientation_index).Subset(map_to_old_simplices);

    if(verbose){
        {std::stringstream ss;ss<<"Total Embedded Nodes: "<<embedded_object.embedded_particles.active_indices.m<<std::endl;LOG::filecout(ss.str());}
        {std::stringstream ss;ss<<"Total Embedded Elements: "<<embedded_object.embedded_mesh.elements.m<<std::endl;LOG::filecout(ss.str());}}
}
//#####################################################################
//  Function Construct_Virtual_Nodes
//#####################################################################
template<class TV,int d> void
Construct_Virtual_Nodes(EMBEDDED_OBJECT<TV,d>& embedded_object,ARRAY<int>& map_to_old_particles,VIRTUAL_NODES& virtual_nodes)
{
    assert(embedded_object.simplicial_object.mesh.number_nodes==embedded_object.particles.array_collection->Size() && embedded_object.embedded_mesh.number_nodes==embedded_object.particles.array_collection->Size()
        && embedded_object.embedded_particles.subset_index_from_point_cloud_index.m==embedded_object.particles.array_collection->Size());
    PARTICLES<TV>& particles=dynamic_cast<PARTICLES<TV>&>(embedded_object.particles);
    SIMPLEX_MESH<d>& mesh=embedded_object.simplicial_object.mesh;

    // construct new virtual nodes
    ARRAY<short> marked(mesh.number_nodes,false);
    for(int node=1;node<=mesh.number_nodes;node++){
        int components=Mark_Disconnected_Components_In_One_Ring(embedded_object,node,marked);
        // make a virtual_node for each component disconnected from the center node
        ARRAY<int> virtual_node_index(components);
        for(int component_index=1;component_index<=components;component_index++) virtual_node_index(component_index)=virtual_nodes.Add_Virtual_Node(node);
        for(int i=1;i<=(*mesh.neighbor_nodes)(node).m;i++){
            int neighbor_node=(*mesh.neighbor_nodes)(node)(i);
            if(marked(neighbor_node)>0) virtual_nodes(virtual_node_index(marked(neighbor_node))).recipients.Append(neighbor_node);}}
    virtual_nodes.Initialize_Replicas();

    // initialize virtual_node indices and update particles
    map_to_old_particles=IDENTITY_ARRAY<>(particles.array_collection->Size());
    for(int p=1;p<=virtual_nodes.replicas.m;p++){if(!virtual_nodes.replicas(p).m) continue;
        int i=1;
        {VIRTUAL_NODE& virtual_node=virtual_nodes(virtual_nodes.replicas(p)(i));
        if(!embedded_object.Node_Near_Material(virtual_node.corresponding_real_node)){i++; // reuse old virtual node
            virtual_node.index=virtual_node.corresponding_real_node;}}
        for(;i<=virtual_nodes.replicas(p).m;i++){VIRTUAL_NODE& virtual_node=virtual_nodes(virtual_nodes.replicas(p)(i));
            virtual_node.index=particles.array_collection->Append(*particles.array_collection,virtual_node.corresponding_real_node); // duplicating particles produces more mass
            map_to_old_particles.Append(virtual_node.corresponding_real_node);}}

    // Update meshes and embedded particles if virtual nodes were added
    embedded_object.embedded_particles.subset_index_from_point_cloud_index.Resize(particles.array_collection->Size());
    embedded_object.simplicial_object.mesh.Set_Number_Nodes(particles.array_collection->Size());
    embedded_object.embedded_mesh.Set_Number_Nodes(particles.array_collection->Size());
}
//#####################################################################
//  Function Mark_Disconnected_Components_In_One_Ring
//#####################################################################
template<class TV> static inline bool
Mark_Disconnected_Components_In_One_Ring_Helper(const EMBEDDED_TRIANGULATED_OBJECT<TV>& embedded_object,const int center_node,ARRAY<short>& marked,const int component)
{
    TRIANGLE_MESH& mesh=embedded_object.simplicial_object.mesh;
    bool changed=false;
    for(int t=1;t<=(*mesh.incident_elements)(center_node).m;t++){
        int triangle=(*mesh.incident_elements)(center_node)(t);
        int a,b;mesh.Other_Two_Nodes(center_node,triangle,a,b);
        if(marked(a)==component && !marked(b) && embedded_object.Nodes_Are_Materially_Connected_In_Simplex(a,b,triangle)){marked(b)=component;changed=true;}
        if(marked(b)==component && !marked(a) && embedded_object.Nodes_Are_Materially_Connected_In_Simplex(b,a,triangle)){marked(a)=component;changed=true;}}
    return changed;
}
template<class T> static inline bool
Mark_Disconnected_Components_In_One_Ring_Helper(const EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& embedded_object,const int center_node,ARRAY<short>& marked,const int component)
{
    TETRAHEDRON_MESH& mesh=embedded_object.simplicial_object.mesh;
    bool changed=false;
    for(int t=1;t<=(*mesh.incident_elements)(center_node).m;t++){
        int tetrahedron=(*mesh.incident_elements)(center_node)(t);
        int a,b,c;mesh.Other_Three_Nodes(center_node,tetrahedron,a,b,c);
        if(marked(a)==component){
            if(!marked(b) && embedded_object.Nodes_Are_Materially_Connected_In_Simplex(a,b,tetrahedron)){marked(b)=component;changed=true;}
            if(!marked(c) && embedded_object.Nodes_Are_Materially_Connected_In_Simplex(a,c,tetrahedron)){marked(c)=component;changed=true;}}
        if(marked(b)==component){
            if(!marked(a) && embedded_object.Nodes_Are_Materially_Connected_In_Simplex(b,a,tetrahedron)){marked(a)=component;changed=true;}
            if(!marked(c) && embedded_object.Nodes_Are_Materially_Connected_In_Simplex(b,c,tetrahedron)){marked(c)=component;changed=true;}}
        if(marked(c)==component){
            if(!marked(a) && embedded_object.Nodes_Are_Materially_Connected_In_Simplex(c,a,tetrahedron)){marked(a)=component;changed=true;}
            if(!marked(b) && embedded_object.Nodes_Are_Materially_Connected_In_Simplex(c,b,tetrahedron)){marked(b)=component;changed=true;}}}
    return changed;
}
template<class TV,int d> int
Mark_Disconnected_Components_In_One_Ring(const EMBEDDED_OBJECT<TV,d>& embedded_object_input,const int center_node,ARRAY<short>& marked)
{
    const typename EMBEDDING_POLICY<TV,d>::EMBEDDED_OBJECT& embedded_object=dynamic_cast<const typename EMBEDDING_POLICY<TV,d>::EMBEDDED_OBJECT&>(embedded_object_input);
    SIMPLEX_MESH<d>& mesh=embedded_object.simplicial_object.mesh;

    // flood fill only nodes that are material in some element in the 1-ring
    INDIRECT_ARRAY<ARRAY<short>,ARRAY<int>&> marked_subset=marked.Subset((*mesh.neighbor_nodes)(center_node));
    ARRAYS_COMPUTATIONS::Fill(marked_subset,(short)-2);
    for(int i=1;i<=(*mesh.incident_elements)(center_node).m;i++){int t=(*mesh.incident_elements)(center_node)(i);
        const VECTOR<int,d+1>& element=mesh.elements(t);
        for(int n=1;n<=element.m;n++) if(embedded_object.node_in_simplex_is_material(t)(n)) marked(element[n])=0;}
    // mark component connected to the center node with -1
    for(int n=1;n<=(*mesh.neighbor_nodes)(center_node).m;n++){
        int node=(*mesh.neighbor_nodes)(center_node)(n);
        if(!embedded_object.Segment_Is_Broken(center_node,node)) marked(node)=-1;}

    while(Mark_Disconnected_Components_In_One_Ring_Helper(embedded_object,center_node,marked,-1));
    // mark other components
    int component=0;
    while(int first_unmarked=Unmarked_Neighbor_Node(mesh,center_node,marked)){
        marked(first_unmarked)=++component;
        while(Mark_Disconnected_Components_In_One_Ring_Helper(embedded_object,center_node,marked,component));}
    return component;
}
//#####################################################################
// Function Unmarked_Neighbor_Node
//#####################################################################
template<int d> static inline int
Unmarked_Neighbor_Node(const SIMPLEX_MESH<d>& mesh,const int center_node,const ARRAY<short>& marked)
{
    for(int n=1;n<=(*mesh.neighbor_nodes)(center_node).m;n++){
        int node=(*mesh.neighbor_nodes)(center_node)(n);
        if(!marked(node)) return node;}
    return 0;
}
//#####################################################################
// Function Add_Embedded_Subelement
//#####################################################################
template<class TV,int d> void
Add_Embedded_Subelement(EMBEDDED_OBJECT<TV,d>& embedded_object,const EMBEDDED_OBJECT<TV,d>& old_embedded_object,const int old_emb_subelement,const int current_element,const int old_element)
{
    VECTOR<int,d> old_emb_nodes=VECTOR<int,d>::Map(old_embedded_object.embedded_particles.subset_index_from_point_cloud_index,
        old_embedded_object.embedded_mesh.elements(old_emb_subelement));
    VECTOR<int,d+1> current_nodes=embedded_object.simplicial_object.mesh.elements(current_element),
        old_nodes=old_embedded_object.simplicial_object.mesh.elements(old_element);

    VECTOR<int,d> emb_nodes;
    for(int i=1;i<=emb_nodes.m;i++){
        for(int j1=1;j1<=current_nodes.m-1;j1++) for(int j2=j1+1;j2<=current_nodes.m;j2++)
            if(old_embedded_object.Are_Parents(VECTOR<int,2>(old_nodes[j1],old_nodes[j2]),old_emb_nodes[i])){
                emb_nodes[i]=embedded_object.Embedded_Particle_On_Segment(current_nodes[j1],current_nodes[j2]);goto NEXT;}
        NEXT:;}

    embedded_object.Add_Embedded_Subelement_If_Not_Already_There(emb_nodes);
}
//#####################################################################
template void Rebuild_Embedded_Object(EMBEDDED_OBJECT<VECTOR<float,2>,2>&,ARRAY<int>&,ARRAY<int>&,ARRAY<int>&,const bool);
template void Rebuild_Embedded_Object(EMBEDDED_OBJECT<VECTOR<float,3>,2>&,ARRAY<int>&,ARRAY<int>&,ARRAY<int>&,const bool);
template void Rebuild_Embedded_Object(EMBEDDED_OBJECT<VECTOR<float,3>,3>&,ARRAY<int>&,ARRAY<int>&,ARRAY<int>&,const bool);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void Rebuild_Embedded_Object(EMBEDDED_OBJECT<VECTOR<double,2>,2>&,ARRAY<int>&,ARRAY<int>&,ARRAY<int>&,const bool);
template void Rebuild_Embedded_Object(EMBEDDED_OBJECT<VECTOR<double,3>,2>&,ARRAY<int>&,ARRAY<int>&,ARRAY<int>&,const bool);
template void Rebuild_Embedded_Object(EMBEDDED_OBJECT<VECTOR<double,3>,3>&,ARRAY<int>&,ARRAY<int>&,ARRAY<int>&,const bool);
#endif
}
}

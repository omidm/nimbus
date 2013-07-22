//#####################################################################
// Copyright 2006-2008, Avi Robinson-Mosher, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RED_GREEN_TRIANGLES
//##################################################################### 
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Dynamics/Meshing/RED_GREEN_TRIANGLES.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void RED_GREEN_TRIANGLES<TV>::
Initialize()
{
    Clean_Memory();
    meshes.Append(new TRIANGLE_MESH(object.mesh));
    meshes(1)->Initialize_Incident_Elements();
    meshes(1)->Initialize_Segment_Mesh(); // this is deleted below
    meshes(1)->Initialize_Element_Edges();
    parent.Append(0); // the first level can't have parents
    children.Append(new ARRAY<VECTOR<int,4> >(meshes(1)->elements.m));
    segment_mesh.number_nodes=meshes(1)->number_nodes;
    element_edges.Resize(1);
    segment_mesh.elements.Exchange(meshes(1)->segment_mesh->elements);
    element_edges(1).Exchange(*meshes(1)->element_edges);
    delete meshes(1)->segment_mesh;meshes(1)->segment_mesh=0; // no longer needed
    delete meshes(1)->element_edges;meshes(1)->element_edges=0; // ditto
    segment_mesh.Initialize_Incident_Elements();
    segment_midpoints.Resize(segment_mesh.elements.m);
    index_in_stack.Append(new ARRAY<int>(meshes(1)->elements.m));
    leaf_levels_and_indices.Resize(meshes(1)->elements.m);
    for(int t=1;t<=meshes(1)->elements.m;t++){leaf_levels_and_indices(t)(1)=1;leaf_levels_and_indices(t)(2)=t;}
    leaf_number.Append(new ARRAY<int>(meshes(1)->elements.m));
    for(int t=1;t<=meshes(1)->elements.m;t++) (*leaf_number(1))(t)=t;
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class TV> void RED_GREEN_TRIANGLES<TV>::
Clean_Memory()
{
    leaf_levels_and_indices.Clean_Memory();leaf_number.Delete_Pointers_And_Clean_Memory();meshes.Delete_Pointers_And_Clean_Memory();
    parent.Delete_Pointers_And_Clean_Memory();children.Delete_Pointers_And_Clean_Memory();segment_mesh.Clean_Memory();
    segment_midpoints.Clean_Memory();stack.Clean_Memory();index_in_stack.Delete_Pointers_And_Clean_Memory();
    element_edges.Clean_Memory();
    if(segment_index_from_midpoint_index){delete segment_index_from_midpoint_index;segment_index_from_midpoint_index=0;}
}
//#####################################################################
// Function Refine_Simplex_List
//#####################################################################
template<class TV> template<class T_ARRAY> void RED_GREEN_TRIANGLES<TV>::
Refine_Simplex_List(const T_ARRAY& triangle_list)
{
    STATIC_ASSERT((IS_SAME<int,typename T_ARRAY::ELEMENT>::value));
    object.particles.array_collection->Preallocate(object.particles.array_collection->Size()+3*triangle_list.Size());
    for(int level=1;level<=index_in_stack.m;level++) ARRAYS_COMPUTATIONS::Fill(*index_in_stack(level),0);
    for(int i=1;i<=triangle_list.Size();i++){
        int level,tri;leaf_levels_and_indices(triangle_list(i)).Get(level,tri);
        if(!Red(level,tri)){tri=(*parent(level))(tri);level-=1;}
        stack.Append(VECTOR<int,2>(level,tri));(*index_in_stack(level))(tri)=stack.m;
        Add_Midpoints_If_Not_Already_There(level,tri);}
    Resolve_Stack();
    Rebuild_Object();
}
//#####################################################################
// Function Resolve_Stack
//#####################################################################
template<class TV> void RED_GREEN_TRIANGLES<TV>::
Resolve_Stack()
{
    while(stack.m){
        int level,tri;stack.Pop().Get(level,tri);
        if(level){assert(tri);(*index_in_stack(level))(tri)=0;if(!Regularly_Refined(level,tri)) Refine_If_Necessary(level,tri);}}
}
//#####################################################################
// Function Refine_If_Necessary
//#####################################################################
template<class TV> void RED_GREEN_TRIANGLES<TV>::
Refine_If_Necessary(const int level,const int tri)
{
    // first check if we're already consistent, and don't need to refine at all.
    ARRAY<int> midpoints(3),subedges(9);Get_Existing_Subindices(level,tri,midpoints,subedges);
    int number_midpoints=3-midpoints.Count_Matches(0);if(number_midpoints==0) return;
    int number_children=0;while(number_children<4 && (*children(level))(tri)(number_children+1)) number_children++;
    if((number_midpoints==1 && number_children==2) || (number_midpoints==3 && number_children==4)) return;

    if(!Red(level,tri)) return Regularly_Refine_Triangle(level-1,(*parent(level))(tri)); //Green case: refine parent
 
    if(number_midpoints>1) return Regularly_Refine_Triangle(level,tri); // too many to go green

    // Green refinement
    Ensure_Level_Exists(level+1);
    ARRAY<int> free_triangle_indices,free_edge_indices;free_triangle_indices.Preallocate(2);free_edge_indices.Preallocate(5);
    Delete_Children(level,tri,free_triangle_indices,free_edge_indices); // get rid of what's here already
    // figure out which green triangles and green interior segment need to be made, and make them
    int i,j,k;meshes(level)->elements(tri).Get(i,j,k);
    // make sure the subedges along triangles's edges exist, and that green bisectors exist
    if(midpoints(1)){ // if edge ij is split
        if(!subedges(1)) Add_Segment(free_edge_indices,i,midpoints(1));if(!subedges(3)) Add_Segment(free_edge_indices,j,midpoints(1));
        if(!midpoints(2) && !midpoints(3) && !segment_mesh.Segment(midpoints(1),k)) Add_Segment(free_edge_indices,midpoints(1),k);}
    if(midpoints(2)){ // if edge jk is split
        if(!subedges(4)) Add_Segment(free_edge_indices,j,midpoints(2));if(!subedges(5)) Add_Segment(free_edge_indices,k,midpoints(2));
        if(!midpoints(1) && !midpoints(3) && !segment_mesh.Segment(midpoints(2),i)) Add_Segment(free_edge_indices,midpoints(2),i);}
    if(midpoints(3)){ // if edge ki is split
        if(!subedges(2)) Add_Segment(free_edge_indices,i,midpoints(3));if(!subedges(6)) Add_Segment(free_edge_indices,k,midpoints(3));
        if(!midpoints(1) && !midpoints(2) && !segment_mesh.Segment(midpoints(3),j)) Add_Segment(free_edge_indices,midpoints(3),j);}
    // go ahead and create any missing edges along with the triangles
    if(midpoints(1)){ // if edge ij is split
        Add_Triangle(free_triangle_indices,level+1,i,midpoints(1),k,tri);Add_Triangle(free_triangle_indices,level+1,midpoints(1),j,k,tri);}
    else if(midpoints(2)){ // if edge jk is split
        Add_Triangle(free_triangle_indices,level+1,i,j,midpoints(2),tri);Add_Triangle(free_triangle_indices,level+1,i,midpoints(2),k,tri);}
    else if(midpoints(3)){ // if edge ki is split
        Add_Triangle(free_triangle_indices,level+1,j,k,midpoints(3),tri);Add_Triangle(free_triangle_indices,level+1,j,midpoints(3),i,tri);}
}
//#####################################################################
// Function Regularly_Refine_Tri
//#####################################################################
template<class TV> void RED_GREEN_TRIANGLES<TV>::
Regularly_Refine_Triangle(const int level,const int tri)
{
    assert(Red(level,tri));
    (*index_in_stack(level))(tri)=-1; // mark tri so we won't consider this tri for resolution in stack
    if(Regularly_Refined(level,tri)) return;
    Ensure_Level_Exists(level+1);
    ARRAY<int> free_triangle_indices,free_edge_indices;
    free_triangle_indices.Preallocate(2);free_edge_indices.Preallocate(5);
    Delete_Children(level,tri,free_triangle_indices,free_edge_indices);
    ARRAY<int> midpoints(3),subedges(9);Get_Existing_Subindices(level,tri,midpoints,subedges);
    for(int a=1;a<=3;a++) if(!midpoints(a)) midpoints(a)=Add_Midpoint(element_edges(level)(tri)(a),level,tri);
    int i,j,k;meshes(level)->elements(tri).Get(i,j,k);
    int ij=midpoints(1),jk=midpoints(2),ki=midpoints(3);
    if(!subedges(1)) Add_Segment(free_edge_indices,i,ij);
    if(!subedges(2)) Add_Segment(free_edge_indices,i,ki);
    if(!subedges(3)) Add_Segment(free_edge_indices,j,ij);
    if(!subedges(4)) Add_Segment(free_edge_indices,j,jk);
    if(!subedges(5)) Add_Segment(free_edge_indices,k,jk);
    if(!subedges(6)) Add_Segment(free_edge_indices,k,ki);
    if(!subedges(7)) Add_Segment(free_edge_indices,ij,jk);
    if(!subedges(8)) Add_Segment(free_edge_indices,jk,ki);
    if(!subedges(9)) Add_Segment(free_edge_indices,ki,ij); // TODO: check that these are always added
    Add_Triangle(free_triangle_indices,level+1,i,ij,ki,tri);Add_Triangle(free_triangle_indices,level+1,ij,j,jk,tri);
    Add_Triangle(free_triangle_indices,level+1,ki,jk,k,tri);Add_Triangle(free_triangle_indices,level+1,ij,jk,ki,tri);
}
//#####################################################################
// Function Get_Existing_Subindices
//#####################################################################
// midpoints (or 0 if no midpoint) listed in the edge order: ij,jk,ki (midpoints has m==3)
// edge-subedges listed in the order: i-ij,i-ki, j-ij,j-jk, k-jk,k-ki, ij-jk,jk-ki,ki-ij (subedges has m==9)
template<class TV> void RED_GREEN_TRIANGLES<TV>::
Get_Existing_Subindices(const int level,const int tri,ARRAY<int>& midpoints,ARRAY<int>& subedges)
{
    assert(midpoints.m==3 && subedges.m==9);
    for(int e=1;e<=3;e++) midpoints(e)=segment_midpoints(element_edges(level)(tri)(e));
    int i,j,k;meshes(level)->elements(tri).Get(i,j,k);
    int ij=midpoints(1),jk=midpoints(2),ki=midpoints(3);
    if(ij){subedges(1)=segment_mesh.Segment(i,ij);subedges(3)=segment_mesh.Segment(j,ij);if(jk) subedges(7)=segment_mesh.Segment(ij,jk);if(ki) subedges(9)=segment_mesh.Segment(ki,ij);}
    if(jk){subedges(4)=segment_mesh.Segment(j,jk);subedges(5)=segment_mesh.Segment(k,jk);if(ki) subedges(8)=segment_mesh.Segment(jk,ki);}
    if(ki){subedges(2)=segment_mesh.Segment(i,ki);subedges(6)=segment_mesh.Segment(k,ki);}
}
//#####################################################################
// Function Ensure_Level_Exists
//#####################################################################
template<class TV> void RED_GREEN_TRIANGLES<TV>::
Ensure_Level_Exists(const int level)
{
    if(meshes.m >= level) return; // it already exists
    meshes.Append(new TRIANGLE_MESH);
    element_edges.Resize(level);
    meshes(level)->number_nodes=object.particles.array_collection->Size();
    meshes(level)->incident_elements=new ARRAY<ARRAY<int> >(object.particles.array_collection->Size());
    parent.Append(new ARRAY<int>(0));children.Append(new ARRAY<VECTOR<int,4> >());
    index_in_stack.Append(new ARRAY<int>(0));leaf_number.Append(new ARRAY<int>(0));
}
//#####################################################################
// Function Delete_Children
//#####################################################################
// This deletes children triangles (should be green) and any green interior edges belonging *only* to those triangles.
// Note that it doesn't touch exterior subedges that fit in red triangles
template<class TV> void RED_GREEN_TRIANGLES<TV>::
Delete_Children(const int level,const int tri,ARRAY<int>& deleted_tri_indices,ARRAY<int>& deleted_edge_indices)
{
    if(Leaf(level,tri)) return;
    // basic information about triangle
    int i,j,k;meshes(level)->elements(tri).Get(i,j,k);
    int ij=segment_midpoints(element_edges(level)(tri)(1)),jk=segment_midpoints(element_edges(level)(tri)(2)),ki=segment_midpoints(element_edges(level)(tri)(3));
     // make the list of deleted triangles (namely, the children) and zero out children
    int p;for(p=1;p<=4&&(*children(level))(tri)(p);p++){deleted_tri_indices.Append((*children(level))(tri)(p));(*children(level))(tri)(p)=0;}
    // get a list of edges to delete (begin by finding all children edges, then filter the red ones out)
    ARRAY<int> children_edges;children_edges.Preallocate(5);
    for(p=1;p<=deleted_tri_indices.m;p++) for(int q=1;q<=3;q++) // get a list of all children edges
        children_edges.Append_Unique(element_edges(level+1)(deleted_tri_indices(p))(q));
    for(p=1;p<=children_edges.m;p++){
        int q,r;segment_mesh.elements(children_edges(p)).Get(q,r);
        if((q==i && r==jk) || (q==j && r==ki) || (q==k && r==ij))
            deleted_edge_indices.Append(children_edges(p));} // delete interior edges
    // now do the actual deletion
    for(p=1;p<=deleted_tri_indices.m;p++) for(int q=1;q<=3;q++){ // remove deleted triangles from incident_triangles
        int node=meshes(level+1)->elements(deleted_tri_indices(p))(q),index=0;
        (*meshes(level+1)->incident_elements)(node).Find(deleted_tri_indices(p),index);assert(index);
        (*meshes(level+1)->incident_elements)(node).Remove_Index_Lazy(index);}
    for(p=1;p<=deleted_edge_indices.m;p++) for(int q=1;q<=2;q++){ // remove deleted edges from incident_segments
        int node=segment_mesh.elements(deleted_edge_indices(p))(q),index=0;
        (*segment_mesh.incident_elements)(node).Find(deleted_edge_indices(p),index);assert(index);
        (*segment_mesh.incident_elements)(node).Remove_Index_Lazy(index);}
    // then zero out occurances in the stack
    for(p=1;p<=deleted_tri_indices.m;p++){
        int t=deleted_tri_indices(p);
        if((*index_in_stack(level+1))(t)){
            stack((*index_in_stack(level+1))(t)).Set(0,0); // can't remove since it would screw up index_in_stack, mark as no longer relevent
            (*index_in_stack(level+1))(t)=0;}}
}
//#####################################################################
// Function Add_Midpoints_If_Not_Already_There
//#####################################################################
template<class TV> void RED_GREEN_TRIANGLES<TV>::
Add_Midpoints_If_Not_Already_There(const int level,const int tri)
{
    for(int i=1;i<=3;i++){int s=element_edges(level)(tri)(i);if(!segment_midpoints(s)) Add_Midpoint(s,level,tri);}
}
//#####################################################################
// Function Add_Midpoint
//#####################################################################
template<class TV> int RED_GREEN_TRIANGLES<TV>::
Add_Midpoint(const int segment,const int level,const int tri)
{
    assert(segment_midpoints(segment)==0);
    int node1,node2,new_node=0;segment_mesh.elements(segment).Get(node1,node2);
    if(free_segment_midpoints && free_segment_midpoints->Get(VECTOR<int,2>(node1,node2).Sorted(),new_node))
        free_segment_midpoints->Delete(VECTOR<int,2>(node1,node2).Sorted());
    else new_node=object.particles.array_collection->Add_Element_From_Deletion_List();
    // define midpoint's position, velocity and acceleration
    object.particles.X(new_node)=(T).5*(object.particles.X(node1)+object.particles.X(node2));
    if(object.particles.store_velocity)
        object.particles.V(new_node)=(T).5*(object.particles.V(node1)+object.particles.V(node2));
    if(rest_position) (*rest_position)(new_node)=(T).5*((*rest_position)(node1)+(*rest_position)(node2));
    // list it
    segment_midpoints(segment)=new_node;
    // add edge neighbors of triangle to the stack if not already there
    int i;for(i=1;i<=(*meshes(level)->incident_elements)(node1).m;i++){
        int s=(*meshes(level)->incident_elements)(node1)(i);
        if(s != tri && !(*index_in_stack(level))(s)){
            int a,b,c;meshes(level)->elements(s).Get(a,b,c);assert(a == node1 || b == node1 || c == node1);
            if(node2 == a || node2 == b || node2 == c){ // not incident on edge
                stack.Append(VECTOR<int,2>(level,s));(*index_in_stack(level))(s)=stack.m;}}}
    // check previous level for triangles incident on the edge
    if(level > 1 && node1 <= meshes(level-1)->number_nodes) for(i=1;i<=(*meshes(level-1)->incident_elements)(node1).m;i++){
        int s=(*meshes(level-1)->incident_elements)(node1)(i);
        if(!(*index_in_stack(level-1))(s)){
            int a,b,c;meshes(level-1)->elements(s).Get(a,b,c);
            if(node2 == a || node2 == b || node2 == c){ // not incident on edge
                stack.Append(VECTOR<int,2>(level-1,s));(*index_in_stack(level-1))(s)=stack.m;}}}
    // check next level for triangles incident on the edge
    if(level < meshes.m && node1 <= meshes(level+1)->number_nodes) for(i=1;i<=(*meshes(level+1)->incident_elements)(node1).m;i++){
        int s=(*meshes(level+1)->incident_elements)(node1)(i);
        if(!(*index_in_stack(level+1))(s)){
            int a,b,c;meshes(level+1)->elements(s).Get(a,b,c);assert(a == node1 || b == node1 || c == node1);
            if(node2 == a || node2 == b || node2 == c){ // not incident on edge
                stack.Append(VECTOR<int,2>(level+1,s));(*index_in_stack(level+1))(s)=stack.m;}}}
    // finally make sure that number_nodes in the triangle meshes and the segment mesh is up to date
    if(new_node > segment_mesh.number_nodes){
        segment_mesh.number_nodes=new_node;segment_mesh.incident_elements->Resize(new_node);}
    for(i=level;i<=meshes.m;i++) if(new_node > meshes(i)->number_nodes){
        meshes(i)->number_nodes=new_node;meshes(i)->incident_elements->Resize(new_node);}
    return new_node;
}
//#####################################################################
// Function Add_Segment
//#####################################################################
template<class TV> void RED_GREEN_TRIANGLES<TV>::
Add_Segment(ARRAY<int>& free_edge_indices,const int node1,const int node2)
{
    assert(!segment_mesh.Segment(node1,node2));
    int index;
    if(free_edge_indices.m == 0){ // no more empty spots
        segment_mesh.elements.Append(VECTOR<int,2>(node1,node2));index=segment_mesh.elements.m;segment_midpoints.Append(0);}
    else{
        index=free_edge_indices.Pop();assert(!segment_midpoints(index));segment_mesh.elements(index).Set(node1,node2);}
    (*segment_mesh.incident_elements)(node1).Append(index);(*segment_mesh.incident_elements)(node2).Append(index);
}
//#####################################################################
// Function Add_Triangle
//#####################################################################
template<class TV> void RED_GREEN_TRIANGLES<TV>::
Add_Triangle(ARRAY<int>& free_triangle_indices,const int level,const int i,const int j,const int k,const int parent_index)
{
    TRIANGLE_MESH& triangle_mesh=*meshes(level);
    int index;
    if(free_triangle_indices.m == 0){
        triangle_mesh.elements.Append(VECTOR<int,3>(i,j,k));
        leaf_number(level)->Append(0);parent(level)->Append(parent_index);index=triangle_mesh.elements.m;
        element_edges(level).Resize(index);children(level)->Resize(index);index_in_stack(level)->Resize(index);}
    else{
        index=free_triangle_indices.Pop();triangle_mesh.elements(index).Set(i,j,k);
        (*parent(level))(index)=parent_index;for(int a=1;a<=4;a++) (*children(level))(index)(a)=0;(*index_in_stack(level))(index)=0;}
    // create triangle edges
    int ij=segment_mesh.Segment(i,j),jk=segment_mesh.Segment(j,k),ki=segment_mesh.Segment(k,i);
    element_edges(level)(index).Set(ij,jk,ki);
    // check if this new tri has a T-junction
    if(segment_midpoints(ij) || segment_midpoints(jk) || segment_midpoints(ki)){stack.Append(VECTOR<int,2>(level,index));(*index_in_stack(level))(index)=stack.m;}
    // add tri to incident_tris
    int a;for(a=1;a<=3;a++) (*triangle_mesh.incident_elements)(triangle_mesh.elements(index)(a)).Append(index);
    // add tri to parent's list of children
    for(a=1;a<=4;a++) if(!(*children(level-1))(parent_index)(a)){(*children(level-1))(parent_index)(a)=index;break;}
}
//#####################################################################
// Function Rebuild_Object
//#####################################################################
template<class TV> void RED_GREEN_TRIANGLES<TV>::
Rebuild_Object()
{
    int number_of_leaves=0;
    int level,tri;for(level=1;level<=meshes.m;level++) for(tri=1;tri<=meshes(level)->elements.m; tri++) if(Leaf(level,tri)) number_of_leaves++;
    object.mesh.elements.Resize(number_of_leaves);
    leaf_levels_and_indices.Exact_Resize(number_of_leaves);
    int index_into_triangles=1;
    for(level=1;level<=meshes.m;level++) for(tri=1;tri<=meshes(level)->elements.m;tri++) 
        if(Leaf(level,tri)){
            object.mesh.elements(index_into_triangles)=meshes(level)->elements(tri);
            leaf_levels_and_indices(index_into_triangles).Set(level,tri);(*leaf_number(level))(tri)=index_into_triangles;index_into_triangles++;}
        else (*leaf_number(level))(tri)=0;
    object.mesh.number_nodes=object.particles.array_collection->Size();
    object.mesh.Refresh_Auxiliary_Structures();
}
//#####################################################################
// Function Unrefined_Parents
//#####################################################################
template<class TV> void RED_GREEN_TRIANGLES<TV>::
Unrefined_Parents(const int node,ARRAY<int>& parents,ARRAY<T>& weights) const
{
    parents.Remove_All();parents.Append(node);weights.Remove_All();weights.Append((T)1);
    for(int i=1;i<=parents.m;){
        if(!(*segment_index_from_midpoint_index)(parents(i))){i++;continue;}
        T old_weight=weights(i);VECTOR<int,2> segment=segment_mesh.elements((*segment_index_from_midpoint_index)(parents(i)));
        parents.Remove_Index_Lazy(i);weights.Remove_Index_Lazy(i);
        for(int j=1;j<=2;j++){
            int new_parent=segment(j);
            int index=parents.Find(new_parent);if(!index){index=parents.Append(new_parent);weights.Append((T)0);}
            weights(index)+=(T).5*old_weight;}}
}
//#####################################################################
// Function Remove_Simplex_List
//#####################################################################
template<class TV> void RED_GREEN_TRIANGLES<TV>::
Remove_Simplex_List(const ARRAY<int>& triangle_list,ARRAY<HASHTABLE<int,int> >* level_simplex_maps)
{
    assert(!stack.m);

    // build a list of all triangles on all levels to be deleted by going down the tree from each input triangle
    ARRAY<ARRAY<int> > level_triangle_list(meshes.m);
    for(int i=1;i<=triangle_list.m;i++){
        int level,tri;leaf_levels_and_indices(triangle_list(i)).Get(level,tri);
        do{level_triangle_list(level).Append_Unique(tri);
        }while(level>1 && (tri=(*parent(level--))(tri)));}
    
    ARRAY<HASHTABLE<int,int> > default_level_simplex_maps;if(!level_simplex_maps) level_simplex_maps=&default_level_simplex_maps;
    (*level_simplex_maps).Resize(meshes.m);

    for(int level=1;level<=meshes.m;level++){
        // remove triangles to be deleted from incident elements
        // TODO: is incident_elements special or not?
        for(int i=1;i<=level_triangle_list(level).m;i++) for(int j=1;j<=3;j++){int node=meshes(level)->elements(level_triangle_list(level)(i))(j);
            int index=0;(*meshes(level)->incident_elements)(node).Find(level_triangle_list(level)(i),index);assert(index);
            (*meshes(level)->incident_elements)(node).Remove_Index_Lazy(index);}

        Sort(level_triangle_list(level));

        HASHTABLE<int,int>& simplex_map=(*level_simplex_maps)(level);simplex_map.Remove_All();
        
        meshes(level)->Delete_Sorted_Elements(level_triangle_list(level),simplex_map);
        if(level<parent.m) for(int i=1;i<=(*parent(level+1)).m;i++) simplex_map.Get((*parent(level+1))(i),(*parent(level+1))(i));
        if(parent(level)) parent(level)->Remove_Sorted_Indices_Lazy(level_triangle_list(level));
        if(level>1) for(int i=1;i<=children(level-1)->m;i++){for(int j=1;j<=4;j++) simplex_map.Get((*children(level-1))(i)(j),(*children(level-1))(i)(j));
            (*children(level-1))(i)=(*children(level-1))(i).Sorted().Reversed();}
        children(level)->Remove_Sorted_Indices_Lazy(level_triangle_list(level));
        for(int i=1;i<=leaf_levels_and_indices.m;i++) if(leaf_levels_and_indices(i)(1)==level) simplex_map.Get(leaf_levels_and_indices(i)(2),leaf_levels_and_indices(i)(2));
        leaf_number(level)->Remove_Sorted_Indices_Lazy(level_triangle_list(level));

        // delete edges which have a node that is no longer referenced by the elements
        ARRAY<int> nodes_referenced(object.particles.array_collection->Size());meshes(level)->Mark_Nodes_Referenced(nodes_referenced,1);
        ARRAY<int> deleted_edge_indices;
        for(int i=1;i<=level_triangle_list(level).m;i++) for(int j=1;j<=3;j++){
            int edge=element_edges(level)(level_triangle_list(level)(i))(j);int node1,node2;segment_mesh.elements(edge).Get(node1,node2);
            if(!nodes_referenced(node1) || !nodes_referenced(node2)){
                int index=0;(*segment_mesh.incident_elements)(node1).Find(edge,index);if(index) (*segment_mesh.incident_elements)(node1).Remove_Index_Lazy(index);
                index=0;(*segment_mesh.incident_elements)(node2).Find(edge,index);if(index) (*segment_mesh.incident_elements)(node2).Remove_Index_Lazy(index);}}}

    Rebuild_Object();
    // TODO: if node has 0 parent move it down to level 1
}
//#####################################################################
// Function Print
//#####################################################################
template<class TV> void RED_GREEN_TRIANGLES<TV>::
Print() const
{
    std::stringstream ss;
    for(int level=1;level<=meshes.m;level++){
        ss << "There are " << meshes(level)->elements.m << " elements at level " << level << std::endl;
        ss << "      and object.particles.array_collection->Size()=" << object.particles.array_collection->Size() << std::endl;
        ss << "      and meshes(level)->incident_elements->m="  << meshes(level)->incident_elements->m << std::endl;}
    LOG::filecout(ss.str());
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template class RED_GREEN_TRIANGLES<VECTOR<T,d> >; \
    template void RED_GREEN_TRIANGLES<VECTOR<T,d> >::Refine_Simplex_List(const ARRAY<int>&); \
    template void RED_GREEN_TRIANGLES<VECTOR<T,d> >::Refine_Simplex_List(const IDENTITY_ARRAY<>&);
INSTANTIATION_HELPER(float,2)
INSTANTIATION_HELPER(float,3)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,2)
INSTANTIATION_HELPER(double,3)
#endif

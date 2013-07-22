//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RED_GREEN_TETRAHEDRA
//##################################################################### 
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Dynamics/Meshing/RED_GREEN_TETRAHEDRA.h>
using namespace PhysBAM; 
template<class T> RED_GREEN_TETRAHEDRA<T>::
RED_GREEN_TETRAHEDRA(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume_input)
    :object(tetrahedralized_volume_input),segment_index_from_midpoint_index(0),rest_position(0)
{
    Initialize();
}
template<class T> RED_GREEN_TETRAHEDRA<T>::
RED_GREEN_TETRAHEDRA(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume_input,ARRAY<VECTOR<T,3> >& rest_position_input)
    :object(tetrahedralized_volume_input),rest_position(&rest_position_input)
{
    PHYSBAM_FATAL_ERROR();Initialize(); // This constructor was broken.
}
template<class T> RED_GREEN_TETRAHEDRA<T>::
~RED_GREEN_TETRAHEDRA()
{
    Clean_Memory();
}
//#####################################################################
// Function Refine_Simplex_List
//#####################################################################
template<class T> void RED_GREEN_TETRAHEDRA<T>::
Refine_Simplex_List(const ARRAY<int>& tetrahedron_list)
{
    object.particles.array_collection->Preallocate(object.particles.array_collection->Size()+6*tetrahedron_list.m);
    for(int level=1;level<=index_in_stack.m;level++) ARRAYS_COMPUTATIONS::Fill(*index_in_stack(level),0);
    for(int i=1;i<=tetrahedron_list.m;i++){
        int level=leaf_levels_and_indices(tetrahedron_list(i))(1),tet=leaf_levels_and_indices(tetrahedron_list(i))(2);
        if(Red(level,tet)){
            stack.Append(VECTOR<int,2>(level,tet));(*index_in_stack(level))(tet)=stack.m;
            for(int j=1;j<=6;j++){int s=(*meshes(level)->element_edges)(tet)(j);if(!segment_midpoints(s)) Add_Midpoint(s,level,tet);}}
        else{ // Green
            int parent_index=(*parent(level))(tet);
            stack.Append(VECTOR<int,2>(level-1,parent_index));(*index_in_stack(level-1))(parent_index)=stack.m;
            for(int j=1;j<=6;j++){int s=(*meshes(level-1)->element_edges)(parent_index)(j);if(!segment_midpoints(s)) Add_Midpoint(s,level-1,parent_index);}}}
    Resolve_Stack();
    Rebuild_Object();
}
//#####################################################################
// Function Subdivide_Segment_List
//#####################################################################
template<class T> void RED_GREEN_TETRAHEDRA<T>::
Subdivide_Segment_List(const ARRAY<VECTOR<int,2> >& segment_list)
{
    object.particles.array_collection->Preallocate(object.particles.array_collection->Size()+4*segment_list.m);
    for(int i=1;i<=index_in_stack.m;i++) ARRAYS_COMPUTATIONS::Fill(*index_in_stack(i),0);
    for(int i=1;i<=segment_list.m;i++){
        int node1,node2;segment_list(i).Get(node1,node2);
        int s,level,tet;
        if(Find_Edge_And_Test_If_Green(node1,node2,&s,&level,&tet)){
            // then we should refine the parent's segments instead - BUT if node1-node2 was a bisector, only subdivide that face!
            int parent_index=(*parent(level))(tet);
            stack.Append(VECTOR<int,2>(level-1,parent_index));(*index_in_stack(level-1))(parent_index)=stack.m;
            int a,b,c,d;meshes(level-1)->elements(parent_index).Get(a,b,c,d);
            int ab=(*meshes(level-1)->element_edges)(parent_index)(1),bc=(*meshes(level-1)->element_edges)(parent_index)(2),
                 ca=(*meshes(level-1)->element_edges)(parent_index)(3),ad=(*meshes(level-1)->element_edges)(parent_index)(4),
                 bd=(*meshes(level-1)->element_edges)(parent_index)(5),cd=(*meshes(level-1)->element_edges)(parent_index)(6);
            if((node1 == a && node2 == segment_midpoints(bc)) || (node2 == a && node1 == segment_midpoints(bc)) || (node1 == b && node2 == segment_midpoints(ca)) || 
                (node2 == b && node1 == segment_midpoints(ca)) || (node1 == c && node2 == segment_midpoints(ab)) || (node2 == c && node1 == segment_midpoints(ab))){ // bisects face abc
                if(!segment_midpoints(ab)) Add_Midpoint(ab,level-1,parent_index);if(!segment_midpoints(bc)) Add_Midpoint(bc,level-1,parent_index);
                if(!segment_midpoints(ca)) Add_Midpoint(ca,level-1,parent_index);}
            else if((node1 == a && node2 == segment_midpoints(bd)) || (node2 == a && node1 == segment_midpoints(bd)) || (node1 == b && node2 == segment_midpoints(ad)) || 
                (node2 == b && node1 == segment_midpoints(ad)) || (node1 == d && node2 == segment_midpoints(ab)) || (node2 == d && node1 == segment_midpoints(ab))){ // bisects face abd
                if(!segment_midpoints(ab)) Add_Midpoint(ab,level-1,parent_index);if(!segment_midpoints(bd)) Add_Midpoint(bd,level-1,parent_index);
                if(!segment_midpoints(ad)) Add_Midpoint(ad,level-1,parent_index);}
            else if((node1 == a && node2 == segment_midpoints(cd)) || (node2 == a && node1 == segment_midpoints(cd)) || (node1 == c && node2 == segment_midpoints(ad)) || 
                (node2 == c && node1 == segment_midpoints(ad)) || (node1 == d && node2 == segment_midpoints(ca)) || (node2 == d && node1 == segment_midpoints(ca))){ // bisects face acd
                if(!segment_midpoints(ca)) Add_Midpoint(ca,level-1,parent_index);if(!segment_midpoints(cd)) Add_Midpoint(cd,level-1,parent_index);
                if(!segment_midpoints(ad)) Add_Midpoint(ad,level-1,parent_index);}
            else if((node1 == b && node2 == segment_midpoints(cd)) || (node2 == b && node1 == segment_midpoints(cd)) || (node1 == c && node2 == segment_midpoints(bd)) || 
                (node2 == c && node1 == segment_midpoints(bd)) || (node1 == d && node2 == segment_midpoints(bc)) || (node2 == d && node1 == segment_midpoints(bc))){ // bisects face bcd
                if(!segment_midpoints(bc)) Add_Midpoint(bc,level-1,parent_index);if(!segment_midpoints(cd)) Add_Midpoint(cd,level-1,parent_index);
                if(!segment_midpoints(bd)) Add_Midpoint(bd,level-1,parent_index);}
            else for(int j=1;j<=6;j++){ // this was an internal green edge, so subdivide all parent edges
                s=(*meshes(level-1)->element_edges)(parent_index)(j);if(!segment_midpoints(s)) Add_Midpoint(s,level-1,parent_index);}}
        else{PHYSBAM_ASSERT(s); // red edge
            stack.Append(VECTOR<int,2>(level,tet));(*index_in_stack(level))(tet)=stack.m;if(!segment_midpoints(s)) Add_Midpoint(s,level,tet);}}

    Resolve_Stack();
    Rebuild_Object();
}
//#####################################################################
// Function Find_Edge_And_Test_If_Green
//#####################################################################
template<class T> bool RED_GREEN_TETRAHEDRA<T>::
Find_Edge_And_Test_If_Green(const int node1,const int node2,int* segment_output,int* level_output,int* tet_output)
{
    int s=segment_mesh.Segment(node1,node2);if(segment_output) *segment_output=s;if(!s) return false;
    
    // now find a tet that contains this segment!
    int level,tet;
    for(int j=1;j<=meshes.m;j++) if(node1 <= meshes(j)->number_nodes && node2 <= meshes(j)->number_nodes){
        ARRAY<int>& incident_tets=(*meshes(j)->incident_elements)(node1);
        for(int k=1;k<=incident_tets.m;k++){
            int a,b,c,d;meshes(j)->elements(incident_tets(k)).Get(a,b,c,d);
            if(a == node2 || b == node2 || c == node2 || d == node2){level=j;tet=incident_tets(k);goto found_tet;}}}
    PHYSBAM_FATAL_ERROR();
    found_tet:
    if(level_output) *level_output=level;if(tet_output) *tet_output=tet;

    if(!Red(level,tet)){
        int parent_index=(*parent(level))(tet);int a,b,c,d;meshes(level-1)->elements(parent_index).Get(a,b,c,d);
        ARRAY<VECTOR<int,6> >& edges=*meshes(level-1)->element_edges;
        if((a == node1 && (segment_midpoints(edges(parent_index)(2)) == node2 || segment_midpoints(edges(parent_index)(5)) == node2 || 
                segment_midpoints(edges(parent_index)(6)) == node2)) 
            || (b == node1 && (segment_midpoints(edges(parent_index)(3)) == node2 || segment_midpoints(edges(parent_index)(4)) == node2 || 
                segment_midpoints(edges(parent_index)(6)) == node2)) 
            || (c == node1 && (segment_midpoints(edges(parent_index)(1)) == node2 || segment_midpoints(edges(parent_index)(4)) == node2 || 
                segment_midpoints(edges(parent_index)(5)) == node2))
            || (d == node1 && (segment_midpoints(edges(parent_index)(1)) == node2 || segment_midpoints(edges(parent_index)(2)) == node2 || 
                segment_midpoints(edges(parent_index)(3)) == node2)) 
            || (a == node2 && (segment_midpoints(edges(parent_index)(2)) == node1 || segment_midpoints(edges(parent_index)(5)) == node1 || 
                segment_midpoints(edges(parent_index)(6)) == node1)) 
            || (b == node2 && (segment_midpoints(edges(parent_index)(3)) == node1 || segment_midpoints(edges(parent_index)(4)) == node1 || 
                segment_midpoints(edges(parent_index)(6)) == node1))
            || (c == node2 && (segment_midpoints(edges(parent_index)(1)) == node1 || segment_midpoints(edges(parent_index)(4)) == node1 || 
                segment_midpoints(edges(parent_index)(5)) == node1)) 
            || (d == node2 && (segment_midpoints(edges(parent_index)(1)) == node1 || segment_midpoints(edges(parent_index)(2)) == node1 || 
                segment_midpoints(edges(parent_index)(3)) == node1)) 
            || (a != node1 && a != node2 && b != node1 && b != node2 && c != node1 && c != node2 && d != node1 && d != node2))
            return true;} // for Green or possibly Green

    return false;
}
//#####################################################################
// Function Resolve_Stack
//#####################################################################
template<class T> void RED_GREEN_TETRAHEDRA<T>::
Resolve_Stack()
{
    while(stack.m){
        int level,tet;stack.Pop().Get(level,tet);
        if(level){PHYSBAM_ASSERT(tet);(*index_in_stack(level))(tet)=0;if(!Regularly_Refined(level,tet)) Refine_If_Necessary(level,tet);}}
}
//#####################################################################
// Function Refine_If_Necessary
//#####################################################################
template<class T> void RED_GREEN_TETRAHEDRA<T>::
Refine_If_Necessary(const int level,const int tet)
{
    // first check if we're already consistent, and don't need to refine at all.
    ARRAY<int> midpoints(6),subedges(24);Get_Existing_Subindices(level,tet,midpoints,subedges);
    int number_midpoints=(midpoints(1)!=0)+(midpoints(2)!=0)+(midpoints(3)!=0)+(midpoints(4)!=0)+(midpoints(5)!=0)+(midpoints(6)!=0);
    if(number_midpoints == 0) return;
    int number_children=0;while(number_children < 8 && (*children(level))(tet)(number_children+1)) number_children++;
    switch(number_midpoints){
        case 1:if(number_children == 2) return;break;
        case 2:if(number_children == 4) return;break;
        case 3:if(number_children == 4 && !(midpoints(1) && midpoints(6)) && !(midpoints(2) && midpoints(4)) && 
                       !(midpoints(3) && midpoints(5))) return;break;
        case 6:if(number_children==8) return;break;
        default: ;}

    if(!Red(level,tet)){ // Green case - most likely refine parent unless we get lucky and can upgrade Green from 2 to 4
        int p=(*parent(level))(tet);
        int number_children=0;while(number_children < 8 && (*children(level-1))(p)(number_children+1)) number_children++;
        if(number_children < 4){ // try to update Green from 2 to 4
            Refine_If_Necessary(level-1,p); // re-refine the parent
            int new_number_children=0;while(new_number_children < 8 && (*children(level-1))(p)(new_number_children+1)) new_number_children++;
            if(new_number_children == number_children) Regularly_Refine_Tet(level-1,p);} // irregular refinement doesn't do anything, try regular
        else Regularly_Refine_Tet(level-1,p); // refine parent
        return;}
 
    if(number_midpoints > 3){Regularly_Refine_Tet(level,tet);return;} // too many to go green

    // add 3rd midpoint to faces with two midpoints
    if(number_midpoints > 1){
        if(!midpoints(1) && ((midpoints(2) && midpoints(3)) || (midpoints(4) && midpoints(5)))){
            midpoints(1)=Add_Midpoint((*meshes(level)->element_edges)(tet)(1),level,tet);number_midpoints++;}
        else if(!midpoints(2) && ((midpoints(3) && midpoints(1)) || (midpoints(5) && midpoints(6)))){
            midpoints(2)=Add_Midpoint((*meshes(level)->element_edges)(tet)(2),level,tet);number_midpoints++;}
        else if(!midpoints(3) && ((midpoints(1) && midpoints(2)) || (midpoints(4) && midpoints(6)))){
            midpoints(3)=Add_Midpoint((*meshes(level)->element_edges)(tet)(3),level,tet);number_midpoints++;}
        else if(!midpoints(4) && ((midpoints(5) && midpoints(1)) || (midpoints(6) && midpoints(3)))){
            midpoints(4)=Add_Midpoint((*meshes(level)->element_edges)(tet)(4),level,tet);number_midpoints++;}
        else if(!midpoints(5) && ((midpoints(1) && midpoints(4)) || (midpoints(6) && midpoints(2)))){
            midpoints(5)=Add_Midpoint((*meshes(level)->element_edges)(tet)(5),level,tet);number_midpoints++;}
        else if(!midpoints(6) && ((midpoints(3) && midpoints(4)) || (midpoints(2) && midpoints(5)))){
            midpoints(6)=Add_Midpoint((*meshes(level)->element_edges)(tet)(6),level,tet);number_midpoints++;}
        // check if we created any more midpoints to make us go red
        if(number_midpoints > 3){Regularly_Refine_Tet(level,tet);return;}}

    // Green refinement
    Ensure_Level_Exists(level+1);
    ARRAY<int> free_tet_indices,free_edge_indices;
    free_tet_indices.Preallocate(8);free_edge_indices.Preallocate(24);
    Delete_Children(level,tet,free_tet_indices,free_edge_indices); // get rid of what's here already
    // figure out which green tets and green or interior segments need to be made, and make them
    int i,j,k,l;meshes(level)->elements(tet).Get(i,j,k,l);
    // make sure the subedges along tet's edges exist, and that green bisectors exist
    if(midpoints(1)){ // if edge ij is split
        if(!subedges(1)) Add_Segment(free_edge_indices,i,midpoints(1));if(!subedges(4)) Add_Segment(free_edge_indices,j,midpoints(1));
        if(!midpoints(2) && !midpoints(3) && !segment_mesh.Segment(midpoints(1),k)) Add_Segment(free_edge_indices,midpoints(1),k);
        if(!midpoints(4) && !midpoints(5) && !segment_mesh.Segment(midpoints(1),l)) Add_Segment(free_edge_indices,midpoints(1),l);}
    if(midpoints(2)){ // if edge jk is split
        if(!subedges(5)) Add_Segment(free_edge_indices,j,midpoints(2));if(!subedges(7)) Add_Segment(free_edge_indices,k,midpoints(2));
        if(!midpoints(1) && !midpoints(3) && !segment_mesh.Segment(midpoints(2),i)) Add_Segment(free_edge_indices,midpoints(2),i);
        if(!midpoints(5) && !midpoints(6) && !segment_mesh.Segment(midpoints(2),l)) Add_Segment(free_edge_indices,midpoints(2),l);}
    if(midpoints(3)){ // if edge ki is split
        if(!subedges(2)) Add_Segment(free_edge_indices,i,midpoints(3));if(!subedges(8)) Add_Segment(free_edge_indices,k,midpoints(3));
        if(!midpoints(1) && !midpoints(2) && !segment_mesh.Segment(midpoints(3),j)) Add_Segment(free_edge_indices,midpoints(3),j);
        if(!midpoints(4) && !midpoints(6) && !segment_mesh.Segment(midpoints(3),l)) Add_Segment(free_edge_indices,midpoints(3),l);}
    if(midpoints(4)){ // if edge il is split
        if(!subedges(3)) Add_Segment(free_edge_indices,i,midpoints(4));if(!subedges(10)) Add_Segment(free_edge_indices,l,midpoints(4));
        if(!midpoints(1) && !midpoints(5) && !segment_mesh.Segment(midpoints(4),j)) Add_Segment(free_edge_indices,midpoints(4),j);
        if(!midpoints(3) && !midpoints(6) && !segment_mesh.Segment(midpoints(4),k)) Add_Segment(free_edge_indices,midpoints(4),k);}
    if(midpoints(5)){ // if edge jl is split
        if(!subedges(6)) Add_Segment(free_edge_indices,j,midpoints(5));if(!subedges(11)) Add_Segment(free_edge_indices,l,midpoints(5));
        if(!midpoints(1) && !midpoints(4) && !segment_mesh.Segment(midpoints(5),i)) Add_Segment(free_edge_indices,midpoints(5),i);
        if(!midpoints(2) && !midpoints(6) && !segment_mesh.Segment(midpoints(5),k)) Add_Segment(free_edge_indices,midpoints(5),k);}
    if(midpoints(6)){ // if edge kl is split
        if(!subedges(9)) Add_Segment(free_edge_indices,k,midpoints(6));if(!subedges(12)) Add_Segment(free_edge_indices,l,midpoints(6));
        if(!midpoints(3) && !midpoints(4) && !segment_mesh.Segment(midpoints(6),i)) Add_Segment(free_edge_indices,midpoints(6),i);
        if(!midpoints(2) && !midpoints(5) && !segment_mesh.Segment(midpoints(6),j)) Add_Segment(free_edge_indices,midpoints(6),j);}
    // go ahead and create any missing edges along with the tets
    if(midpoints(1) && midpoints(2) && midpoints(3)){ // if face ijk is subdivided but the other faces aren't
        if(!subedges(13)) Add_Segment(free_edge_indices,midpoints(1),midpoints(2));if(!subedges(14)) Add_Segment(free_edge_indices,midpoints(2),midpoints(3));
        if(!subedges(15)) Add_Segment(free_edge_indices,midpoints(3),midpoints(1));
        Add_Tetrahedron(free_tet_indices,level+1,i,midpoints(1),midpoints(3),l,tet);Add_Tetrahedron(free_tet_indices,level+1,j,midpoints(2),midpoints(1),l,tet);
        Add_Tetrahedron(free_tet_indices,level+1,k,midpoints(3),midpoints(2),l,tet);Add_Tetrahedron(free_tet_indices,level+1,midpoints(1),midpoints(2),midpoints(3),l,tet);}
    else if(midpoints(1) && midpoints(4) && midpoints(5)){ // if face ilj is subdivided but the other faces aren't
        if(!subedges(16)) Add_Segment(free_edge_indices,midpoints(1),midpoints(5));if(!subedges(17)) Add_Segment(free_edge_indices,midpoints(5),midpoints(4));
        if(!subedges(18)) Add_Segment(free_edge_indices,midpoints(4),midpoints(1));
        Add_Tetrahedron(free_tet_indices,level+1,i,midpoints(4),midpoints(1),k,tet);Add_Tetrahedron(free_tet_indices,level+1,l,midpoints(5),midpoints(4),k,tet);
        Add_Tetrahedron(free_tet_indices,level+1,j,midpoints(1),midpoints(5),k,tet);Add_Tetrahedron(free_tet_indices,level+1,midpoints(1),midpoints(4),midpoints(5),k,tet);}
    else if(midpoints(3) && midpoints(4) && midpoints(6)){ // if face ikl is subdivided but the other faces aren't
        if(!subedges(19)) Add_Segment(free_edge_indices,midpoints(3),midpoints(6));if(!subedges(20)) Add_Segment(free_edge_indices,midpoints(6),midpoints(4));
        if(!subedges(21)) Add_Segment(free_edge_indices,midpoints(4),midpoints(3));
        Add_Tetrahedron(free_tet_indices,level+1,i,midpoints(3),midpoints(4),j,tet);Add_Tetrahedron(free_tet_indices,level+1,midpoints(4),midpoints(3),midpoints(6),j,tet);
        Add_Tetrahedron(free_tet_indices,level+1,midpoints(6),midpoints(3),k,j,tet);Add_Tetrahedron(free_tet_indices,level+1,l,midpoints(4),midpoints(6),j,tet);}
    else if(midpoints(2) && midpoints(5) && midpoints(6)){ // if face jlk is subdivided but the other faces aren't
        if(!subedges(22)) Add_Segment(free_edge_indices,midpoints(2),midpoints(6));if(!subedges(23)) Add_Segment(free_edge_indices,midpoints(6),midpoints(5));
        if(!subedges(24)) Add_Segment(free_edge_indices,midpoints(5),midpoints(2));
        Add_Tetrahedron(free_tet_indices,level+1,midpoints(2),j,midpoints(5),i,tet);Add_Tetrahedron(free_tet_indices,level+1,midpoints(2),midpoints(6),k,i,tet);
        Add_Tetrahedron(free_tet_indices,level+1,midpoints(2),midpoints(5),midpoints(6),i,tet);Add_Tetrahedron(free_tet_indices,level+1,midpoints(6),midpoints(5),l,i,tet);}
    else if(midpoints(1) && midpoints(6)){ // if edges ij and kl are split
        Add_Segment(free_edge_indices,midpoints(1),midpoints(6));
        Add_Tetrahedron(free_tet_indices,level+1,i,midpoints(1),k,midpoints(6),tet);Add_Tetrahedron(free_tet_indices,level+1,k,midpoints(1),j,midpoints(6),tet);
        Add_Tetrahedron(free_tet_indices,level+1,i,l,midpoints(1),midpoints(6),tet);Add_Tetrahedron(free_tet_indices,level+1,midpoints(1),l,j,midpoints(6),tet);}
    else if(midpoints(2) && midpoints(4)){ // if edges jk and il are split
        Add_Segment(free_edge_indices,midpoints(2),midpoints(4));
        Add_Tetrahedron(free_tet_indices,level+1,i,j,midpoints(2),midpoints(4),tet);Add_Tetrahedron(free_tet_indices,level+1,i,midpoints(2),k,midpoints(4),tet);
        Add_Tetrahedron(free_tet_indices,level+1,j,l,midpoints(2),midpoints(4),tet);Add_Tetrahedron(free_tet_indices,level+1,midpoints(2),l,k,midpoints(4),tet);}
    else if(midpoints(3) && midpoints(5)){ // if edges ki and jl are split
        Add_Segment(free_edge_indices,midpoints(3),midpoints(5));
        Add_Tetrahedron(free_tet_indices,level+1,j,i,midpoints(5),midpoints(3),tet);Add_Tetrahedron(free_tet_indices,level+1,midpoints(5),i,l,midpoints(3),tet);
        Add_Tetrahedron(free_tet_indices,level+1,k,j,midpoints(5),midpoints(3),tet);Add_Tetrahedron(free_tet_indices,level+1,k,midpoints(5),l,midpoints(3),tet);}
    else if(midpoints(1)){ // if just edge ij is split
        Add_Tetrahedron(free_tet_indices,level+1,i,midpoints(1),k,l,tet);Add_Tetrahedron(free_tet_indices,level+1,midpoints(1),j,k,l,tet);}
    else if(midpoints(2)){ // if just edge jk is split
        Add_Tetrahedron(free_tet_indices,level+1,i,j,midpoints(2),l,tet);Add_Tetrahedron(free_tet_indices,level+1,i,midpoints(2),k,l,tet);}
    else if(midpoints(3)){ // if just edge ki is split
        Add_Tetrahedron(free_tet_indices,level+1,j,k,midpoints(3),l,tet);Add_Tetrahedron(free_tet_indices,level+1,j,midpoints(3),i,l,tet);}
    else if(midpoints(4)){ // if just edge il is split
        Add_Tetrahedron(free_tet_indices,level+1,j,i,midpoints(4),k,tet);Add_Tetrahedron(free_tet_indices,level+1,j,midpoints(4),l,k,tet);}
    else if(midpoints(5)){ // if just edge jl is split
        Add_Tetrahedron(free_tet_indices,level+1,i,midpoints(5),j,k,tet);Add_Tetrahedron(free_tet_indices,level+1,i,l,midpoints(5),k,tet);}
    else if(midpoints(6)){ // if just edge kl is split
        Add_Tetrahedron(free_tet_indices,level+1,i,k,midpoints(6),j,tet);Add_Tetrahedron(free_tet_indices,level+1,i,midpoints(6),l,j,tet);}
}
//#####################################################################
// Function Regularly_Refine_Tet
//#####################################################################
template<class T> void RED_GREEN_TETRAHEDRA<T>::
Regularly_Refine_Tet(const int level,const int tet)
{
    PHYSBAM_ASSERT(Red(level,tet));
    (*index_in_stack(level))(tet)=-1; // mark tet so we won't consider this tet for resolution in stack
    if(Regularly_Refined(level,tet)) return;
    Ensure_Level_Exists(level+1);
    ARRAY<int> free_tet_indices,free_edge_indices;
    free_tet_indices.Preallocate(8);free_edge_indices.Preallocate(24);
    Delete_Children(level,tet,free_tet_indices,free_edge_indices);
    ARRAY<int> midpoints(6),subedges(24);Get_Existing_Subindices(level,tet,midpoints,subedges);
    // there are three pairs of valid midpoints for the interior edge.  if possible, we would like to pick a pair that exists already.
    int first_midpoint_index,second_midpoint_index;
    ARRAY<VECTOR<int,6> >& element_edges=*meshes(level)->element_edges;
    if(refinement_target_segments.Contains(element_edges(tet)(3)) && refinement_target_segments.Contains(element_edges(tet)(5))){first_midpoint_index=3;second_midpoint_index=5;}
    else if(refinement_target_segments.Contains(element_edges(tet)(2)) && refinement_target_segments.Contains(element_edges(tet)(4))){first_midpoint_index=2;second_midpoint_index=4;}
    else{first_midpoint_index=1;second_midpoint_index=6;}

    for(int a=1;a<=6;a++) if(!midpoints(a)) midpoints(a)=Add_Midpoint((*meshes(level)->element_edges)(tet)(a),level,tet);
    int i,j,k,l;meshes(level)->elements(tet).Get(i,j,k,l);
    int ij=midpoints(1),jk=midpoints(2),ki=midpoints(3),il=midpoints(4),jl=midpoints(5),kl=midpoints(6);
    if(!subedges(1)) Add_Segment(free_edge_indices,i,ij);if(!subedges(2)) Add_Segment(free_edge_indices,i,ki);
    if(!subedges(3)) Add_Segment(free_edge_indices,i,il);if(!subedges(4)) Add_Segment(free_edge_indices,j,ij);
    if(!subedges(5)) Add_Segment(free_edge_indices,j,jk);if(!subedges(6)) Add_Segment(free_edge_indices,j,jl);
    if(!subedges(7)) Add_Segment(free_edge_indices,k,jk);if(!subedges(8)) Add_Segment(free_edge_indices,k,ki);
    if(!subedges(9)) Add_Segment(free_edge_indices,k,kl);if(!subedges(10)) Add_Segment(free_edge_indices,l,il);
    if(!subedges(11)) Add_Segment(free_edge_indices,l,jl);if(!subedges(12)) Add_Segment(free_edge_indices,l,kl);
    if(!subedges(13)) Add_Segment(free_edge_indices,ij,jk);if(!subedges(14)) Add_Segment(free_edge_indices,jk,ki);
    if(!subedges(15)) Add_Segment(free_edge_indices,ki,ij);if(!subedges(16)) Add_Segment(free_edge_indices,ij,jl);
    if(!subedges(17)) Add_Segment(free_edge_indices,jl,il);if(!subedges(18)) Add_Segment(free_edge_indices,il,ij);
    if(!subedges(19)) Add_Segment(free_edge_indices,ki,kl);if(!subedges(20)) Add_Segment(free_edge_indices,kl,il);
    if(!subedges(21)) Add_Segment(free_edge_indices,il,ki);if(!subedges(22)) Add_Segment(free_edge_indices,jk,kl);
    if(!subedges(23)) Add_Segment(free_edge_indices,kl,jl);if(!subedges(24)) Add_Segment(free_edge_indices,jl,jk);
    Add_Segment(free_edge_indices,midpoints(first_midpoint_index),midpoints(second_midpoint_index)); // interior edge always needs to be added
    Add_Tetrahedron(free_tet_indices,level+1,i,ij,ki,il,tet);Add_Tetrahedron(free_tet_indices,level+1,ij,j,jk,jl,tet);
    Add_Tetrahedron(free_tet_indices,level+1,ki,jk,k,kl,tet);Add_Tetrahedron(free_tet_indices,level+1,il,jl,kl,l,tet);

    if(first_midpoint_index==1){
        Add_Tetrahedron(free_tet_indices,level+1,kl,ij,il,ki,tet);Add_Tetrahedron(free_tet_indices,level+1,kl,ij,jk,jl,tet);
        Add_Tetrahedron(free_tet_indices,level+1,kl,ij,ki,jk,tet);Add_Tetrahedron(free_tet_indices,level+1,kl,ij,jl,il,tet);}
    else if(first_midpoint_index==2){
        Add_Tetrahedron(free_tet_indices,level+1,jk,il,ij,jl,tet);Add_Tetrahedron(free_tet_indices,level+1,jk,il,ki,ij,tet);
        Add_Tetrahedron(free_tet_indices,level+1,jk,il,kl,ki,tet);Add_Tetrahedron(free_tet_indices,level+1,jk,il,jl,kl,tet);}
    else{
        Add_Tetrahedron(free_tet_indices,level+1,jk,ki,jl,kl,tet);Add_Tetrahedron(free_tet_indices,level+1,ij,ki,jl,jk,tet);
        Add_Tetrahedron(free_tet_indices,level+1,il,ki,jl,ij,tet);Add_Tetrahedron(free_tet_indices,level+1,kl,ki,jl,il,tet);}
}
//#####################################################################
// Function Get_Existing_Subindices
//#####################################################################
// midpoints (or 0 if no midpoint) listed in the edge order: ij,jk,ki,il,jl,kl (midpoints has m==6)
// face-subedges listed in the order: i-ij,i-ki,i-il, j-ij,j-jk,j-jl, k-jk,k-ki,k-kl, l-il,l-jl,l-kl, ij-jk,jk-ki,ki-ij, ij-jl,jl-il,il-ij, ki-kl,kl-il,il-ki, jk-kl,kl-jl,jl-jk
template<class T> void RED_GREEN_TETRAHEDRA<T>::
Get_Existing_Subindices(const int level,const int tet,ARRAY<int>& midpoints,ARRAY<int>& subedges)
{
    PHYSBAM_ASSERT(midpoints.m == 6 && subedges.m == 24);
    for(int e=1;e<=6;e++) midpoints(e)=segment_midpoints((*meshes(level)->element_edges)(tet)(e));
    int i,j,k,l;meshes(level)->elements(tet).Get(i,j,k,l);
    int ij=midpoints(1),jk=midpoints(2),ki=midpoints(3),il=midpoints(4),jl=midpoints(5),kl=midpoints(6);
    if(ij){
        subedges(1)=segment_mesh.Segment(i,ij);subedges(4)=segment_mesh.Segment(j,ij);
        if(jk) subedges(13)=segment_mesh.Segment(ij,jk);if(ki) subedges(15)=segment_mesh.Segment(ki,ij);
        if(il) subedges(18)=segment_mesh.Segment(il,ij);if(jl) subedges(16)=segment_mesh.Segment(ij,jl);}
    if(jk){
        subedges(5)=segment_mesh.Segment(j,jk);subedges(7)=segment_mesh.Segment(k,jk);
        if(ki) subedges(14)=segment_mesh.Segment(jk,ki);if(jl) subedges(24)=segment_mesh.Segment(jl,jk);if(kl) subedges(22)=segment_mesh.Segment(jk,kl);}
    if(ki){
        subedges(2)=segment_mesh.Segment(i,ki);subedges(8)=segment_mesh.Segment(k,ki);
        if(il) subedges(21)=segment_mesh.Segment(il,ki);if(kl) subedges(19)=segment_mesh.Segment(ki,kl);}
    if(il){
        subedges(3)=segment_mesh.Segment(i,il);subedges(10)=segment_mesh.Segment(l,il);
        if(jl) subedges(17)=segment_mesh.Segment(jl,il);if(kl) subedges(20)=segment_mesh.Segment(kl,il);}
    if(jl){
        subedges(6)=segment_mesh.Segment(j,jl);subedges(11)=segment_mesh.Segment(l,jl);
        if(kl) subedges(23)=segment_mesh.Segment(kl,jl);}
    if(kl){subedges(9)=segment_mesh.Segment(k,kl);subedges(12)=segment_mesh.Segment(l,kl);}    
}
//#####################################################################
// Function Ensure_Level_Exists
//#####################################################################
template<class T> void RED_GREEN_TETRAHEDRA<T>::
Ensure_Level_Exists(const int level)
{
    if(meshes.m >= level) return; // it already exists
    meshes.Append(new TETRAHEDRON_MESH);
    meshes(level)->number_nodes=object.particles.array_collection->Size();
    meshes(level)->incident_elements=new ARRAY<ARRAY<int> >(object.particles.array_collection->Size());
    meshes(level)->element_edges=new ARRAY<VECTOR<int,6> >();
    parent.Append(new ARRAY<int>(0));children.Append(new ARRAY<VECTOR<int,8> >());
    index_in_stack.Append(new ARRAY<int>(0));leaf_number.Append(new ARRAY<int>(0));
}
//#####################################################################
// Function Delete_Children
//#####################################################################
// This deletes children tets (should be green) and any green or interior edges belonging *only* to those tets.
// Note that it doesn't touch exterior subedges that fit in red tetrahedra.
template<class T> void RED_GREEN_TETRAHEDRA<T>::
Delete_Children(const int level,const int tet,ARRAY<int>& deleted_tet_indices,ARRAY<int>& deleted_edge_indices)
{
    if(Leaf(level,tet)) return;
    // basic information about tet
    int i,j,k,l;meshes(level)->elements(tet).Get(i,j,k,l);
    int ij=segment_midpoints((*meshes(level)->element_edges)(tet)(1)),jk=segment_midpoints((*meshes(level)->element_edges)(tet)(2)),
         ki=segment_midpoints((*meshes(level)->element_edges)(tet)(3)),il=segment_midpoints((*meshes(level)->element_edges)(tet)(4)),
         jl=segment_midpoints((*meshes(level)->element_edges)(tet)(5)),kl=segment_midpoints((*meshes(level)->element_edges)(tet)(6));
    // make the list of deleted tets (namely, the children) and zero out children
    int p;for(p=1;p<=8&&(*children(level))(tet)(p);p++){deleted_tet_indices.Append((*children(level))(tet)(p));(*children(level))(tet)(p)=0;}
    // get a list of edges to delete (begin by finding all children edges, then filter the red ones out)
    ARRAY<int> children_edges;children_edges.Preallocate(25);
    for(p=1;p<=deleted_tet_indices.m;p++) for(int q=1;q<=6;q++) // get a list of all children edges
        children_edges.Append_Unique((*meshes(level+1)->element_edges)(deleted_tet_indices(p))(q));
    for(p=1;p<=children_edges.m;p++){
        int q,r;segment_mesh.elements(children_edges(p)).Get(q,r);
        if((q==ij&&r==kl) || (q==jk&&r==il) || (q==ki&&r==jl) || (q==il&&r==jk) || (q==jl&&r==ki) || (q==kl&&r==ij)) 
            deleted_edge_indices.Append(children_edges(p)); // delete interior edges
        else if((q==i && (r==jk || r==jl || r==kl)) || (q==j && (r==ki || r==il || r==kl)) || (q==k && (r==ij || r==il || r==jl)) || (q==l && (r==ij || r==jk || r==ki)) || 
            (r==i && (q==jk || q==jl || q==kl)) || (r==j && (q==ki || q==il || q==kl)) || (r==k && (q==ij || q==il || q==jl)) || (r==l && (q==ij || q==jk || q==ki))){ // delete face bisectors
            for(int s=1;s<=(*meshes(level+1)->incident_elements)(q).m;s++){ // check incident tets for this edge
                int t=(*meshes(level+1)->incident_elements)(q)(s);
                if((*parent(level+1))(t)!=tet){ // only look at tets with a different parent (i.e. not children that should be deleted)
                    int a,b,c,d;meshes(level+1)->elements(t).Get(a,b,c,d);
                    if(a==r||b==r||c==r||d==r) goto edge_owned_by_other_tet_so_break_out_of_loop_without_appending;}}
            deleted_edge_indices.Append(children_edges(p));
            edge_owned_by_other_tet_so_break_out_of_loop_without_appending:;}}
    // now do the actual deletion
    for(p=1;p<=deleted_tet_indices.m;p++) for(int q=1;q<=4;q++){ // remove deleted tets from incident_tets
        int node=meshes(level+1)->elements(deleted_tet_indices(p))(q),index=0;
        (*meshes(level+1)->incident_elements)(node).Find(deleted_tet_indices(p),index);PHYSBAM_ASSERT(index);
        (*meshes(level+1)->incident_elements)(node).Remove_Index_Lazy(index);}
    for(p=1;p<=deleted_edge_indices.m;p++) for(int q=1;q<=2;q++){ // remove deleted edges from incident_segments
        int node=segment_mesh.elements(deleted_edge_indices(p))(q),index=0;
        (*segment_mesh.incident_elements)(node).Find(deleted_edge_indices(p),index);PHYSBAM_ASSERT(index);
        (*segment_mesh.incident_elements)(node).Remove_Index_Lazy(index);}
    // then zero out occurances in the stack
    for(p=1;p<=deleted_tet_indices.m;p++){
        int t=deleted_tet_indices(p);
        if((*index_in_stack(level+1))(t)){
            stack((*index_in_stack(level+1))(t)).Set(0,0); // can't remove since it would screw up index_in_stack, mark as no longer relevent
            (*index_in_stack(level+1))(t)=0;}}
}
//#####################################################################
// Function Add_Midpoint
//#####################################################################
template<class T> int RED_GREEN_TETRAHEDRA<T>::
Add_Midpoint(const int segment,const int level,const int tet)
{
    PHYSBAM_ASSERT(segment_midpoints(segment)==0);
    int node1,node2;segment_mesh.elements(segment).Get(node1,node2);int new_node=object.particles.array_collection->Add_Element();
    // define midpoint's position, velocity and acceleration
    object.particles.X(new_node)=(T).5*(object.particles.X(node1)+object.particles.X(node2));
    if(object.particles.store_velocity)
        object.particles.V(new_node)=(T).5*(object.particles.V(node1)+object.particles.V(node2));
    if(rest_position) (*rest_position)(new_node)=(T).5*((*rest_position)(node1)+(*rest_position)(node2));
    // list it
    segment_midpoints(segment)=new_node;
    // add edge neighbors of tet to the stack if not already there
    int i;for(i=1;i<=(*meshes(level)->incident_elements)(node1).m;i++){
        int s=(*meshes(level)->incident_elements)(node1)(i);
        if(s != tet && !(*index_in_stack(level))(s)){
            int a,b,c,d;meshes(level)->elements(s).Get(a,b,c,d);assert(a == node1 || b == node1 || c == node1 || d == node1);
            if(node2 == a || node2 == b || node2 == c || node2 == d){ // not incident on edge
                stack.Append(VECTOR<int,2>(level,s));(*index_in_stack(level))(s)=stack.m;}}}
    // check previous level for tets incident on the edge
    if(level > 1 && node1 <= meshes(level-1)->number_nodes) for(i=1;i<=(*meshes(level-1)->incident_elements)(node1).m;i++){
        int s=(*meshes(level-1)->incident_elements)(node1)(i);
        if(!(*index_in_stack(level-1))(s)){
            int a,b,c,d;meshes(level-1)->elements(s).Get(a,b,c,d);
            if(node2 == a || node2 == b || node2 == c || node2 == d){ // not incident on edge
                stack.Append(VECTOR<int,2>(level-1,s));(*index_in_stack(level-1))(s)=stack.m;}}}
    // check next level for tets incident on the edge
    if(level < meshes.m && node1 <= meshes(level+1)->number_nodes) for(i=1;i<=(*meshes(level+1)->incident_elements)(node1).m;i++){
        int s=(*meshes(level+1)->incident_elements)(node1)(i);
        if(!(*index_in_stack(level+1))(s)){
            int a,b,c,d;meshes(level+1)->elements(s).Get(a,b,c,d);assert(a == node1 || b == node1 || c == node1 || d == node1);
            if(node2 == a || node2 == b || node2 == c || node2 == d){ // not incident on edge
                stack.Append(VECTOR<int,2>(level+1,s));(*index_in_stack(level+1))(s)=stack.m;}}}
    // finally make sure that number_nodes in the tet meshes and the segment mesh is up to date
    if(new_node > segment_mesh.number_nodes){
        segment_mesh.number_nodes=new_node;segment_mesh.incident_elements->Resize(new_node);}
    for(i=level;i<=meshes.m;i++) if(new_node > meshes(i)->number_nodes){
        meshes(i)->number_nodes=new_node;meshes(i)->incident_elements->Resize(new_node);}
    return new_node;
}
//#####################################################################
// Function Add_Segment
//#####################################################################
template<class T> void RED_GREEN_TETRAHEDRA<T>::
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
// Function Add_Tetrahedron
//#####################################################################
template<class T> void RED_GREEN_TETRAHEDRA<T>::
Add_Tetrahedron(ARRAY<int>& free_tet_indices,const int level,const int i,const int j,const int k,const int l,const int parent_index)
{
    TETRAHEDRON_MESH& tet_mesh=*meshes(level);
    int index;
    if(free_tet_indices.m == 0){
        tet_mesh.elements.Append(VECTOR<int,4>(i,j,k,l));
        leaf_number(level)->Append(0);parent(level)->Append(parent_index);index=tet_mesh.elements.m;
        tet_mesh.element_edges->Resize(index);children(level)->Resize(index);index_in_stack(level)->Resize(index);}
    else{
        index=free_tet_indices.Pop();tet_mesh.elements(index).Set(i,j,k,l);
        (*parent(level))(index)=parent_index;for(int a=1;a<=8;a++) (*children(level))(index)(a)=0;(*index_in_stack(level))(index)=0;}
    // create tet edges
    int ij=segment_mesh.Segment(i,j),jk=segment_mesh.Segment(j,k),ki=segment_mesh.Segment(k,i),il=segment_mesh.Segment(i,l),
         jl=segment_mesh.Segment(j,l),kl=segment_mesh.Segment(k,l);
    (*meshes(level)->element_edges)(index).Set(ij,jk,ki,il,jl,kl);
    // check if this new tet has a T-junction
    if(segment_midpoints(ij) || segment_midpoints(jk) || segment_midpoints(ki) || segment_midpoints(il) || segment_midpoints(jl) || 
        segment_midpoints(kl)){stack.Append(VECTOR<int,2>(level,index));(*index_in_stack(level))(index)=stack.m;}
    // add tet to incident_tets
    int a;for(a=1;a<=4;a++) (*tet_mesh.incident_elements)(tet_mesh.elements(index)(a)).Append(index);
    // add tet to parent's list of children
    for(a=1;a<=8;a++) if(!(*children(level-1))(parent_index)(a)){(*children(level-1))(parent_index)(a)=index;break;}
}
//#####################################################################
// Function Rebuild_Object
//#####################################################################
template<class T> void RED_GREEN_TETRAHEDRA<T>::
Rebuild_Object()
{
    int number_of_leaves=0;
    int level,tet;for(level=1;level<=meshes.m;level++) for(tet=1;tet<=meshes(level)->elements.m;tet++) if(Leaf(level,tet)) number_of_leaves++;
    object.mesh.elements.Resize(number_of_leaves);
    leaf_levels_and_indices.Exact_Resize(number_of_leaves);
    int index_into_tets=1;
    for(level=1;level<=meshes.m;level++) for(tet=1;tet<=meshes(level)->elements.m;tet++) 
        if(Leaf(level,tet)){
            for(int i=1;i<=4;i++) object.mesh.elements(index_into_tets)(i)=meshes(level)->elements(tet)(i);
            leaf_levels_and_indices(index_into_tets).Set(level,tet);(*leaf_number(level))(tet)=index_into_tets;index_into_tets++;}
        else (*leaf_number(level))(tet)=0;
    object.mesh.number_nodes=object.particles.array_collection->Size();
    object.Refresh_Auxiliary_Structures();
}
//#####################################################################
// Function Coarsen_Green_Refinements
//#####################################################################
// Assumes all children of any green refined tetrahedron are present in the mesh
template<class T> void RED_GREEN_TETRAHEDRA<T>::
Coarsen_Green_Refinements(TETRAHEDRON_MESH& final_mesh,ARRAY<int>& t_junctions,ARRAY<VECTOR<int,2> >& t_junction_parents) const
{
    ARRAY<bool> parent_already_coarsened(object.mesh.elements.m);
    ARRAY<VECTOR<int,2> > coarse_parents;
    final_mesh.Clean_Memory();
    for(int t=1;t<=object.mesh.elements.m;t++){
        int level,tet;leaf_levels_and_indices(t).Get(level,tet);
        if(Red(level,tet)) final_mesh.elements.Append(object.mesh.elements(t));
        else if(!parent_already_coarsened(t)){
            const int parent_tet=(*parent(level))(tet);
            final_mesh.elements.Append(meshes(level-1)->elements(parent_tet));
            coarse_parents.Append(VECTOR<int,2>(level-1,parent_tet));
            assert(!(*children(level-1))(parent_tet)(5));
            for(int i=1;i<=4;i++){
                int child=(*children(level-1))(parent_tet)(i);if(!child) break;
                parent_already_coarsened((*leaf_number(level))(child))=true;}}}
    ARRAY<bool> occurs_in_final_mesh(object.particles.array_collection->Size());final_mesh.Mark_Nodes_Referenced(occurs_in_final_mesh,true);
    ARRAY<bool> is_t_junction(object.particles.array_collection->Size());
    t_junctions.Clean_Memory();t_junction_parents.Clean_Memory();
    for(int t=1;t<=coarse_parents.m;t++){
        int level,tet;coarse_parents(t).Get(level,tet);
        for(int e=1;e<=6;e++){
            int s=(*meshes(level)->element_edges)(tet)(e),midpoint=segment_midpoints(s);
            if(midpoint && occurs_in_final_mesh(midpoint) && !is_t_junction(midpoint)){
                is_t_junction(midpoint)=true;
                t_junctions.Append(midpoint);t_junction_parents.Append(segment_mesh.elements(s));}}}
    final_mesh.number_nodes=object.mesh.number_nodes;
}
//#####################################################################
// Function Coarsen_Complete_Refinements_Of_Subset
//#####################################################################
template<class T> void RED_GREEN_TETRAHEDRA<T>::
Coarsen_Complete_Refinements_Of_Subset(TETRAHEDRON_MESH& final_mesh,ARRAY<bool>& subset,ARRAY<int>& t_junctions,ARRAY<VECTOR<int,2> >& t_junction_parents,
    bool allow_red_coarsening,ARRAY<bool>* node_is_uncoarsenable) const
{
    // Coarsen refinement steps for which all children are present in the subset
    TETRAHEDRON_MESH& mesh=object.mesh;
    ARRAY<ARRAY<bool> > tetrahedron_kept(meshes.m);
    for(int level=1;level<=meshes.m;level++){tetrahedron_kept(level).Resize(meshes(level)->elements.m);ARRAYS_COMPUTATIONS::Fill(tetrahedron_kept(level),true);}
    for(int t=1;t<=mesh.elements.m;t++) if(!subset(t)){
        int tet,level;leaf_levels_and_indices(t).Get(level,tet);
        tetrahedron_kept(level)(tet)=false;Unmark_Parents(tetrahedron_kept,level,tet);}
    if(!allow_red_coarsening) for(int level=meshes.m;level>=2;level--) for(int tet=1;tet<=meshes(level)->elements.m;tet++)
        if(tetrahedron_kept(level)(tet) && Red(level,tet)) Unmark_Parents(tetrahedron_kept,level,tet);
    // Prevent turning uncoarsenable nodes into T-junctions
    if(node_is_uncoarsenable) for(int level=meshes.m;level>=2;level--) for(int tet=1;tet<=meshes(level)->elements.m;tet++)
        if(tetrahedron_kept(level)(tet)) for(int i=1;i<=4;i++){int node=meshes(level)->elements(tet)(i);
            if((*node_is_uncoarsenable)(node) && !meshes(level-1)->elements((*parent(level))(tet)).Contains(node)) Unmark_Parents(tetrahedron_kept,level,tet);}
    final_mesh.Clean_Memory();
    for(int level=1;level<=meshes.m;level++) for(int tet=1;tet<=meshes(level)->elements.m;tet++) if(tetrahedron_kept(level)(tet)){
        Unmark_Children(tetrahedron_kept,level,tet);final_mesh.elements.Append(meshes(level)->elements(tet));}
    final_mesh.number_nodes=mesh.number_nodes;

    // T-junctions are those nodes whose topological 1-ring subtends different solid angles after coarsening
    ARRAY<bool> is_t_junction(mesh.number_nodes);
    PHYSBAM_ASSERT(object.mesh.incident_elements);final_mesh.Initialize_Incident_Elements();
    OPERATION_HASH<> hash(mesh.elements.m);
    for(int node=1;node<=mesh.number_nodes;node++){
        int adjusted_incident_elements=0;
        for(int t=1;t<=(*mesh.incident_elements)(node).m;t++){
            int incident_tet=(*mesh.incident_elements)(node)(t);
            if(!subset(incident_tet) || hash.Is_Marked_Current(incident_tet)) continue;
            int tet,level;leaf_levels_and_indices(incident_tet).Get(level,tet);
            adjusted_incident_elements++;if(Red(level,tet)) continue;
            int parent_tet=(*parent(level))(tet);
            if(!meshes(level-1)->elements(parent_tet).Contains(node)) continue;
            VECTOR<int,8>& siblings=(*children(level-1))(parent_tet);
            bool collapse_siblings=true;
            for(int i=1;siblings(i);i++) if(!subset((*leaf_number(level))(siblings(i)))){collapse_siblings=false;break;}
            if(collapse_siblings) for(int i=1;siblings(i);i++) hash.Mark((*leaf_number(level))(siblings(i)));}
        if(adjusted_incident_elements!=(*final_mesh.incident_elements)(node).m) is_t_junction(node)=true;
        PHYSBAM_ASSERT(!node_is_uncoarsenable || !is_t_junction(node) || !(*node_is_uncoarsenable)(node));
        hash.Next_Operation();}
    t_junctions.Clean_Memory();t_junction_parents.Clean_Memory();
    for(int s=1;s<=segment_midpoints.m;s++) if(segment_midpoints(s) && is_t_junction(segment_midpoints(s))){
        t_junctions.Append(segment_midpoints(s));t_junction_parents.Append(segment_mesh.elements(s));}
}
//#####################################################################
// Function Unmark_Parents
//#####################################################################
template<class T> void RED_GREEN_TETRAHEDRA<T>::
Unmark_Parents(ARRAY<ARRAY<bool> >& tetrahedron_kept,int level,int tet) const
{
    while(parent(level)){tet=(*parent(level--))(tet);tetrahedron_kept(level)(tet)=false;}
}
//#####################################################################
// Function Unmark_Children
//#####################################################################
template<class T> void RED_GREEN_TETRAHEDRA<T>::
Unmark_Children(ARRAY<ARRAY<bool> >& tetrahedron_kept,int level,int tet) const
{
    if(children(level)) for(int i=1;i<=8;i++){
        int child=(*children(level))(tet)(i);if(!child) break;
        tetrahedron_kept(level+1)(child)=false;Unmark_Children(tetrahedron_kept,level+1,child);}
}
//#####################################################################
// Function Remove_Simplex_List
//#####################################################################
template<class T> void RED_GREEN_TETRAHEDRA<T>::
Remove_Simplex_List(const ARRAY<int>& tetrahedron_list,ARRAY<HASHTABLE<int,int> >* level_simplex_maps)
{
    PHYSBAM_ASSERT(!stack.m);

    // build a list of all tetrahedrons on all levels to be deleted by going down the tree from each input tetrahedron
    ARRAY<ARRAY<int> > level_tetrahedron_list(meshes.m);
    for(int i=1;i<=tetrahedron_list.m;i++){
        int level,tet;leaf_levels_and_indices(tetrahedron_list(i)).Get(level,tet);
        do{level_tetrahedron_list(level).Append_Unique(tet);
        }while(level>1 && (tet=(*parent(level--))(tet)));}
    
    ARRAY<HASHTABLE<int,int> > default_level_simplex_maps;if(!level_simplex_maps) level_simplex_maps=&default_level_simplex_maps;
    (*level_simplex_maps).Resize(meshes.m);

    for(int level=1;level<=meshes.m;level++){
        // remove tetrahedrons to be deleted from incident elements
        // TODO: is incident_elements special or not?
        for(int i=1;i<=level_tetrahedron_list(level).m;i++) for(int j=1;j<=3;j++){int node=meshes(level)->elements(level_tetrahedron_list(level)(i))(j);
            int index=0;(*meshes(level)->incident_elements)(node).Find(level_tetrahedron_list(level)(i),index);PHYSBAM_ASSERT(index);
            (*meshes(level)->incident_elements)(node).Remove_Index_Lazy(index);}

        Sort(level_tetrahedron_list(level));

        HASHTABLE<int,int>& simplex_map=(*level_simplex_maps)(level);simplex_map.Remove_All();
        
        meshes(level)->Delete_Sorted_Elements(level_tetrahedron_list(level),simplex_map);
        if(level<parent.m) for(int i=1;i<=(*parent(level+1)).m;i++) simplex_map.Get((*parent(level+1))(i),(*parent(level+1))(i));
        if(parent(level)) parent(level)->Remove_Sorted_Indices_Lazy(level_tetrahedron_list(level));
        if(level>1) for(int i=1;i<=children(level-1)->m;i++){for(int j=1;j<=4;j++) simplex_map.Get((*children(level-1))(i)(j),(*children(level-1))(i)(j));
            // really we just need the zeros at the end...bubble it?
            VECTOR<int,number_of_red_children>& child_list=(*children(level-1))(i);
            for(int i=1;i<=number_of_red_children;i++) if(child_list(i)==0) for(int j=i+1;j<=number_of_red_children;j++) if(child_list(j)>0){exchange(child_list(j),child_list(i));break;}}
        children(level)->Remove_Sorted_Indices_Lazy(level_tetrahedron_list(level));
        for(int i=1;i<=leaf_levels_and_indices.m;i++) if(leaf_levels_and_indices(i)(1)==level) simplex_map.Get(leaf_levels_and_indices(i)(2),leaf_levels_and_indices(i)(2));
        leaf_number(level)->Remove_Sorted_Indices_Lazy(level_tetrahedron_list(level));}

    Rebuild_Object();
}
template<class T> void RED_GREEN_TETRAHEDRA<T>::
Initialize()
{
    Clean_Memory();
    meshes.Append(new TETRAHEDRON_MESH(object.mesh));
    meshes(1)->Initialize_Incident_Elements();
    meshes(1)->Initialize_Segment_Mesh(); // this is deleted below
    meshes(1)->Initialize_Element_Edges();
    parent.Append(0); // the first level can't have parents
    children.Append(new ARRAY<VECTOR<int,8> >(meshes(1)->elements.m));
    segment_mesh.number_nodes=meshes(1)->number_nodes;
    segment_mesh.elements.Exchange(meshes(1)->segment_mesh->elements);
    delete meshes(1)->segment_mesh;meshes(1)->segment_mesh=0; // no longer needed
    segment_mesh.Initialize_Incident_Elements();
    segment_midpoints.Resize(segment_mesh.elements.m);
    index_in_stack.Append(new ARRAY<int>(meshes(1)->elements.m));
    leaf_levels_and_indices.Resize(meshes(1)->elements.m);
    int t;for(t=1;t<=meshes(1)->elements.m;t++){leaf_levels_and_indices(t)(1)=1;leaf_levels_and_indices(t)(2)=t;}
    leaf_number.Append(new ARRAY<int>(meshes(1)->elements.m));
    for(t=1;t<=meshes(1)->elements.m;t++) (*leaf_number(1))(t)=t;
}
template<class T> void RED_GREEN_TETRAHEDRA<T>::
Clean_Memory()
{
    leaf_levels_and_indices.Clean_Memory();
    leaf_number.Delete_Pointers_And_Clean_Memory();
    meshes.Delete_Pointers_And_Clean_Memory();
    parent.Delete_Pointers_And_Clean_Memory();
    children.Delete_Pointers_And_Clean_Memory();
    segment_mesh.Clean_Memory();
    segment_midpoints.Clean_Memory();
    stack.Clean_Memory();
    index_in_stack.Delete_Pointers_And_Clean_Memory();
}
template<class T> void RED_GREEN_TETRAHEDRA<T>::
Initialize_Segment_Index_From_Midpoint_Index() // TODO: check that this is correct
{
    if(segment_index_from_midpoint_index) delete segment_index_from_midpoint_index;
    segment_index_from_midpoint_index=new ARRAY<int>(object.particles.array_collection->Size());
    for(int s=1;s<=segment_midpoints.m;s++) if(segment_midpoints(s)) (*segment_index_from_midpoint_index)(segment_midpoints(s))=s;
}
template<class T> void RED_GREEN_TETRAHEDRA<T>::
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
template class RED_GREEN_TETRAHEDRA<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RED_GREEN_TETRAHEDRA<double>;
#endif

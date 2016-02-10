//#####################################################################
// Copyright 2006-2007, Kevin Der, Geoffrey Irving, Andrew Selle, Eftychios Sifakis, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUTTING_TETRAHEDRALIZED_VOLUME
//#####################################################################
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_OP.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/EXACT_FLOAT.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/EXACT_RATIONAL.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/EXPANSION_ARITHMETIC.h>
#include <PhysBAM_Tools/Arrays_Computations/HEAPIFY.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/DIRECTED_GRAPH.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/clamp.h>
#include <PhysBAM_Tools/Math_Tools/exchange.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_HASHTABLE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Adaptive_Geometry/Barycentric_Coordinates.h>
#include <PhysBAM_Geometry/Adaptive_Geometry/Intersection_Coordinates.h>
#include <PhysBAM_Geometry/Adaptive_Geometry/Intersects.h>
#include <PhysBAM_Geometry/Adaptive_Geometry/Intersects_And_Intersection_Coordinates.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Fracture/CUTTING_SIMPLEX.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Read_Write/Topology/READ_WRITE_POLYGON_MESH.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/ROBUST_SIMPLEX_INTERACTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/SIMPLEX_INTERACTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Dynamics/Fracture/CUTTING_POLYGON.h>
#include <PhysBAM_Dynamics/Fracture/CUTTING_SIMPLICES.h>
#include <PhysBAM_Dynamics/Fracture/CUTTING_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Dynamics/Meshing/POLYGONAL_TRIANGULATION.h>
#include <cmath>
#include <limits>
using namespace PhysBAM;
//#####################################################################
// Function CUTTING_TETRAHEDRALIZED_VOLUME
//#####################################################################
template<class T> CUTTING_TETRAHEDRALIZED_VOLUME<T>::
CUTTING_TETRAHEDRALIZED_VOLUME(const bool verbose)
    : current_tetrahedralized_volume(0),next_tetrahedralized_volume(0),cutting_simplices(new CUTTING_SIMPLICES<T,3>),intersection_registry(new INTERSECTION_REGISTRY<T,3>(*cutting_simplices)),
      intersection_counter(0),verbose(verbose),abs_tol_for_level1_barycentric_coordinates(std::numeric_limits<T>::denorm_min()),abs_tol_for_level2_barycentric_coordinates(std::numeric_limits<T>::denorm_min())
{
    EXPANSION_ARITHMETIC<T>::Check_IEEE_Compliance();
    Set_Intersection_Thickness();
}
//#####################################################################
// Function ~CUTTING_TETRAHEDRALIZED_VOLUME
//#####################################################################
template<class T> CUTTING_TETRAHEDRALIZED_VOLUME<T>::
~CUTTING_TETRAHEDRALIZED_VOLUME()
{
    delete cutting_simplices;
    delete intersection_registry;
    delete current_tetrahedralized_volume;
    delete next_tetrahedralized_volume;
}
//#####################################################################
// Function Initialize_Original_Embedding
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Initialize_Original_Embedding(const TETRAHEDRALIZED_VOLUME<T>& original_tetrahedralized_volume)
{
    // copy original tet volume to current and initialize accel structures
    current_tetrahedralized_volume=TETRAHEDRALIZED_VOLUME<T>::Create();
    current_tetrahedralized_volume->particles.array_collection->Initialize(*original_tetrahedralized_volume.particles.array_collection);
    current_tetrahedralized_volume->mesh.Initialize_Mesh(original_tetrahedralized_volume.mesh);
    current_tetrahedralized_volume->mesh.Initialize_Incident_Elements();
    current_tetrahedralized_volume->mesh.Initialize_Triangle_Mesh();
    current_tetrahedralized_volume->mesh.triangle_mesh->Initialize_Incident_Elements();

    const ARRAY<VECTOR<int,4> >& original_embedding_tetrahedrons=current_tetrahedralized_volume->mesh.elements;
    int original_embedding_faces_count=current_tetrahedralized_volume->mesh.triangle_mesh->elements.m;
    const ARRAY<VECTOR<int,4> >& embedding_tetrahedrons=current_tetrahedralized_volume->mesh.elements;
    const TRIANGLE_MESH& embedding_faces=*current_tetrahedralized_volume->mesh.triangle_mesh;
    const GEOMETRY_PARTICLES<TV>& embedding_vertices=current_tetrahedralized_volume->particles;

    simplices_per_current_tet.Resize(original_embedding_tetrahedrons.m); regions_per_tet.Resize(original_embedding_tetrahedrons.m);
    intersection_registry->Preallocate_Intersections(embedding_vertices.array_collection->Size());
    intersection_registry->Preallocate_Simplices(original_embedding_faces_count);
    cutting_particles.Preallocate(embedding_vertices.array_collection->Size());
    cutting_polygons_per_cutting_simplex.Preallocate(original_embedding_faces_count);
    polygon_mesh.elements.Preallocate(original_embedding_faces_count);

    // cutting simplices, fill polygon mesh and generate lists of cutting polygons
    HASHTABLE<int,int> embedding_face_to_simplex;
    for(int face=1;face<=embedding_faces.elements.m;face++){const VECTOR<int,3>& face_nodes=embedding_faces.elements(face);
        int parent=cutting_simplices->Add_Simplex(face_nodes,CUTTING_SIMPLEX<T,3>::GLOBAL_EMBEDDING_FACE); // TODO: call this a boundary type???
        sorted_to_original_boundary_nodes.Insert(face_nodes.Sorted(),face_nodes);
        cutting_polygons_per_cutting_simplex.Append(ARRAY<int>());
        embedding_face_to_simplex.Insert(face,parent);
        ARRAY<int> tets_on_face;
        current_tetrahedralized_volume->mesh.Tetrahedrons_On_Face(face_nodes,tets_on_face);
        CUTTING_POLYGON::POLYGON_TYPE polygon_type=(tets_on_face.m==1?CUTTING_POLYGON::FACE_BOUNDARY:CUTTING_POLYGON::FACE_INTERIOR);
        for(int i=1;i<=tets_on_face.m;i++){int tet=tets_on_face(i);
            VECTOR<TV,3> weights=Get_Face_Weights_From_Tet(face_nodes,embedding_tetrahedrons(tet));
            int index=cutting_simplices->Add_Simplex(face_nodes,CUTTING_SIMPLEX<T,3>::LOCAL_EMBEDDING_FACE,weights,parent,tet);
            simplices_per_current_tet(tet).Append(index);
            // make polygon
            polygon_mesh.elements.Append(ARRAY<ARRAY<int> >());ARRAY<ARRAY<int> >& new_polygon=polygon_mesh.elements.Last();
            bool face_reversed=Face_Reversed_In_Tetrahedron(embedding_tetrahedrons(tet),face_nodes,embedding_vertices);
            VECTOR<int,3> face_nodes_ordered=face_nodes;
            if(face_reversed) exchange(face_nodes_ordered[2],face_nodes_ordered[3]);
            new_polygon.Append(ARRAY<int>(face_nodes_ordered));
            int cutting_polygon_index=current_cutting_polygons.Append(CUTTING_POLYGON(polygon_mesh.elements.m,index,face_reversed,polygon_type));
            cutting_polygons_per_cutting_simplex.Append(ARRAY<int>(VECTOR<int,1>(cutting_polygon_index)));
            if(!regions_per_tet(tet).m) regions_per_tet(tet).Append(ARRAY<int>());
            regions_per_tet(tet)(1).Append(current_cutting_polygons.m);}}
    // intersection registry (vertices)
    const ARRAY<ARRAY<int> >& faces_on_vertices=*current_tetrahedralized_volume->mesh.triangle_mesh->incident_elements;
    for(int node=1;node<=current_tetrahedralized_volume->mesh.number_nodes;node++){ARRAY<int> simplices;
        ARRAY<VECTOR<T,2> > all_weights(faces_on_vertices(node).m);
        for(int i=1;i<=faces_on_vertices(node).m;i++){int face=faces_on_vertices(node)(i);
            int parent;bool found=embedding_face_to_simplex.Get(face,parent);if(!found) PHYSBAM_FATAL_ERROR();
            simplices.Append(parent);
            all_weights(i)=Get_Node_Weights_From_Triangle(node,embedding_faces.elements(face));}
        intersection_registry->Register_Intersection(simplices,all_weights,++intersection_counter);
        cutting_particles.Add_Tet_Node_And_Intersection_Id(node,intersection_counter);}
}
//#####################################################################
// Function Cut_Material
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Cut_Material(const TRIANGULATED_SURFACE<T>& cutting_triangulated_surface_input,const int first_new)
{
    LOG::SCOPE scope("CUT_MATERIAL","Cut Material");

    if(next_tetrahedralized_volume) delete next_tetrahedralized_volume;
    next_tetrahedralized_volume=TETRAHEDRALIZED_VOLUME<T>::Create();
    next_tetrahedralized_volume->particles.Store_Velocity();
    dynamic_cast<PARTICLES<TV>&>(next_tetrahedralized_volume->particles).Store_Mass();
    cutting_triangulated_surface=&cutting_triangulated_surface_input;

    LOG::Time("Initializing acceleration structures");
    TRIANGLE_MESH& cutting_tri_mesh=const_cast<TRIANGLE_MESH&>(cutting_triangulated_surface->mesh);
    cutting_tri_mesh.Initialize_Element_Edges();
    cutting_tri_mesh.Initialize_Edge_Triangles();
    current_tetrahedralized_volume->mesh.Initialize_Incident_Elements();
    current_tetrahedralized_volume->mesh.Initialize_Triangle_Mesh();
    current_tetrahedralized_volume->mesh.triangle_mesh->Initialize_Incident_Elements(); // used in TRIANGLE_MESH::Initialize_Element_Edges and TETRAHEDRON_MESH::Initialize_Element_Edges
    current_tetrahedralized_volume->mesh.triangle_mesh->Initialize_Element_Edges();
    current_tetrahedralized_volume->mesh.triangle_mesh->Initialize_Edge_Triangles();
    current_tetrahedralized_volume->mesh.Initialize_Segment_Mesh_From_Triangle_Mesh(); // make TETRAHEDRON_MESH::segment_mesh and TETRAHEDRON_MESH::triangle_mesh->segment_mesh consistent
    current_tetrahedralized_volume->mesh.Initialize_Element_Edges();
    current_tetrahedralized_volume->Initialize_Hierarchy();
    current_tetrahedralized_volume->particles.Store_Velocity();
    dynamic_cast<PARTICLES<TV>&>(current_tetrahedralized_volume->particles).Store_Mass();

    LOG::Time("Finding triangle tet intersections");
    Find_Triangle_Tet_Intersections(first_new);

    LOG::Time("Finding triangle triangle intersections");
    Find_Triangle_Triangle_Intersections();

    LOG::Time("Add cuts to intersection registry");
    Fill_Intersection_Registry_With_Cuts();

    LOG::Time("Splitting existing polygon edges");
    Split_Existing_Polygon_Edges();

    LOG::Time("Splitting existing polygons");
    Split_Existing_Polygons();

    LOG::Time("Adding polygons on new simplices");
    Add_Polygons_On_New_Simplices();

    LOG::Time("Finding material regions");
    Find_Material_Regions();

    LOG::Time("Finding duplicate tets and duplicate particles");
    Determine_Duplicate_Tets_And_Duplicate_Particles();

    LOG::Time("Merging tets");
    Duplicate_And_Merge_Elements();

    LOG::Time("Creating boundary surface");
    Create_Boundary_Surface();

    LOG::Time("Rebuild state");
    Build_New_Cutting_Simplices_And_Intersection_Registry();
}
//#####################################################################
// Function Find_Triangle_Tet_Intersections
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Find_Triangle_Tet_Intersections(const int first_new)
{
    const ARRAY<VECTOR<int,4> >& current_embedding_tetrahedrons=current_tetrahedralized_volume->mesh.elements;
    const GEOMETRY_PARTICLES<TV>& current_embedding_particles=current_tetrahedralized_volume->particles;
    const ARRAY<VECTOR<int,3> >& cutting_triangles=cutting_triangulated_surface->mesh.elements;
    const ARRAY<VECTOR<int,2> >& cutting_segments=cutting_triangulated_surface->mesh.segment_mesh->elements;
    const GEOMETRY_PARTICLES<TV>& cutting_particles=cutting_triangulated_surface->particles;
    const ARRAY<ARRAY<int> >& cutting_edge_triangles=*cutting_triangulated_surface->mesh.edge_triangles;
    const TETRAHEDRON_HIERARCHY<T>& current_embedding_tetrahedron_hierarchy=*current_tetrahedralized_volume->hierarchy;

    cutting_simplices->Set_Index_For_Last_Old_Cutting_Simplex();
    last_old_simplex_index_per_current_tet.Exact_Resize(current_embedding_tetrahedrons.m);
    assert(simplices_per_current_tet.m==current_embedding_tetrahedrons.m);
    for(int tet=1;tet<=current_embedding_tetrahedrons.m;++tet)
        last_old_simplex_index_per_current_tet(tet)=simplices_per_current_tet(tet).m;
    assert(first_new>=1);
    // real newly added triangles
    ARRAY<int> tet_intersection_list;
    for(int cutting_triangle_index=first_new;cutting_triangle_index<=cutting_triangles.m;cutting_triangle_index++){
        // insert global cutting face
        const VECTOR<int,3>& cut_nodes=cutting_triangles(cutting_triangle_index);
        int parent=cutting_simplices->Add_Simplex(-cut_nodes,CUTTING_SIMPLEX<T,3>::GLOBAL_CUT_FACE);
        // find intersections
        VECTOR<TV,3> cut_X(cutting_particles.X.Subset(cut_nodes));
        tet_intersection_list.Remove_All();
        current_embedding_tetrahedron_hierarchy.Intersection_List(
            RANGE<TV>::Bounding_Box(cut_X),tet_intersection_list,intersection_thickness);
        for(int i=1;i<=tet_intersection_list.m;++i){
            int etet=tet_intersection_list(i);
            VECTOR<TV,4> embedded_tet_X(current_embedding_particles.X.Subset(current_embedding_tetrahedrons(etet)));
            bool is_degenerate;
            bool intersects=Intersects<EXACT_FLOAT<T> >(embedded_tet_X,cut_X,&is_degenerate);
            if(is_degenerate){
                LOG::cerr<<"Find_Triangle_Tet_Intersections(...) line "<<__LINE__<<":  degeneracy; cutting_triangle_index="<<cutting_triangle_index<<", etet="<<etet<<std::endl
                    <<"    cut_nodes="<<cut_nodes<<std::endl<<"    current_embedding_tetrahedrons(etet)="<<current_embedding_tetrahedrons(etet)<<std::endl
                    <<"    cut_X="<<cut_X<<std::endl<<"    embedded_tet_X="<<embedded_tet_X<<std::endl;
                PHYSBAM_FATAL_ERROR();}
            if(intersects){
                // create local copy of cutting face
                int index=cutting_simplices->Add_Simplex(-cut_nodes,CUTTING_SIMPLEX<T,3>::LOCAL_CUT_FACE,abs_tol_for_level1_barycentric_coordinates,parent,etet,embedded_tet_X,cut_X);
                simplices_per_current_tet(etet).Append(index);}}}
    // fake triangles on boundary edges of newly added triangles
    for(int cseg=1;cseg<=cutting_segments.m;++cseg) if(cutting_edge_triangles(cseg).m==1 && cutting_edge_triangles(cseg)(1)>=first_new){
        // insert global fake cutting triangle
        const VECTOR<int,3>& real_bordering_triangle_nodes=cutting_triangles(cutting_edge_triangles(cseg)(1));
        const VECTOR<int,2>& cut_nodes=cutting_segments(cseg);
        VECTOR<int,2> ordered_cut_nodes=cut_nodes;
        if(Segment_Reversed_In_Triangle(real_bordering_triangle_nodes, cut_nodes))
            exchange(ordered_cut_nodes(1),ordered_cut_nodes(2));
        VECTOR<int,3> fake_triangle_nodes=ordered_cut_nodes.Append(0);
        int parent=cutting_simplices->Add_Simplex(-fake_triangle_nodes,CUTTING_SIMPLEX<T,3>::GLOBAL_CUT_FACE);
        // find intersections
        VECTOR<TV,2> cut_X(cutting_particles.X.Subset(ordered_cut_nodes));
        tet_intersection_list.Remove_All();
        current_embedding_tetrahedron_hierarchy.Intersection_List(
            RANGE<TV>::Bounding_Box(cut_X), tet_intersection_list, intersection_thickness);
        for(int i=1;i<=tet_intersection_list.m;++i){
            int etet=tet_intersection_list(i);
            VECTOR<TV,4> embedded_tet_X(current_embedding_particles.X.Subset(current_embedding_tetrahedrons(etet)));
            bool is_degenerate;
            bool intersects=Intersects<EXACT_FLOAT<T> >(embedded_tet_X,cut_X,&is_degenerate);
            if(is_degenerate){
                LOG::cerr<<"Find_Triangle_Tet_Intersections(...) line "<<__LINE__<<":  degeneracy; cseg="<<cseg<<", etet="<<etet<<std::endl
                    <<"    ordered_cut_nodes="<<ordered_cut_nodes<<std::endl<<"    current_embedding_tetrahedrons(etet)="<<current_embedding_tetrahedrons(etet)<<std::endl
                    <<"    cut_X="<<cut_X<<std::endl<<"    embedded_tet_X="<<embedded_tet_X<<std::endl;
                PHYSBAM_FATAL_ERROR();}
            if(intersects){
                // create local copy of fake cutting tri
                int index=cutting_simplices->Add_Simplex(-fake_triangle_nodes,CUTTING_SIMPLEX<T,3>::LOCAL_CUT_FACE,abs_tol_for_level1_barycentric_coordinates,parent,etet,embedded_tet_X,cut_X.Append(TV()));
                simplices_per_current_tet(etet).Append(index);}}}
    if(verbose) cutting_simplices->Print();
}
//#####################################################################
// Function Find_Triangle_Triangle_Intersections
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Find_Triangle_Triangle_Intersections()
{
    intersecting_simplices_per_simplex.Clean_Memory();
    intersecting_simplices_per_simplex.Exact_Resize(cutting_simplices->simplices.m);
    VECTOR<int,2> one_two(1,2);
    for(int tet=1;tet<=current_tetrahedralized_volume->mesh.elements.m;++tet){
        const ARRAY<int>& simplices=simplices_per_current_tet(tet);
        // TODO: this only needs to create a box hierarchy.
        TRIANGULATED_SURFACE<T>& triangulated_surface=*TRIANGULATED_SURFACE<T>::Create();
        triangulated_surface.particles.array_collection->Add_Elements(3*simplices.m);
        triangulated_surface.mesh.elements.Preallocate(simplices.m);
        for(int i=1;i<=simplices.m;++i){
            assert(!cutting_simplices->Simplex_Is_Parent(simplices(i)));
            const VECTOR<TV,3>& triangle_weights=cutting_simplices->simplices(simplices(i)).weights;
            int p=3*(i-1);
            for(int j=1;j<=3;++j) triangulated_surface.particles.X(p+j)=triangle_weights(j);
            triangulated_surface.mesh.elements.Append(VECTOR<int,3>(p+1,p+2,p+3));}
        triangulated_surface.Update_Number_Nodes();
        triangulated_surface.Initialize_Hierarchy();
        triangulated_surface.Update_Triangle_List();
        for(int i1=1;i1<=simplices.m;++i1){
            int simplex1=simplices(i1);
            RANGE<TV> bounding_box=(*triangulated_surface.triangle_list)(i1).Bounding_Box();
            ARRAY<int> intersection_list;
            triangulated_surface.hierarchy->Intersection_List(bounding_box,intersection_list);
            if(verbose) {std::stringstream ss;ss<<"for simplex index="<<simplex1<<", simplices.Subset(intersection_list)="<<simplices.Subset(intersection_list);LOG::filecout(ss.str());}
            for(int j=1;j<=intersection_list.m;++j){
                int i2=intersection_list(j);
                int simplex2=simplices(i2);
                if(simplex1>=simplex2) continue;
                bool intersects,is_degenerate;
                ARRAY<int> shared_nodes;
                cutting_simplices->Shared_Nodes_On_Simplices(VECTOR<int,2>(simplex1,simplex2),shared_nodes);

                typedef typename CUTTING_SIMPLEX<T,3>::GET_ADAPTIVE_WEIGHTS_RESULT_TYPE WEIGHT_TYPE;
                if(shared_nodes.m>0){
                    intersects=true;is_degenerate=false;}
                else if(cutting_simplices->Simplex_Is_Fake(simplex1)){
                    if(cutting_simplices->Simplex_Is_Fake(simplex2)) continue; // assume no edge-edge intersections with no shared nodes
                    VECTOR<VECTOR<WEIGHT_TYPE,3>,2> adaptive_weights_for_simplex1;
                    VECTOR<VECTOR<WEIGHT_TYPE,3>,3> adaptive_weights_for_simplex2;
                    cutting_simplices->simplices(simplex1).Get_Adaptive_Weights(adaptive_weights_for_simplex1,one_two);
                    cutting_simplices->simplices(simplex2).Get_Adaptive_Weights(adaptive_weights_for_simplex2);
                    intersects=Intersects<void>(adaptive_weights_for_simplex1,adaptive_weights_for_simplex2,&is_degenerate);}
                else if(cutting_simplices->Simplex_Is_Fake(simplex2)){
                    VECTOR<VECTOR<WEIGHT_TYPE,3>,3> adaptive_weights_for_simplex1;
                    VECTOR<VECTOR<WEIGHT_TYPE,3>,2> adaptive_weights_for_simplex2;
                    cutting_simplices->simplices(simplex1).Get_Adaptive_Weights(adaptive_weights_for_simplex1);
                    cutting_simplices->simplices(simplex2).Get_Adaptive_Weights(adaptive_weights_for_simplex2,one_two);
                    intersects=Intersects<void>(adaptive_weights_for_simplex1,adaptive_weights_for_simplex2,&is_degenerate);}
                else{
                    VECTOR<VECTOR<WEIGHT_TYPE,3>,3> adaptive_weights_for_simplex1;
                    VECTOR<VECTOR<WEIGHT_TYPE,3>,3> adaptive_weights_for_simplex2;
                    cutting_simplices->simplices(simplex1).Get_Adaptive_Weights(adaptive_weights_for_simplex1);
                    cutting_simplices->simplices(simplex2).Get_Adaptive_Weights(adaptive_weights_for_simplex2);
                    intersects=Intersects<void>(adaptive_weights_for_simplex1,adaptive_weights_for_simplex2,&is_degenerate);}

                if(is_degenerate){
                    LOG::cerr<<"Find_Triangle_Triangle_Intersections() line "<<__LINE__<<":  degeneracy"<<"; tet="<<tet<<", simplex1="<<simplex1<<", simplex2="<<simplex2<<std::endl;
                    LOG::cerr<<"simplex1 coordinates="<<cutting_simplices->simplices(simplex1).simplex_original_coordinates<<std::endl;
                    LOG::cerr<<"simplex2 coordinates="<<cutting_simplices->simplices(simplex2).simplex_original_coordinates<<std::endl;
                    PHYSBAM_FATAL_ERROR();}
                if(intersects){intersecting_simplices_per_simplex(simplex1).Append(simplex2);intersecting_simplices_per_simplex(simplex2).Append(simplex1);}}
            if(verbose){std::stringstream ss;ss<<"    intersecting simplices: "<<intersecting_simplices_per_simplex(simplex1);LOG::filecout(ss.str());}}
        delete &triangulated_surface;}
}
//#####################################################################
// Function Fill_Intersection_Registry_With_Cuts
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Fill_Intersection_Registry_With_Cuts()
{
    intersection_counter=intersection_registry->simplices_on_intersection.m;
    intersection_registry->Resize_Simplices(cutting_simplices->simplices.m);
    intersection_registry->Set_Index_For_Last_Old_Intersection();
    const ARRAY<VECTOR<int,4> >& current_embedding_tetrahedrons=current_tetrahedralized_volume->mesh.elements;
    for(int tet=1;tet<=current_embedding_tetrahedrons.m;tet++){
        const ARRAY<int>& simplices=simplices_per_current_tet(tet);
        int last_old_simplex_index=last_old_simplex_index_per_current_tet(tet);
        for(int i=last_old_simplex_index+1;i<=simplices.m;++i){
            int new_simplex=simplices(i);
            const ARRAY<int>& intersecting_simplices=intersecting_simplices_per_simplex(new_simplex);
            for(int j1=1;j1<intersecting_simplices.m;++j1){
                int old_simplex1=intersecting_simplices(j1);
                if(old_simplex1>=new_simplex) continue; // assume simplices are ordered old first, new last
                for(int j2=j1+1;j2<=intersecting_simplices.m;++j2){
                    int old_simplex2=intersecting_simplices(j2);
                    if(old_simplex2>=new_simplex) continue; // assume simplices are ordered old first, new last
                    // TODO: ensure that old_simplex2 is in intersecting_simplices_per_simplex(old_simplex1) to improve efficiency ?
                    if(verbose) {std::stringstream ss;ss<<" Checking old simplices "<<old_simplex1<<" and "<<old_simplex2<<" with new simplex "<<new_simplex<<std::endl;LOG::filecout(ss.str());}
                    Perform_Smart_Simplex_Intersection(VECTOR<int,3>(old_simplex1,old_simplex2,new_simplex));}}}}
    if(verbose) cutting_particles.Print();
    // fix polygon mesh structures
    polygon_mesh.Set_Number_Nodes(intersection_registry->Number_Of_Intersections());
    polygon_mesh.Initialize_Segment_Mesh();
    polygon_mesh.segment_mesh->Initialize_Incident_Elements();
    polygon_mesh.Initialize_Element_Oriented_Edges();
    polygon_mesh.Initialize_Edge_Elements();
}
//#####################################################################
// Function Add_Non_Tet_Node_Intersection_To_Registry
//#####################################################################
template<class T> template<class T_ARRAY> int CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Add_Non_Tet_Node_Intersection_To_Registry(const T_ARRAY& simplices,const typename REBIND<T_ARRAY,VECTOR<T,2> >::TYPE& weights) // TODO: Rename this function to be different than next
{
    Add_Non_Tet_Node_Intersection_To_Registry(simplices,weights,++intersection_counter);
    cutting_particles.Add_Intersection_Id(intersection_counter);
    return intersection_counter;
}
//#####################################################################
// Function Add_Non_Tet_Node_Intersection_To_Registry
//#####################################################################
template<class T> template<class T_ARRAY> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Add_Non_Tet_Node_Intersection_To_Registry(const T_ARRAY& simplices,const typename REBIND<T_ARRAY,VECTOR<T,2> >::TYPE& weights,const int index)
{
    assert(simplices.Size()==weights.Size());
    intersection_registry->Register_Intersection(simplices,weights,index);
    if(verbose){std::stringstream ss;ss<<"Intersection: "<<index<<std::endl;
        for(int i=1;i<=simplices.m;i++) ss<<"simplex "<<simplices(i)<<" has weights "<<weights(i)<<std::endl;ss<<std::endl;LOG::filecout(ss.str());}
}
//#####################################################################
// Function Perform_Smart_Simplex_Intersection
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Perform_Smart_Simplex_Intersection(const VECTOR<int,3>& simplices) // simplices=(old1,old2,new)
{
    if(cutting_simplices->Simplex_Is_Fake(simplices[1]) && cutting_simplices->Simplex_Is_Fake(simplices[2]) && cutting_simplices->Simplex_Is_Fake(simplices[3])) return;
    bool converted;VECTOR<int,3> converted_simplices=cutting_simplices->Convert_Simplex_Indices_To_Global_Indices(simplices,converted);
    // already exists?
    {ARRAY<int> intersection_list;if(intersection_registry->Intersection_List(converted_simplices,intersection_list)) return;}
    int new_simplex=simplices[3];
    // 1 or more point is shared by all simplices?
    ARRAY<int> shared_nodes;
    cutting_simplices->Shared_Nodes_On_Simplices(converted_simplices,shared_nodes);
    if(shared_nodes.m>=2) return; // they are all at least on the same edge
    else if(shared_nodes.m==1){
        int shared_node=shared_nodes(1); // they share a single point
        // Intersection point must be interior to a tet, i.e. all simplices are cutting simplices
        assert(simplices==converted_simplices);
        VECTOR<VECTOR<T,2>,3> all_weights;
        for(int i=1;i<=3;i++) all_weights[i]=Get_Node_Weights_From_Triangle(shared_node,cutting_simplices->simplices(simplices[i]).nodes);
        // check if the node is inside the tet
        const CUTTING_SIMPLEX<T,3>& simplex1=cutting_simplices->simplices(simplices[1]);
        assert((simplex1.type==CUTTING_SIMPLEX<T,3>::LOCAL_CUT_FACE));
        if(!simplex1.node_in_embedded_simplex[simplex1.nodes.Find(shared_node)]) return;
        // it's inside, so add it if it doesn't already exist (if it already exists, it must be a shared node on child simplices in this tet, so we check for that first)
        ARRAY<int> simplices_containing_shared_node;
        const ARRAY<int>& intersecting_simplices=intersecting_simplices_per_simplex(new_simplex);
        for(int i=1;i<=intersecting_simplices.m;++i){
            int old_simplex=intersecting_simplices(i);
            if(old_simplex>=new_simplex) continue;
            if(cutting_simplices->simplices(old_simplex).nodes.Contains(shared_node)) simplices_containing_shared_node.Append(old_simplex);}
        // add the new simplex to the list, as it may have intersected with 2 prior old simplices
        simplices_containing_shared_node.Append(new_simplex);
        for(int i=1;i<simplices_containing_shared_node.m;++i){
            int simplex_i=simplices_containing_shared_node(i);
            for(int j=i+1;j<simplices_containing_shared_node.m;++j){
                int simplex_j=simplices_containing_shared_node(j);
                for(int k=j+1;k<=simplices_containing_shared_node.m;++k){
                    int simplex_k=simplices_containing_shared_node(k);
                    // since new_simplex was added on the end, it will only appear at simplex_k
                    if(simplex_k==new_simplex && ((simplex_i==simplices[1] && simplex_j==simplices[2]) || (simplex_i==simplices[2] && simplex_j==simplices[1])))
                        // don't look for an intersection among the input simplices! that will certainly not show up in the intersection registry
                        continue;
                    VECTOR<int,3> triplet(simplex_i,simplex_j,simplex_k);
                    shared_nodes.Resize(0);
                    cutting_simplices->Shared_Nodes_On_Simplices(triplet,shared_nodes);
                    if(shared_nodes.m>1) continue; // The intersection does not define a single point
                    assert(shared_nodes.m==1 && shared_nodes(1)==shared_node);
                    int intersection_index=intersection_registry->Intersection(triplet);
                    if(intersection_index!=0){Add_Non_Tet_Node_Intersection_To_Registry(simplices,all_weights,intersection_index);return;}}}}
        Add_Non_Tet_Node_Intersection_To_Registry(simplices,all_weights);
        return;}
    // If any pair of simplices share precisely one node (not on the 3rd simplex), terminate early.
    for(int i=1;i<=2;++i) for(int j=i+1;j<=3;++j){
        shared_nodes.Resize(0);
        cutting_simplices->Shared_Nodes_On_Simplices(VECTOR<int,2>(converted_simplices[i],converted_simplices[j]),shared_nodes);
        if(shared_nodes.m==1) return;}
    // Does the new simplex share an edge with (exactly) one of the old simplices?
    for(int i=1;i<=2;i++){
        shared_nodes.Resize(0);
        cutting_simplices->Shared_Nodes_On_Simplices(VECTOR<int,2>(converted_simplices[i],converted_simplices[3]),shared_nodes);
        if(shared_nodes.m==2 && !cutting_simplices->Simplex_Is_Fake(converted_simplices[3-i])){ // the simplex intersecting the edge shouldn't be fake (2 edges never intersect)
            VECTOR<int,2> shared_edge(shared_nodes(1),shared_nodes(2));
            // Find if there is another set of simplices that shares the same edge. If there is an intersection, it already has it recorded, if not then no intersection exists.
            ARRAY<int> simplices_containing_shared_edge;
            const ARRAY<int>& intersecting_simplices=intersecting_simplices_per_simplex(new_simplex);
            for(int j=1;j<=intersecting_simplices.m;j++){
                int old_simplex=intersecting_simplices(j);
                if(old_simplex>=new_simplex) continue;
                if(cutting_simplices->simplices(old_simplex).nodes.Contains_All(shared_edge)) simplices_containing_shared_edge.Append(old_simplex);}
            for(int j1=1;j1<simplices_containing_shared_edge.m;++j1){
                int local_simplex1=simplices_containing_shared_edge(j1);
                for(int j2=j1+1;j2<=simplices_containing_shared_edge.m;++j2){
                    int local_simplex2=simplices_containing_shared_edge(j2);
                    if(cutting_simplices->simplices(local_simplex1).nodes==cutting_simplices->simplices(local_simplex2).nodes) continue;
                    int old_simplex1=converted?cutting_simplices->simplices(local_simplex1).parent:local_simplex1;
                    int old_simplex2=converted?cutting_simplices->simplices(local_simplex2).parent:local_simplex2;
                    int intersection_index=intersection_registry->Intersection(VECTOR<int,3>(converted_simplices[3-i],old_simplex1,old_simplex2));
                    if(intersection_index!=0){
                        VECTOR<VECTOR<T,2>,3> all_weights;
                        Get_Simplex_Weights_For_Edge_Triangle_Intersection(simplices,3-i,shared_edge,all_weights);
                        Add_Non_Tet_Node_Intersection_To_Registry(VECTOR<int,1>(converted_simplices[3]),VECTOR<VECTOR<T,2>,1>(all_weights[3]),intersection_index);}
                    return;}}
            // None already found, check if there is an intersection
            VECTOR<VECTOR<T,2>,3> all_weights;
            bool intersects=Intersects_And_Get_Simplex_Weights_For_Edge_Triangle_Intersection(simplices,3-i,shared_edge,all_weights);
            if(intersects) Add_Non_Tet_Node_Intersection_To_Registry(converted_simplices,all_weights);
            return;}}
    // old simplices share edge?
    shared_nodes.Resize(0);
    cutting_simplices->Shared_Nodes_On_Simplices(VECTOR<int,2>(converted_simplices[1],converted_simplices[2]),shared_nodes);
    if(shared_nodes.m==2&&!cutting_simplices->Simplex_Is_Fake(converted_simplices[3])){
        // the simplex intersecting the edge shouldn't be fake
        VECTOR<int,2> shared_edge(shared_nodes(1),shared_nodes(2));
        // check if it already exists
        if(intersection_registry->intersections_on_simplex.m>=converted_simplices[3]){
            const ARRAY<int>& intersections_on_new_simplex=intersection_registry->intersections_on_simplex(converted_simplices[3]);
            for(int i=1;i<=intersections_on_new_simplex.m;++i){
                int old_intersection=intersections_on_new_simplex(i);
                const ARRAY<int>& simplices_on_old_intersection=intersection_registry->simplices_on_intersection(old_intersection);
                int count=0;
                for(int j=1;j<=simplices_on_old_intersection.m;++j){
                    const VECTOR<int,3>& simplex_nodes=cutting_simplices->simplices(simplices_on_old_intersection(j)).nodes;
                    if(simplex_nodes.Contains_All(shared_edge)) ++count;
                    if(count>=2){
                        VECTOR<VECTOR<T,2>,3> all_weights;
                        Get_Simplex_Weights_For_Edge_Triangle_Intersection(simplices,3,shared_edge,all_weights);
                        Add_Non_Tet_Node_Intersection_To_Registry(converted_simplices,all_weights,old_intersection);
                        return;}}}}
        // it doesn't exist, check intersection
        VECTOR<VECTOR<T,2>,3> all_weights;
        bool intersects=Intersects_And_Get_Simplex_Weights_For_Edge_Triangle_Intersection(simplices,3,shared_edge,all_weights);
        if(intersects) Add_Non_Tet_Node_Intersection_To_Registry(converted_simplices,all_weights);
        return;}
    // normal intersection?
    for(int i=1;i<=3;i++) if(cutting_simplices->Simplex_Is_Fake(simplices[i])) return; // we know not all are fake from above, but if any are fake we cannot intersect
    typedef typename CUTTING_SIMPLEX<T,3>::GET_ADAPTIVE_WEIGHTS_RESULT_TYPE WEIGHT_TYPE;
    // 3 tris, 3 nodes per tri, 3 coordinates per node
    VECTOR<VECTOR<VECTOR<WEIGHT_TYPE,3>,3>,3> adaptive_weights_for_simplices;
    for(int i=1;i<=3;++i) cutting_simplices->simplices(simplices[i]).Get_Adaptive_Weights(adaptive_weights_for_simplices[i]);
    VECTOR<TV,3> all_weights;
    bool is_degenerate;
    bool intersects=Intersects_And_Intersection_Coordinates<void>(adaptive_weights_for_simplices[1],adaptive_weights_for_simplices[2],adaptive_weights_for_simplices[3],
        all_weights[1],all_weights[2],all_weights[3],abs_tol_for_level2_barycentric_coordinates,&is_degenerate);
    if(is_degenerate){LOG::cerr<<"Perform_Smart_Simplex_Intersection(...) line "<<__LINE__<<":  degeneracy; simplices = "<<simplices<<std::endl;PHYSBAM_FATAL_ERROR();}
    if(!intersects) return;
    VECTOR<VECTOR<T,2>,3> final_weights;
    for(int i=1;i<=3;++i) for(int j=1;j<=2;++j) final_weights[i][j]=all_weights[i][j];
    Add_Non_Tet_Node_Intersection_To_Registry(converted_simplices,final_weights);
}
//#####################################################################
// Function Split_Existing_Polygon_Edges
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Split_Existing_Polygon_Edges()
{
    // TODO: maybe do this per polygon for parallelization...
    for(int i=intersection_registry->index_for_last_old_intersection+1;i<=intersection_registry->Number_Of_Intersections();i++){
        const ARRAY<int>& simplices=intersection_registry->simplices_on_intersection(i);if(!simplices.m) PHYSBAM_FATAL_ERROR("TODO: Kevin & Efty, Why?");
        HASHTABLE<VECTOR<int,2> > edges_seen;
        // loop over every pair of simplices
        for(int j=1;j<=simplices.m;j++) for(int k=j+1;k<=simplices.m;k++){
            int simplex_1=simplices(j),simplex_2=simplices(k);
            // TODO: why quit if both are global cut faces
            if(cutting_simplices->simplices(simplex_1).type==CUTTING_SIMPLEX<T,3>::GLOBAL_CUT_FACE && cutting_simplices->simplices(simplex_2).type==CUTTING_SIMPLEX<T,3>::GLOBAL_CUT_FACE) continue;
            ARRAY<int> all_particles;Get_Particles_On_Simplices(VECTOR<int,2>(simplex_1,simplex_2),all_particles);
            if(all_particles.m<=2) continue; // we have an edge which has not been cut
            // get 2D weights on any simplex
            const VECTOR<T,2>& new_particle_weights_on_simplex=
                intersection_registry->Get_Simplex_Weights_Of_Intersection(i,simplex_1);
            int closest_particle_left=0,closest_particle_right=0;
            T closest_left_val=-std::numeric_limits<T>::max(),closest_right_val=std::numeric_limits<T>::max();
            // find best conditioned axis to sort particle locations
            VECTOR<T,2> delta_weights=intersection_registry->Get_Simplex_Weights_Of_Intersection(all_particles(1),simplex_1)-
                intersection_registry->Get_Simplex_Weights_Of_Intersection(all_particles(2),simplex_1);
            int dominant_axis=delta_weights.Dominant_Axis();
            // find left and right particle that we wish to insert the new ith particle between
            for(int p=1;p<=all_particles.m;p++) if(all_particles(p)<i){ // only consider particles that existed before or have already split polygon edges
                const VECTOR<T,2>& particle_weights_on_simplex=intersection_registry->Get_Simplex_Weights_Of_Intersection(all_particles(p),simplex_1);
                if(particle_weights_on_simplex[dominant_axis]<new_particle_weights_on_simplex[dominant_axis] && particle_weights_on_simplex[dominant_axis]>closest_left_val){
                    closest_particle_left=all_particles(p);closest_left_val=particle_weights_on_simplex[dominant_axis];}
                else if(particle_weights_on_simplex[dominant_axis]>new_particle_weights_on_simplex[dominant_axis] && particle_weights_on_simplex[dominant_axis]<closest_right_val){
                    closest_particle_right=all_particles(p);closest_right_val=particle_weights_on_simplex[dominant_axis];}}
            if(!closest_particle_left || !closest_particle_right) continue; // there is no edge to split
            // split the edge
            if(!polygon_mesh.segment_mesh->Simplex(VECTOR<int,2>(closest_particle_left,closest_particle_right))) continue; // no real edge found
            ARRAY<int> nodes_shared;cutting_simplices->Shared_Nodes_On_Simplices(VECTOR<int,2>(simplex_1,simplex_2),nodes_shared);
            // look at a real edge only once
            if(nodes_shared.m==2 && !edges_seen.Set(VECTOR<int,2>(nodes_shared(1),nodes_shared(2)).Sorted())) continue;
            polygon_mesh.Split_Polygon_Edge(closest_particle_left,closest_particle_right,i);
            if(verbose) {std::stringstream ss;ss<<"Split edge "<<closest_particle_left<<", "<<closest_particle_right<<" with particle "<<i<<std::endl;LOG::filecout(ss.str());}}}
}
//#####################################################################
// Function Split_Existing_Polygons
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Split_Existing_Polygons()
{
    // TODO: What's the purpose of polygons_per_element???
    polygons_per_element=CONSTANT_ARRAY<ARRAY<int> >(current_tetrahedralized_volume->mesh.elements.m,ARRAY<int>());
    for(int tet=1;tet<=current_tetrahedralized_volume->mesh.elements.m;tet++){
        const ARRAY<int>& simplices=simplices_per_current_tet(tet);
        int last_old_simplex_index=last_old_simplex_index_per_current_tet(tet);
        for(int i=1;i<=last_old_simplex_index;++i){
            int old_simplex=simplices(i);
            ARRAY<VECTOR<int,2> > new_unoriented_segments_on_simplex;
            const ARRAY<int>& intersecting_simplices=intersecting_simplices_per_simplex(old_simplex);
            for(int j=1;j<=intersecting_simplices.m;++j){
                int new_simplex=intersecting_simplices(j);
                if(new_simplex<=cutting_simplices->index_for_last_old_cutting_simplex) continue;// assume simplices are ordered old first, new last
                Get_Unoriented_Segments_On_Two_Simplices(old_simplex,new_simplex,new_unoriented_segments_on_simplex);}
            const ARRAY<int>& cutting_polygons_on_simplex=cutting_polygons_per_cutting_simplex(old_simplex);
            for(int j=1;j<=cutting_polygons_on_simplex.m;++j){
                int cutting_polygon_index=cutting_polygons_on_simplex(j);
                if(!new_unoriented_segments_on_simplex.m) polygons_per_element(tet).Append(cutting_polygon_index);
                else{
                    const CUTTING_POLYGON& cutting_polygon=current_cutting_polygons(cutting_polygon_index);
                    if(cutting_polygon.simplex_owner!=old_simplex) PHYSBAM_FATAL_ERROR("TODO: remove this check");
                    ARRAY<VECTOR<int,2> > all_polygonal_segments_on_simplex;
                    Get_Polygon_Edges(cutting_polygon.polygon_index,all_polygonal_segments_on_simplex);
                    ARRAY<ARRAY<ARRAY<int> > > final_polygon_element_particles;
                    Divide_Polygon_Particles_With_New_Segments(all_polygonal_segments_on_simplex,new_unoriented_segments_on_simplex,
                        polygon_mesh.elements(cutting_polygon.polygon_index),old_simplex,cutting_polygon.flipped,final_polygon_element_particles);
                    // replace old one with first new
                    polygon_mesh.elements(cutting_polygon.polygon_index)=final_polygon_element_particles(1);
                    polygons_per_element(tet).Append(cutting_polygon_index);
                    // make all the rest
                    for(int k=2;k<=final_polygon_element_particles.m;++k){
                        int polygon_element_index=polygon_mesh.elements.Append(final_polygon_element_particles(k));
                        int new_cutting_polygon_index=current_cutting_polygons.Append(
                            CUTTING_POLYGON(polygon_element_index,old_simplex,cutting_polygon.flipped,cutting_polygon.polygon_type));
                        polygons_per_element(tet).Append(new_cutting_polygon_index);}}}}}
    if(verbose){std::stringstream ss;ss<<std::endl;for(int i=1;i<=current_cutting_polygons.m;++i)
        ss<<"Cutting polygon "<<i<<" has polygon element "<<current_cutting_polygons(i).polygon_index<<": "<<polygon_mesh.elements(current_cutting_polygons(i).polygon_index)
                <<std::endl;LOG::filecout(ss.str());}
}
//#####################################################################
// Function Add_Polygons_On_New_Simplices
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Add_Polygons_On_New_Simplices()
{
    for(int tet=1;tet<=current_tetrahedralized_volume->mesh.elements.m;tet++){
        if(verbose) {std::stringstream ss;ss<<"polygons_per_element("<<tet<<")="<<polygons_per_element(tet)<<std::endl;LOG::filecout(ss.str());}
        const ARRAY<int>& simplices=simplices_per_current_tet(tet);
        int last_old_simplex_index=last_old_simplex_index_per_current_tet(tet);
        for(int i=last_old_simplex_index+1;i<=simplices.m;++i){
            int new_simplex=simplices(i);
            if(cutting_simplices->Simplex_Is_Fake(new_simplex)) continue;
            if(verbose) {std::stringstream ss;ss<<"new_simplex="<<new_simplex<<std::endl;LOG::filecout(ss.str());}
            ARRAY<VECTOR<int,2> > unoriented_existing_segments;
            const ARRAY<int>& intersecting_simplices=intersecting_simplices_per_simplex(new_simplex);
            for(int j=1;j<=intersecting_simplices.m;++j){
                int old_simplex=intersecting_simplices(j);
                if(old_simplex>cutting_simplices->index_for_last_old_cutting_simplex) continue; // assume simplices are ordered old first, new last
                if(cutting_simplices->Simplex_Is_Fake(old_simplex)) continue;
                Get_Unoriented_Segments_On_Two_Simplices(new_simplex,old_simplex,unoriented_existing_segments);}
            if(verbose) {std::stringstream ss;ss<<"unoriented_existing_segments="<<unoriented_existing_segments<<std::endl;LOG::filecout(ss.str());}
            // orient them or add segments both ways
            ARRAY<VECTOR<int,2> > oriented_segments_on_material_and_simplex_boundaries;
            ARRAY<VECTOR<int,2> > killer_segments;
            TV simplex_normal=Get_Child_Cutting_Simplex_Local_Normal(new_simplex);
            for(int j=1;j<=unoriented_existing_segments.m;++j){
                const VECTOR<int,2>& unoriented_segment=unoriented_existing_segments(j);
                VECTOR<int,2> reversed_segment(unoriented_segment(2),unoriented_segment(1));
                int any_valid_cutting_polygon=-1;
                int count=0;
                for(int k=1;k<=polygons_per_element(tet).m;++k){
                    int cutting_polygon_index=polygons_per_element(tet)(k);
                    int polygon_element_index=current_cutting_polygons(cutting_polygon_index).polygon_index;
                    ARRAY<VECTOR<int,2> > polygon_segments;
                    Get_Polygon_Edges(polygon_element_index,polygon_segments);
                    if(polygon_segments.Contains(unoriented_segment)){++count;any_valid_cutting_polygon=cutting_polygon_index;}
                    if(polygon_segments.Contains(reversed_segment)){++count;any_valid_cutting_polygon=cutting_polygon_index;}}
                if(count==0) continue;
                assert(count%2==0&&any_valid_cutting_polygon>0);
                // the two simplices are topologically connected, so any segments on boundary will be contributed by new simplex down below
                int cutting_simplex=current_cutting_polygons(any_valid_cutting_polygon).simplex_owner;
                ARRAY<int> nodes_shared;
                cutting_simplices->Shared_Nodes_On_Simplices(VECTOR<int,2>(cutting_simplex,new_simplex),nodes_shared);
                if(nodes_shared.m==2) continue;
                // they actually intersect
                if(count==2){
                    TV polygon_normal=Get_Child_Cutting_Simplex_Local_Normal(current_cutting_polygons(any_valid_cutting_polygon).simplex_owner);
                    if(current_cutting_polygons(any_valid_cutting_polygon).flipped) polygon_normal*=(T)-1;
                    TV proposed_edge_dir=VECTOR<T,3>::Cross_Product(simplex_normal,polygon_normal);
                    TV actual_edge_dir=(Get_Particle_Position_In_Local_Tet_Space(unoriented_segment(2),new_simplex)-
                        Get_Particle_Position_In_Local_Tet_Space(unoriented_segment(1),new_simplex)).Normalized();
                    oriented_segments_on_material_and_simplex_boundaries.Append_Unique(unoriented_segment);
                    oriented_segments_on_material_and_simplex_boundaries.Append_Unique(reversed_segment);
                    if(VECTOR<T,3>::Dot_Product(proposed_edge_dir,actual_edge_dir)<0) killer_segments.Append(unoriented_segment);
                    else killer_segments.Append(reversed_segment);}
                else{
                    oriented_segments_on_material_and_simplex_boundaries.Append_Unique(unoriented_segment);
                    oriented_segments_on_material_and_simplex_boundaries.Append_Unique(reversed_segment);}}
            // add non-existent segments on triangle boundary
            ARRAY<VECTOR<int,2> > new_unoriented_segments_to_examine_later;
            for(int j=1;j<=intersecting_simplices.m;++j){
                int other_simplex=intersecting_simplices(j);
                if(other_simplex==new_simplex) continue;
                ARRAY<int> shared_nodes;
                cutting_simplices->Shared_Nodes_On_Simplices(VECTOR<int,2>(new_simplex,other_simplex),shared_nodes);
                if(shared_nodes.m==2){
                    ARRAY<VECTOR<int,2> > new_oriented_segments_on_simplex;
                    Get_Oriented_Segments_On_Edge_Of_Two_Simplices(new_simplex,other_simplex,shared_nodes(1),shared_nodes(2),new_oriented_segments_on_simplex);
                    oriented_segments_on_material_and_simplex_boundaries.Append_Unique_Elements(new_oriented_segments_on_simplex);}
                else if(shared_nodes.m<=1) Get_Unoriented_Segments_On_Two_Simplices(new_simplex,other_simplex,new_unoriented_segments_to_examine_later);}
            for(int j=new_unoriented_segments_to_examine_later.m;j>=1;--j)
                if(oriented_segments_on_material_and_simplex_boundaries.Contains(new_unoriented_segments_to_examine_later(j)))
                    new_unoriented_segments_to_examine_later.Remove_Index(j);
            if(verbose){
                std::stringstream ss;
                ss<<"Killer segments: "<<killer_segments<<std::endl;
                ss<<"Segments to examine later: "<<new_unoriented_segments_to_examine_later<<std::endl;
                ss<<"Performing 2D region finding on segments: "<<std::endl<<oriented_segments_on_material_and_simplex_boundaries<<std::endl;
                LOG::filecout(ss.str());}
            // 2D region finding for prior material boundary
            ARRAY<ARRAY<VECTOR<int,2> > > unconnected_polygonal_regions;
            Two_Dimensional_Region_Finding_On_Cutting_Simplex(new_simplex,false,oriented_segments_on_material_and_simplex_boundaries,unconnected_polygonal_regions,true);
            if(verbose){std::stringstream ss;ss<<"**** Original material boundary within simplex "<<new_simplex<<std::endl<<unconnected_polygonal_regions<<std::endl;LOG::filecout(ss.str());}
            // prune away regions with killer segments
            for(int j=1;j<=killer_segments.m;++j) for(int k=unconnected_polygonal_regions.m;k>=1;--k)
                if(unconnected_polygonal_regions(k).Contains(killer_segments(j))) unconnected_polygonal_regions.Remove_Index(k);
            // inside outside classification to form polygons
            ARRAY<ARRAY<ARRAY<int> > > final_polygon_element_particles;
            Inside_Outside_Determination_For_Unconnected_Polygonal_Regions(new_simplex,false,unconnected_polygonal_regions,final_polygon_element_particles);
            if(verbose){std::stringstream ss;ss<<"**** Polygon before adding new edges: "<<std::endl<<final_polygon_element_particles<<std::endl;LOG::filecout(ss.str());}
            // add segments from new triangles that don't share an edge, cut polygons, and find final regions
            ARRAY<ARRAY<ARRAY<int> > > final_final_polygon_element_particles;
            for(int j=1;j<=final_polygon_element_particles.m;++j){
                ARRAY<VECTOR<int,2> > all_final_segments;
                for(int k=1;k<=final_polygon_element_particles(j).m;++k) for(int ell=1;ell<=final_polygon_element_particles(j)(k).m;++ell){
                    int next_ell=ell%final_polygon_element_particles(j)(k).m+1;
                    all_final_segments.Append(VECTOR<int,2>(final_polygon_element_particles(j)(k)(ell),final_polygon_element_particles(j)(k)(next_ell)));}
                ARRAY<ARRAY<ARRAY<int> > > new_polygon_element_particles_for_this_polygon;
                Divide_Polygon_Particles_With_New_Segments(all_final_segments,new_unoriented_segments_to_examine_later,final_polygon_element_particles(j),
                    new_simplex,false,new_polygon_element_particles_for_this_polygon);
                final_final_polygon_element_particles.Append_Elements(new_polygon_element_particles_for_this_polygon);}
            // make new polygons
            for(int k=1;k<=final_final_polygon_element_particles.m;++k){
                ARRAY<ARRAY<int> > reversed_polygon_particles(final_final_polygon_element_particles(k).m);
                for(int p=1;p<=final_final_polygon_element_particles(k).m;++p)
                    for(int q=final_final_polygon_element_particles(k)(p).m;q>=1;--q)
                        reversed_polygon_particles(p).Append(final_final_polygon_element_particles(k)(p)(q));
                int polygon_element_index_1=polygon_mesh.elements.Append(final_final_polygon_element_particles(k));
                int new_cutting_polygon_index_1=current_cutting_polygons.Append(
                    CUTTING_POLYGON(polygon_element_index_1,new_simplex,false,CUTTING_POLYGON::TRIANGLE_CLIPPED));
                polygons_per_element(tet).Append(new_cutting_polygon_index_1);
                int polygon_element_index_2=polygon_mesh.elements.Append(reversed_polygon_particles);
                int new_cutting_polygon_index_2=current_cutting_polygons.Append(
                    CUTTING_POLYGON(polygon_element_index_2,new_simplex,true,CUTTING_POLYGON::TRIANGLE_CLIPPED));
                polygons_per_element(tet).Append(new_cutting_polygon_index_2);
                if(verbose){
                    std::stringstream ss;
                    ss<<"***** Made final polygon "<<final_final_polygon_element_particles(k)<<std::endl;
                    ss<<"***** Made final polygon "<<reversed_polygon_particles<<std::endl;
                    LOG::filecout(ss.str());}}}
        if(verbose) {std::stringstream ss;ss<<"polygons_per_element("<<tet<<")="<<polygons_per_element(tet)<<std::endl;LOG::filecout(ss.str());}}
    polygon_mesh.Refresh_Auxiliary_Structures();
}
//#####################################################################
// Function Divide_Polygon_Particles_With_New_Segments
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Divide_Polygon_Particles_With_New_Segments(ARRAY<VECTOR<int,2> >& all_segments,const ARRAY<VECTOR<int,2> >& possible_segments_to_add,const ARRAY<ARRAY<int> >& polygon_particles,
    const int cutting_simplex,const bool flipped,ARRAY<ARRAY<ARRAY<int > > >& final_polygon_element_particles) const
{
    if(verbose){std::stringstream ss;ss<<"All segments on simplex "<<cutting_simplex<<" before adding new: ";for(int j=1;j<=all_segments.m;j++) ss<<all_segments(j)<<"; ";ss<<std::endl;LOG::filecout(ss.str());}
    // choose subset of new segments to add to this polygon
    bool added_any=false;
    for(int k=1;k<=possible_segments_to_add.m;k++){const VECTOR<int,2> possible_segment_to_add=possible_segments_to_add(k);
        if(all_segments.Contains(possible_segment_to_add) || all_segments.Contains(possible_segment_to_add.Reversed())) continue;
        if(Potential_Segment_Should_Be_Added_To_Polygon(polygon_particles,flipped,cutting_simplex,possible_segment_to_add)){
            added_any=true;all_segments.Append(possible_segment_to_add);
            all_segments.Append(possible_segment_to_add.Reversed());}}
    if(verbose){std::stringstream ss;ss<<"All segments on simplex "<<cutting_simplex<<" after adding new: ";for(int j=1;j<=all_segments.m;j++) ss<<all_segments(j)<<"; ";ss<<std::endl;LOG::filecout(ss.str());}
    if(!added_any) // do nothing if no segments were added
        final_polygon_element_particles.Append(polygon_particles);
    else{ // 2D region finding
        ARRAY<ARRAY<VECTOR<int,2> > > unconnected_polygonal_regions;
        Two_Dimensional_Region_Finding_On_Cutting_Simplex(cutting_simplex,flipped,all_segments,unconnected_polygonal_regions);
        if(verbose) {std::stringstream ss;ss<<"Simplex "<<cutting_simplex<<" has unconnected polygon regions "<<std::endl<<unconnected_polygonal_regions<<std::endl;LOG::filecout(ss.str());}
        Inside_Outside_Determination_For_Unconnected_Polygonal_Regions(cutting_simplex,flipped,unconnected_polygonal_regions,final_polygon_element_particles);
        if(verbose){std::stringstream ss;ss<<"***** Final polygon: "<<std::endl<<final_polygon_element_particles<<std::endl;LOG::filecout(ss.str());}}
}
//#####################################################################
// Function Get_Unoriented_Segments_On_Two_Simplices
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Get_Unoriented_Segments_On_Two_Simplices(const int simplex_1,const int simplex_2,ARRAY<VECTOR<int,2> >& new_unoriented_segments_on_simplex) const
{
    ARRAY<int> particles;Get_Particles_On_Simplices(VECTOR<int,2>(simplex_1,simplex_2),particles);
    if(particles.m<2) return;
    // sort the particles
    ARRAY<VECTOR<T,2> > particle_weights;
    for(int k=1;k<=particles.m;k++) particle_weights.Append(intersection_registry->Get_Simplex_Weights_Of_Intersection(particles(k),simplex_1));
    int dominant_axis=(particle_weights(2)-particle_weights(1)).Dominant_Axis();
    ARRAY<int> permutation(IDENTITY_ARRAY<>(particles.m));
    Sort(permutation,Indirect_Comparison(particle_weights.Project(dominant_axis)));
    for(int i=1;i<permutation.m;i++) new_unoriented_segments_on_simplex.Append_Unique(VECTOR<int,2>(particles(permutation(i)),particles(permutation(i+1))));
    if(verbose) {std::stringstream ss;ss<<"Simplex "<<simplex_1<<" and "<<simplex_2<<" have particles "<<particles<<" segments "<<new_unoriented_segments_on_simplex<<std::endl;LOG::filecout(ss.str());}
}
//#####################################################################
// Function Get_Oriented_Segments_On_Edge_Of_Two_Simplices
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Get_Oriented_Segments_On_Edge_Of_Two_Simplices(const int simplex_1,const int simplex_2,const int shared_node_1,const int shared_node_2,ARRAY<VECTOR<int,2> >& new_oriented_segments_on_simplex) const
{
    ARRAY<int> particles;Get_Particles_On_Simplices(VECTOR<int,2>(simplex_1,simplex_2),particles);
    if(particles.m<2) return;
    // sort the particles
    ARRAY<VECTOR<T,2> > particle_weights;
    for(int k=1;k<=particles.m;k++) particle_weights.Append(intersection_registry->Get_Simplex_Weights_Of_Intersection(particles(k),simplex_1));
    int edge_number=0;
    if(particle_weights(1)(1)==0 && particle_weights(2)(1)==0) edge_number=2;
    else if(particle_weights(1)(2)==0 && particle_weights(2)(2)==0) edge_number=3;
    else edge_number=1;
    int axis=(edge_number==1 || edge_number==3)?1:2;
    for(int j=2;j<=particle_weights.m;j++){VECTOR<T,2> object=particle_weights(j);int k=j;int particle=particles(j);
        while(k>1 && particle_weights(k-1)[axis]>object[axis]){particle_weights(k)=particle_weights(k-1);particles(k)=particles(k-1);k--;}
        particle_weights(k)=object;particles(k)=particle;}
    // make the right way on simplex_1
    for(int k=1;k<particles.m;k++){
        if(edge_number==1 || edge_number==2)
            new_oriented_segments_on_simplex.Append_Unique(VECTOR<int,2>(particles(k+1),particles(k)));
        else new_oriented_segments_on_simplex.Append_Unique(VECTOR<int,2>(particles(k),particles(k+1)));}
}
//#####################################################################
// Function Potential_Segment_Should_Be_Added_To_Polygon
//#####################################################################
template<class T> bool CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Potential_Segment_Should_Be_Added_To_Polygon(const ARRAY<ARRAY<int > >& particles_for_polygon,const bool flipped,const int simplex,const VECTOR<int,2>& nodes) const
{
    // check inside/outside rings
    VECTOR<bool,2> nodes_on_boundary;for(int i=1;i<=2;i++) for(int j=1;j<=particles_for_polygon.m;j++) if(particles_for_polygon(j).Contains(nodes[i])) nodes_on_boundary[i]=true;
    if(nodes_on_boundary.Contains(false)){
        for(int i=1;i<=2;i++) if(!nodes_on_boundary[i]) for(int j=1;j<=particles_for_polygon.m;j++)
            if((j==1)^Point_Is_Inside_Unoriented_Polygon(particles_for_polygon(j),simplex,nodes[i])) return false; // first loop is outer ring and others are negative holes
        return true;}
    else{
        // both endpoints of segment are on existing points, so operate on nodes[1] arbitrarily
        for(int i=1;i<=particles_for_polygon.m;i++){const ARRAY<int>& particles=particles_for_polygon(i);
            if(int j=particles.Find(nodes[1])){
                int p1=particles(j-1<1?particles.m:j-1),p2=particles(j),p3=particles(j+1>particles.m?1:j+1);
                VECTOR<T,2> x1=intersection_registry->Get_Simplex_Weights_Of_Intersection(p1,simplex);
                VECTOR<T,2> x2=intersection_registry->Get_Simplex_Weights_Of_Intersection(p2,simplex);
                VECTOR<T,2> x3=intersection_registry->Get_Simplex_Weights_Of_Intersection(p3,simplex);
                VECTOR<T,2> x0=intersection_registry->Get_Simplex_Weights_Of_Intersection(nodes[2],simplex);
                VECTOR<T,2> v12=x1-x2,v32=x3-x2,v02=x0-x2;
                T angle31=VECTOR<T,2>::Oriented_Angle_Between(v32,v12),angle30=VECTOR<T,2>::Oriented_Angle_Between(v32,v02);
                const T two_pi=2*(T)pi;
                if(angle31<0) angle31+=two_pi;
                if(angle30<0) angle30+=two_pi;
                if(angle30==0||angle31==angle30) PHYSBAM_FATAL_ERROR("nodes[2] should not be collinear with either v12 or v32; probable roundoff error?");
                return angle31==0||((i>1)^flipped^(angle31>angle30));}}
        PHYSBAM_FATAL_ERROR("Point should not be reached");}
}
//#####################################################################
// Function Point_Is_Inside_Unoriented_Polygon
//#####################################################################
template<class T> bool CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Point_Is_Inside_Unoriented_Polygon(const ARRAY<int>& polygon_particles,const int simplex_owner,const int point) const
{
    T sum=0;VECTOR<T,2> point_position=intersection_registry->Get_Simplex_Weights_Of_Intersection(point,simplex_owner);
    for(int k=1;k<=polygon_particles.m;k++){
        int node_1=polygon_particles(k);int node_2=polygon_particles(k+1>polygon_particles.m?1:k+1);
        VECTOR<T,2> node_1_position=intersection_registry->Get_Simplex_Weights_Of_Intersection(node_1,simplex_owner);
        VECTOR<T,2> node_2_position=intersection_registry->Get_Simplex_Weights_Of_Intersection(node_2,simplex_owner);
        VECTOR<T,2> direction_1=node_1_position-point_position;VECTOR<T,2> direction_2=node_2_position-point_position;
        sum+=VECTOR<T,2>::Oriented_Angle_Between(direction_2,direction_1);}
    return abs(sum)>(T)pi;
}
//#####################################################################
// Function Two_Dimensional_Region_Finding_On_Cutting_Simplex
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Two_Dimensional_Region_Finding_On_Cutting_Simplex(const int cutting_simplex,const bool flipped,const ARRAY<VECTOR<int,2> >& segments,
    ARRAY<ARRAY<VECTOR<int,2> > >& unconnected_polygonal_regions,const bool survive_incomplete_regions) const
{
    // TODO: figure out why we need to survive incomplete regions
    ARRAY<VECTOR<int,2> > unused_segments=segments;
    if(!survive_incomplete_regions){
        while(unused_segments.m){
            VECTOR<int,2> nodes=unused_segments.Pop();
            ARRAY<VECTOR<int,2> > segment_loop;segment_loop.Append(nodes);
            int start=nodes[1];
            for(;;){
                T min_angle=100;int next_segment=0;
                for(int i=1;i<=segments.m;i++){VECTOR<int,2> next_nodes=segments(i);
                    if(segments(i)[1]==nodes[2]){
                        T angle;
                        if(segments(i)[2]!=nodes[1]){
                            VECTOR<T,2> X1=intersection_registry->Get_Simplex_Weights_Of_Intersection(nodes[1],cutting_simplex);
                            VECTOR<T,2> X2=intersection_registry->Get_Simplex_Weights_Of_Intersection(nodes[2],cutting_simplex);
                            VECTOR<T,2> next_X1=intersection_registry->Get_Simplex_Weights_Of_Intersection(next_nodes[1],cutting_simplex);
                            VECTOR<T,2> next_X2=intersection_registry->Get_Simplex_Weights_Of_Intersection(next_nodes[2],cutting_simplex);
                            VECTOR<T,2> direction=X2-X1,next_direction=next_X2-next_X1;
                            //if(direction.Magnitude()<1e-12) PHYSBAM_FATAL_ERROR(); TODO: remove debugging
                            //if(next_direction.Magnitude()<1e-12) PHYSBAM_FATAL_ERROR();
                            angle=VECTOR<T,2>::Oriented_Angle_Between(next_direction,direction);
                            if(flipped) angle=-angle;}
                        else angle=50;
                        if(angle<min_angle){min_angle=angle;next_segment=i;}}}
                if(!next_segment) break;
                int next_unused_segment_index=unused_segments.Find(segments(next_segment));
                if(!next_unused_segment_index) break;
                nodes=unused_segments(next_unused_segment_index);
                unused_segments.Remove_Index_Lazy(next_unused_segment_index);
                segment_loop.Append(nodes);}
            if(nodes[2]!=start) PHYSBAM_FATAL_ERROR("did not arrive back at start node");
            unconnected_polygonal_regions.Append(segment_loop);}}
    else{
        while(unused_segments.m){
            int starting_index=0;
            for(int j=1;j<=unused_segments.m;j++){
                if(!unused_segments.Contains(unused_segments(j).Reversed())){
                    int count=0;for(int k=1;k<=unused_segments.m;k++) if(unused_segments(k)[2]==unused_segments(j)[1]) count++;
                    if(count<=1){starting_index=j;break;}}}
            if(!starting_index) for(int j=1;j<=unused_segments.m;j++){
                int count=0;for(int k=1;k<=unused_segments.m;k++) if(unused_segments(k)[2]==unused_segments(j)[1]) count++;
                if(count<=1){starting_index=j;break;}}
            if(!starting_index) starting_index=1;
            ARRAY<VECTOR<int,2> > segments_for_this_fragment;segments_for_this_fragment.Append(unused_segments(starting_index));
            unused_segments.Remove_Index(starting_index);
            const int node_to_end_at=segments_for_this_fragment(1)[1];
            for(;;){
                const VECTOR<int,2> seg_to_examine=segments_for_this_fragment.Last();
                // best segment going the right way
                T min_value_found_correct=1e2;int next_seg_in_whole_list_correct=0;
                for(int i=1;i<=segments.m;i++) if(segments(i)[1]==seg_to_examine[2]){T angle=0;
                    if(segments(i)[2]!=seg_to_examine[1]){
                        VECTOR<T,2> cus1=intersection_registry->Get_Simplex_Weights_Of_Intersection(seg_to_examine[1],cutting_simplex);
                        VECTOR<T,2> cus2=intersection_registry->Get_Simplex_Weights_Of_Intersection(seg_to_examine[2],cutting_simplex);
                        VECTOR<T,2> cas1=intersection_registry->Get_Simplex_Weights_Of_Intersection(segments(i)[1],cutting_simplex);
                        VECTOR<T,2> cas2=intersection_registry->Get_Simplex_Weights_Of_Intersection(segments(i)[2],cutting_simplex);
                        VECTOR<T,2> seg_to_examine_direction=cus2-cus1,candidate_seg_direction=cas2-cas1;
                        //if(seg_to_examine_direction.Magnitude()<1e-12) PHYSBAM_FATAL_ERROR(); TODO: remove debugging
                        //if(candidate_seg_direction.Magnitude()<1e-12) PHYSBAM_FATAL_ERROR();
                        angle=VECTOR<T,2>::Oriented_Angle_Between(candidate_seg_direction,seg_to_examine_direction);
                        if(flipped) angle=-angle;}
                    else angle=99;
                    if(angle<min_value_found_correct){min_value_found_correct=angle;next_seg_in_whole_list_correct=i;}}
                // best segment going the wrong way
                T min_value_found_reversed=1e2;int next_seg_in_whole_list_reversed=0;
                for(int i=1;i<=segments.m;i++) if(segments(i)[2]==seg_to_examine[2]){T angle=0;
                    if(segments(i)[1]!=seg_to_examine[1]){
                        VECTOR<T,2> cus1=intersection_registry->Get_Simplex_Weights_Of_Intersection(seg_to_examine[1],cutting_simplex);
                        VECTOR<T,2> cus2=intersection_registry->Get_Simplex_Weights_Of_Intersection(seg_to_examine[2],cutting_simplex);
                        VECTOR<T,2> cas1=intersection_registry->Get_Simplex_Weights_Of_Intersection(segments(i)[2],cutting_simplex);
                        VECTOR<T,2> cas2=intersection_registry->Get_Simplex_Weights_Of_Intersection(segments(i)[1],cutting_simplex);
                        VECTOR<T,2> seg_to_examine_direction=cus2-cus1,candidate_seg_direction=cas2-cas1;
                        //if(seg_to_examine_direction.Magnitude()<1e-12) PHYSBAM_FATAL_ERROR(); // TODO: remove debugging
                        //if(candidate_seg_direction.Magnitude()<1e-12) PHYSBAM_FATAL_ERROR();
                        angle=VECTOR<T,2>::Oriented_Angle_Between(candidate_seg_direction,seg_to_examine_direction);
                        if(flipped) angle=-angle;}
                    else angle=99;
                    if(angle<min_value_found_reversed){min_value_found_reversed=angle;next_seg_in_whole_list_reversed=i;}}
                // invalid if reversed is sharper than non-reversed
                if(!next_seg_in_whole_list_correct){
                    if(survive_incomplete_regions) goto Next_Region;
                    else PHYSBAM_FATAL_ERROR();}
                bool reversed_is_sharper=next_seg_in_whole_list_reversed && min_value_found_reversed<min_value_found_correct
                    && segments(next_seg_in_whole_list_correct)!=segments(next_seg_in_whole_list_reversed).Reversed();
                if(reversed_is_sharper){
                    if(survive_incomplete_regions) goto Next_Region;
                    else PHYSBAM_FATAL_ERROR();}
                int index_of_best_segment_in_unused_list=unused_segments.Find(segments(next_seg_in_whole_list_correct));
                // add segment?
                if(index_of_best_segment_in_unused_list){
                    segments_for_this_fragment.Append(unused_segments(index_of_best_segment_in_unused_list));
                    unused_segments.Remove_Index(index_of_best_segment_in_unused_list);}
                else if(segments_for_this_fragment.Last()[2]==node_to_end_at) break;
                else if(survive_incomplete_regions) goto Next_Region;
                else PHYSBAM_FATAL_ERROR();}
            if(segments_for_this_fragment.m==1){if(!survive_incomplete_regions) PHYSBAM_FATAL_ERROR();}
            else unconnected_polygonal_regions.Append(segments_for_this_fragment);
            Next_Region:;}}
}
//#####################################################################
// Function Inside_Outside_Determination_For_Unconnected_Polygonal_Regions
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Inside_Outside_Determination_For_Unconnected_Polygonal_Regions(const int simplex,const bool flipped,const ARRAY<ARRAY<VECTOR<int,2> > >& unconnected_polygonal_regions,
    ARRAY<ARRAY<ARRAY<int > > >& final_polygon_element_particles) const
{
    if(verbose){
        std::stringstream ss;ss<<"Starting inside/outside on unconnected regions: \n"<<unconnected_polygonal_regions;LOG::filecout(ss.str());
        Draw_Polygon(simplex,flipped,unconnected_polygonal_regions);}
    // get all inside/outside flags
    DIRECTED_GRAPH<int> in_out_graph(unconnected_polygonal_regions.m); // TODO: possibly rename graph to reflect correct direction
    bool something_is_contained=false;
    for(int i=1;i<=unconnected_polygonal_regions.m;i++) for(int j=1;j<=unconnected_polygonal_regions.m;j++) if(i!=j){
        int particle=0;
        for(int a=1;a<=unconnected_polygonal_regions(i).m;a++){
            int p=unconnected_polygonal_regions(i)(a)[1];bool ok=true;
            for(int k=1;k<=unconnected_polygonal_regions(j).m;k++){
                if(unconnected_polygonal_regions(j)(k).Contains(p)){ok=false;break;}}
            if(ok){particle=p;break;}}
        if(particle){
            ARRAY<int> particles_on_j_segments(unconnected_polygonal_regions(j).Project(1));
            if(Point_Is_Inside_Unoriented_Polygon(particles_on_j_segments,simplex,particle)){
                if(verbose) PHYSBAM_DEBUG_PRINT("Adding j,i",particles_on_j_segments,simplex,particle,j,i);
                in_out_graph.Add_Edge(j,i);something_is_contained=true;}}}
    if(!something_is_contained) // nothing is inside anything else...
        for(int i=1;i<=unconnected_polygonal_regions.m;i++){
            ARRAY<int> particles(unconnected_polygonal_regions(i).Project(1)); // if nothing is contained, it's a closed region
            if(particles.m<3) PHYSBAM_FATAL_ERROR("degenerate polygons should always be contained in other ones");
            final_polygon_element_particles.Append(ARRAY<ARRAY<int> >());
            final_polygon_element_particles.Last().Append(particles);}
    else{ // something is inside...
        ARRAY<int> depths;in_out_graph.Maximal_Depth_On_Acyclic_Graph(depths);
        ARRAY<bool> unjoined_region_orientations(unconnected_polygonal_regions.m); // true means positively oriented
        ARRAY<bool> unjoined_region_has_no_area(unconnected_polygonal_regions.m);
        ARRAY<int> regions_in_order;
        {ARRAY<int> finish_times;in_out_graph.Topological_Sort_Assuming_Cycle_Free(finish_times,regions_in_order);}
        for(int i=1;i<=unconnected_polygonal_regions.m;i++){
            if(unconnected_polygonal_regions(i).m==2) continue; // degenerate polygons are always negatively oriented, and have no area (TODO: set has_no_area?)
            // determine if region has no area
            bool all_segments_appear_both_ways=true;
            for(int j=1;j<=unconnected_polygonal_regions(i).m;j++){
                if(!unconnected_polygonal_regions(i).Contains(unconnected_polygonal_regions(i)(j).Reversed())){
                    all_segments_appear_both_ways=false;break;}}
            unjoined_region_has_no_area(i)=all_segments_appear_both_ways;
            // figure out orientation
            T sum=0;
            for(int j=1;j<=unconnected_polygonal_regions(i).m;j++){
                const VECTOR<int,2>& seg_1=unconnected_polygonal_regions(i)(j);
                const VECTOR<int,2>& seg_2=unconnected_polygonal_regions(i)(j+1>unconnected_polygonal_regions(i).m?1:j+1);
                VECTOR<T,2> seg_11=intersection_registry->Get_Simplex_Weights_Of_Intersection(seg_1(1),simplex);
                VECTOR<T,2> seg_12=intersection_registry->Get_Simplex_Weights_Of_Intersection(seg_1(2),simplex);
                VECTOR<T,2> seg_21=intersection_registry->Get_Simplex_Weights_Of_Intersection(seg_2(1),simplex);
                VECTOR<T,2> seg_22=intersection_registry->Get_Simplex_Weights_Of_Intersection(seg_2(2),simplex);
                VECTOR<T,2> direction_1=seg_12-seg_11,direction_2=seg_22-seg_21;
                sum+=VECTOR<T,2>::Oriented_Angle_Between(direction_2,direction_1);}
            unjoined_region_orientations(i)=flipped^(sum<0);}
        ARRAY<bool> used(unconnected_polygonal_regions.m);
        for(int i=regions_in_order.m;i>=1;i--){int region_id=regions_in_order(i);if(!used(region_id)){
            used(region_id)=true; // TODO: Continue from here
            final_polygon_element_particles.Append(ARRAY<ARRAY<int> >());
            ARRAY<ARRAY<int> >& particles_on_polygon=final_polygon_element_particles.Last();
            // add particles for parent
            particles_on_polygon.Append(ARRAY<int>(unconnected_polygonal_regions(region_id).Project(1)));
            // add particles for each child
            for(int j=1;j<=in_out_graph.Children(region_id).m;j++){int child_region=in_out_graph.Children(region_id)(j);
                if(used(child_region)) continue;
                if(depths(region_id)+1==depths(child_region) &&
                    (unjoined_region_orientations(region_id)!=unjoined_region_orientations(child_region) || unjoined_region_has_no_area(child_region))){
                        used(child_region)=true;
                        particles_on_polygon.Append(ARRAY<int>(unconnected_polygonal_regions(child_region).Project(1)));}}}}}
}
//#####################################################################
// Function Find_Material_Regions
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Find_Material_Regions()
{
    if(verbose) for(int i=1;i<=current_cutting_polygons.m;i++){
        const CUTTING_POLYGON& cutting_polygon=current_cutting_polygons(i);
        std::stringstream ss;
        ss<<"Cutting polygon "<<i<<" has polygon element index "<<cutting_polygon.polygon_index<<" and is flipped? "<<(cutting_polygon.flipped?"yes":"no")
            <<" and simplex owner "<<cutting_polygon.simplex_owner<<" with particles "<<polygon_mesh.elements(cutting_polygon.polygon_index)
            <<" and untransformed normal "<<Get_Child_Cutting_Simplex_Local_Normal(cutting_polygon.simplex_owner)
            <<" and transformed normal "<<Get_Child_Cutting_Simplex_Local_Normal(cutting_polygon.simplex_owner) * (T)(cutting_polygon.flipped?-1:1)
            <<std::endl<<std::endl;LOG::filecout(ss.str());}
    // polygons have already been duplicated
    const SEGMENT_MESH& polygon_segment_mesh=*polygon_mesh.segment_mesh;
    const ARRAY<ARRAY<ARRAY<PAIR<int,bool> > > >& element_oriented_edges=*polygon_mesh.element_oriented_edges;
    regions_per_tet.Clean_Memory();
    regions_per_tet.Resize(current_tetrahedralized_volume->mesh.elements.m);
    for(int tet=1;tet<=current_tetrahedralized_volume->mesh.elements.m;tet++){const ARRAY<int>& polygons=polygons_per_element(tet);
        HASHTABLE<int,int> polygon_element_to_cutting_polygon;
        for(int i=1;i<=polygons.m;i++) polygon_element_to_cutting_polygon.Insert(current_cutting_polygons(polygons(i)).polygon_index,polygons(i));
        if(verbose) {std::stringstream ss;ss<<"Tet "<<tet<<" has all polygons: "<<polygons<<std::endl;LOG::filecout(ss.str());}
        ARRAY<ARRAY<int> >& regions_for_tet=regions_per_tet(tet);
        int count=0;HASHTABLE<int> polygons_used; // bool does nothing
        while(count<polygons.m){regions_for_tet.Append(ARRAY<int>());ARRAY<int>& polygons_for_new_region=regions_for_tet.Last();
            // find starting polygon
            bool found_starting_index=false;int starting_index=-1;
            for(int i=1;i<=polygons.m;i++) if(!polygons_used.Contains(polygons(i))){found_starting_index=true;starting_index=i;break;}
            assert(found_starting_index);int cutting_polygon_index=polygons(starting_index);
            //g++-4.6 fix - void cast used to avoid set but unused warning in release mode
            (void)found_starting_index;
            STACK<int> polygons_to_examine;polygons_to_examine.Push(cutting_polygon_index);
            // repeat until the list of faces-to-examine is exhausted
            while(!polygons_to_examine.Empty()){
                int polygon_to_examine=polygons_to_examine.Pop();
                if(!polygons_used.Set(polygon_to_examine)) continue;
                polygons_for_new_region.Append(polygon_to_examine);count++;
                const CUTTING_POLYGON& cutting_polygon=current_cutting_polygons(polygon_to_examine);
                int cutting_simplex_owner=cutting_polygon.simplex_owner,polygon_element_index=cutting_polygon.polygon_index;
                TV polygon_to_examine_normal=Get_Child_Cutting_Simplex_Local_Normal(cutting_simplex_owner);
                if(cutting_polygon.flipped) polygon_to_examine_normal*=(T)-1;
                const ARRAY<ARRAY<PAIR<int,bool> > >& edges_for_polygon=element_oriented_edges(polygon_element_index);
                for(int i=1;i<=edges_for_polygon.m;i++) for(int j=1;j<=edges_for_polygon(i).m;j++){
                    const PAIR<int,bool>& oriented_edge=edges_for_polygon(i)(j);int segment_index=oriented_edge.x;bool segment_flipped_in_mesh=oriented_edge.y;
                    // obtain transform for this polygon
                    int node_1=(segment_flipped_in_mesh?polygon_segment_mesh.elements(segment_index)(2):polygon_segment_mesh.elements(segment_index)(1));
                    int node_2=(segment_flipped_in_mesh?polygon_segment_mesh.elements(segment_index)(1):polygon_segment_mesh.elements(segment_index)(2));
                    TV edge_direction=(Get_Particle_Position_In_Local_Tet_Space(node_1,cutting_simplex_owner)-
                                       Get_Particle_Position_In_Local_Tet_Space(node_2,cutting_simplex_owner)).Normalized();
                    assert(edge_direction.Magnitude_Squared()>0);
                    TV e_cross_n=VECTOR<T,3>::Cross_Product(edge_direction,polygon_to_examine_normal);
                    MATRIX<T,3> transform(e_cross_n,-polygon_to_examine_normal,edge_direction);transform.Invert();
                    // find polygon with sharpest angle
                    ARRAY<int> polygons_on_this_edge;T min_value_found=(T)1e2;int best_polygon_so_far=-1;
                    int number_of_incident_polygons=polygon_mesh.Elements_On_Oriented_Edge(node_1,node_2,&polygons_on_this_edge);assert(number_of_incident_polygons>0);
                    for(int k=1;k<=number_of_incident_polygons;k++){int candidate_polygon_element_index=polygons_on_this_edge(k);int candidate_cutting_polygon_index;
                        T angle=0;bool found=polygon_element_to_cutting_polygon.Get(candidate_polygon_element_index,candidate_cutting_polygon_index);if(!found) continue;
                        if(polygon_element_index==polygon_mesh.Opposite_Oriented_Element(candidate_polygon_element_index)){angle=99.;}
                        else{
                            int candidate_simplex_owner=current_cutting_polygons(candidate_cutting_polygon_index).simplex_owner;
                            TV candidate_normal=Get_Child_Cutting_Simplex_Local_Normal(candidate_simplex_owner);
                            if(current_cutting_polygons(candidate_cutting_polygon_index).flipped) candidate_normal*=(T)-1;
                            TV transformed_normal=transform*candidate_normal;
                            angle=VECTOR<T,2>::Oriented_Angle_Between(VECTOR<T,2>(transformed_normal.x,transformed_normal.y),VECTOR<T,2>(0,-1));}
                    if(angle<min_value_found){min_value_found=angle;best_polygon_so_far=candidate_cutting_polygon_index;}}assert(best_polygon_so_far>0);
                    // found the best polygon; if new, add it to the region, enqueue polygon to visit it next
                    if(polygons_used.Contains(best_polygon_so_far)) continue;
                    polygons_to_examine.Push(best_polygon_so_far);}}
            if(polygons_for_new_region.m<=3){
                std::stringstream ss;ss<<" Count is "<<count<<" of "<<polygons.m<<", wants to add polygons "<<polygons_for_new_region<<std::endl;LOG::filecout(ss.str());
                PHYSBAM_FATAL_ERROR();}
            if(verbose){std::stringstream ss;ss<<"Tet "<<tet<<" added a region: "<<polygons_for_new_region<<std::endl;LOG::filecout(ss.str());}}}
}
//#####################################################################
// Function Determine_Duplicate_Tets_And_Duplicate_Particles
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Determine_Duplicate_Tets_And_Duplicate_Particles()
{
    const ARRAY<VECTOR<int,4> >& current_embedding_tetrahedrons=current_tetrahedralized_volume->mesh.elements;
    uncollapsed_new_particles_per_all_current_particle_ids.Clean_Memory();
    uncollapsed_new_particles_per_all_current_particle_ids.Resize(cutting_particles.Number());
    new_tets_per_current.Clean_Memory();
    new_tets_per_current.Resize(current_tetrahedralized_volume->mesh.elements.m);num_new_tets=0;num_new_particles=0;
    new_parents_per_new_particle.Clean_Memory();
    new_parent_weights_per_new_particle.Clean_Memory();
    current_particle_id_per_uncollapsed_new_particle.Remove_All();
    hash_all_new_uncollapsed_particles.Remove_All();
    // duplicate particles in each tet
    for(int otet=1;otet<=regions_per_tet.m;otet++){const VECTOR<int,4>& vertices_for_tet=current_embedding_tetrahedrons(otet);
        for(int j=1;j<=regions_per_tet(otet).m;j++){
            HASHTABLE<int> particle_ids_checked;
            new_tets_per_current(otet).Append(++num_new_tets);
            int dup_tet=num_new_tets;const ARRAY<int> cutting_polygons_in_this_region=regions_per_tet(otet)(j);
            // duplicate tet nodes
            VECTOR<int,4> duplicate_particles_for_nodes;
            for(int k=1;k<=4;k++){int dup_particle_index=++num_new_particles;duplicate_particles_for_nodes(k)=dup_particle_index;
                int particle_id=cutting_particles.Particle_Id_From_Tet_Node(vertices_for_tet(k));
                uncollapsed_new_particles_per_all_current_particle_ids(particle_id).Append(dup_particle_index);
                particle_ids_checked.Insert(particle_id);
                hash_all_new_uncollapsed_particles.Insert(PAIR<int,int>(dup_tet,particle_id),dup_particle_index);
                current_particle_id_per_uncollapsed_new_particle.Append(particle_id);
                ARRAY<int> parents;parents.Append(num_new_particles);new_parents_per_new_particle.Append(parents);
                ARRAY<T> parent_weights;parent_weights.Append((T)1);new_parent_weights_per_new_particle.Append(parent_weights);}
            // look at all original nodes in this region and create a duplicate copy
            for(int k=1;k<=cutting_polygons_in_this_region.m;k++){int polygon_element_index=Cutting_Polygon_To_Element(cutting_polygons_in_this_region(k));
                const ARRAY<ARRAY<int > >& particles_for_polygon=polygon_mesh.elements(polygon_element_index);
                for(int p=1;p<=particles_for_polygon.m;p++) for(int q=1;q<=particles_for_polygon(p).m;q++){
                    int orig_intersection=particles_for_polygon(p)(q);int particle_id=cutting_particles.Particle_Id_From_Intersection(orig_intersection);
                    if(!particle_ids_checked.Set(particle_id)) continue;
                    // find parents for non-node particles
                    int uncollapsed_particle=++num_new_particles;
                    hash_all_new_uncollapsed_particles.Insert(PAIR<int,int>(dup_tet,particle_id),uncollapsed_particle);
                    uncollapsed_new_particles_per_all_current_particle_ids(particle_id).Append(uncollapsed_particle);
                    current_particle_id_per_uncollapsed_new_particle.Append(particle_id);
                    ARRAY<int> new_parents;ARRAY<T> parent_weights;
                    Compute_Parents_And_Weights_For_Nodes(orig_intersection,vertices_for_tet,duplicate_particles_for_nodes,new_parents,parent_weights);
                    new_parents_per_new_particle.Append(new_parents);
                    new_parent_weights_per_new_particle.Append(parent_weights);}}}}
    if(verbose){std::stringstream ss;ss<<"number of original tets: "<<current_tetrahedralized_volume->mesh.elements.m<<std::endl<<"number of total tets: "<<num_new_tets<<std::endl;LOG::filecout(ss.str());}
}
//#####################################################################
// Function Compute_Parents_And_Weights_For_Nodess
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Compute_Parents_And_Weights_For_Nodes(const int orig_intersection,const VECTOR<int,4>& vertices_for_tet,const VECTOR<int,4>& duplicate_particles_for_nodes,
    ARRAY<int>& new_parents,ARRAY<T>& parent_weights) const
{
    ARRAY<int> boundaries;
    for(int i=1;i<=intersection_registry->simplices_on_intersection(orig_intersection).m;i++){
        int simplex=intersection_registry->simplices_on_intersection(orig_intersection)(i);
        if(cutting_simplices->simplices(simplex).type==CUTTING_SIMPLEX<T,3>::GLOBAL_EMBEDDING_FACE)
            boundaries.Append(simplex);}
    // interior node
    if(boundaries.m==0){
        new_parents=duplicate_particles_for_nodes;
        int cutting_triangle_simplex=intersection_registry->simplices_on_intersection(orig_intersection)(1);
        TV orig_particle_weights_wrt_tet=Get_Particle_Position_In_Local_Tet_Space(orig_intersection,cutting_triangle_simplex);
        for(int i=1;i<=3;i++) parent_weights.Append(orig_particle_weights_wrt_tet(i));
        parent_weights.Append((T)1-orig_particle_weights_wrt_tet.Sum());}
    // face node
    else if(boundaries.m==1){
        int simplex=boundaries(1);const VECTOR<int,3> boundary_nodes=cutting_simplices->simplices(simplex).nodes;
        for(int i=1;i<=3;i++) new_parents.Append(duplicate_particles_for_nodes(vertices_for_tet.Find(boundary_nodes(i))));
        const VECTOR<T,2>& weights_on_simplex=intersection_registry->Get_Simplex_Weights_Of_Intersection(orig_intersection,simplex);
        for(int i=1;i<=2;i++) parent_weights.Append(weights_on_simplex(i));
        parent_weights.Append((T)1-weights_on_simplex.Sum());}
    // edge node
    else{
        ARRAY<int> shared_nodes;cutting_simplices->Shared_Nodes_On_Simplices(VECTOR<int,2>(boundaries(1),boundaries(2)),shared_nodes);assert(shared_nodes.m==2);
        for(int i=1;i<=2;i++) new_parents.Append(duplicate_particles_for_nodes(vertices_for_tet.Find(shared_nodes(i))));
        const VECTOR<T,2>& weights_on_simplex=intersection_registry->Get_Simplex_Weights_Of_Intersection(orig_intersection,boundaries(1));
        const VECTOR<int,3> nodes_on_simplex=cutting_simplices->simplices(boundaries(1)).nodes;
        for(int i=1;i<=2;i++){int index=nodes_on_simplex.Find(shared_nodes(i));assert(index>0);
            if(index<=2) parent_weights.Append(weights_on_simplex(index));
            else parent_weights.Append((T)1-weights_on_simplex.Sum());}}
}
//#####################################################################
// Function Duplicate_And_Merge_Elements
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Duplicate_And_Merge_Elements()
{
    // new structures
    union_vertices=UNION_FIND<>(num_new_particles);ARRAYS_COMPUTATIONS::Fill(union_vertices.parents,0);
    PARTICLES<TV>& new_particles=dynamic_cast<PARTICLES<TV>&>(next_tetrahedralized_volume->particles);
    new_particle_indices.Resize(num_new_particles,true,false);
    // original embedding structures
    const ARRAY<VECTOR<int,4> >& current_embedding_tetrahedrons=current_tetrahedralized_volume->mesh.elements;
    const TRIANGLE_MESH& current_embedding_faces=*current_tetrahedralized_volume->mesh.triangle_mesh;
    const PARTICLES<TV>& current_embedding_vertices=dynamic_cast<PARTICLES<TV>&>(current_tetrahedralized_volume->particles);
    // collapse faces
    for(int eface=1;eface<=current_embedding_faces.elements.m;eface++){const VECTOR<int,3>& particles_on_current_face=current_embedding_faces.elements(eface);
        ARRAY<int> tet_indices_for_this_face;
        current_tetrahedralized_volume->mesh.Tetrahedrons_On_Face(particles_on_current_face,tet_indices_for_this_face);
        // look at each possible pair of tets
        for(int p=1;p<=tet_indices_for_this_face.m;p++) for(int q=p+1;q<=tet_indices_for_this_face.m;q++){
            int tet_index_1=tet_indices_for_this_face(p);int tet_index_2=tet_indices_for_this_face(q);
            // look at each possible pair of regions
            for(int j=1;j<=regions_per_tet(tet_index_1).m;j++) for(int k=1;k<=regions_per_tet(tet_index_2).m;k++){bool join=false;
                const ARRAY<int>& region_1=regions_per_tet(tet_index_1)(j);const ARRAY<int>& region_2=regions_per_tet(tet_index_2)(k);
                for(int v=1;v<=region_1.m;v++){int polygon_element_index_1=Cutting_Polygon_To_Element(region_1(v));
                    int element_opposite=polygon_mesh.Opposite_Oriented_Element(polygon_element_index_1);if(!element_opposite) continue;
                    for(int w=1;w<=region_2.m;w++){int polygon_element_index_2=Cutting_Polygon_To_Element(region_2(w));
                        if(polygon_element_index_2==element_opposite){join=true;goto Know_Join;}}}Know_Join:;
                if(join){int dup_tet_1=new_tets_per_current(tet_index_1)(j);int dup_tet_2=new_tets_per_current(tet_index_2)(k);
                    for(int y=1;y<=3;y++){int orig_tet_node=particles_on_current_face(y);int particle_id=cutting_particles.Particle_Id_From_Tet_Node(orig_tet_node);
                        int particle_1=-1;hash_all_new_uncollapsed_particles.Get(PAIR<int,int>(dup_tet_1,particle_id),particle_1);
                        int particle_2=-1;hash_all_new_uncollapsed_particles.Get(PAIR<int,int>(dup_tet_2,particle_id),particle_2);
                        assert(particle_1>0 && particle_2>0);union_vertices.Union(particle_1,particle_2);}
                    if(verbose){std::stringstream ss;ss<<"Joining original tets "<<tet_index_1<<" and "<<tet_index_2<<" with copies "<<j<<" and "<<k<<std::endl;LOG::filecout(ss.str());}}}}}
    // collapse dup particles per original particle
    for(int opar=1;opar<=uncollapsed_new_particles_per_all_current_particle_ids.m;opar++){
        for(int i=1;i<=uncollapsed_new_particles_per_all_current_particle_ids(opar).m;i++) for(int j=i+1;j<=uncollapsed_new_particles_per_all_current_particle_ids(opar).m;j++){
            int dup_particle_1=uncollapsed_new_particles_per_all_current_particle_ids(opar)(i);int dup_particle_2=uncollapsed_new_particles_per_all_current_particle_ids(opar)(j);
            const ARRAY<int>& parents_1=new_parents_per_new_particle(dup_particle_1);
            const ARRAY<int>& parents_2=new_parents_per_new_particle(dup_particle_2);
            assert(parents_1.m==parents_2.m);bool all_parents_unioned=true;
            for(int k=1;k<=parents_1.m;k++) if(union_vertices.Find(parents_1(k))!=union_vertices.Find(parents_2(k))){
                all_parents_unioned=false;break;}
            if(all_parents_unioned) union_vertices.Union(dup_particle_1,dup_particle_2);}}
    // make the particles, which are now all collapsed, and set their positions
    // assign original cutting particle index per final duplicate particles (-1 if it the duplicate is not a cutting particle)
    int dup_particle_count=0;for(int i=1;i<=union_vertices.parents.m;i++) if(union_vertices.Is_Root(i)) dup_particle_count++;
    new_particles.array_collection->Preallocate(new_particles.array_collection->Size()+dup_particle_count);
    previous_particle_index_per_new_particle_index.Resize(new_particles.array_collection->Size()+dup_particle_count,false,false);
    ARRAYS_COMPUTATIONS::Fill(previous_particle_index_per_new_particle_index,0);
    current_particle_id_per_collapsed_new_particle.Remove_All();
    current_particle_id_per_collapsed_new_particle.Preallocate(dup_particle_count);
    old_particle_per_new_collapsed_particle.Remove_All();
    old_particle_per_new_collapsed_particle.Preallocate(dup_particle_count);
    new_collapsed_tet_particle_per_old_tet_particle.Clean_Memory();
    new_collapsed_tet_particle_per_old_tet_particle.Resize(current_tetrahedralized_volume->particles.array_collection->Size());
    for(int i=1;i<=union_vertices.parents.m;i++) if(union_vertices.Is_Root(i)){int particle_index=new_particle_indices(i)=new_particles.array_collection->Add_Element();
        int current_particle_id_for_this_new_particle=current_particle_id_per_uncollapsed_new_particle(i);
        new_particles.X(particle_index)=Compute_World_Space_Position_Of_Uncollapsed_Particle(i);
        current_particle_id_per_collapsed_new_particle.Append(current_particle_id_for_this_new_particle);
        // assign old particle per new particle, if the old partiele existed (some point ids were just made, but old point ids correspond to particles)
        if(current_particle_id_for_this_new_particle<=current_tetrahedralized_volume->particles.array_collection->Size()){
            // the current particle index is either the tet node or the intersection, or both
            int old_tet_node=cutting_particles.tet_node_indices(current_particle_id_for_this_new_particle);
            old_particle_per_new_collapsed_particle.Append(old_tet_node);
            if(old_tet_node)
                new_collapsed_tet_particle_per_old_tet_particle(old_tet_node).Append(particle_index);}
        else old_particle_per_new_collapsed_particle.Append(0);
        // set mass/velocity of tet nodes
        if(cutting_particles.particle_ids_types(current_particle_id_for_this_new_particle)!=CUTTING_PARTICLES::INTERSECTION_ID){
            int tet_node=cutting_particles.tet_node_indices(current_particle_id_for_this_new_particle);if(!tet_node) PHYSBAM_FATAL_ERROR();
            previous_particle_index_per_new_particle_index(particle_index)=tet_node;
            new_particles.V(particle_index)=current_embedding_vertices.V(tet_node);
            new_particles.mass(particle_index)=current_embedding_vertices.mass(tet_node);}}
    // make new triangulated area and create unique elements, deleting identical ones from duplicated_triangles_per_original as we go
    next_tetrahedralized_volume->Update_Number_Nodes();
    next_tetrahedralized_volume->mesh.Initialize_Incident_Elements();
    ARRAY<bool> remove_flags;
    for(int otet=1;otet<=new_tets_per_current.m;otet++){const VECTOR<int,4>& original_node=current_embedding_tetrahedrons(otet);
        ARRAY<VECTOR<int,4> > tets_to_add;
        for(int i=new_tets_per_current(otet).m;i>=1;i--){int dtet_1=new_tets_per_current(otet)(i);
            // make mesh, merge duplicates
            VECTOR<int,4> collapsed_dup_particles_1;for(int k=1;k<=4;k++){int particle_index=-1;
                int particle_id=cutting_particles.Particle_Id_From_Tet_Node(original_node(k));
                bool found=hash_all_new_uncollapsed_particles.Get(PAIR<int,int>(dtet_1,particle_id),particle_index);if(!found) PHYSBAM_FATAL_ERROR();
                collapsed_dup_particles_1(k)=new_particle_indices(union_vertices.Find(particle_index));
                assert(collapsed_dup_particles_1(k)>0 && collapsed_dup_particles_1(k)<=new_particles.array_collection->Size());}
            bool create=true;
            for(int j=i-1;j>=1;j--){int dtet_2=new_tets_per_current(otet)(j);VECTOR<int,4> collapsed_dup_particles_2;
                for(int k=1;k<=4;k++){int particle_index=-1;int poind_id=cutting_particles.Particle_Id_From_Tet_Node(original_node(k));
                    bool found=hash_all_new_uncollapsed_particles.Get(PAIR<int,int>(dtet_2,poind_id),particle_index);if(!found) PHYSBAM_FATAL_ERROR();
                    collapsed_dup_particles_2(k)=new_particle_indices(union_vertices.Find(particle_index));
                    assert(collapsed_dup_particles_2(k)>0 && collapsed_dup_particles_2(k)<=new_particles.array_collection->Size());}
                if(collapsed_dup_particles_1==collapsed_dup_particles_2){create=false;
                    for(int p=1;p<=regions_per_tet(otet)(i).m;p++){int polygon_element_index=Cutting_Polygon_To_Element(regions_per_tet(otet)(i)(p));
                        const ARRAY<ARRAY<int> >& polygon_nodes=polygon_mesh.elements(polygon_element_index);
                        for(int q=1;q<=polygon_nodes.m;q++) for(int r=1;r<=polygon_nodes(q).m;r++){
                            int original_intersection=polygon_nodes(q)(r);int particle_id=cutting_particles.Particle_Id_From_Intersection(original_intersection);
                            int particle_index=-1;
                            if(!hash_all_new_uncollapsed_particles.Get(PAIR<int,int>(dtet_1,particle_id),particle_index)) continue;
                            hash_all_new_uncollapsed_particles.Delete(PAIR<int,int>(dtet_1,particle_id));
                            if(!hash_all_new_uncollapsed_particles.Contains(PAIR<int,int>(dtet_2,particle_id))) // TODO: check if we can use Set here
                                hash_all_new_uncollapsed_particles.Insert(PAIR<int,int>(dtet_2,particle_id),particle_index);}}
                    new_tets_per_current(otet).Remove_Index(i);
                    regions_per_tet(otet)(j).Append_Elements(regions_per_tet(otet)(i));
                    regions_per_tet(otet).Remove_Index(i);
                    if(verbose){std::stringstream ss;ss<<"Deleting duplicate tet "<<otet<<" with copy ids "<<i<<" and "<<j<<std::endl;LOG::filecout(ss.str());}break;}}
            if(create){assert(!next_tetrahedralized_volume->mesh.Simplex(collapsed_dup_particles_1));
                tets_to_add.Append(collapsed_dup_particles_1);}}
        ARRAYS_COMPUTATIONS::Reverse_In_Place(tets_to_add);
        next_tetrahedralized_volume->mesh.elements.Append_Elements(tets_to_add);}
    // fix embedding map
    final_parent_weights_per_new_particle.Remove_All();final_parents_per_new_particle.Remove_All();
    for(int pcdup=1;pcdup<=new_parents_per_new_particle.m;pcdup++){
        if(union_vertices.Is_Root(pcdup)){ARRAY<int>& parents=new_parents_per_new_particle(pcdup);ARRAY<int> final_parents(parents.m);
            for(int i=1;i<=parents.m;i++){
                final_parents(i)=new_particle_indices(union_vertices.Find(parents(i)));assert(final_parents(i)<=next_tetrahedralized_volume->particles.array_collection->Size());}
            final_parents_per_new_particle.Append(final_parents);
            final_parent_weights_per_new_particle.Append(new_parent_weights_per_new_particle(pcdup));}}
    // fix dup tet numbers
    dup_tet_before_to_after_collapse.Remove_All();dup_tet_after_to_before_collapse.Remove_All();int count=0;
    for(int i=1;i<=new_tets_per_current.m;i++) for(int j=1;j<=new_tets_per_current(i).m;j++){int num=++count;
        dup_tet_before_to_after_collapse.Insert(new_tets_per_current(i)(j),num);
        dup_tet_after_to_before_collapse.Insert(num,new_tets_per_current(i)(j));}
}
//#####################################################################
// Function Compute_World_Space_Position_Of_Uncollapsed_Particle
//#####################################################################
template<class T> VECTOR<T,3> CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Compute_World_Space_Position_Of_Uncollapsed_Particle(const int uncollapsed_particle) const
{
    const GEOMETRY_PARTICLES<TV>& current_embedding_particles=current_tetrahedralized_volume->particles;
    const ARRAY<int>& uncollapsed_parents=new_parents_per_new_particle(uncollapsed_particle);
    const ARRAY<T>& uncollapsed_parent_weights=new_parent_weights_per_new_particle(uncollapsed_particle);
    assert(uncollapsed_parents.m==uncollapsed_parent_weights.m);TV position=TV();
    for(int i=1;i<=uncollapsed_parents.m;i++){
        int original_particle_id=current_particle_id_per_uncollapsed_new_particle(uncollapsed_parents(i));
        assert(cutting_particles.particle_ids_types(original_particle_id)!=CUTTING_PARTICLES::INTERSECTION_ID);
        int original_tet_node=cutting_particles.tet_node_indices(original_particle_id);
        assert(original_tet_node>=0 && original_tet_node<=current_embedding_particles.array_collection->Size());
        position+=current_embedding_particles.X(original_tet_node)*uncollapsed_parent_weights(i);}
    return position;
}
//#####################################################################
// Function Create_Boundary_Surface
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Create_Boundary_Surface()
{
    const GEOMETRY_PARTICLES<TV>& duplicate_particles=next_tetrahedralized_volume->particles;
    POLYGON_MESH boundary_polygons;
    ARRAY<int> original_cutting_polygon_per_new_boundary_polygon;
    final_duplicated_boundary_mesh.elements.Remove_All();
    final_duplicated_boundary_mesh.Set_Number_Nodes(duplicate_particles.array_collection->Size());
    // make duplicate polygons
    for(int otet=1;otet<=new_tets_per_current.m;otet++){
        for(int i=1;i<=regions_per_tet(otet).m;i++){
            int dtet=new_tets_per_current(otet)(i);
            const ARRAY<int>& cutting_polygons_for_region=regions_per_tet(otet)(i);
            for(int j=1;j<=cutting_polygons_for_region.m;j++){
                const CUTTING_POLYGON& cutting_polygon=current_cutting_polygons(cutting_polygons_for_region(j));
                CUTTING_POLYGON::POLYGON_TYPE polygon_type=cutting_polygon.polygon_type;
                if(polygon_type==CUTTING_POLYGON::FACE_INTERIOR) continue;
                original_cutting_polygon_per_new_boundary_polygon.Append(cutting_polygons_for_region(j));
                const ARRAY<ARRAY<int> >& polygon_particles=polygon_mesh.elements(cutting_polygon.polygon_index);
                boundary_polygons.elements.Append(ARRAY<ARRAY<int> >());
                ARRAY<ARRAY<int> >& new_polygon_particles=boundary_polygons.elements.Last();
                for(int k=1;k<=polygon_particles.m;k++){
                    new_polygon_particles.Append(ARRAY<int>());
                    ARRAY<int>& new_polygon=new_polygon_particles.Last();
                    for(int p=1;p<=polygon_particles(k).m;p++){
                        int original_intersection=polygon_particles(k)(p);int particle_id=cutting_particles.Particle_Id_From_Intersection(original_intersection);
                        int real_dup_particle=-1;
                        if(!hash_all_new_uncollapsed_particles.Get(PAIR<int,int>(dtet,particle_id),real_dup_particle)) PHYSBAM_FATAL_ERROR();
                        real_dup_particle=new_particle_indices(union_vertices.Find(real_dup_particle));
                        new_polygon.Append(real_dup_particle);}}}}}
    // set up mesh
    boundary_polygons.Set_Number_Nodes(duplicate_particles.array_collection->Size());
    boundary_polygons.Initialize_Segment_Mesh();
    boundary_polygons.segment_mesh->Initialize_Incident_Elements();
    boundary_polygons.Initialize_Element_Oriented_Edges();
    boundary_polygons.Initialize_Edge_Elements();
    // triangulate each polygon, adding the result to the final boundary triangle mesh
    for(int new_polygon=1;new_polygon<=boundary_polygons.elements.m;new_polygon++){
        if(boundary_polygons.Opposite_Oriented_Element(new_polygon)) continue;
        int original_cutting_polygon=original_cutting_polygon_per_new_boundary_polygon(new_polygon);
        const CUTTING_POLYGON& cutting_polygon=current_cutting_polygons(original_cutting_polygon);
        int child_simplex_for_polygon=cutting_polygon.simplex_owner;
        // work on the polygon particles
        ARRAY<int> real_particle_per_fake_particle;ARRAY<VECTOR<T,2> > fake_positions;ARRAY<ARRAY<int> > fake_polygon;
        const ARRAY<ARRAY<int> >& polygon_particles=boundary_polygons.elements(new_polygon);
        for(int k=1;k<=polygon_particles.m;k++) if(polygon_particles(k).m>=3){
            const ARRAY<int>& particles_for_piece=polygon_particles(k);
            for(int p=1;p<=particles_for_piece.m;p++){
                // this particle index is a real new particle (after duplication) with an actual position into the particles array
                int dup_particle=particles_for_piece(p);real_particle_per_fake_particle.Append(dup_particle);
                int particle_id=current_particle_id_per_collapsed_new_particle(dup_particle);
                assert(cutting_particles.particle_ids_types(particle_id)!=CUTTING_PARTICLES::TET_NODE_ID);
                int original_particle_in_intersection_registry=cutting_particles.intersection_indices(particle_id);
                VECTOR<T,2> bary=intersection_registry->Get_Simplex_Weights_Of_Intersection(
                    original_particle_in_intersection_registry,child_simplex_for_polygon);
                if(cutting_polygon.flipped) exchange(bary[1],bary[2]);
                fake_positions.Append(VECTOR<T,2>(clamp<T>(bary[1],0,1),clamp<T>(bary[2],0,1)));}
            fake_polygon.Append(ARRAY<int>());
            ARRAY<int>& fake_loop=fake_polygon.Last();
            fake_loop.Preallocate(polygon_particles(k).m);
            for(int p=1;p<=polygon_particles(k).m;p++)
                fake_loop.Append(real_particle_per_fake_particle.Find(polygon_particles(k)(p)));}
        ARRAY<VECTOR<int,3> > triangles;
        POLYGONAL_TRIANGULATION<T>::Triangulate_Nonconvex_Nonsimple_Polygon(fake_positions,fake_polygon,triangles);
        for(int k=1;k<=triangles.m;k++){
            final_duplicated_boundary_mesh.elements.Append(VECTOR<int,3>(
                real_particle_per_fake_particle(triangles(k)[1]),real_particle_per_fake_particle(triangles(k)[3]),real_particle_per_fake_particle(triangles(k)[2])));}}
}
//#####################################################################
// Function Build_New_Cutting_Simplices_And_Intersection_Registry
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Build_New_Cutting_Simplices_And_Intersection_Registry()
{
    LOG::SCOPE scope("BUILD_NEW_CUTTING","Build New Cutting Simplices and Intersection Registry");
    // presumably regions per tet has already been updated in duplicate tets
    // initialize new volume
    next_tetrahedralized_volume->mesh.Delete_Auxiliary_Structures();
    next_tetrahedralized_volume->mesh.Initialize_Triangle_Mesh();
    // rebuild new structures
    CUTTING_SIMPLICES<T,3>* new_cutting_simplices=new CUTTING_SIMPLICES<T,3>();
    INTERSECTION_REGISTRY<T,3>* new_intersection_registry=new INTERSECTION_REGISTRY<T,3>(*new_cutting_simplices);
    ARRAY<ARRAY<int> > simplices_per_new_tet(next_tetrahedralized_volume->mesh.elements.m);
    ARRAY<CUTTING_POLYGON> new_cutting_polygons;
    ARRAY<ARRAY<ARRAY<int> > > new_regions_per_tet(next_tetrahedralized_volume->mesh.elements.m);
    POLYGON_MESH new_polygon_mesh;
    new_polygon_mesh.elements.Preallocate(polygon_mesh.elements.m);
    // take care of faces and boundaries
    HASHTABLE<VECTOR<int,3>,int> volume_triangle_mesh_element_to_new_boundary;
    HASHTABLE<VECTOR<int,3>,VECTOR<int,3> > new_sorted_to_original_boundary_nodes;
    LOG::Time(" Adding global embedding faces");
    Add_Boundaries_To_New_Cutting_Simplices(*new_cutting_simplices,volume_triangle_mesh_element_to_new_boundary,new_sorted_to_original_boundary_nodes);
    HASHTABLE<int,int> old_cut_to_new_cut;
    LOG::Time(" Adding global cut faces");
    Add_Cuts_To_New_Cutting_Simplices(*new_cutting_simplices,old_cut_to_new_cut);
    // loop over each dup tet and fill structures 
    HASHTABLE<PAIR<int,int>,int> old_child_simplex_to_new_child_simplex;
    HASHTABLE<int,int> new_child_simplex_to_old_child_simplex;
    ARRAY<ARRAY<int> > new_cutting_polygons_per_cutting_simplex;
    LOG::Time("Adding local faces");
    Add_New_Child_Simplices_And_Create_New_Polygon_Mesh(new_polygon_mesh,*new_cutting_simplices,simplices_per_new_tet,new_cutting_polygons,new_regions_per_tet,
        old_cut_to_new_cut,volume_triangle_mesh_element_to_new_boundary,old_child_simplex_to_new_child_simplex,
        new_child_simplex_to_old_child_simplex,new_sorted_to_original_boundary_nodes,new_cutting_polygons_per_cutting_simplex);
    // fix intersection registry
    LOG::Time("Building new intersection registry");
    CUTTING_PARTICLES new_cutting_particles;
    Build_New_Intersection_Registry(new_polygon_mesh,new_cutting_polygons,old_cut_to_new_cut,volume_triangle_mesh_element_to_new_boundary,
        old_child_simplex_to_new_child_simplex,*new_cutting_simplices,*new_intersection_registry,new_child_simplex_to_old_child_simplex,new_cutting_particles);
    if(verbose){
        std::stringstream ss;
        ss<<"Old cuttings simplices"<<std::endl;cutting_simplices->Print();
        ss<<"New cuttings simplices"<<std::endl;new_cutting_simplices->Print();
        ss<<"Intersection registry before: "<<std::endl;intersection_registry->Print();
        ss<<"Intersection registry after: "<<std::endl;new_intersection_registry->Print();
        ss<<"Old cutting polygons: "<<std::endl;
        for(int i=1;i<=current_cutting_polygons.m;i++)
            ss<<"Polygon index = "<<current_cutting_polygons(i).polygon_index
                <<", simplex owner = "<<current_cutting_polygons(i).simplex_owner
                <<", flipped = "<<(current_cutting_polygons(i).flipped?"true":"false")
                <<", particles = "<<polygon_mesh.elements(current_cutting_polygons(i).polygon_index)<<std::endl;
        ss<<"New cutting polygons: "<<std::endl;
        for(int i=1;i<=current_cutting_polygons.m;i++)
            ss<<"Polygon index = "<<new_cutting_polygons(i).polygon_index
                <<", simplex owner = "<<new_cutting_polygons(i).simplex_owner
                <<", flipped = "<<(new_cutting_polygons(i).flipped?"true":"false")
                <<", particles = "<<new_polygon_mesh.elements(new_cutting_polygons(i).polygon_index)<<std::endl;
        ss<<"Regions per tet before: "<<regions_per_tet<<std::endl;
        ss<<"Regions per tet after: "<<new_regions_per_tet<<std::endl;
        for(int i=1;i<=cutting_polygons_per_cutting_simplex.m;i++)
            ss<<"Cutting polygons for simplex "<<i<<": "<<cutting_polygons_per_cutting_simplex(i)<<std::endl;
        LOG::filecout(ss.str());}
    // recopy over all data for next iteration
    LOG::Time("Replacing old structures with new ones");
    exchange(current_tetrahedralized_volume,next_tetrahedralized_volume);
    exchange(cutting_simplices,new_cutting_simplices);delete new_cutting_simplices;
    exchange(intersection_registry,new_intersection_registry);delete new_intersection_registry;
    simplices_per_current_tet.Exchange(simplices_per_new_tet);
    current_cutting_polygons.Exchange(new_cutting_polygons);
    regions_per_tet.Exchange(new_regions_per_tet);
    sorted_to_original_boundary_nodes=new_sorted_to_original_boundary_nodes;
    cutting_particles=new_cutting_particles;
    cutting_polygons_per_cutting_simplex.Exchange(new_cutting_polygons_per_cutting_simplex);
    polygon_mesh.Initialize_Mesh(new_polygon_mesh);
}
//#####################################################################
// Function Add_Cuts_To_New_Cutting_Simplices
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Add_Cuts_To_New_Cutting_Simplices(CUTTING_SIMPLICES<T,3>& new_cutting_simplices,HASHTABLE<int,int>& old_cut_to_new_cut) const
{
    for(int i=1;i<=cutting_simplices->simplices.m;i++) if(cutting_simplices->simplices(i).type==CUTTING_SIMPLEX<T,3>::GLOBAL_CUT_FACE){
        int new_index=new_cutting_simplices.Add_Simplex(cutting_simplices->simplices(i).nodes,CUTTING_SIMPLEX<T,3>::GLOBAL_CUT_FACE);
        old_cut_to_new_cut.Insert(i,new_index);}
}
//#####################################################################
// Function Add_Boundaries_To_New_Cutting_Simplices
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Add_Boundaries_To_New_Cutting_Simplices(CUTTING_SIMPLICES<T,3>& new_cutting_simplices,HASHTABLE<VECTOR<int,3>,int>& volume_triangle_mesh_element_to_new_boundary,
    HASHTABLE<VECTOR<int,3>,VECTOR<int,3> >& new_sorted_to_original_boundary_nodes) const
{
    const TRIANGLE_MESH& embedding_faces=*next_tetrahedralized_volume->mesh.triangle_mesh;
    for(int face=1;face<=embedding_faces.elements.m;face++){const VECTOR<int,3>& face_vertices=embedding_faces.elements(face);
        // unscramble face vertices to match old boundary that this one came from
        VECTOR<int,3> previous_face_vertices;
        for(int i=1;i<=3;i++){int particle_id=current_particle_id_per_collapsed_new_particle(face_vertices(i));
            assert(cutting_particles.particle_ids_types(particle_id)!=CUTTING_PARTICLES::INTERSECTION_ID);
            previous_face_vertices(i)=cutting_particles.tet_node_indices(particle_id);}
        VECTOR<int,3> correct_previous_face_vertices,correct_new_face_vertices;
        if(!sorted_to_original_boundary_nodes.Get(previous_face_vertices.Sorted(),correct_previous_face_vertices)) PHYSBAM_FATAL_ERROR();
        for(int i=1;i<=3;i++){
            int index_into_unsorted=previous_face_vertices.Find(correct_previous_face_vertices(i));
            correct_new_face_vertices(i)=face_vertices(index_into_unsorted);}
        // store info for next iteration
        if(!new_sorted_to_original_boundary_nodes.Contains(correct_new_face_vertices.Sorted())) // TODO: do we really need the unsorted version stored in the hashtable?
            new_sorted_to_original_boundary_nodes.Insert(correct_new_face_vertices.Sorted(),correct_new_face_vertices);
        // make simplex
        int parent=new_cutting_simplices.Add_Simplex(correct_new_face_vertices,CUTTING_SIMPLEX<T,3>::GLOBAL_EMBEDDING_FACE);
        if(!volume_triangle_mesh_element_to_new_boundary.Contains(correct_new_face_vertices)) // TODO: can we use Set here?
            volume_triangle_mesh_element_to_new_boundary.Insert(correct_new_face_vertices,parent);}
}
//#####################################################################
// Function Add_New_Child_Simplices_And_Create_New_Polygon_Mesh
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Add_New_Child_Simplices_And_Create_New_Polygon_Mesh(POLYGON_MESH& new_polygon_mesh,CUTTING_SIMPLICES<T,3>& new_cutting_simplices,ARRAY<ARRAY<int> >& simplices_per_new_tet,
    ARRAY<CUTTING_POLYGON>& new_cutting_polygons,ARRAY<ARRAY<ARRAY<int> > >& new_regions_per_tet,HASHTABLE<int,int>& old_cut_to_new_cut,
    HASHTABLE<VECTOR<int,3>,int>& volume_triangle_mesh_element_to_new_boundary,HASHTABLE<PAIR<int,int>,int>& old_child_simplex_to_new_child_simplex,
    HASHTABLE<int,int>& new_child_simplex_to_old_child_simplex,HASHTABLE<VECTOR<int,3>,VECTOR<int,3> >& new_sorted_to_original_boundary_nodes,
    ARRAY<ARRAY<int> >& new_cutting_polygons_per_cutting_simplex) const
{
    const ARRAY<VECTOR<int,4> >& new_embedding_tetrahedrons=next_tetrahedralized_volume->mesh.elements;
    for(int otet=1;otet<=regions_per_tet.m;otet++){
        const ARRAY<int>& unique_child_simplices=simplices_per_current_tet(otet);
        // for each child simplex in this region, add it as a new simplex
        for(int region_index=1;region_index<=regions_per_tet(otet).m;region_index++){
            const ARRAY<int>& cutting_polygons_for_region=regions_per_tet(otet)(region_index);
            int dtet=new_tets_per_current(otet)(region_index);
            int actual_dtet;
            if(!dup_tet_before_to_after_collapse.Get(dtet,actual_dtet)) PHYSBAM_FATAL_ERROR();
            assert(actual_dtet<=new_embedding_tetrahedrons.m);
            ARRAY<int>& simplices_in_actual_dtet=simplices_per_new_tet(actual_dtet);
            simplices_in_actual_dtet.Preallocate(simplices_in_actual_dtet.m+unique_child_simplices.m);
            for(int j=1;j<=unique_child_simplices.m;j++){
                int simplex_index=unique_child_simplices(j);
                const CUTTING_SIMPLEX<T,3>& cutting_simplex=cutting_simplices->simplices(simplex_index);
                int new_child;
                if(cutting_simplex.type==CUTTING_SIMPLEX<T,3>::LOCAL_EMBEDDING_FACE){
                    VECTOR<int,3> new_nodes;
                    const VECTOR<int,3>& old_nodes=cutting_simplex.nodes;
                    for(int k=1;k<=3;k++){
                        int particle_id=cutting_particles.Particle_Id_From_Tet_Node(old_nodes(k));
                        assert(cutting_particles.particle_ids_types(particle_id)!=CUTTING_PARTICLES::INTERSECTION_ID);
                        if(!hash_all_new_uncollapsed_particles.Get(PAIR<int,int>(dtet,particle_id),new_nodes(k))) PHYSBAM_FATAL_ERROR();
                        new_nodes(k)=new_particle_indices(union_vertices.Find(new_nodes(k)));}
                    int new_parent;
                    if(!volume_triangle_mesh_element_to_new_boundary.Get(new_nodes,new_parent)) PHYSBAM_FATAL_ERROR();
                    VECTOR<TV,3> weights=Get_Face_Weights_From_Tet(new_nodes,new_embedding_tetrahedrons(actual_dtet));
                    new_child=new_cutting_simplices.Add_Simplex(new_nodes,CUTTING_SIMPLEX<T,3>::LOCAL_EMBEDDING_FACE,weights,new_parent,actual_dtet);}
                else if(cutting_simplex.type==CUTTING_SIMPLEX<T,3>::LOCAL_CUT_FACE){
                    int new_parent;
                    if(!old_cut_to_new_cut.Get(cutting_simplices->simplices(simplex_index).parent, new_parent)) PHYSBAM_FATAL_ERROR();
                    new_child=new_cutting_simplices.simplices.Append(cutting_simplex);
                    CUTTING_SIMPLEX<T,3>& new_child_simplex=new_cutting_simplices.simplices.Last();
                    new_child_simplex.nodes=new_cutting_simplices.simplices(new_parent).nodes;
                    new_child_simplex.parent=new_parent;
                    new_child_simplex.element_owner=actual_dtet;}
                else PHYSBAM_FATAL_ERROR();
                assert(!simplices_in_actual_dtet.Contains(new_child));
                simplices_in_actual_dtet.Append(new_child);
                old_child_simplex_to_new_child_simplex.Insert(PAIR<int,int>(actual_dtet,simplex_index),new_child);
                new_child_simplex_to_old_child_simplex.Insert(new_child,simplex_index);}
            // add cutting polygons
            new_cutting_polygons_per_cutting_simplex.Resize(new_cutting_simplices.simplices.m);
            new_regions_per_tet(actual_dtet).Append(ARRAY<int>());
            ARRAY<int>& new_cutting_polygons_in_region=new_regions_per_tet(actual_dtet)(1);
            for(int j=1;j<=cutting_polygons_for_region.m;j++){
                const CUTTING_POLYGON& cutting_polygon=current_cutting_polygons(cutting_polygons_for_region(j));
                const ARRAY<ARRAY<int> >& polygon_particles=polygon_mesh.elements(cutting_polygon.polygon_index);
                int new_polygon_index=new_polygon_mesh.elements.Append(ARRAY<ARRAY<int> >());
                ARRAY<ARRAY<int> >& new_polygon_particles=new_polygon_mesh.elements.Last();
                for(int k=1;k<=polygon_particles.m;k++){
                    new_polygon_particles.Append(ARRAY<int>());
                    ARRAY<int>& new_polygon=new_polygon_particles.Last();
                    for(int p=1;p<=polygon_particles(k).m;p++){
                        int original_intersection=polygon_particles(k)(p),particle_id=cutting_particles.Particle_Id_From_Intersection(original_intersection);
                        assert(cutting_particles.particle_ids_types(particle_id)!=CUTTING_PARTICLES::TET_NODE_ID);
                        int real_dup_particle=-1;
                        if(!hash_all_new_uncollapsed_particles.Get(PAIR<int,int>(dtet,particle_id),real_dup_particle)) PHYSBAM_FATAL_ERROR();
                        real_dup_particle=new_particle_indices(union_vertices.Find(real_dup_particle));
                        new_polygon.Append(real_dup_particle);}}
                // make CUTTING_POLYGON
                int new_child_simplex=old_child_simplex_to_new_child_simplex.Get(Tuple(actual_dtet,cutting_polygon.simplex_owner));
                int new_cutting_polygon_index=new_cutting_polygons.Append(
                    CUTTING_POLYGON(new_polygon_index,new_child_simplex,cutting_polygon.flipped,cutting_polygon.polygon_type));
                new_cutting_polygons_per_cutting_simplex(new_child_simplex).Append(new_cutting_polygon_index);
                new_cutting_polygons_in_region.Append(new_cutting_polygon_index);}}}
}
//#####################################################################
// Function Build_New_Intersection_Registry
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Build_New_Intersection_Registry(POLYGON_MESH& new_polygon_mesh,ARRAY<CUTTING_POLYGON>& new_cutting_polygons,HASHTABLE<int,int>& old_cut_to_new_cut,
    HASHTABLE<VECTOR<int,3>,int>& volume_triangle_mesh_element_to_new_boundary,HASHTABLE<PAIR<int,int>,int>& old_child_simplex_to_new_child_simplex,CUTTING_SIMPLICES<T,3>& new_cutting_simplices,
    INTERSECTION_REGISTRY<T,3>& new_intersection_registry,const HASHTABLE<int,int>& new_child_simplex_to_old_child_simplex,CUTTING_PARTICLES& new_cutting_particles) const
{
    // get all tet nodes that will be added (possibly as intersections)
    const ARRAY<VECTOR<int,4> >& tet_elements=next_tetrahedralized_volume->mesh.elements;
    ARRAY<int> all_tet_nodes;tet_elements.Flattened().Get_Unique(all_tet_nodes);
    if(verbose){std::stringstream ss;ss<<"All new tet nodes: "<<all_tet_nodes<<std::endl;LOG::filecout(ss.str());}
    ARRAY<bool> tet_particle_flags(next_tetrahedralized_volume->particles.array_collection->Size());
    ARRAYS_COMPUTATIONS::Fill(tet_particle_flags,false);
    INDIRECT_ARRAY<ARRAY<bool>,ARRAY<int>&> subset=tet_particle_flags.Subset(all_tet_nodes);ARRAYS_COMPUTATIONS::Fill(subset,true);
    // build intersections
    for(int cp=1;cp<=new_cutting_polygons.m;cp++){const CUTTING_POLYGON& cutting_polygon=new_cutting_polygons(cp);
        int new_child_simplex=cutting_polygon.simplex_owner;int old_child_simplex;
        int actual_dtet=new_cutting_simplices.simplices(new_child_simplex).element_owner;
        int dtet=-1;bool found_dtet=dup_tet_after_to_before_collapse.Get(actual_dtet,dtet);if(!found_dtet) PHYSBAM_FATAL_ERROR();
        bool found=new_child_simplex_to_old_child_simplex.Get(new_child_simplex,old_child_simplex);if(!found) PHYSBAM_FATAL_ERROR();
        ARRAY<ARRAY<int> > particles_on_polygon=new_polygon_mesh.elements(cutting_polygon.polygon_index);
        for(int i=1;i<=particles_on_polygon.m;i++) for(int j=1;j<=particles_on_polygon(i).m;j++){
            int new_particle=particles_on_polygon(i)(j);
            int particle_id=current_particle_id_per_collapsed_new_particle(new_particle);
            assert(cutting_particles.particle_ids_types(particle_id)!=CUTTING_PARTICLES::TET_NODE_ID);
            int old_intersection=cutting_particles.intersection_indices(particle_id);
            // get list of new simplices corresponding to all old ones (except for boundaries which might be in a different tet)
            ARRAY<VECTOR<T,2> > all_weights_for_new_simplices;ARRAY<int> new_simplices;
            const ARRAY<int>& old_simplices_on_old_particle=intersection_registry->simplices_on_intersection(old_intersection);
            for(int k=1;k<=old_simplices_on_old_particle.m;k++){int old_simplex=old_simplices_on_old_particle(k);
                int new_simplex=-1;bool found=false;
                if(cutting_simplices->simplices(old_simplex).type==CUTTING_SIMPLEX<T,3>::GLOBAL_CUT_FACE){
                    found=old_cut_to_new_cut.Get(old_simplex,new_simplex);if(!found) PHYSBAM_FATAL_ERROR();}
                else if(cutting_simplices->simplices(old_simplex).type==CUTTING_SIMPLEX<T,3>::LOCAL_CUT_FACE){
                    found=old_child_simplex_to_new_child_simplex.Get(PAIR<int,int>(actual_dtet,old_simplex),new_simplex);if(!found) PHYSBAM_FATAL_ERROR();}
                else if(cutting_simplices->simplices(old_simplex).type==CUTTING_SIMPLEX<T,3>::GLOBAL_EMBEDDING_FACE){
                    VECTOR<int,3> new_nodes;VECTOR<int,3> old_nodes=cutting_simplices->simplices(old_simplex).nodes;bool found_all_nodes=true;
                    for(int l=1;l<=3;l++){int particle_id=cutting_particles.Particle_Id_From_Tet_Node(old_nodes(l));
                        assert(cutting_particles.particle_ids_types(particle_id)!=CUTTING_PARTICLES::INTERSECTION_ID);
                        found_all_nodes&=hash_all_new_uncollapsed_particles.Get(PAIR<int,int>(dtet,particle_id),new_nodes(l));
                        if(!found_all_nodes) break;
                        else new_nodes(l)=new_particle_indices(union_vertices.Find(new_nodes(l)));}
                    if(found_all_nodes){
                        found=volume_triangle_mesh_element_to_new_boundary.Get(new_nodes,new_simplex);if(!found) PHYSBAM_FATAL_ERROR();}}
                if(new_simplex<0) continue;
                VECTOR<T,2> weights_for_old_particle_on_old_simplex=intersection_registry->Get_Simplex_Weights_Of_Intersection(old_intersection,old_simplex);
                all_weights_for_new_simplices.Append(weights_for_old_particle_on_old_simplex);
                new_simplices.Append(new_simplex);}
            // register and add to point id structure
            new_intersection_registry.Register_Intersection(new_simplices,all_weights_for_new_simplices,new_particle);
            bool is_tet_node=(cutting_particles.particle_ids_types(particle_id)!=CUTTING_PARTICLES::INTERSECTION_ID);
            if(is_tet_node){
                const VECTOR<int,4>& new_vertics_in_tet=tet_elements(actual_dtet);bool found=false;int new_tet_node=0;
                for(int k=1;k<=4;k++) if(current_particle_id_per_collapsed_new_particle(new_vertics_in_tet(k))==particle_id){found=true;new_tet_node=new_vertics_in_tet(k);}
                if(!found) PHYSBAM_FATAL_ERROR();
                if(tet_particle_flags(new_tet_node)){
                    new_cutting_particles.Add_Tet_Node_And_Intersection_Id(new_tet_node,new_particle);
                    tet_particle_flags(new_tet_node)=false;}}
            else new_cutting_particles.Add_Intersection_Id(new_particle);}}
    // build boundary point ids
    for(int i=1;i<=tet_particle_flags.m;i++) if(tet_particle_flags(i)) new_cutting_particles.Add_Tet_Node_Id(i);
    // print
    for(int i=1;i<=new_intersection_registry.simplices_on_intersection.m;i++){
        if(verbose){std::stringstream ss;ss<<"Intersection "<<i<<" has simplices "<<new_intersection_registry.simplices_on_intersection(i)<<std::endl;LOG::filecout(ss.str());}
        if(new_intersection_registry.simplices_on_intersection(i).m==1 || new_intersection_registry.simplices_on_intersection(i).m==2) PHYSBAM_FATAL_ERROR();}
}
//#####################################################################
// Function Cutting_Polygon_To_Element
//#####################################################################
template<class T> int CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Cutting_Polygon_To_Element(const int index) const
{
    return current_cutting_polygons(index).polygon_index;
}
//#####################################################################
// Function Get_Particles_On_Simplices
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Get_Particles_On_Simplices(const VECTOR<int,2>& simplices,ARRAY<int>& particles) const
{
    bool is_all_cutting_simplices=true;
    for(int i=1;i<=simplices.m;i++) if(cutting_simplices->simplices(simplices(i)).type!=CUTTING_SIMPLEX<T,3>::LOCAL_CUT_FACE){is_all_cutting_simplices=false;break;}
    intersection_registry->Intersection_List(simplices,particles);
    // also grab particles on cuts if all cutting simplices
    if(is_all_cutting_simplices){VECTOR<int,2> converted_simplices;ARRAY<int> particles_on_cuts;
        for(int i=1;i<=simplices.m;i++){converted_simplices(i)=cutting_simplices->simplices(simplices(i)).parent;
        assert((cutting_simplices->simplices(converted_simplices(i)).type==CUTTING_SIMPLEX<T,3>::GLOBAL_CUT_FACE));}
        int tet=cutting_simplices->simplices(simplices(1)).element_owner;
        const VECTOR<int,4>& tet_nodes=current_tetrahedralized_volume->mesh.elements(tet);
        intersection_registry->Intersection_List_For_Cuts(converted_simplices,tet_nodes,particles_on_cuts); // only keep particle if simplices contained in our tet
        particles.Append_Elements(particles_on_cuts);}
}
//#####################################################################
// Function Get_Node_Weights_From_Triangle
//#####################################################################
template<class T> VECTOR<T,2> CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Get_Node_Weights_From_Triangle(const int vertex,const VECTOR<int,3>& triangle_vertices) const
{
    VECTOR<T,2> weights;int index=triangle_vertices.Find(vertex);if(index<=2) weights(index)=(T)1;
    return weights;
}
//#####################################################################
// Function Get_Face_Weights_From_Tet
//#####################################################################
template<class T> VECTOR<VECTOR<T,3>,3> CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Get_Face_Weights_From_Tet(const VECTOR<int,3>& face_vertices,const VECTOR<int,4>& tet_vertices) const
{
    VECTOR<VECTOR<T,3>,3> weights;
    for(int i=1;i<=3;i++){int index=tet_vertices.Find(face_vertices(i));weights(i)=VECTOR<T,3>();if(index<=3) weights(i)(index)=(T)1;}
    return weights;
}
//#####################################################################
// Function Get_Weight_From_Segment_On_Face
//#####################################################################
template<class T> VECTOR<T,2> CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Get_Weight_From_Segment_On_Face(T bary,const int node_1,const int node_2,VECTOR<int,3> triangle) const
{
    int permutation=1;for(;permutation<=6;permutation++) if(permute_three(triangle,permutation).Remove_Index(3)==VECTOR<int,2>(node_1,node_2)) break; // find correct permutationd
    assert(permutation<=6);
    return permute_three_inverse(VECTOR<T,3>(bary,1-bary,0),permutation).Remove_Index(3);
}
//#####################################################################
// Function Get_Triangle_Weights_From_Tet
//#####################################################################
template<class T> VECTOR<VECTOR<T,3>,3> CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Get_Triangle_Weights_From_Tet(const VECTOR<TV,3>& triangle_vertex_positions,const VECTOR<int,4>& tet_vertices) const
{
    const GEOMETRY_PARTICLES<TV>& current_embedding_particles=current_tetrahedralized_volume->particles;
    VECTOR<VECTOR<T,3>,3> all_weights;
    for(int i=1;i<=3;i++){
        all_weights(i)=TETRAHEDRON<T>::First_Three_Barycentric_Coordinates(triangle_vertex_positions(i),current_embedding_particles.X(tet_vertices(1)),
            current_embedding_particles.X(tet_vertices(2)),current_embedding_particles.X(tet_vertices(3)),current_embedding_particles.X(tet_vertices(4)));}
    return all_weights;
}
//#####################################################################
// Function Get_Simplex_Weights_For_Edge_Triangle_Intersection
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Get_Simplex_Weights_For_Edge_Triangle_Intersection(const VECTOR<int,3>& simplices,const int triangle_array_index,const VECTOR<int,2>& shared_edge,VECTOR<VECTOR<T,2>,3>& all_weights) const
{
    Intersects_And_Get_Simplex_Weights_For_Edge_Triangle_Intersection_Helper(simplices,triangle_array_index,shared_edge,all_weights,false);
}
//#####################################################################
// Function Intersects_And_Get_Simplex_Weights_For_Edge_Triangle_Intersection
//#####################################################################
template<class T> bool CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Intersects_And_Get_Simplex_Weights_For_Edge_Triangle_Intersection(const VECTOR<int,3>& simplices,const int triangle_array_index,const VECTOR<int,2>& shared_edge,
    VECTOR<VECTOR<T,2>,3>& all_weights) const
{
    return Intersects_And_Get_Simplex_Weights_For_Edge_Triangle_Intersection_Helper(simplices,triangle_array_index,shared_edge,all_weights,true);
}
//#####################################################################
// Function Intersects_And_Get_Simplex_Weights_For_Edge_Triangle_Intersection_Helper
//#####################################################################
template<class T> bool CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Intersects_And_Get_Simplex_Weights_For_Edge_Triangle_Intersection_Helper(const VECTOR<int,3>& simplices,const int triangle_array_index,
    const VECTOR<int,2>& shared_edge,VECTOR<VECTOR<T,2>,3>& all_weights,bool check_if_intersects) const
{
    int temp_triangle_index=(triangle_array_index==1?2:1);
    const CUTTING_SIMPLEX<T,3>& temp_triangle=cutting_simplices->simplices(simplices[temp_triangle_index]);
    VECTOR<T,3> triangle_coordinates;
    VECTOR<T,2> segment_coordinates;
    typedef typename CUTTING_SIMPLEX<T,3>::GET_ADAPTIVE_WEIGHTS_RESULT_TYPE WEIGHT_TYPE;
    VECTOR<VECTOR<WEIGHT_TYPE,3>,3> adaptive_triangle;
    cutting_simplices->simplices(simplices[triangle_array_index]).Get_Adaptive_Weights(adaptive_triangle);
    VECTOR<VECTOR<WEIGHT_TYPE,3>,2> adaptive_segment;
    VECTOR<int,2> edge_node_indices;
    edge_node_indices[1]=temp_triangle.nodes.Find(shared_edge[1]);
    edge_node_indices[2]=temp_triangle.nodes.Find(shared_edge[2]);
    temp_triangle.Get_Adaptive_Weights(adaptive_segment,edge_node_indices);
    if(check_if_intersects){
        bool is_degenerate;
        bool intersects=Intersects<void>(adaptive_triangle,adaptive_segment,&is_degenerate);
        if(is_degenerate) PHYSBAM_FATAL_ERROR("degeneracy");
        if(!intersects) return false;}
    Intersection_Coordinates<void>(adaptive_triangle,adaptive_segment,triangle_coordinates,segment_coordinates,abs_tol_for_level2_barycentric_coordinates);
    all_weights[triangle_array_index][1]=triangle_coordinates[1];
    all_weights[triangle_array_index][2]=triangle_coordinates[2];
    T segment_weight=segment_coordinates[1];
    for(int i=1;i<=3;++i){
        if(i==triangle_array_index) continue;
        all_weights[i]=Get_Weight_From_Segment_On_Face(segment_weight,shared_edge[1],shared_edge[2],cutting_simplices->simplices(simplices[i]).nodes);}
    return true;
}
//#####################################################################
// Function Segment_Reversed_In_Triangle
//#####################################################################
template<class T> bool CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Segment_Reversed_In_Triangle(const VECTOR<int,3>& triangle_nodes,const VECTOR<int,2>& segment_nodes) const
{
    int i=triangle_nodes.Find(segment_nodes[1]);if(triangle_nodes[i%3+1]==segment_nodes[2]) return false;
    assert(triangle_nodes[(i+1)%3+1]==segment_nodes[2]);return true;
}
//#####################################################################
// Function Face_Reversed_In_Tetrahedron
//#####################################################################
template<class T> bool CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Face_Reversed_In_Tetrahedron(const VECTOR<int,4>& tetrahedron_nodes,const VECTOR<int,3>& triangle_nodes,const GEOMETRY_PARTICLES<TV>& particles) const
{
    int other_node=TETRAHEDRON_MESH::Other_Node(tetrahedron_nodes,triangle_nodes);
    return TETRAHEDRON<T>::Signed_Volume(particles.X(triangle_nodes(1)),particles.X(triangle_nodes(2)),particles.X(triangle_nodes(3)),particles.X(other_node))<0;
}
//#####################################################################
// Function Get_Polygon_Edges
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Get_Polygon_Edges(const int polygon_element_index,ARRAY<VECTOR<int,2> >& polygonal_segments) const
{
    const ARRAY<ARRAY<int> >& polygon_particles=polygon_mesh.elements(polygon_element_index);
    for(int i=1;i<=polygon_particles.m;i++){
        const ARRAY<int>& component=polygon_particles(i);
        for(int j=1;j<component.m;j++) polygonal_segments.Append(VECTOR<int,2>(component(j),component(j+1)));
        if(component.m) polygonal_segments.Append(VECTOR<int,2>(component.Last(),component(1)));}
}
//#####################################################################
// Function Get_Particle_Position_In_Local_Tet_Space
//#####################################################################
template<class T> VECTOR<T,3> CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Get_Particle_Position_In_Local_Tet_Space(const int particle,const int child_simplex) const
{
    assert(!cutting_simplices->Simplex_Is_Parent(child_simplex));
    const VECTOR<T,2>& particle_weights_on_simplex=intersection_registry->Get_Simplex_Weights_Of_Intersection(particle,child_simplex);
    const VECTOR<TV,3>& simplex_weights_wrt_tet=cutting_simplices->simplices(child_simplex).weights;
    TV particle_weights_wrt_tet=particle_weights_on_simplex(1)*simplex_weights_wrt_tet(1)+
        particle_weights_on_simplex(2)*simplex_weights_wrt_tet(2)+
        ((T)1-particle_weights_on_simplex.Sum())*simplex_weights_wrt_tet(3);
    return particle_weights_wrt_tet;
}
//#####################################################################
// Function Get_Child_Cutting_Simplex_Local_Normal
//#####################################################################
template<class T> VECTOR<T,3> CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Get_Child_Cutting_Simplex_Local_Normal(const int simplex) const
{
    assert(!cutting_simplices->Simplex_Is_Parent(simplex));
    const VECTOR<TV,3>& simplex_nodes_in_local_tet_space=cutting_simplices->simplices(simplex).weights;
    return TRIANGLE_3D<T>::Normal(simplex_nodes_in_local_tet_space(1),simplex_nodes_in_local_tet_space(2),simplex_nodes_in_local_tet_space(3));
}
//#####################################################################
// Function Local_Planar_Coordinates_For_Intersection
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Local_Planar_Coordinates_For_Intersection(int intersection_index,int cutting_simplex_index,VECTOR<T,2>& planar_coordinates) const
{
    const ARRAY<int>& simplices_on_intersection=intersection_registry->simplices_on_intersection(intersection_index);
    const ARRAY<VECTOR<T,2> >& simplex_weights_on_intersection=intersection_registry->simplex_weights_on_intersection(intersection_index);
    int s=simplices_on_intersection.Find(cutting_simplex_index);
    if(s!=0){planar_coordinates=simplex_weights_on_intersection(s);return;}
    int parent_simplex_index=cutting_simplices->simplices(cutting_simplex_index).parent;
    assert(parent_simplex_index!=0);
    s=simplices_on_intersection.Find(parent_simplex_index);
    if(s!=0){planar_coordinates=simplex_weights_on_intersection(s);return;}
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function World_Coordinates_For_Intersection
//#####################################################################
template<class T> bool CUTTING_TETRAHEDRALIZED_VOLUME<T>::
World_Coordinates_For_Intersection(int intersection_index,TV& world_coordinates) const
{
    // If only it were this easy...
    //return current_tetrahedralized_volume->particles.X(
    //    cutting_particles.Particle_Id_From_Intersection(intersection_index));
    const ARRAY<int>& simplex_indices=intersection_registry->simplices_on_intersection(intersection_index);
    bool has_GLOBAL_EMBEDDING_FACE=false,has_LOCAL_EMBEDDING_FACE=false,has_GLOBAL_CUT_FACE=false,has_LOCAL_CUT_FACE=false;
    for(int s=1;s<=simplex_indices.m;++s){
        switch(cutting_simplices->simplices(simplex_indices(s)).type){
        case CUTTING_SIMPLEX<T,3>::GLOBAL_EMBEDDING_FACE: has_GLOBAL_EMBEDDING_FACE=true;break;
        case CUTTING_SIMPLEX<T,3>::LOCAL_EMBEDDING_FACE: has_LOCAL_EMBEDDING_FACE=true;break;
        case CUTTING_SIMPLEX<T,3>::GLOBAL_CUT_FACE: has_GLOBAL_CUT_FACE=true;break;
        case CUTTING_SIMPLEX<T,3>::LOCAL_CUT_FACE: has_LOCAL_CUT_FACE=true;break;
        default: PHYSBAM_FATAL_ERROR();}}
    // TODO: determine if these don't ever happen
    if(has_LOCAL_EMBEDDING_FACE) PHYSBAM_FATAL_ERROR();
    if(!has_GLOBAL_EMBEDDING_FACE && !has_GLOBAL_CUT_FACE && !has_LOCAL_CUT_FACE) PHYSBAM_FATAL_ERROR();
    if(has_GLOBAL_CUT_FACE && !has_GLOBAL_EMBEDDING_FACE) PHYSBAM_FATAL_ERROR();
    if(has_LOCAL_CUT_FACE && (has_GLOBAL_EMBEDDING_FACE||has_GLOBAL_CUT_FACE)) PHYSBAM_FATAL_ERROR();
    for(int s=1;s<=simplex_indices.m;++s){
        const CUTTING_SIMPLEX<T,3>& cutting_simplex=cutting_simplices->simplices(simplex_indices(s));
        typename CUTTING_SIMPLEX<T,3>::SIMPLEX_TYPE simplex_type=cutting_simplex.type;
        const VECTOR<T,2>& simplex_weights=intersection_registry->simplex_weights_on_intersection(intersection_index)(s);
        const GEOMETRY_PARTICLES<TV>& particles=current_tetrahedralized_volume->particles;
        if(simplex_type==CUTTING_SIMPLEX<T,3>::GLOBAL_EMBEDDING_FACE){
            world_coordinates=particles.X(cutting_simplex.nodes[1])*simplex_weights(1)+particles.X(cutting_simplex.nodes[2])*simplex_weights(2)+
                particles.X(cutting_simplex.nodes[3])*(1-simplex_weights.Sum());
            return has_GLOBAL_CUT_FACE;}
        if(simplex_type==CUTTING_SIMPLEX<T,3>::LOCAL_CUT_FACE){
            TV local_coordinates=cutting_simplex.weights(1)*simplex_weights(1)+cutting_simplex.weights(2)*simplex_weights(2)+
                cutting_simplex.weights(3)*(1-simplex_weights.Sum());
            int element_index=cutting_simplex.element_owner;
            const VECTOR<int,4>& element=current_tetrahedralized_volume->mesh.elements(element_index);
            world_coordinates=particles.X(element(1))*local_coordinates(1)+particles.X(element(2))*local_coordinates(2)+
                particles.X(element(3))*local_coordinates(3)+particles.X(element(4))*(1-local_coordinates.Sum());
            return true;}}
    PHYSBAM_FATAL_ERROR();
    return false;
}
//#####################################################################
// Function World_Coordinates_For_Intersection
//#####################################################################
template<class T> bool CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Is_Cutting_Particle(int intersection_index) const
{
    const ARRAY<int>& simplex_indices=intersection_registry->simplices_on_intersection(intersection_index);
    for(int s=1;s<=simplex_indices.m;++s){
        typename CUTTING_SIMPLEX<T,3>::SIMPLEX_TYPE type=cutting_simplices->simplices(simplex_indices(s)).type;
        if(type==CUTTING_SIMPLEX<T,3>::GLOBAL_CUT_FACE || type==CUTTING_SIMPLEX<T,3>::LOCAL_CUT_FACE) return true; }
    return false;
}
//#####################################################################
// Function Write_Consistent_State
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Write_Consistent_State(const STREAM_TYPE stream_type,const std::string& output_directory,const int suffix) const
{
    std::string snum=FILE_UTILITIES::Number_To_String(suffix);
    std::string prefix=output_directory+"/cutting_state_";
    // for more cutting
    FILE_UTILITIES::Write_To_File(stream_type,prefix+"tet_volume."+snum,*current_tetrahedralized_volume);
    FILE_UTILITIES::Write_To_File(stream_type,prefix+"cutting_simplices."+snum,*cutting_simplices);
    FILE_UTILITIES::Write_To_File(stream_type,prefix+"intersection_registry."+snum,*intersection_registry);
    FILE_UTILITIES::Write_To_File(stream_type,prefix+"simplices_per_tet."+snum,simplices_per_current_tet);
    FILE_UTILITIES::Write_To_File(stream_type,prefix+"cutting_polygons."+snum,current_cutting_polygons);
    FILE_UTILITIES::Write_To_File(stream_type,prefix+"regions_per_tet."+snum,regions_per_tet);
    FILE_UTILITIES::Write_To_File(stream_type,prefix+"boundary_nodes."+snum,sorted_to_original_boundary_nodes);
    FILE_UTILITIES::Write_To_File(stream_type,prefix+"particle_ids."+snum,cutting_particles);
    FILE_UTILITIES::Write_To_File(stream_type,prefix+"polygons_per_simplex."+snum,cutting_polygons_per_cutting_simplex);
    FILE_UTILITIES::Write_To_File(stream_type,prefix+"polygon_mesh."+snum,polygon_mesh);
    // for simulation (in addition to tet volume)
    if(suffix){
        FILE_UTILITIES::Write_To_File(stream_type,prefix+"final_cutting_surface."+snum,*cutting_triangulated_surface);
        FILE_UTILITIES::Write_To_File(stream_type,prefix+"final_boundary_mesh."+snum,final_duplicated_boundary_mesh);
        FILE_UTILITIES::Write_To_File(stream_type,prefix+"final_parents."+snum,final_parents_per_new_particle);
        FILE_UTILITIES::Write_To_File(stream_type,prefix+"final_weights."+snum,final_parent_weights_per_new_particle);
        FILE_UTILITIES::Write_To_File(stream_type,prefix+"final_map_to_previous_particles."+snum,previous_particle_index_per_new_particle_index);}
}
//#####################################################################
// Function Read_Consistent_State
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Read_Consistent_State(const STREAM_TYPE stream_type,const std::string& input_directory,const int suffix)
{
    // TODO: Implement restart!
    if(current_tetrahedralized_volume) delete current_tetrahedralized_volume;
    current_tetrahedralized_volume=TETRAHEDRALIZED_VOLUME<T>::Create();
    std::string snum=FILE_UTILITIES::Number_To_String(suffix);
    std::string prefix=input_directory+"/cutting_state_";
    FILE_UTILITIES::Read_From_File(stream_type,prefix+"tet_volume."+snum,*current_tetrahedralized_volume);
    FILE_UTILITIES::Read_From_File(stream_type,prefix+"cutting_simplices."+snum,*cutting_simplices);
    FILE_UTILITIES::Read_From_File(stream_type,prefix+"intersection_registry."+snum,*intersection_registry);
    FILE_UTILITIES::Read_From_File(stream_type,prefix+"simplices_per_tet."+snum,simplices_per_current_tet);
    FILE_UTILITIES::Read_From_File(stream_type,prefix+"cutting_polygons."+snum,current_cutting_polygons);
    FILE_UTILITIES::Read_From_File(stream_type,prefix+"regions_per_tet."+snum,regions_per_tet);
    FILE_UTILITIES::Read_From_File(stream_type,prefix+"boundary_nodes."+snum,sorted_to_original_boundary_nodes);
    FILE_UTILITIES::Read_From_File(stream_type,prefix+"particle_ids."+snum,cutting_particles);
    FILE_UTILITIES::Read_From_File(stream_type,prefix+"polygons_per_simplex."+snum,cutting_polygons_per_cutting_simplex);
    FILE_UTILITIES::Read_From_File(stream_type,prefix+"polygon_mesh."+snum,polygon_mesh);
}
//#####################################################################
// Function Draw_Polygon
//#####################################################################
template<class T> void CUTTING_TETRAHEDRALIZED_VOLUME<T>::
Draw_Polygon(const int simplex,const bool flipped,const ARRAY<ARRAY<VECTOR<int,2> > >& unconnected_polygonal_regions) const
{
    for(int i=1;i<=unconnected_polygonal_regions.m;i++){
        T hsb=(T)(.3*(i%10)/9);
        std::stringstream ss;
        ss<<"% Run "<<i<<std::endl;
        ss<<" newpath "<<std::endl;
        for(int j=1;j<=unconnected_polygonal_regions(i).m;j++){
            const int p=unconnected_polygonal_regions(i)(j)[1];
            VECTOR<T,2> X1=intersection_registry->Get_Simplex_Weights_Of_Intersection(p,simplex);
            ss<<"  "<<X1<<(j==1?" moveto ":" lineto ");}
        ss<<"\n"<<hsb<<" .5 1 sethsbcolor fill \n newpath "<<std::endl;
        for(int j=1;j<=unconnected_polygonal_regions(i).m;j++){
            const int p=unconnected_polygonal_regions(i)(j)[1];
            VECTOR<T,2> X1=intersection_registry->Get_Simplex_Weights_Of_Intersection(p,simplex);
            ss<<"  "<<X1<<" ("<<p<<") point ";}
        ss<<"\n newpath "<<std::endl;
        for(int j=1;j<=unconnected_polygonal_regions(i).m;j++){
            const VECTOR<int,2>& p=unconnected_polygonal_regions(i)(j);
            VECTOR<T,2> X1=intersection_registry->Get_Simplex_Weights_Of_Intersection(p[1],simplex);
            VECTOR<T,2> X2=intersection_registry->Get_Simplex_Weights_Of_Intersection(p[2],simplex);
            ss<<"  "<<X1<<" "<<X2<<" stem head headl cut "<<(flipped?-1:1)<<" arrow "<<std::endl;}
        ss<<hsb<<" 1 1 sethsbcolor fill"<<"% END RUN"<<std::endl;
        LOG::filecout(ss.str());}
}
//#####################################################################
template class CUTTING_TETRAHEDRALIZED_VOLUME<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class CUTTING_TETRAHEDRALIZED_VOLUME<double>;
#endif

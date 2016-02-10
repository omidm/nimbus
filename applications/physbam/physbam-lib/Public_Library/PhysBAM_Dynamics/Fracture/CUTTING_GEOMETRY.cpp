//#####################################################################
// Copyright 2007, Kevin Der, Geoffrey Irving, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/SIMPLEX_INTERACTIONS.h>
#include <PhysBAM_Dynamics/Fracture/CUTTING_GEOMETRY.h>
using namespace PhysBAM;
//#####################################################################
// Function CUTTING_GEOMETRY
//#####################################################################
template<class TV,int d_input> CUTTING_GEOMETRY<TV,d_input>::
CUTTING_GEOMETRY(const bool verbose)
    :cutting(*new T_CUT_MESH,*new PARTICLES<TV>()),cutting_simplices(new T_CUTTING_SIMPLICES),intersection_registry(new T_INTERSECTION_REGISTRY(*cutting_simplices)),verbose(verbose),
    first_new_cut_element(1)
{
    Set_Intersection_Thickness();
}
//#####################################################################
// Function ~CUTTING_GEOMETRY
//#####################################################################
template<class TV,int d_input> CUTTING_GEOMETRY<TV,d_input>::
~CUTTING_GEOMETRY()
{
    delete &cutting.mesh;
    delete &cutting.particles;
}
//#####################################################################
// Function Initialize_Cutting_Acceleration_Structures
//#####################################################################
namespace{
    template<class T> void Initialize_Cutting_Acceleration_Structures(TRIANGULATED_SURFACE<T>& cutting,TETRAHEDRALIZED_VOLUME<T>& embedding)
    {
        cutting.mesh.Initialize_Element_Edges();
        cutting.mesh.Initialize_Edge_Triangles();
        embedding.mesh.Initialize_Incident_Elements();
        embedding.mesh.Initialize_Triangle_Mesh();
        embedding.mesh.triangle_mesh->Initialize_Incident_Elements();
        embedding.mesh.triangle_mesh->Initialize_Element_Edges();
        embedding.mesh.triangle_mesh->Initialize_Edge_Triangles();
        embedding.mesh.Initialize_Segment_Mesh_From_Triangle_Mesh(); // IMPORTANT make TET_MESH::segment_mesh consistent with TET_MESH::triangle_mesh->segment_mesh
        embedding.mesh.Initialize_Element_Edges();
        embedding.Initialize_Hierarchy();
        embedding.particles.Store_Velocity();
        dynamic_cast<PARTICLES<VECTOR<T,3> >&>(embedding.particles).Store_Mass();
    }
    template<class T> void Initialize_Cutting_Acceleration_Structures(SEGMENTED_CURVE_2D<T>& cutting,TRIANGULATED_AREA<T>& embedding)
    {
        embedding.mesh.Initialize_Incident_Elements();
        embedding.mesh.Initialize_Segment_Mesh();
        embedding.mesh.Initialize_Element_Edges();
        embedding.mesh.Initialize_Edge_Triangles();
        embedding.Initialize_Hierarchy();
        embedding.particles.Store_Velocity();
        dynamic_cast<PARTICLES<VECTOR<T,2> >&>(embedding.particles).Store_Mass();
        cutting.mesh.Initialize_Incident_Elements();
    }
    template<class T> void Initialize_Cutting_Acceleration_Structures(SEGMENTED_CURVE<VECTOR<T,3> >& cutting,TRIANGULATED_SURFACE<T>& embedding)
    {
        embedding.mesh.Initialize_Incident_Elements();
        embedding.mesh.Initialize_Segment_Mesh();
        embedding.mesh.Initialize_Element_Edges();
        embedding.mesh.Initialize_Edge_Triangles();
        embedding.Initialize_Hierarchy();
        embedding.particles.Store_Velocity();
        dynamic_cast<PARTICLES<VECTOR<T,3> >&>(embedding.particles).Store_Mass();
    }
}
//#####################################################################
// Function Barycentric_Coordinates
//#####################################################################
namespace{
    template<class T,class T_ARRAY> VECTOR<T,3> Barycentric_Coordinates(const VECTOR<T,3>& X,const T_ARRAY& tetrahedron_X)
    {
        return TETRAHEDRON<T>::First_Three_Barycentric_Coordinates(X,tetrahedron_X);
    }
    template<class T,class T_ARRAY> VECTOR<T,2> Barycentric_Coordinates(const VECTOR<T,2>& X,const T_ARRAY& triangle_X)
    {
        return TRIANGLE_2D<T>::Barycentric_Coordinates(X,triangle_X).Remove_Index(3);
    }
}
//#####################################################################
// Function Get_Cutting_Simplex_Weights_From_Embedding_Simplex
//#####################################################################
namespace{
    template<class T,int d> VECTOR<VECTOR<T,d>,d>
    Get_Cutting_Simplex_Weights_From_Embedding_Simplex(const VECTOR<VECTOR<T,d>,d>& cutting_simplex_X,const VECTOR<VECTOR<T,d>,d+1>& embedding_simplex_X)
    {
        VECTOR<VECTOR<T,d>,d> all_weights;
        for(int i=1;i<=d;i++) all_weights(i)=Barycentric_Coordinates(cutting_simplex_X(i),embedding_simplex_X);
        return all_weights;
    }
}
//#####################################################################
// Function Fake_Triangle_Intersections
//#####################################################################
namespace{
    template<class T> void Fake_Triangle_Intersections(CUTTING_GEOMETRY<VECTOR<T,2>,2>& cutting,ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,2> >& robust_interactions)
    {
        typedef VECTOR<T,2> TV; typedef CUTTING_GEOMETRY<TV,2> T_CUTTING_GEOMETRY;
        typedef typename T_CUTTING_GEOMETRY::T_CUTTING_SIMPLICES T_CUTTING_SIMPLICES;
        typedef typename T_CUTTING_GEOMETRY::T_CUTTING_SIMPLEX T_CUTTING_SIMPLEX;

        const SEGMENT_MESH& cutting_mesh=cutting.cutting.mesh;
        ARRAY<int> tri_intersection_list;
        for(int i=cutting.first_new_cut_element;i<=cutting_mesh.elements.m;i++) for(int j=1;j<=2;j++){
            int node=cutting_mesh.elements(i)(j);if((*cutting_mesh.incident_elements)(node).m==1){ // vertex has only the new cutting segment incident on it so add fake segment
                // find intersections
                ARRAY<int> intersecting_triangles;
                TV cut_X(cutting.cutting.particles.X(node));
                tri_intersection_list.Remove_All();
                cutting.current_embedding->hierarchy->Intersection_List(cut_X,tri_intersection_list,cutting.intersection_thickness);
                for(int i=1;i<=tri_intersection_list.m;i++){int embedding_tri=tri_intersection_list(i);
                    VECTOR<TV,3> embedding_tri_X(cutting.current_embedding->particles.X.Subset(cutting.current_embedding->mesh.elements(embedding_tri)));
                    bool is_robust,intersects=robust_interactions.Intersection(embedding_tri_X,VECTOR<TV,1>(cut_X),&is_robust);if(!is_robust) PHYSBAM_FATAL_ERROR();
                    if(intersects) intersecting_triangles.Append(embedding_tri);}
                // build simplices 
                VECTOR<int,2> fake_segment_nodes(node,0); // TODO: ordering shouldn't matter here as we are combining with a point not with a segment like in 3D
                int parent=cutting.cutting_simplices->Add_Simplex(-fake_segment_nodes,T_CUTTING_SIMPLEX::GLOBAL_CUT_FACE);
                for(int i=1;i<=intersecting_triangles.m;i++){int tri=intersecting_triangles(i);
                    VECTOR<TV,2> fake_segment_X(cut_X,TV());
                    VECTOR<TV,3> embedding_triangle_X(cutting.current_embedding->particles.X.Subset(cutting.current_embedding->mesh.elements(tri)));
                    VECTOR<VECTOR<T,2>,2> weights=Get_Cutting_Simplex_Weights_From_Embedding_Simplex(fake_segment_X,embedding_triangle_X);
                    int index=cutting.cutting_simplices->Add_Simplex(-fake_segment_nodes,T_CUTTING_SIMPLEX::LOCAL_CUT_FACE,weights,parent,tri);
                    cutting.simplex_stack_per_current_embedding_simplex(tri).Push(index);}}}
        if(cutting.verbose) cutting.cutting_simplices->Print();
    }
    template<class T> void Fake_Triangle_Intersections(CUTTING_GEOMETRY<VECTOR<T,3>,3>& cutting,ROBUST_SIMPLEX_INTERACTIONS<VECTOR<T,3> >& robust_interactions)
    {
        typedef VECTOR<T,3> TV; typedef CUTTING_GEOMETRY<TV,3> T_CUTTING_GEOMETRY;
        typedef typename T_CUTTING_GEOMETRY::T_FACE_MESH T_CUTTING_MESH;typedef typename T_CUTTING_GEOMETRY::T_EMBEDDING_MESH T_EMBEDDING_MESH;
        typedef typename T_CUTTING_GEOMETRY::T_CUTTING_SIMPLICES T_CUTTING_SIMPLICES;
        typedef typename T_CUTTING_GEOMETRY::T_CUTTING_SIMPLEX T_CUTTING_SIMPLEX;
        const ARRAY<ARRAY<int> >& cutting_edge_triangles=*cutting.cutting.mesh.edge_triangles;
        const ARRAY<VECTOR<int,2> >& cutting_segments=cutting.cutting.mesh.segment_mesh->elements;

        ARRAY<int> tet_intersection_list;
        for(int cseg=1;cseg<=cutting_segments.m;cseg++) if(cutting_edge_triangles(cseg).m==1 && cutting_edge_triangles(cseg)(1)>=cutting.first_new_cut_element){
            // find intersections
            ARRAY<int> intersecting_tets;
            const VECTOR<int,2>& cut_nodes=cutting_segments(cseg);VECTOR<TV,2> cut_X(cutting.cutting.particles.X.Subset(cut_nodes));
            tet_intersection_list.Remove_All();
            cutting.current_embedding->hierarchy->Intersection_List(RANGE<VECTOR<T,3> >::Bounding_Box(cut_X),tet_intersection_list,cutting.intersection_thickness);
            for(int i=1;i<=tet_intersection_list.m;i++){int embedding_tet=tet_intersection_list(i);
                VECTOR<TV,4> embedding_tet_X(cutting.current_embedding->particles.X.Subset(cutting.current_embedding->mesh.elements(embedding_tet)));
                bool is_robust,intersects=robust_interactions.Intersection(embedding_tet_X,cut_X,&is_robust);if(!is_robust) PHYSBAM_FATAL_ERROR();
                if(intersects) intersecting_tets.Append(embedding_tet);}
            // build simplices
            const VECTOR<int,3>& real_bordering_triangle_nodes=cutting.cutting.mesh.elements(cutting_edge_triangles(cseg)(1));
            VECTOR<int,2> ordered_cut_nodes=T_CUTTING_MESH::Face_Reversed_In_Simplex(cut_nodes,real_bordering_triangle_nodes)?cut_nodes.Reversed():cut_nodes;
            VECTOR<int,3> fake_triangle_nodes=ordered_cut_nodes.Append(0);
            int parent=cutting.cutting_simplices->Add_Simplex(-fake_triangle_nodes,T_CUTTING_SIMPLEX::GLOBAL_CUT_FACE);
            for(int i=1;i<=intersecting_tets.m;i++){int tet=intersecting_tets(i);
                VECTOR<TV,3> fake_triangle_X(VECTOR<TV,2>(cutting.cutting.particles.X.Subset(ordered_cut_nodes)).Append(TV()));
                VECTOR<TV,4> embedding_tetrahedron_X(cutting.current_embedding->particles.X.Subset(cutting.current_embedding->mesh.elements(tet)));
                VECTOR<VECTOR<T,3>,3> weights=Get_Cutting_Simplex_Weights_From_Embedding_Simplex(fake_triangle_X,embedding_tetrahedron_X);
                int index=cutting.cutting_simplices->Add_Simplex(-fake_triangle_nodes,T_CUTTING_SIMPLEX::LOCAL_CUT_FACE,weights,parent,tet);
                cutting.simplex_stack_per_current_embedding_simplex(tet).Push(index);}}
        if(cutting.verbose) cutting.cutting_simplices->Print();}
}
//#####################################################################
// Function Initialize_Original_Embedding
//#####################################################################
template<class TV,int d_input> void CUTTING_GEOMETRY<TV,d_input>::
Initialize_Original_Embedding(const T_EMBEDDING_OBJECT& original_embedding)
{  
    // copy original tet volume to current and initialize accel structures
    current_embedding=T_EMBEDDING_OBJECT::Create();
    current_embedding->particles.array_collection->Initialize(*original_embedding.particles.array_collection);
    current_embedding->mesh.Initialize_Mesh(original_embedding.mesh);
    current_embedding->mesh.Initialize_Incident_Elements();
    Initialize_Face_Mesh(current_embedding->mesh);

    const ARRAY<VECTOR<int,d_embed+1> >& original_embedding_simplices=current_embedding->mesh.elements;
    const ARRAY<VECTOR<int,d_embed+1> >& embedding_simplices=current_embedding->mesh.elements;
    const T_FACE_MESH& embedding_faces=Face_Mesh(current_embedding->mesh);

    simplices_per_current_embedding_simplex.Resize(original_embedding_simplices.m);regions_per_embedding_simplex.Resize(original_embedding_simplices.m);
    // cutting simplices, fill polygon mesh and generate lists of cutting polygons
    HASHTABLE<int,int> embedding_face_to_simplex;
    ARRAY<int> embedding_simplices_on_face;
    for(int face=1;face<=embedding_faces.elements.m;face++){const VECTOR<int,d_cut+1>& face_nodes=embedding_faces.elements(face); 
        int parent=cutting_simplices->Add_Simplex(face_nodes,T_CUTTING_SIMPLEX::GLOBAL_EMBEDDING_FACE); // TODO: call this a boundary type???
        sorted_to_original_boundary_nodes.Insert(face_nodes.Sorted(),face_nodes);
        cutting_polygons_per_cutting_simplex.Append(ARRAY<int>());
        embedding_face_to_simplex.Insert(face,parent);
        embedding_simplices_on_face.Remove_All();
        current_embedding->mesh.Simplices_On_Subsimplex(face_nodes,embedding_simplices_on_face);
        CUTTING_POLYGON::POLYGON_TYPE polygon_type=(embedding_simplices_on_face.m==1?CUTTING_POLYGON::FACE_BOUNDARY:CUTTING_POLYGON::FACE_INTERIOR);
        for(int i=1;i<=embedding_simplices_on_face.m;i++){int embedding_simplex=embedding_simplices_on_face(i);
            T_CUTTING_SIMPLEX_WEIGHTS weights=Get_Face_Weights_From_Embedding_Simplex(face_nodes,embedding_simplices(embedding_simplex));
            int simplex_index=cutting_simplices->Add_Simplex(face_nodes,T_CUTTING_SIMPLEX::LOCAL_EMBEDDING_FACE,weights,parent,embedding_simplex);
            simplices_per_current_embedding_simplex(embedding_simplex).Append(simplex_index);
            // make polygon
            polygon_mesh.elements.Append(ARRAY<ARRAY<int> >());ARRAY<ARRAY<int> >& new_polygon=polygon_mesh.elements.Last();
            bool face_reversed=T_EMBEDDING_MESH::Face_Reversed_In_Simplex(face_nodes,embedding_simplices(embedding_simplex));
            VECTOR<int,d_cut+1> face_nodes_ordered=face_nodes;
            // TODO: reverse face
            if(face_reversed) exchange(face_nodes_ordered[d_cut],face_nodes_ordered[d_cut+1]);
            new_polygon.Append(ARRAY<int>(face_nodes_ordered));
            int cutting_polygon_index=current_cutting_polygons.Append(CUTTING_POLYGON(polygon_mesh.elements.m,simplex_index,face_reversed,polygon_type));
            cutting_polygons_per_cutting_simplex.Append(ARRAY<int>(VECTOR<int,1>(cutting_polygon_index)));
            if(!regions_per_embedding_simplex(embedding_simplex).m) regions_per_embedding_simplex(embedding_simplex).Append(ARRAY<int>());
            regions_per_embedding_simplex(embedding_simplex)(1).Append(current_cutting_polygons.m);}}
    // intersection registry (vertices)
    const ARRAY<ARRAY<int> >& faces_on_vertices=*Face_Mesh(current_embedding->mesh).incident_elements;
    ARRAY<int> embedding_nodes;
    current_embedding->mesh.elements.Flattened().Get_Unique(embedding_nodes);
    for(int i=1;i<=embedding_nodes.m;i++){int node=embedding_nodes(i);ARRAY<int> simplices;
        ARRAY<T_CUT_WEIGHTS> all_weights(faces_on_vertices(node).m);
        for(int i=1;i<=faces_on_vertices(node).m;i++){int face=faces_on_vertices(node)(i);
            int parent;bool found=embedding_face_to_simplex.Get(face,parent);if(!found) PHYSBAM_FATAL_ERROR("Could not find simplex");
            simplices.Append(parent);
            all_weights(i)=Get_Node_Weights_From_Face(node,embedding_faces.elements(face));}
        int intersection=intersection_registry->Number_Of_Intersections()+1;
        intersection_registry->Register_Intersection(simplices,all_weights,intersection);
        cutting_particles.Add_Tet_Node_And_Intersection_Id(node,intersection);}
    if(verbose){cutting_simplices->Print();cutting_particles.Print();intersection_registry->Print();}
}
//#####################################################################
// Function Cutting_Intersections_Helper
//#####################################################################
namespace{
    template<class T_CUTTING_GEOMETRY> void Cutting_Intersections_Helper(T_CUTTING_GEOMETRY& cutting)
    {
        typedef typename T_CUTTING_GEOMETRY::VECTOR_T TV;typedef typename TV::SCALAR T;const int d_cut=T_CUTTING_GEOMETRY::d_cut;const int d_embed=T_CUTTING_GEOMETRY::d_embed;
        typedef typename T_CUTTING_GEOMETRY::T_CUTTING_SIMPLICES T_CUTTING_SIMPLICES;
        typedef typename T_CUTTING_GEOMETRY::T_CUTTING_SIMPLEX T_CUTTING_SIMPLEX;
        typedef typename T_CUTTING_GEOMETRY::T_CUTTING_SIMPLEX_WEIGHTS T_CUTTING_SIMPLEX_WEIGHTS;
    
        ROBUST_SIMPLEX_INTERACTIONS<TV> robust_interactions;
        HASHTABLE<int,int> cut_to_simplex;
        cutting.cutting_simplices->Set_Index_For_Last_Old_Cutting_Simplex();
        cutting.simplex_stack_per_current_embedding_simplex.Clean_Memory();cutting.simplex_stack_per_current_embedding_simplex.Resize(cutting.current_embedding->mesh.elements.Size());
        // iterate over every new triangle
        ARRAY<int> intersection_list;
        for(int cutting_simplex_index=cutting.first_new_cut_element;cutting_simplex_index<=cutting.cutting.mesh.elements.m;cutting_simplex_index++){
            ARRAY<int> intersecting_embedding_simplices;
            // find intersections
            const VECTOR<int,d_cut+1>& cut_nodes=cutting.cutting.mesh.elements(cutting_simplex_index);
            VECTOR<TV,d_cut+1> cut_X(cutting.cutting.particles.X.Subset(cut_nodes));
            intersection_list.Remove_All();cutting.current_embedding->hierarchy->Intersection_List(
                RANGE<TV>::Bounding_Box(cut_X),intersection_list,cutting.intersection_thickness);
            for(int i=1;i<=intersection_list.m;i++){int embedding_simplex_index=intersection_list(i);
                VECTOR<TV,d_embed+1> embedding_simplex_X(cutting.current_embedding->particles.X.Subset(cutting.current_embedding->mesh.elements(embedding_simplex_index)));
                bool is_robust,intersects=robust_interactions.Intersection(embedding_simplex_X,cut_X,&is_robust);
                if(!is_robust) PHYSBAM_FATAL_ERROR();
                if(intersects) intersecting_embedding_simplices.Append(embedding_simplex_index);}
            // build simplices
            int parent=cutting.cutting_simplices->Add_Simplex(-cut_nodes,T_CUTTING_SIMPLEX::GLOBAL_CUT_FACE);
            cut_to_simplex.Insert(cutting_simplex_index,parent);
            for(int i=1;i<=intersecting_embedding_simplices.m;i++){int embedding_simplex=intersecting_embedding_simplices(i);
                VECTOR<TV,d_embed+1> embedding_simplex_X(cutting.current_embedding->particles.X.Subset(cutting.current_embedding->mesh.elements(embedding_simplex)));
                T_CUTTING_SIMPLEX_WEIGHTS weights=Get_Cutting_Simplex_Weights_From_Embedding_Simplex(cut_X,embedding_simplex_X);
                int index=cutting.cutting_simplices->Add_Simplex(-cut_nodes,T_CUTTING_SIMPLEX::LOCAL_CUT_FACE,weights,parent,embedding_simplex);
                cutting.simplex_stack_per_current_embedding_simplex(embedding_simplex).Push(index);}}
        // do the fake triangles
        Fake_Triangle_Intersections<T>(cutting,robust_interactions);
    }
    template<class T> void Cutting_Intersections_Helper(CUTTING_GEOMETRY<VECTOR<T,3>,2>& cutting)
    {
        PHYSBAM_NOT_IMPLEMENTED(); // implement for S3D to process from new cutting triangles to make cutting segments
    }
}
//#####################################################################
// Function Fill_Intersection_Registry_With_Cuts
//#####################################################################
template<class TV,int d_input> void CUTTING_GEOMETRY<TV,d_input>::
Fill_Intersection_Registry_With_Cuts()
{
    intersection_registry->Resize_Simplices(cutting_simplices->simplices.m);
    intersection_registry->Set_Index_For_Last_Old_Intersection();
    const ARRAY<VECTOR<int,d_embed+1> >& current_embedding_simplices=current_embedding->mesh.elements;
    for(int tet=1;tet<=current_embedding_simplices.m;tet++){
        STACK<int>& new_simplices=simplex_stack_per_current_embedding_simplex(tet);
        while(!new_simplices.Empty()){
            int new_simplex=new_simplices.Pop();
            Intersect_Simplex_With_Old_Simplices_In_Embedding(tet,new_simplex);
            simplices_per_current_embedding_simplex(tet).Append(new_simplex);}}
    if(verbose) cutting_particles.Print();
    // fix polygon mesh structures
    polygon_mesh.Set_Number_Nodes(intersection_registry->Number_Of_Intersections());
    polygon_mesh.Initialize_Segment_Mesh();
    polygon_mesh.segment_mesh->Initialize_Incident_Elements();
    polygon_mesh.Initialize_Element_Oriented_Edges();
    polygon_mesh.Initialize_Edge_Elements();
}
//#####################################################################
// Function Cut_Material
//#####################################################################
template<class TV,int d_input> void CUTTING_GEOMETRY<TV,d_input>::
Cut_Material(T_EMBEDDING_OBJECT& next_embedding_input)
{
    next_embedding=&next_embedding_input;
    if(next_embedding->mesh.elements.m || next_embedding->particles.array_collection->Size()) PHYSBAM_FATAL_ERROR("Non empty inputted next embedding");
    next_embedding->particles.Store_Velocity();
    dynamic_cast<PARTICLES<TV>&>(next_embedding->particles).Store_Mass();
    
    LOG::Time("Initializing acceleration structures");Initialize_Cutting_Acceleration_Structures<T>(cutting,*current_embedding);
    LOG::Time("Finding cutting intersections");Cutting_Intersections_Helper(*this);
    LOG::Time("Fill_Intersection_Registry_With_Cuts");Fill_Intersection_Registry_With_Cuts();
    LOG::Time("Split Existing Polygon Edges");Split_Existing_Polygon_Edges();
    LOG::Time("Split Existing Polygons");Split_Existing_Polygons();
    // TODO: Continue here

    first_new_cut_element=cutting.mesh.elements.m+1;
}
//#####################################################################
// Function Get_Node_Weights_From_Face
//#####################################################################
template<class TV,int d_input> typename CUTTING_GEOMETRY<TV,d_input>::T_CUT_WEIGHTS CUTTING_GEOMETRY<TV,d_input>::
Get_Node_Weights_From_Face(const int node,const VECTOR<int,d_cut+1>& face_nodes) const
{
    T_CUT_WEIGHTS weights;int index=face_nodes.Find(node);if(index<=T_CUT_WEIGHTS::dimension) weights(index)=(T)1;
    return weights;
}
//#####################################################################
// Function Get_Face_Weights_From_Embedding_Simplex
//#####################################################################
template<class TV,int d_input> typename CUTTING_GEOMETRY<TV,d_input>::T_CUTTING_SIMPLEX_WEIGHTS CUTTING_GEOMETRY<TV,d_input>::
Get_Face_Weights_From_Embedding_Simplex(const VECTOR<int,d_cut+1>& face_nodes,const VECTOR<int,d_embed+1>& embedding_simplex_nodes) const
{
    T_CUTTING_SIMPLEX_WEIGHTS weights;
    for(int i=1;i<=T_CUTTING_SIMPLEX_WEIGHTS::dimension;i++){
        int index=embedding_simplex_nodes.Find(face_nodes(i));weights(i)=T_EMBEDDING_WEIGHTS();
        if(index<=T_EMBEDDING_WEIGHTS::dimension) weights(i)(index)=(T)1;}
    return weights;
}
//#####################################################################
// Function Get_Polygon_Edges
//#####################################################################
template<class TV,int d_input> void CUTTING_GEOMETRY<TV,d_input>::
Get_Polygon_Edges(const int polygon_element_index,ARRAY<VECTOR<int,2> >& polygonal_segments) const
{
    const ARRAY<ARRAY<int> >& polygon_particles=polygon_mesh.elements(polygon_element_index);
    for(int i=1;i<=polygon_particles.m;i++){
        const ARRAY<int>& component=polygon_particles(i);
        for(int j=1;j<component.m;j++) polygonal_segments.Append(VECTOR<int,2>(component(j),component(j+1)));
        if(component.m) polygonal_segments.Append(VECTOR<int,2>(component.Last(),component(1)));}
}
//#####################################################################
// Function Initialize_Face_Mesh
//#####################################################################
namespace{
    void Initialize_Face_Mesh(TETRAHEDRON_MESH& mesh)
    {
        mesh.Initialize_Triangle_Mesh();mesh.triangle_mesh->Initialize_Incident_Elements();
    }
    void Initialize_Face_Mesh(TRIANGLE_MESH& mesh)
    {
        mesh.Initialize_Segment_Mesh();mesh.segment_mesh->Initialize_Incident_Elements();
    }
}
//#####################################################################
// Function Face_Mesh
//#####################################################################
namespace{
    TRIANGLE_MESH& Face_Mesh(TETRAHEDRON_MESH& mesh)
    {
        return *mesh.triangle_mesh;
    }
    SEGMENT_MESH& Face_Mesh(TRIANGLE_MESH& mesh)
    {
        return *mesh.segment_mesh;
    }
}
//#####################################################################
template class CUTTING_GEOMETRY<VECTOR<float,2>,2>;
template class CUTTING_GEOMETRY<VECTOR<float,3>,2>;
template class CUTTING_GEOMETRY<VECTOR<float,3>,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class CUTTING_GEOMETRY<VECTOR<double,2>,2>;
template class CUTTING_GEOMETRY<VECTOR<double,3>,2>;
template class CUTTING_GEOMETRY<VECTOR<double,3>,3>;
#endif

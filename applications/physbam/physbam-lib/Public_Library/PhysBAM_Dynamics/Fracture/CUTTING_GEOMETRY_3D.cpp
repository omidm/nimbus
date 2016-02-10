//#####################################################################
// Copyright 2007, Kevin Der, Geoffrey Irving, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/DIRECTED_GRAPH.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
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
#include <PhysBAM_Dynamics/Fracture/CUTTING_GEOMETRY_3D.h>
using namespace PhysBAM;
//#####################################################################
// Function Split_Existing_Polygons
//#####################################################################
template<class TV,int d_input> void CUTTING_GEOMETRY_3D<TV,d_input>::
Split_Existing_Polygon_Edges()
{
    //typedef typename T_CUTTING_GEOMETRY::T T;typedef typename T_CUTTING_GEOMETRY::VECTOR_T TV;typedef typename T_CUTTING_GEOMETRY::T_CUTTING_SIMPLICES T_CUTTING_SIMPLICES;
    for(int new_intersection=intersection_registry->index_for_last_old_intersection+1;new_intersection<=intersection_registry->Number_Of_Intersections();new_intersection++){
        const ARRAY<int>& simplices=intersection_registry->simplices_on_intersection(new_intersection);if(!simplices.m) PHYSBAM_FATAL_ERROR("TODO: This should be impossible");
        HASHTABLE<VECTOR<int,2> > edges_seen;
        // loop over every pair of simplices
        for(int j=1;j<=simplices.m;j++) for(int k=j+1;k<=simplices.m;k++){
            int simplex_1=simplices(j),simplex_2=simplices(k);
            // TODO: why quit if both are global cut faces
            if(cutting_simplices->simplices(simplex_1).type==T_CUTTING_SIMPLEX::GLOBAL_CUT_FACE
                && cutting_simplices->simplices(simplex_2).type==T_CUTTING_SIMPLEX::GLOBAL_CUT_FACE) continue;
            ARRAY<int> all_particles;Get_Particles_On_Simplices(VECTOR<int,2>(simplex_1,simplex_2),all_particles);
            if(all_particles.m<=2) continue; // we have an edge which has not been cut // TODO: This is 2 in 3D, should this be 1 in 2D
            // get 2D weights on any simplex
            const VECTOR<T,2>& new_particle_weights_on_simplex=intersection_registry->Get_Simplex_Weights_Of_Intersection(new_intersection,simplex_1);
            int closest_particle_left=0,closest_particle_right=0;T closest_left_val=-FLT_MAX,closest_right_val=FLT_MAX;
            // find best conditioned axis to sort particle locations
            VECTOR<T,2> delta_weights=intersection_registry->Get_Simplex_Weights_Of_Intersection(all_particles(1),simplex_1)-
                intersection_registry->Get_Simplex_Weights_Of_Intersection(all_particles(2),simplex_1);
            int dominant_axis=delta_weights.Dominant_Axis();
            // find left and right particle that we wish to insert the new ith particle between
            for(int p=1;p<=all_particles.m;p++) if(all_particles(p)<new_intersection){ // only consider particles that existed before or have already split polygon edges
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
            polygon_mesh.Split_Polygon_Edge(closest_particle_left,closest_particle_right,new_intersection);
            if(verbose) {std::stringstream ss;ss<<"Split edge "<<closest_particle_left<<", "<<closest_particle_right<<" with particle "<<new_intersection<<std::endl;LOG::filecout(ss.str());}}}
}
//#####################################################################
// Function Split_Existing_Polygons_Edges
//#####################################################################
template<class TV,int d_input> void CUTTING_GEOMETRY_3D<TV,d_input>::
Split_Existing_Polygons()
{
    polygons_per_element=CONSTANT_ARRAY<ARRAY<int> >(current_embedding->mesh.elements.m,ARRAY<int>());
    for(int i=1;i<=cutting_simplices->index_for_last_old_cutting_simplex;i++) if(!cutting_simplices->Simplex_Is_Parent(i)){
        int tet_owner=cutting_simplices->simplices(i).element_owner;
        // obtain all new segments on simplex
        ARRAY<VECTOR<int,2> > new_unoriented_segments_on_simplex;
        // TODO: this next for makes quadratic time
        for(int j=cutting_simplices->index_for_last_old_cutting_simplex+1;j<=cutting_simplices->simplices.m;j++) if(!cutting_simplices->Simplex_Is_Fake(j)){
            Get_Unoriented_Segments_On_Two_Simplices(i,j,new_unoriented_segments_on_simplex);}
        // operate on each polygon on the simplex
        const ARRAY<int>& cutting_polygons_on_simplex=cutting_polygons_per_cutting_simplex(i);
        for(int j=1;j<=cutting_polygons_on_simplex.m;j++){
            int cutting_polygon_index=cutting_polygons_on_simplex(j);
            if(!new_unoriented_segments_on_simplex.m){
                polygons_per_element(tet_owner).Append(cutting_polygon_index);}
            else{
                const CUTTING_POLYGON& cutting_polygon=current_cutting_polygons(cutting_polygon_index);
                if(cutting_polygon.simplex_owner!=i) PHYSBAM_FATAL_ERROR("TODO: remove this check");
                ARRAY<VECTOR<int,2> > all_polygonal_segments_on_simplex;
                Get_Polygon_Edges(cutting_polygon.polygon_index,all_polygonal_segments_on_simplex);
                ARRAY<ARRAY<ARRAY<int > > > final_polygon_element_particles;
                Divide_Polygon_Particles_With_New_Segments(all_polygonal_segments_on_simplex,new_unoriented_segments_on_simplex,
                    polygon_mesh.elements(cutting_polygon.polygon_index),i,cutting_polygon.flipped,final_polygon_element_particles);
                // replace old one with first new
                polygon_mesh.elements(cutting_polygon.polygon_index)=final_polygon_element_particles(1);
                polygons_per_element(tet_owner).Append(cutting_polygon_index);
                // make all the rest
                for(int k=2;k<=final_polygon_element_particles.m;k++){
                    int polygon_element_index=polygon_mesh.elements.Append(final_polygon_element_particles(k));
                    int new_cutting_polygon_index=current_cutting_polygons.Append(CUTTING_POLYGON(polygon_element_index,i,cutting_polygon.flipped,cutting_polygon.polygon_type));
                    polygons_per_element(tet_owner).Append(new_cutting_polygon_index);}}}}
    if(verbose){std::stringstream ss;ss<<std::endl;for(int i=1;i<=current_cutting_polygons.m;i++)
        ss<<"Cutting polygon "<<i<<" has polygon element "
                 <<current_cutting_polygons(i).polygon_index<<": "<<polygon_mesh.elements(current_cutting_polygons(i).polygon_index)
                 <<std::endl;LOG::filecout(ss.str());}
}
//#####################################################################
// Function Get_Simplex_Weights_For_Edge_Triangle_Intersection
//#####################################################################
template<class TV,int d_input> void CUTTING_GEOMETRY_3D<TV,d_input>::
Get_Simplex_Weights_For_Edge_Triangle_Intersection(const VECTOR<int,3>& simplices,const int triangle_array_index,const VECTOR<int,2>& shared_edge,
    VECTOR<VECTOR<typename TV::SCALAR,2>,3>& all_weights)
{    
    ROBUST_SIMPLEX_INTERACTIONS<TV> robust_interactions;T segment_weight;
    robust_interactions.Triangle_Segment_Intersection_Weights(cutting_simplices->simplices(simplices(triangle_array_index)).weights,
        cutting_simplices->Weight_For_Nodes_In_Simplex(simplices(1==triangle_array_index?2:1),shared_edge),all_weights(triangle_array_index),segment_weight);
    
    for(int i=1;i<=3;i++) if(i!= triangle_array_index){
        const VECTOR<int,3>& triangle_nodes=cutting_simplices->simplices(simplices(i)).nodes;
        int permutation=1;for(;permutation<=6;permutation++) if(permute_three(triangle_nodes,permutation).Remove_Index(3)==shared_edge) break;
        assert(permutation<=6);
        all_weights(i)=permute_three(VECTOR<T,3>(segment_weight,1-segment_weight,0),permutation).Remove_Index(3);}
}
//#####################################################################
// Function Intersect_Simplex_With_Old_Simplices_in_Embedding
//#####################################################################
template<class TV,int d_input> void CUTTING_GEOMETRY_3D<TV,d_input>::
Intersect_Simplex_With_Old_Simplices_In_Embedding(const int tet,const int new_simplex)
{
    ARRAY<int>& old_simplices=simplices_per_current_embedding_simplex(tet);
    for(int i=1;i<=old_simplices.m;i++) for(int j=i+1;j<=old_simplices.m;j++){
        if(verbose) {std::stringstream ss;ss<<" Checking old simplices "<<old_simplices(i)<<" and "<<old_simplices(j)<<" with new simplex "<<new_simplex<<std::endl;LOG::filecout(ss.str());}
        VECTOR<int,3> simplices(old_simplices(i),old_simplices(j),new_simplex);
        bool converted;VECTOR<int,3> converted_simplices=cutting_simplices->Convert_Simplex_Indices_To_Global_Indices(simplices,converted);

        // skip if all simplices are fake or intersection already found
        if(cutting_simplices->Simplex_Is_Fake(simplices(1)) && cutting_simplices->Simplex_Is_Fake(simplices(2)) && cutting_simplices->Simplex_Is_Fake(simplices(3)))
            continue;
        {ARRAY<int> intersection_list;if(intersection_registry->Intersection_List(converted_simplices,intersection_list)) continue;}

        // Case 1 - Sharing of one or more points
        {ARRAY<int> points_all_shared;cutting_simplices->Shared_Nodes_On_Simplices(converted_simplices,points_all_shared);
        if(points_all_shared.m>=2) continue; // they all share at least an edge
        else if(points_all_shared.m==1){int shared_node=points_all_shared(1); // they all share 1 point
            VECTOR<VECTOR<T,2>,3> all_weights;
            assert(simplices==converted_simplices); // intersection point should be all internal to a tet, i.e. so all simplices are cutting simplices
            for(int i=1;i<=3;i++) all_weights(i)=Get_Node_Weights_From_Face(shared_node,cutting_simplices->simplices(simplices(i)).nodes);
            // ensure intersection point is inside the tet
            const VECTOR<TV,3>& simplex_weights_wrt_tet=cutting_simplices->simplices(simplices(1)).weights;
            TV particle_weights_wrt_tet=all_weights(1)(1)*simplex_weights_wrt_tet(1)+all_weights(1)(2)*simplex_weights_wrt_tet(2)+
                ((T)1-all_weights(1).Sum())*simplex_weights_wrt_tet(3);
            if(particle_weights_wrt_tet.Min()<0 || particle_weights_wrt_tet.Sum()>(T)1) goto NEXT_OLD_SIMPLEX;
            // it's not outside, so add it if it doesn't already exist (if it already exists, it must be a shared node on child simplices in this tet)
            for(int i=1;i<=old_simplices.m;i++) if(cutting_simplices->simplices(old_simplices(i)).nodes.Contains(shared_node))
                for(int j=i+1;j<=old_simplices.m;j++) if(cutting_simplices->simplices(old_simplices(j)).nodes.Contains(shared_node))
                    for(int k=j+1;k<=old_simplices.m;k++) if(cutting_simplices->simplices(old_simplices(k)).nodes.Contains(shared_node)){
                        ARRAY<int> nodes_shared_on_triplet;
                        cutting_simplices->Shared_Nodes_On_Simplices(VECTOR<int,3>(old_simplices(i),old_simplices(j),old_simplices(k)),nodes_shared_on_triplet);
                        if(nodes_shared_on_triplet.m>1) continue; // The intersection does not define a single point
                        if(int intersection_index=intersection_registry->Intersection(VECTOR<int,3>(old_simplices(i),old_simplices(j),old_simplices(k))))
                            Register_Cut_Intersection(converted_simplices,all_weights,intersection_index);
                        goto NEXT_OLD_SIMPLEX;}
            Register_Cut_Intersection(converted_simplices,all_weights,0);goto NEXT_OLD_SIMPLEX;}}
        // Case 2: new simplex shares an edge with (exactly) one of the old simplices
        for(int i=1;i<=2;i++){
            ARRAY<int> shared_nodes;cutting_simplices->Shared_Nodes_On_Simplices(VECTOR<int,2>(converted_simplices(i),converted_simplices(3)),shared_nodes);
            if(shared_nodes.m==2 && !cutting_simplices->Simplex_Is_Fake(converted_simplices(3-i))){ // the simplex intersecting the edge shouldn't be fake (2 edges never intersect)
                VECTOR<int,2> shared_edge(shared_nodes(1),shared_nodes(2));
                // Find if there is another set of simplices that shares the same edge. If there is an intersection, it already has it recorded, if not then no intersection exists.
                for(int j=1;j<=old_simplices.m;j++) if(old_simplices(j)!=converted_simplices(3) && old_simplices(j)!=simplices(3)) // TODO: why could new simplex be in old simplices
                    if(cutting_simplices->simplices(old_simplices(j)).nodes.Contains_All(shared_edge)){
                        for(int k=j+1;k<=old_simplices.m;k++) if(old_simplices(k)!=converted_simplices(3) && old_simplices(k)!=simplices(3)){
                            if(cutting_simplices->simplices(old_simplices(j)).nodes==cutting_simplices->simplices(old_simplices(k)).nodes) continue;
                            if(cutting_simplices->simplices(old_simplices(k)).nodes.Contains_All(shared_edge)){
                                int old_simplex_index_1=converted?cutting_simplices->simplices(old_simplices(j)).parent:old_simplices(j);
                                int old_simplex_index_2=converted?cutting_simplices->simplices(old_simplices(k)).parent:old_simplices(k);
                                if(int particle=intersection_registry->Intersection(VECTOR<int,3>(converted_simplices(3-i),old_simplex_index_1,old_simplex_index_2))){
                                    VECTOR<VECTOR<T,2>,3> all_weights;
                                    Get_Simplex_Weights_For_Edge_Triangle_Intersection(VECTOR<int,3>(simplices(3),simplices(i),simplices(3-i)),3,shared_edge,all_weights);
                                    Register_Cut_Intersection(VECTOR<int,1>(converted_simplices(3)),VECTOR<VECTOR<T,2>,1>(all_weights(1)),particle);}
                                goto NEXT_OLD_SIMPLEX;}}}
                // None already found, check if therea is an intersection
                bool intersects=SIMPLEX_INTERACTIONS<T>::Intersection(cutting_simplices->simplices(simplices(3-i)).weights,
                    cutting_simplices->Weight_For_Nodes_In_Simplex(simplices(i),shared_edge));
                if(intersects){
                    VECTOR<VECTOR<T,2>,3> all_weights;
                    Get_Simplex_Weights_For_Edge_Triangle_Intersection(simplices,3-i,shared_edge,all_weights);
                    Register_Cut_Intersection(converted_simplices,all_weights,0);}
                goto NEXT_OLD_SIMPLEX;}}
        // Case 3: old simplices share an edge
        {ARRAY<int> shared_nodes;cutting_simplices->Shared_Nodes_On_Simplices(VECTOR<int,2>(converted_simplices(1),converted_simplices(2)),shared_nodes);
        if(shared_nodes.m==2 && !cutting_simplices->Simplex_Is_Fake(converted_simplices(3))){ // can't intersected shared edge on old simplices with fake new simplex
            VECTOR<int,2> shared_edge(shared_nodes(1),shared_nodes(2));
            // check through other pairs of simplices that share the same edge to see if it already exists
            if(intersection_registry->intersections_on_simplex.m>=converted_simplices(3)){ // TODO: why can this use intersections_on_simplex on the previous two don't?
                const ARRAY<int>& intersections_on_new_simplex=intersection_registry->intersections_on_simplex(converted_simplices(3));
                for(int i=1;i<=intersections_on_new_simplex.m;i++){int old_intersection=intersections_on_new_simplex(i);
                    const ARRAY<int>& simplices_on_old_intersection=intersection_registry->simplices_on_intersection(old_intersection);int count=0;
                    for(int j=1;j<=simplices_on_old_intersection.m;j++){
                        const VECTOR<int,3>& simplex_nodes=cutting_simplices->simplices(simplices_on_old_intersection(j)).nodes;
                        if(simplex_nodes.Contains_All(shared_edge)) count++;
                        if(count>=2){
                            VECTOR<VECTOR<T,2>,3> all_weights;
                            Get_Simplex_Weights_For_Edge_Triangle_Intersection(simplices,3,shared_edge,all_weights);
                            Register_Cut_Intersection(converted_simplices,all_weights,old_intersection);goto NEXT_OLD_SIMPLEX;}}}}
            // it doesn't exist, check intersection
            bool intersects=SIMPLEX_INTERACTIONS<T>::Intersection(cutting_simplices->simplices(simplices(3)).weights,
                cutting_simplices->Weight_For_Nodes_In_Simplex(simplices(1),shared_edge));
            if(intersects){
                VECTOR<VECTOR<T,2>,3> all_weights;
                Get_Simplex_Weights_For_Edge_Triangle_Intersection(simplices,3,shared_edge,all_weights);
                Register_Cut_Intersection(converted_simplices,all_weights,0);}
            goto NEXT_OLD_SIMPLEX;}}
        // Case 4: normal intersection
        {for(int i=1;i<=3;i++) if(cutting_simplices->Simplex_Is_Fake(simplices(i))) return; // we know not all are fake from above, but if any are fake we cannot normal intersect
        VECTOR<VECTOR<VECTOR<T,3>,3>,3> weights_for_simplices;
        for(int i=1;i<=3;i++) weights_for_simplices[i]=cutting_simplices->simplices(simplices[i]).weights;
        VECTOR<VECTOR<T,2>,3> all_weights;
        bool intersects=SIMPLEX_INTERACTIONS<T>::Three_Triangle_Intersection_Barycentric_Coordinates(weights_for_simplices[1],weights_for_simplices[2],weights_for_simplices[3],
            all_weights[1],all_weights[2],all_weights[3]);
        if(intersects) Register_Cut_Intersection(converted_simplices,all_weights,0);}
      NEXT_OLD_SIMPLEX:;   
    }
}
//#####################################################################
// Function Get_Unoriented_Segments_On_Two_Simplices
//#####################################################################
template<class TV,int d_input> void CUTTING_GEOMETRY_3D<TV,d_input>::
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
// Function Get_Unoriented_Segments_On_Two_Simplices
//#####################################################################
template<class TV,int d_input> void CUTTING_GEOMETRY_3D<TV,d_input>::
Divide_Polygon_Particles_With_New_Segments(ARRAY<VECTOR<int,2> >& all_segments,const ARRAY<VECTOR<int,2> >& possible_segments_to_add,const ARRAY<ARRAY<int> >& polygon_particles,
    const int cutting_simplex,const bool flipped,ARRAY<ARRAY<ARRAY<int > > >& final_polygon_element_particles) const
{
    // TODO: Continue here
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
// Function Potential_Segment_Should_Be_Added_To_Polygon
//#####################################################################
template<class TV,int d_input> bool CUTTING_GEOMETRY_3D<TV,d_input>::
Potential_Segment_Should_Be_Added_To_Polygon(const ARRAY<ARRAY<int > >& particles_for_polygon,const bool flipped,const int simplex,const VECTOR<int,2>& nodes) const
{
    // TODO: Continue here
    // check inside/outside rings
    VECTOR<bool,2> nodes_on_boundary;for(int i=1;i<=2;i++) for(int j=1;j<=particles_for_polygon.m;j++) if(particles_for_polygon(j).Contains(nodes(i))) nodes_on_boundary[i]=true;
    if(nodes_on_boundary.Contains(false)){
        for(int i=1;i<=2;i++) if(!nodes_on_boundary[i]) for(int j=1;j<=particles_for_polygon.m;j++)
            if((j==1)^Point_Is_Inside_Unoriented_Polygon(particles_for_polygon(j),simplex,nodes[i])) return false; // first loop is outer ring and others are negative holes
        return true;}
    else{
        // both endpoints of segment are on existing points, so operate on nodes[1] arbitrarily
        for(int i=1;i<=particles_for_polygon.m;i++){const ARRAY<int>& particles=particles_for_polygon(i);
            if(int j=particles.Find(nodes[1])){ // TODO: this seems to be just a call to Point_Is_Inside_Unoriented_Polygon
                int p1=particles(j-1<1?particles.m:j-1),p2=particles(j),p3=particles(j+1>particles.m?1:j+1);
                VECTOR<T,2> X1=intersection_registry->Get_Simplex_Weights_Of_Intersection(p1,simplex),
                    X2=intersection_registry->Get_Simplex_Weights_Of_Intersection(p2,simplex),
                    X3=intersection_registry->Get_Simplex_Weights_Of_Intersection(p3,simplex),
                    X=intersection_registry->Get_Simplex_Weights_Of_Intersection(nodes[2],simplex),
                    v1=X2-X1,v2=X3-X2,v3=X-X2;
                T angle_1=VECTOR<T,2>::Oriented_Angle_Between(v2,v1),angle_2=VECTOR<T,2>::Oriented_Angle_Between(v3,v1);
                return (i>1)^flipped^(angle_1>angle_2);}} // first loop is outer ring and others are negative holes
        PHYSBAM_FATAL_ERROR("Point should not be reached");}
    return false;
}
//#####################################################################
// Function Point_Is_Inside_Unoriented_Polygon
//#####################################################################
template<class TV,int d_input> bool CUTTING_GEOMETRY_3D<TV,d_input>::
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
template<class TV,int d_input> void CUTTING_GEOMETRY_3D<TV,d_input>::
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
                            if(direction.Magnitude()<1e-12) PHYSBAM_FATAL_ERROR();
                            if(next_direction.Magnitude()<1e-12) PHYSBAM_FATAL_ERROR();
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
                        if(seg_to_examine_direction.Magnitude()<1e-12) PHYSBAM_FATAL_ERROR();
                        if(candidate_seg_direction.Magnitude()<1e-12) PHYSBAM_FATAL_ERROR();
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
                        if(seg_to_examine_direction.Magnitude()<1e-12) PHYSBAM_FATAL_ERROR();
                        if(candidate_seg_direction.Magnitude()<1e-12) PHYSBAM_FATAL_ERROR();
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
template<class TV,int d_input> void CUTTING_GEOMETRY_3D<TV,d_input>::
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
// Function Draw_Polygon
//#####################################################################
template<class TV,int d_input> void CUTTING_GEOMETRY_3D<TV,d_input>::
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
        ss<<hsb<<" 1 1 sethsbcolor fill"<<"% END RUN"<<std::endl;LOG::filecout(ss.str());}
}
//#####################################################################
template class CUTTING_GEOMETRY_3D<VECTOR<float,3>,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class CUTTING_GEOMETRY_3D<VECTOR<double,3>,3>;
#endif

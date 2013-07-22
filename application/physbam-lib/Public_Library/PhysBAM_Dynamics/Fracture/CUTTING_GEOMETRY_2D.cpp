//#####################################################################
// Copyright 2007, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
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
#include <PhysBAM_Dynamics/Fracture/CUTTING_GEOMETRY_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Intersect_Simplex_With_Old_Simplices_in_Embedding
//#####################################################################
template<class TV,int d_input> void CUTTING_GEOMETRY_2D<TV,d_input>::
Intersect_Simplex_With_Old_Simplices_In_Embedding(const int tri,const int new_simplex)
{
    ARRAY<int>& old_simplices=simplices_per_current_embedding_simplex(tri);
    for(int dummy_i=1;dummy_i<=old_simplices.m;dummy_i++) {
        int old_simplex=old_simplices(dummy_i);VECTOR<int,2> simplices(old_simplex,new_simplex);
        if(verbose) {std::stringstream ss;ss<<" Checking old simplices "<<old_simplex<<" with new simplex "<<new_simplex<<std::endl;LOG::filecout(ss.str());}
        bool converted;VECTOR<int,2> converted_simplices=cutting_simplices->Convert_Simplex_Indices_To_Global_Indices(simplices,converted);

        // skip intersection if both simplices are fake or already found in intersection list
        if(cutting_simplices->Simplex_Is_Fake(old_simplex) && cutting_simplices->Simplex_Is_Fake(new_simplex)) continue;
        {ARRAY<int> intersection_list;if(intersection_registry->Intersection_List(converted_simplices,intersection_list)) continue;}
        
        // Case 1 - Sharing of 1 or more nodes
        ARRAY<int> points_all_shared;cutting_simplices->Shared_Nodes_On_Simplices(converted_simplices,points_all_shared);
        if(points_all_shared.m>=2) continue; // multiple nodes means no intersection at a point
        else if(points_all_shared.m==1){int shared_node=points_all_shared(1);
            assert(simplices==converted_simplices);
            VECTOR<VECTOR<T,1>,2> all_weights;for(int i=1;i<=2;i++) all_weights(i)=Get_Node_Weights_From_Face(shared_node,cutting_simplices->simplices(simplices(i)).nodes);
            // skip intersection if it is outside the tet
            const VECTOR<VECTOR<T,2>,2>& simplex_weights_in_embedding=cutting_simplices->simplices(simplices(1)).weights;
            VECTOR<T,2> shared_node_weights_in_embedding=simplex_weights_in_embedding(1)*all_weights(1)[1]+simplex_weights_in_embedding(2)*(1-all_weights(1)[1]);
            if(shared_node_weights_in_embedding.Min()<0 || shared_node_weights_in_embedding.Sum()>(T)1) goto NEXT_OLD_SIMPLEX;
            // skip intersection if already added under a different name
            for(int i=1;i<=old_simplices.m;i++) if(cutting_simplices->simplices(old_simplices(i)).nodes.Contains(shared_node))
                for(int j=i+1;j<=old_simplices.m;j++) if(cutting_simplices->simplices(old_simplices(j)).nodes.Contains(shared_node)){
                    ARRAY<int> nodes_shared_on_pair;VECTOR<int,2> alternative_simplex=VECTOR<int,2>(old_simplices(i),old_simplices(j));
                    cutting_simplices->Shared_Nodes_On_Simplices(alternative_simplex,nodes_shared_on_pair);
                    if(nodes_shared_on_pair.m>1) continue; // intersection is not a point so skip
                    if(int intersection_index=intersection_registry->Intersection(alternative_simplex))
                        Register_Cut_Intersection(converted_simplices,all_weights,intersection_index);
                    goto NEXT_OLD_SIMPLEX;}
            Register_Cut_Intersection(converted_simplices,all_weights,0);goto NEXT_OLD_SIMPLEX;}
        // Case 2 - Normal Intersection
        if(cutting_simplices->Simplex_Is_Fake(old_simplex) || cutting_simplices->Simplex_Is_Fake(new_simplex)) goto NEXT_OLD_SIMPLEX;
        {VECTOR<VECTOR<VECTOR<T,2>,2>,2> weights_for_simplices;
        for(int i=1;i<=2;i++) weights_for_simplices[i]=cutting_simplices->simplices(simplices[i]).weights;
        VECTOR<VECTOR<T,1>,2> intersection_weights_in_segments;
        if(SIMPLEX_INTERACTIONS<T>::Two_Segment_Intersection_Barycentric_Coordinates(weights_for_simplices[1],weights_for_simplices[2],
                intersection_weights_in_segments[1],intersection_weights_in_segments[2]))
            Register_Cut_Intersection(converted_simplices,intersection_weights_in_segments,0);}
      NEXT_OLD_SIMPLEX:;   
    }
}
//#####################################################################
// Function Split_Existing_Polygons
//#####################################################################
template<class TV,int d_input> void CUTTING_GEOMETRY_2D<TV,d_input>::
Split_Existing_Polygon_Edges()
{

}
//#####################################################################
// Function Split_Existing_Polygons_Edges
//#####################################################################
template<class TV,int d_input> void CUTTING_GEOMETRY_2D<TV,d_input>::
Split_Existing_Polygons()
{
    polygons_per_element=CONSTANT_ARRAY<ARRAY<int> >(current_embedding->mesh.elements.m,ARRAY<int>());
    for(int i=1;i<=cutting_simplices->index_for_last_old_cutting_simplex;i++) if(!cutting_simplices->Simplex_Is_Parent(i)){
        int tet_owner=cutting_simplices->simplices(i).element_owner;
        // obtain all new segments on simplex
        ARRAY<int> new_particles_on_simplex;
        // TODO: this next for makes quadratic time
        ARRAY<int> temp_particles;
        for(int j=cutting_simplices->index_for_last_old_cutting_simplex+1;j<=cutting_simplices->simplices.m;j++) if(!cutting_simplices->Simplex_Is_Fake(j)){
            temp_particles.Remove_All();Get_Particles_On_Simplices(VECTOR<int,2>(i,j),temp_particles);
            new_particles_on_simplex.Append_Elements(temp_particles);}
        std::stringstream ss;ss<<"new particles on simplex "<<new_particles_on_simplex<<std::endl;LOG::filecout(ss.str());
        // operate on each polygon on the simplex
        const ARRAY<int>& cutting_polygons_on_simplex=cutting_polygons_per_cutting_simplex(i);
        for(int j=1;j<=cutting_polygons_on_simplex.m;j++){
            int cutting_polygon_index=cutting_polygons_on_simplex(j);
            if(!new_particles_on_simplex.m){
                polygons_per_element(tet_owner).Append(cutting_polygon_index);}
            else{
                const CUTTING_POLYGON& cutting_polygon=current_cutting_polygons(cutting_polygon_index);
                if(cutting_polygon.simplex_owner!=i) PHYSBAM_FATAL_ERROR("TODO: remove this check");
                ARRAY<VECTOR<int,2> > all_polygonal_segments_on_simplex;
                for(int run=1;run<=polygon_mesh.elements(cutting_polygon.polygon_index).m;run++){
                    ARRAY<int>& pair=polygon_mesh.elements(cutting_polygon.polygon_index)(run);
                    assert(pair.m==2);
                    all_polygonal_segments_on_simplex.Append(VECTOR<int,2>(pair(1),pair(2)));}
                ARRAY<ARRAY<ARRAY<int > > > final_polygon_element_particles;
                Divide_Polygon_Particles_With_New_Segments(all_polygonal_segments_on_simplex,new_particles_on_simplex,
                    polygon_mesh.elements(cutting_polygon.polygon_index),i,cutting_polygon.flipped,final_polygon_element_particles);
                // replace old one with first new
                polygon_mesh.elements(cutting_polygon.polygon_index)=final_polygon_element_particles(1);
                polygons_per_element(tet_owner).Append(cutting_polygon_index);
                // make all the rest
                for(int k=2;k<=final_polygon_element_particles.m;k++){
                    int polygon_element_index=polygon_mesh.elements.Append(final_polygon_element_particles(k));
                    int new_cutting_polygon_index=current_cutting_polygons.Append(CUTTING_POLYGON(polygon_element_index,i,cutting_polygon.flipped,cutting_polygon.polygon_type));
                    polygons_per_element(tet_owner).Append(new_cutting_polygon_index);}}}}
    if(verbose){std::stringstream ss1;ss1<<std::endl;for(int i=1;i<=current_cutting_polygons.m;i++)
        ss1<<"Cutting polygon "<<i<<" has polygon element "
                 <<current_cutting_polygons(i).polygon_index<<": "<<polygon_mesh.elements(current_cutting_polygons(i).polygon_index)
                 <<std::endl;LOG::filecout(ss1.str());}
}
//#####################################################################
// Function Divide_Polygon_Particles_With_New_Segments
//#####################################################################
template<class TV,int d_input> void CUTTING_GEOMETRY_2D<TV,d_input>::
Divide_Polygon_Particles_With_New_Segments(ARRAY<VECTOR<int,2> >& all_segments,const ARRAY<int>& possible_particles_to_add,const ARRAY<ARRAY<int> >& polygon_particles,
    const int cutting_simplex,const bool flipped,ARRAY<ARRAY<ARRAY<int > > >& final_polygon_element_particles) const
{
    if(verbose){std::stringstream ss;ss<<"All particles on simplex "<<cutting_simplex<<" before adding new ";for(int j=1;j<=all_segments.m;j++) ss<<all_segments(j)<<"; ";ss<<std::endl;ss<<"     new potential particles are "<<possible_particles_to_add<<std::endl;LOG::filecout(ss.str());}
    // choose subset of particles to add to this polygon
    for(int k=1;k<=possible_particles_to_add.m;k++){const int possible_particle_to_add=possible_particles_to_add(k);
        for(int i=1;i<=all_segments.m;i++){
            if(all_segments(i).Contains(possible_particle_to_add)) goto NEXT_POSSIBLE_PARTICLES;
            int node1=all_segments(i)(1),node2=all_segments(i)(2);
            T X=intersection_registry->Get_Simplex_Weights_Of_Intersection(possible_particle_to_add,cutting_simplex).x,
                X1=intersection_registry->Get_Simplex_Weights_Of_Intersection(node1,cutting_simplex).x,
                X2=intersection_registry->Get_Simplex_Weights_Of_Intersection(node2,cutting_simplex).x;
            if(verbose) PHYSBAM_DEBUG_PRINT("   trying to split ",flipped,possible_particle_to_add,node1,node2,X,X1,X2);
            if(flipped ? X>X1 && X<X2 : X<X1 && X>X2){
                all_segments(i)(2)=possible_particle_to_add;all_segments.Append(VECTOR<int,2>(possible_particle_to_add,node2));}}
      NEXT_POSSIBLE_PARTICLES:;
    }
    if(verbose){std::stringstream ss;ss<<"All particles on simplex "<<cutting_simplex<<" after adding new ";for(int j=1;j<=all_segments.m;j++) ss<<all_segments(j)<<"; ";ss<<std::endl;LOG::filecout(ss.str());}
    for(int i=1;i<=all_segments.m;i++){
        int polygon=final_polygon_element_particles.Append(ARRAY<ARRAY<int> >());
        int component=final_polygon_element_particles(polygon).Append(ARRAY<int>());
        final_polygon_element_particles(polygon)(component).Append(all_segments(i)(1));
        final_polygon_element_particles(polygon)(component).Append(all_segments(i)(2));}
}
//#####################################################################
template class CUTTING_GEOMETRY_2D<VECTOR<float,2>,2>;
template class CUTTING_GEOMETRY_2D<VECTOR<float,3>,2>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class CUTTING_GEOMETRY_2D<VECTOR<double,2>,2>;
template class CUTTING_GEOMETRY_2D<VECTOR<double,3>,2>;
#endif

//#####################################################################
// Copyright 2007, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PAINTED_GRAIN_BOUNDARIES
//##################################################################### 
#ifndef __PAINTED_GRAIN_BOUNDARIES__
#define __PAINTED_GRAIN_BOUNDARIES__

#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_CALLBACKS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_GRAIN_BOUNDARIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_OBJECT.h>
namespace PhysBAM{

template<class TV,int d>
class PAINTED_GRAIN_BOUNDARIES:public FRACTURE_GRAIN_BOUNDARIES<TV,d>
{
    typedef typename TV::SCALAR T;

public:
    int total_number_of_regions;
    HASHTABLE<VECTOR<int,2>,int> edge_hashtable;
    ARRAY<VECTOR<T,2> > edge_ratios;
    ARRAY<VECTOR<T,3> > particle_positions;

    using FRACTURE_GRAIN_BOUNDARIES<TV,d>::mesh;using FRACTURE_GRAIN_BOUNDARIES<TV,d>::seed_positions;using FRACTURE_GRAIN_BOUNDARIES<TV,d>::seed_weakness_multipliers;
    using FRACTURE_GRAIN_BOUNDARIES<TV,d>::seed_weakness_multipliers_callback;using FRACTURE_GRAIN_BOUNDARIES<TV,d>::node_smallest_distance;using FRACTURE_GRAIN_BOUNDARIES<TV,d>::is_breakable;
    using FRACTURE_GRAIN_BOUNDARIES<TV,d>::node_region;using FRACTURE_GRAIN_BOUNDARIES<TV,d>::fracture_callbacks;

    using FRACTURE_GRAIN_BOUNDARIES<TV,d>::Number_Of_Nodes_In_Region;

    PAINTED_GRAIN_BOUNDARIES(const PARTICLES<TV>& particles,const SIMPLEX_MESH<d>& mesh_input,const ARRAY<TV>& seed_positions_input,const ARRAY<T>& seed_weakness_multipliers_input,const FRAME<TV> frame,
                              const FRACTURE_CALLBACKS<TV>* fracture_callbacks_input=0)
        :FRACTURE_GRAIN_BOUNDARIES<TV,d>(particles,mesh_input,seed_positions_input,seed_weakness_multipliers_input,frame,fracture_callbacks_input),total_number_of_regions(0)
    {
        particle_positions.Resize(particles.array_collection->Size());
        for(int i=1;i<=particles.array_collection->Size();i++) particle_positions(i)=frame*particles.X(i);
        Calculate_Grain_Boundaries(particles.X,frame,node_smallest_distance,node_region,total_number_of_regions,fracture_callbacks);
    }

    void Calculate_Grain_Boundaries(ARRAY_VIEW<const TV> positions,const FRAME<TV> frame,ARRAY<T>& node_smallest_distance,ARRAY<int>& node_region,int& total_number_of_regions,
                                           const FRACTURE_CALLBACKS<TV>* fracture_callbacks=0)
    {
        assert(fracture_callbacks); // we have to have fracture callbacks for this one to access the field
        
        int number_of_positions=positions.Size();
        node_smallest_distance.Resize(number_of_positions);node_region.Resize(number_of_positions);
        fracture_callbacks->Node_Regions(positions,node_region,frame);
            
        int max_node_region=ARRAYS_COMPUTATIONS::Max(node_region);
        ARRAY<int> regions_selected(max_node_region);
        for(int node=1;node<=number_of_positions;node++){
            node_smallest_distance(node) = 1; // this should be fixed
            regions_selected(node_region(node))=1;}

        // Consolidate the regions
        ARRAY<int> region_mapping;
        for(int region=1;region<=max_node_region;region++) 
            if(regions_selected(region)) {
                region_mapping.Append(region);
                regions_selected(region)=region_mapping.m;}
        total_number_of_regions=region_mapping.m;
        for(int node=1;node<=number_of_positions;node++) node_region(node)=regions_selected(node_region(node));

        // refine the node_smallest_distance
        int number_of_edge_divisions=10;
        T divisor=(T)1/(number_of_edge_divisions+1);
        
        int number_of_edges=0;
        for(int tet=1;tet<=mesh.elements.m;tet++){
            for(int v1=1;v1<=d;v1++) for(int v2=v1+1;v2<=d+1;v2++){
                int vertex1=mesh.elements(tet)[v1],vertex2=mesh.elements(tet)[v2];
                if(node_region(vertex1)!=node_region(vertex2)){ // we have a crossover edge
                    VECTOR<int,2> edge=VECTOR<int,2>(vertex1,vertex2).Sorted();
                    if(!edge_hashtable.Contains(edge)) edge_hashtable.Insert(edge,++number_of_edges);}}}
        ARRAY<VECTOR<T,3> > edge_divisions(number_of_edges*number_of_edge_divisions);
        ARRAY<int> edge_division_regions(number_of_edges*number_of_edge_divisions);
        for(HASHTABLE_ITERATOR<VECTOR<int,2>,int> iterator(edge_hashtable);iterator.Valid();iterator.Next()){
            int current_index=(iterator.Data()-1)*number_of_edge_divisions;
            VECTOR<int,2> edge_vertices=iterator.Key();
            for(int division=1;division<=number_of_edge_divisions;division++){
                T t=division*divisor;
                edge_divisions(current_index+division)=positions(edge_vertices.x)*t+positions(edge_vertices.y)*(1-t);}}

        fracture_callbacks->Node_Regions(edge_divisions,edge_division_regions,frame);

        edge_ratios.Resize(number_of_edges);
        for(HASHTABLE_ITERATOR<VECTOR<int,2>,int> iterator(edge_hashtable);iterator.Valid();iterator.Next()){
            VECTOR<int,2> edge_vertices=iterator.Key();
            VECTOR<T,3> vertex1=positions(edge_vertices.x),vertex2=positions(edge_vertices.y);
            int region1=node_region(edge_vertices.x),region2=node_region(edge_vertices.y);
            int region1count=0,region2count=0;
            int current_index=(iterator.Data()-1)*number_of_edge_divisions;
            for(int division=1;division<=number_of_edge_divisions;division++){
                if(regions_selected(edge_division_regions(current_index+division))==region1) region1count++;
                if(regions_selected(edge_division_regions(current_index+division))==region2) region2count++;}
            T vertex1dist=(T)(region1count+(T)1)/(number_of_edge_divisions+(T)2);
            T vertex2dist=(T)(region2count+(T)1)/(number_of_edge_divisions+(T)2);
            //{std::stringstream ss;ss<<"vertex1dist: "<<vertex1dist<<", vertex2dist: "<<vertex2dist<<std::endl;LOG::filecout(ss.str());}
            edge_ratios(iterator.Data())=VECTOR<T,2>(vertex1dist,vertex2dist);
            node_smallest_distance(edge_vertices.x)=min(node_smallest_distance(edge_vertices.x),vertex1dist);
            node_smallest_distance(edge_vertices.y)=min(node_smallest_distance(edge_vertices.y),vertex2dist);}
    }

    void Phi_For_Region_In_Element(const int element,const int region,VECTOR<T,d+1>& phi)
    {
        ARRAY<VECTOR<T,3> > points;
        for(int v1=1;v1<=d;v1++){// loop through each node
            int current_node=mesh.elements(element)[v1];
            for(int v2=v1+1;v2<=d+1;v2++){ // loop through the opposite nodes
                if(v2!=v1){ 
                    int opposite_node=mesh.elements(element)[v2];
                    VECTOR<int,2> edge=VECTOR<int,2>(current_node,opposite_node).Sorted();
                    int edge_index;
                    if(edge_hashtable.Get(edge,edge_index)){ // if the edge is a cut edge, then we want to calculate the point
                        T normalizer=(edge_ratios(edge_index).x+edge_ratios(edge_index).y);
                        VECTOR<T,3> blend_point=edge_ratios(edge_index).y/normalizer*particle_positions(edge.x)+edge_ratios(edge_index).x/normalizer*particle_positions(edge.y);
                        points.Append(blend_point);}}}}
        assert(points.m>=3);
        PLANE<T> cut_plane(points(1),points(2),points(3));
        for(int v=1;v<=d+1;v++){
            int current_node=mesh.elements(element)[v];
            phi[v]=fabs(cut_plane.Signed_Distance(particle_positions(current_node)));
            if(node_region(current_node)==region) phi[v]=-phi[v];}
    }

    int Number_Of_Regions()
    {return total_number_of_regions;}

//#####################################################################
};  
}
#endif

//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Geoffrey Irving, Neil Molino, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRAIN_BOUNDARIES
//##################################################################### 
#ifndef __GRAIN_BOUNDARIES__
#define __GRAIN_BOUNDARIES__

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_OBJECT.h>
namespace PhysBAM{

template<class TV,int d>
class GRAIN_BOUNDARIES:public NONCOPYABLE
{
private:
    typedef typename TV::SCALAR T;

    FRACTURE_OBJECT<TV,d>& fracture_object;
    bool adjacent_elements_defined,neighbor_nodes_defined;
    RANDOM_NUMBERS<T> random_numbers;
    SIMPLEX_MESH<d>& mesh;
    ARRAY<int> regions;
    ARRAY<int> node_regions;
public:
    ARRAY<int>* seed_element_index;
    ARRAY<int>* seed_node_index;

    GRAIN_BOUNDARIES(FRACTURE_OBJECT<TV,d>& fracture_object_input,SIMPLEX_MESH<d>& mesh_input)
        :fracture_object(fracture_object_input),mesh(mesh_input),seed_element_index(0),seed_node_index(0)
    {
        adjacent_elements_defined=mesh.adjacent_elements!=0;if(!adjacent_elements_defined) mesh.Initialize_Adjacent_Elements();
        neighbor_nodes_defined=mesh.neighbor_nodes!=0;if(!neighbor_nodes_defined) mesh.Initialize_Neighbor_Nodes();
        random_numbers.Set_Seed(9765);
    }

    ~GRAIN_BOUNDARIES()
    {delete seed_element_index;delete seed_node_index;
    if(!adjacent_elements_defined){delete mesh.adjacent_elements;mesh.adjacent_elements=0;}
    if(!neighbor_nodes_defined){delete mesh.neighbor_nodes;mesh.neighbor_nodes=0;}}

    void Set_Random_Fracture_Bias_Directions()
    {for(int t=1;t<=mesh.elements.m;t++) fracture_object.fracture_bias_direction(t)=random_numbers.template Get_Direction<VECTOR<T,d> >();}

    void Smooth_Fracture_Bias_Directions(const int passes)
    {for(int iteration=1;iteration<=passes;iteration++) for(int t=1;t<=mesh.elements.m;t++){
        VECTOR<T,d> average=fracture_object.fracture_bias_direction(t);
        for(int i=1;i<=(*mesh.adjacent_elements)(t).m;i++) average+=fracture_object.fracture_bias_direction((*mesh.adjacent_elements)(t)(i));
        fracture_object.fracture_bias_direction(t)=average.Normalized();}}

    void Seed_Regions(const int number_of_regions)
    {regions.Resize(mesh.elements.m);
    for(int r=1;r<=number_of_regions;r++){
        int element=seed_element_index?(*seed_element_index)(r):random_numbers.Get_Uniform_Integer(1,mesh.elements.m);
        while(regions(element))element=random_numbers.Get_Uniform_Integer(1,mesh.elements.m);
        regions(element)=r;}}

    void Seed_Node_Regions(const int number_of_node_regions)
    {node_regions.Resize(mesh.number_nodes);
    for(int r=1;r<=number_of_node_regions;r++){
        int node=seed_node_index?(*seed_node_index)(r):random_numbers.Get_Uniform_Integer(1,mesh.number_nodes);
        while(node_regions(node))node=random_numbers.Get_Uniform_Integer(1,mesh.number_nodes);
        node_regions(node)=r;}}

    int Mark_All_Neigbors_Of_Region(const int region,const T coin_flip_threshold=(T)0)
    {assert(regions.m==mesh.elements.m);
    ARRAY<bool> to_be_marked(mesh.elements.m);
    for(int t=1;t<=mesh.elements.m;t++) if(regions(t)==region) for(int a=1;a<=(*mesh.adjacent_elements)(t).m;a++){ // mark neighbors of t
        int adj_t=(*mesh.adjacent_elements)(t)(a);
        if(!coin_flip_threshold) to_be_marked(adj_t)=true;
        else if(random_numbers.Get_Uniform_Number((T)0,(T)1) > coin_flip_threshold)to_be_marked(adj_t)=true;}
    int number_marked=0;
    for(int t=1;t<=mesh.elements.m;t++) if(!regions(t) && to_be_marked(t)){regions(t)=region;number_marked++;}
    return number_marked;}

    int Mark_All_Node_Neighbors_Of_Region(const int node_region,const T coin_flip_threshold=(T)0)
    {assert(node_regions.m==mesh.number_nodes);
    ARRAY<bool> to_be_marked(mesh.number_nodes);
    for(int n=1;n<=mesh.number_nodes;n++) if(node_regions(n)==node_region) for(int nbr=1;nbr<=(*mesh.neighbor_nodes)(n).m;nbr++){
        int nbr_node=(*mesh.neighbor_nodes)(n)(nbr);
        if(!coin_flip_threshold)to_be_marked(nbr_node)=true;
        else if(random_numbers.Get_Uniform_Number((T)0,(T)1) > coin_flip_threshold)to_be_marked(nbr_node)=true;}
    int number_marked=0;
    for(int n=1;n<=mesh.number_nodes;n++) if(!node_regions(n) && to_be_marked(n)){node_regions(n)=node_region;number_marked++;}
    return number_marked;}

    int Unmarked_Element()
    {return regions.Find(0);}

    int Unmarked_Node()
    {return node_regions.Find(0);}

    void Floodfill_Regions(const int number_of_regions)
    {Seed_Regions(number_of_regions);
    while(Unmarked_Element()) for(int r=1;r<=number_of_regions;r++) Mark_All_Neigbors_Of_Region(r);}

    void Floodfill_Node_Regions(const int number_of_node_regions)
    {Seed_Node_Regions(number_of_node_regions);
    while(Unmarked_Node()) for(int r=1;r<=number_of_node_regions;r++) Mark_All_Node_Neighbors_Of_Region(r);}

    void Fill_Regions_With_Uniform_Vectors(const int number_of_regions)
    {Floodfill_Regions(number_of_regions);
    ARRAY<VECTOR<T,d> > random_vectors(number_of_regions);
    for(int r=1;r<=number_of_regions;r++) random_vectors(r)=random_numbers.template Get_Direction<VECTOR<T,d> >();
    for(int r=1;r<=regions.m;r++) fracture_object.fracture_bias_direction(r)=random_vectors(regions(r));}

    void Fill_Regions_With_Uniform_Spatial_Fracture_Bias_Directions(ARRAY<TV>& spatial_fracture_bias_direction,const int number_of_regions)
    {assert(spatial_fracture_bias_direction.m==mesh.elements.m);
    Floodfill_Regions(number_of_regions);
    ARRAY<TV> random_vectors(number_of_regions);
    for(int r=1;r<=number_of_regions;r++) random_vectors(r)=random_numbers.template Get_Direction<TV>();
    for(int r=1;r<=regions.m;r++) spatial_fracture_bias_direction(r)=random_vectors(regions(r));}

    void Fill_Node_Regions_With_Uniform_Vectors(const int number_of_node_regions)
    {Floodfill_Node_Regions(number_of_node_regions);
    ARRAY<VECTOR<T,d> > random_vectors(number_of_node_regions);
    for(int r=1;r<=number_of_node_regions;r++) random_vectors(r)=random_numbers.template Get_Direction<VECTOR<T,d> >();
    for(int t=1;t<=mesh.elements.m;t++) 
        if(!Element_Contains_Nodes_From_Different_Regions(t)) fracture_object.fracture_bias_direction(t)=random_vectors(node_regions(mesh.elements(t)(1)));
        else fracture_object.fracture_bias_direction(t)=VECTOR<T,d>();}

    bool Is_On_Border_With_Another_Region(const int element)
    {int region=regions(element);
    for(int a=1;a<=(*mesh.adjacent_elements)(element).m;a++) if(regions((*mesh.adjacent_elements)(element)(a)) != region) return true;
    return false;}

    bool Element_Contains_Nodes_From_Different_Regions(const int element)
    {VECTOR<int,d+1> element_node_regions(node_regions.Subset(mesh.elements(element)));
    for(int i=2;i<=element_node_regions.m;i++) if(element_node_regions[1]!=element_node_regions[i]) return true;
    return false;}

    int Number_In_Region(const int region)
    {int count=0;for(int t=1;t<=mesh.elements.m;t++) if(regions(t)==region) count++;return count;}

    void Print_Number_In_Regions(const int number_of_regions)
    {for(int r=1;r<=number_of_regions;r++) {std::stringstream ss;ss<<"there are "<<Number_In_Region(r)<<" in region "<<r<<std::endl;LOG::filecout(ss.str());}}

    void Print_Regions()
    {for(int t=1;t<=mesh.elements.m;t++) {std::stringstream ss;ss<<"region("<<t<<")="<<regions(t)<<std::endl;LOG::filecout(ss.str());}}

//#####################################################################
};  
}
#endif

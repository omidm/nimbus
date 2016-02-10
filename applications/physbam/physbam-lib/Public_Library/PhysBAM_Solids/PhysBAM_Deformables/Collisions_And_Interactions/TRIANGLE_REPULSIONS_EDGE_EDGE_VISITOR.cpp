//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Neil Molino, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR
//##################################################################### 
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_3D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY_DEFINITION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/STRUCTURE_INTERACTION_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<TV>::
TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR(ARRAY<EDGE_EDGE_REPULSION_PAIR<TV> >& pairs,const STRUCTURE_INTERACTION_GEOMETRY<TV>& edge_structure1,
    const STRUCTURE_INTERACTION_GEOMETRY<TV>& edge_structure2,ARRAY_VIEW<const TV> X_other,const TRIANGLE_REPULSIONS<TV>& repulsions,int& pruned)
    :pairs(pairs),edges1(edge_structure1.Edges()),edges2(edge_structure2.Edges()),
    X_other(X_other),X_self_collision_free(repulsions.geometry.X_self_collision_free),repulsion_thickness(repulsions.repulsion_thickness),
    thickness_multiplier(repulsions.repulsion_thickness_detection_multiplier*repulsions.hierarchy_repulsion_thickness_multiplier),pruned(pruned),
    perform_attractions(repulsions.perform_attractions),intersecting_edge_edge_pairs(repulsions.geometry.intersecting_edge_edge_pairs),
    omit_edge_edge_repulsion_pairs(repulsions.omit_edge_edge_repulsion_pairs)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<TV>::
~TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR()
{}
//#####################################################################
// Function Store_Helper
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<TV>::
Store_Helper(const int point1_index,const int point2_index,const VECTOR<T,2>&)
{
    int p1=edges1(point1_index),p2=edges2(point2_index);
    EDGE_EDGE_REPULSION_PAIR<TV> pair;pair.nodes=VECTOR<int,2>(p1,p2);
    POINT_2D<T> point1(X_other(p1)),point2(X_other(p2));
    T total_repulsion_thickness=thickness_multiplier*pair.Total_Repulsion_Thickness(repulsion_thickness);
    if(!point1.Edge_Edge_Interaction(point2,total_repulsion_thickness,pair.distance,pair.normal)) pruned++;
    else{
        if(omit_edge_edge_repulsion_pairs.Size() && omit_edge_edge_repulsion_pairs.Contains(pair.nodes)) return;
        if(intersecting_edge_edge_pairs.Size() && intersecting_edge_edge_pairs.Contains(pair.nodes)) return;
        if(&X_other!=&X_self_collision_free && perform_attractions){
            TV point1_collision_free(X_self_collision_free(p1)),point2_collision_free(X_self_collision_free(p2));
            pair.collision_free_normal=point1_collision_free-point2_collision_free;}
        else pair.collision_free_normal=pair.normal;
        pair.collision_free_normal.Normalize();
        pairs.Append(pair);}
}
//#####################################################################
// Function Store_Helper_Helper
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<TV>::
Store_Helper_Helper(const int segment1_index,const int segment2_index)
{
    const VECTOR<int,2> &segment1_nodes=edges1(segment1_index),&segment2_nodes=edges2(segment2_index);
    EDGE_EDGE_REPULSION_PAIR<TV> pair;pair.nodes=VECTOR<int,4>(segment1_nodes[1],segment1_nodes[2],segment2_nodes[1],segment2_nodes[2]);
    SEGMENT_3D<T> segment1(X_other.Subset(segment1_nodes)),segment2(X_other.Subset(segment2_nodes));
    T total_repulsion_thickness=thickness_multiplier*pair.Total_Repulsion_Thickness(repulsion_thickness);
    if(!segment1.Edge_Edge_Interaction(segment2,total_repulsion_thickness,pair.distance,pair.normal,pair.weights,false)) pruned++;
    else{
        if(omit_edge_edge_repulsion_pairs.Size() && omit_edge_edge_repulsion_pairs.Contains(pair.nodes)) return;
        if(intersecting_edge_edge_pairs.Size() && intersecting_edge_edge_pairs.Contains(pair.nodes)) return;
        if(&X_other!=&X_self_collision_free && perform_attractions){
            VECTOR<T,2> collision_free_weights;
            SEGMENT_3D<T> segment1_collision_free(X_self_collision_free.Subset(segment1_nodes)),segment2_collision_free(X_self_collision_free.Subset(segment2_nodes));
            pair.collision_free_normal=segment1_collision_free.Shortest_Vector_Between_Segments(segment2_collision_free,collision_free_weights);}
        else pair.collision_free_normal=pair.normal;
        pair.collision_free_normal.Normalize();
        pairs.Append(pair);}
}
//####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template void BOX_HIERARCHY<VECTOR<T,d> >::Intersection_List<TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<VECTOR<T,d> >,ZERO>(BOX_HIERARCHY<VECTOR<T,d> > const&, \
        TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<VECTOR<T,d> >&,ZERO) const; \
    template void BOX_HIERARCHY<VECTOR<T,d> >::Intersection_List<TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<VECTOR<T,d> >,T>(BOX_HIERARCHY<VECTOR<T,d> > const&, \
        TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<VECTOR<T,d> >&,T) const; \
    template void BOX_HIERARCHY<VECTOR<T,d> >::Intersection_List<BOX_VISITOR_MPI<TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<VECTOR<T,d> > >,ZERO>(BOX_HIERARCHY<VECTOR<T,d> > const&, \
        BOX_VISITOR_MPI<TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<VECTOR<T,d> > >&,ZERO) const; \
    template void BOX_HIERARCHY<VECTOR<T,d> >::Intersection_List<BOX_VISITOR_MPI<TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<VECTOR<T,d> > >,T>(BOX_HIERARCHY<VECTOR<T,d> > const&, \
        BOX_VISITOR_MPI<TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<VECTOR<T,d> > >&,T) const;

template void TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<VECTOR<float,2> >::Store_Helper(int,int,VECTOR<float,2> const&);
template TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<VECTOR<float,2> >::TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR(ARRAY<EDGE_EDGE_REPULSION_PAIR<VECTOR<float,2> >,int>&,
    STRUCTURE_INTERACTION_GEOMETRY<VECTOR<float,2> > const&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<float,2> > const&,ARRAY_VIEW<VECTOR<float,2> const,int>,
    TRIANGLE_REPULSIONS<VECTOR<float,2> > const&,int&);
template TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<VECTOR<float,2> >::~TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR();
template void TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<VECTOR<float,3> >::Store_Helper_Helper(int,int);
template TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<VECTOR<float,3> >::TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR(ARRAY<EDGE_EDGE_REPULSION_PAIR<VECTOR<float,3> >,int>&,
    STRUCTURE_INTERACTION_GEOMETRY<VECTOR<float,3> > const&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<float,3> > const&,ARRAY_VIEW<VECTOR<float,3> const,int>,
    TRIANGLE_REPULSIONS<VECTOR<float,3> > const&,int&);
template TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<VECTOR<float,3> >::~TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR();

INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<VECTOR<double,2> >::Store_Helper(int,int,VECTOR<double,2> const&);
template TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<VECTOR<double,2> >::TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR(ARRAY<EDGE_EDGE_REPULSION_PAIR<VECTOR<double,2> >,int>&,
    STRUCTURE_INTERACTION_GEOMETRY<VECTOR<double,2> > const&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<double,2> > const&,ARRAY_VIEW<VECTOR<double,2> const,int>,
    TRIANGLE_REPULSIONS<VECTOR<double,2> > const&,int&);
template TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<VECTOR<double,2> >::~TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR();
template void TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<VECTOR<double,3> >::Store_Helper_Helper(int,int);
template TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<VECTOR<double,3> >::TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR(ARRAY<EDGE_EDGE_REPULSION_PAIR<VECTOR<double,3> >,int>&,
    STRUCTURE_INTERACTION_GEOMETRY<VECTOR<double,3> > const&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<double,3> > const&,ARRAY_VIEW<VECTOR<double,3> const,int>,
    TRIANGLE_REPULSIONS<VECTOR<double,3> > const&,int&);
template TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<VECTOR<double,3> >::~TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR();

INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif

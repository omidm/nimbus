//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Neil Molino, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR
//##################################################################### 
#include <PhysBAM_Tools/Arrays_Computations/SUMMATIONS.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY_DEFINITION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/STRUCTURE_INTERACTION_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<TV>::
TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR(ARRAY<VECTOR<int,2*d-2> >& pairs_internal, ARRAY<VECTOR<int,2*d-2> >& pairs_external,
    const STRUCTURE_INTERACTION_GEOMETRY<TV>& edge_structure1,const STRUCTURE_INTERACTION_GEOMETRY<TV>& edge_structure2,
    const TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry,const T collision_thickness,MPI_SOLIDS<TV>* mpi_solids)
    :pairs_internal(pairs_internal),pairs_external(pairs_external),edges1(edge_structure1.Edges()),edges2(edge_structure2.Edges()),
    X(geometry.deformable_body_collection.particles.X),
    X_self_collision_free(geometry.X_self_collision_free),collision_thickness(collision_thickness),edge1_box_modified(edge_structure1.Edge_Modified()),
    edge2_box_modified(edge_structure2.Edge_Modified()),intersecting_edge_edge_pairs(geometry.intersecting_edge_edge_pairs),mpi_solids(mpi_solids)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<TV>::
~TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR()
{}
//#####################################################################
// Function Store_Helper
//#####################################################################
template<class TV> void TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<TV>::
Store_Helper(const int point1,const int point2,VECTOR<T,2>)
{
    int p1=edges1(point1),p2=edges2(point2);
    TV dX_average=(X(p1)-X_self_collision_free(p1))+(X(p2)-X_self_collision_free(p2));
    dX_average*=(T).5; // TODO: This is separated because of compiler bugs
    RANGE<TV> box1=RANGE<TV>(X_self_collision_free(p1))+dX_average;box1.Enlarge_Nonempty_Box_To_Include_Point(X(p1));
    RANGE<TV> box2=RANGE<TV>(X_self_collision_free(p2))+dX_average;box2.Enlarge_Nonempty_Box_To_Include_Point(X(p2));
    if(!box1.Intersection(box2,collision_thickness)) return;
    VECTOR<int,2> nodes(p1,p2);
    if(intersecting_edge_edge_pairs.Size() && intersecting_edge_edge_pairs.Contains(nodes)) return;
    if (mpi_solids){
        VECTOR<PARTITION_ID,2> processors(mpi_solids->partition_id_from_particle_index.Subset(nodes));
        if (processors(1)!=processors(2)) pairs_external.Append(nodes);
        else pairs_internal.Append(nodes);}
    else pairs_internal.Append(nodes);
}
//#####################################################################
// Function Store_Helper_Helper
//#####################################################################
template<class TV> void TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<TV>::
Store_Helper_Helper(const int segment1,const int segment2)
{
    const VECTOR<int,2> &segment1_nodes=edges1(segment1),&segment2_nodes=edges2(segment2);
    TV dX_average=(ARRAYS_COMPUTATIONS::Average(X.Subset(segment1_nodes)) // TODO: optimize scalar multiplications
                   -ARRAYS_COMPUTATIONS::Average(X_self_collision_free.Subset(segment1_nodes))
                   +ARRAYS_COMPUTATIONS::Average(X.Subset(segment2_nodes))
                   -ARRAYS_COMPUTATIONS::Average(X_self_collision_free.Subset(segment2_nodes)));
    dX_average*=(T).5; // TODO: This is separated because of compiler bugs
    RANGE<TV> box1=RANGE<TV>::Bounding_Box(X_self_collision_free.Subset(segment1_nodes))+dX_average;box1.Enlarge_Nonempty_Box_To_Include_Points(X.Subset(segment1_nodes));
    RANGE<TV> box2=RANGE<TV>::Bounding_Box(X_self_collision_free.Subset(segment2_nodes))+dX_average;box2.Enlarge_Nonempty_Box_To_Include_Points(X.Subset(segment2_nodes));
    if(!box1.Intersection(box2,collision_thickness)) return;
    VECTOR<int,4> nodes(segment1_nodes[1],segment1_nodes[2],segment2_nodes[1],segment2_nodes[2]);
    if(intersecting_edge_edge_pairs.Size() && intersecting_edge_edge_pairs.Contains(nodes)) return;
    if (mpi_solids){
        VECTOR<PARTITION_ID,4> processors(mpi_solids->partition_id_from_particle_index.Subset(nodes));
        for(int i=1;i<=3;i++) if (processors(i)!=processors(4)) {pairs_external.Append(nodes);return;}
        pairs_internal.Append(nodes);}
    else pairs_internal.Append(nodes);
}
//####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template void BOX_HIERARCHY<VECTOR<T,d> >::Intersection_List<BOX_VISITOR_MPI<TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<VECTOR<T,d> > >,T>(BOX_HIERARCHY<VECTOR<T,d> > const&, \
        BOX_VISITOR_MPI<TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<VECTOR<T,d> > >&,T) const; \
    template void BOX_HIERARCHY<VECTOR<T,d> >::Intersection_List<TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<VECTOR<T,d> >,T>(BOX_HIERARCHY<VECTOR<T,d> > const&, \
        TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<VECTOR<T,d> >&,T) const;

template TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<VECTOR<float,2> >::TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR(ARRAY<VECTOR<int,2>,int>&,ARRAY<VECTOR<int,2>,int>&,
    STRUCTURE_INTERACTION_GEOMETRY<VECTOR<float,2> > const&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<float,2> > const&,
    TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<VECTOR<float,2> > const&,float,MPI_SOLIDS<VECTOR<float,2> >*);
template TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<VECTOR<float,2> >::~TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR();
template TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<VECTOR<float,3> >::TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR(ARRAY<VECTOR<int,4>,int>&,ARRAY<VECTOR<int,4>,int>&,
    STRUCTURE_INTERACTION_GEOMETRY<VECTOR<float,3> > const&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<float,3> > const&,
    TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<VECTOR<float,3> > const&,float,MPI_SOLIDS<VECTOR<float,3> >*);
template TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<VECTOR<float,3> >::~TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR();
template void TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<VECTOR<float,2> >::Store_Helper(int,int,VECTOR<float,2>);
template void TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<VECTOR<float,3> >::Store_Helper_Helper(int,int);

INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<VECTOR<double,2> >::TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR(ARRAY<VECTOR<int,2>,int>&,ARRAY<VECTOR<int,2>,int>&,
    STRUCTURE_INTERACTION_GEOMETRY<VECTOR<double,2> > const&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<double,2> > const&,
    TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<VECTOR<double,2> > const&,double,MPI_SOLIDS<VECTOR<double,2> >*);
template TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<VECTOR<double,2> >::~TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR();
template TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<VECTOR<double,3> >::TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR(ARRAY<VECTOR<int,4>,int>&,ARRAY<VECTOR<int,4>,int>&,
    STRUCTURE_INTERACTION_GEOMETRY<VECTOR<double,3> > const&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<double,3> > const&,
    TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<VECTOR<double,3> > const&,double,MPI_SOLIDS<VECTOR<double,3> >*);
template TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<VECTOR<double,3> >::~TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR();
template void TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<VECTOR<double,2> >::Store_Helper(int,int,VECTOR<double,2>);
template void TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<VECTOR<double,3> >::Store_Helper_Helper(int,int);

INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif

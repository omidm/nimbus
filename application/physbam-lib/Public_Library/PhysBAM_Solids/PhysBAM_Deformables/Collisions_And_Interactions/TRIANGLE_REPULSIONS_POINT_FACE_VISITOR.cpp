//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Neil Molino, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_REPULSIONS_POINT_FACE_VISITOR
//##################################################################### 
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY_DEFINITION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/STRUCTURE_INTERACTION_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_POINT_FACE_VISITOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> TRIANGLE_REPULSIONS_POINT_FACE_VISITOR<TV>::
TRIANGLE_REPULSIONS_POINT_FACE_VISITOR(ARRAY<POINT_FACE_REPULSION_PAIR<TV> >& pairs,const STRUCTURE_INTERACTION_GEOMETRY<TV>& particle_structure,
    const STRUCTURE_INTERACTION_GEOMETRY<TV>& face_structure,ARRAY_VIEW<const TV> X_other,const TRIANGLE_REPULSIONS<TV>& repulsions,int& pruned)
    :pairs(pairs),particle_active_indices(particle_structure.collision_particles.active_indices),faces(face_structure.Face_Mesh_Object()->mesh.elements),
    X_other(X_other),X_self_collision_free(repulsions.geometry.X_self_collision_free),repulsion_thickness(repulsions.repulsion_thickness),
    thickness_multiplier(repulsions.repulsion_thickness_detection_multiplier*repulsions.hierarchy_repulsion_thickness_multiplier),pruned(pruned),
    perform_attractions(repulsions.perform_attractions),intersecting_point_face_pairs(repulsions.geometry.intersecting_point_face_pairs),
    omit_point_face_repulsion_pairs(repulsions.omit_point_face_repulsion_pairs)
{}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> TRIANGLE_REPULSIONS_POINT_FACE_VISITOR<TV>::
~TRIANGLE_REPULSIONS_POINT_FACE_VISITOR()
{}
//#####################################################################
// Function Store
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS_POINT_FACE_VISITOR<TV>::
Store(const int point_index,const int face_index)
{
    const VECTOR<int,d>& face_nodes=faces(face_index);int p=particle_active_indices(point_index);
    if(face_nodes.Contains(p)) return;
    POINT_FACE_REPULSION_PAIR<TV> pair;pair.nodes=face_nodes.Insert(p,1);
    T total_repulsion_thickness=thickness_multiplier*pair.Total_Repulsion_Thickness(repulsion_thickness);
    T_FACE face(X_other.Subset(face_nodes));
    if(!face.Point_Face_Interaction(X_other(pair.nodes[1]),total_repulsion_thickness,false,pair.distance)) pruned++;
    else{
        if(omit_point_face_repulsion_pairs.Size() && omit_point_face_repulsion_pairs.Contains(pair.nodes)) return;
        if(intersecting_point_face_pairs.Size() && intersecting_point_face_pairs.Contains(pair.nodes)) return;
        T distance;
        if(&X_other!=&X_self_collision_free && perform_attractions){
            T_FACE face(X_self_collision_free.Subset(face_nodes));
            face.Point_Face_Interaction(X_self_collision_free(pair.nodes[1]),total_repulsion_thickness,false,distance);}
        else distance=pair.distance;
        if(distance<0){exchange(pair.nodes[d],pair.nodes[d+1]);pair.distance=-pair.distance;}
        pairs.Append(pair);}
}
//####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template struct TRIANGLE_REPULSIONS_POINT_FACE_VISITOR<VECTOR<T,d> >; \
    template void BOX_HIERARCHY<VECTOR<T,d> >::Intersection_List<TRIANGLE_REPULSIONS_POINT_FACE_VISITOR<VECTOR<T,d> >,ZERO>(BOX_HIERARCHY<VECTOR<T,d> > const&, \
        TRIANGLE_REPULSIONS_POINT_FACE_VISITOR<VECTOR<T,d> >&,ZERO) const; \
    template void BOX_HIERARCHY<VECTOR<T,d> >::Intersection_List<BOX_VISITOR_MPI<TRIANGLE_REPULSIONS_POINT_FACE_VISITOR<VECTOR<T,d> > >,ZERO>( \
        BOX_HIERARCHY<VECTOR<T,d> > const&,BOX_VISITOR_MPI<TRIANGLE_REPULSIONS_POINT_FACE_VISITOR<VECTOR<T,d> > >&,ZERO) const; \
    template void BOX_HIERARCHY<VECTOR<T,d> >::Intersection_List<TRIANGLE_REPULSIONS_POINT_FACE_VISITOR<VECTOR<T,d> >,T>(BOX_HIERARCHY<VECTOR<T,d> > const&, \
        TRIANGLE_REPULSIONS_POINT_FACE_VISITOR<VECTOR<T,d> >&,T) const; \
    template void BOX_HIERARCHY<VECTOR<T,d> >::Intersection_List<BOX_VISITOR_MPI<TRIANGLE_REPULSIONS_POINT_FACE_VISITOR<VECTOR<T,d> > >,T>(BOX_HIERARCHY<VECTOR<T,d> > const&, \
        BOX_VISITOR_MPI<TRIANGLE_REPULSIONS_POINT_FACE_VISITOR<VECTOR<T,d> > >&,T) const;

INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif

//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR
//#####################################################################
#ifndef __TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR__
#define __TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
namespace PhysBAM{
template<class TV> struct EDGE_EDGE_REPULSION_PAIR;
template<class TV> class STRUCTURE_INTERACTION_GEOMETRY;
template<class TV> class TRIANGLE_REPULSIONS;
template<class TV>
struct TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
    typedef typename IF<d==2,int,VECTOR<int,2> >::TYPE T_EDGE;

    ARRAY<EDGE_EDGE_REPULSION_PAIR<TV> >& pairs;
    ARRAY_VIEW<const T_EDGE> edges1,edges2;
    ARRAY_VIEW<const TV> X_other,X_self_collision_free;
    ARRAY_VIEW<const T> repulsion_thickness;
    const T thickness_multiplier;
    int& pruned;
    const bool perform_attractions;
    const HASHTABLE<VECTOR<int,2*d-2> > &intersecting_edge_edge_pairs,&omit_edge_edge_repulsion_pairs;

    TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR(ARRAY<EDGE_EDGE_REPULSION_PAIR<TV> >& pairs,const STRUCTURE_INTERACTION_GEOMETRY<TV>& edge_structure1,
        const STRUCTURE_INTERACTION_GEOMETRY<TV>& edge_structure2,ARRAY_VIEW<const TV> X_other,const TRIANGLE_REPULSIONS<TV>& repulsions,int& pruned);

    ~TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR();

    bool Cull_Self(const int edge) const
    {return false;}

    bool Cull(const int edge1,const int edge2) const
    {return false;}

    void Store_Helper(const int segment1_index,const int segment2_index,const VECTOR<T,3>&) PHYSBAM_ALWAYS_INLINE
    {const VECTOR<int,2> &segment1_nodes=edges1(segment1_index),&segment2_nodes=edges2(segment2_index);
    if(!segment1_nodes.Contains_Any(segment2_nodes))
        Store_Helper_Helper(segment1_index,segment2_index);}

    void Store(const int edge1,const int edge2) PHYSBAM_ALWAYS_INLINE
    {Store_Helper(edge1,edge2,TV());}

//#####################################################################
    void Store_Helper(const int point1_index,const int point2_index,const VECTOR<T,2>&);
    void Store_Helper_Helper(const int segment1_index,const int segment2_index);
//#####################################################################
};
}
#endif

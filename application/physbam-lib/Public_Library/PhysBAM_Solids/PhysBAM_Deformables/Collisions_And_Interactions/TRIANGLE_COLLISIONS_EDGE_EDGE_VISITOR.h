//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR
//#####################################################################
#ifndef __TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR__
#define __TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>

namespace PhysBAM{
template<class TV> class STRUCTURE_INTERACTION_GEOMETRY;
template<class TV> class TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY;
template<class TV>
struct TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
    typedef typename IF<d==2,int,VECTOR<int,2> >::TYPE T_EDGE;

    ARRAY<VECTOR<int,2*d-2> > &pairs_internal,&pairs_external;
    ARRAY_VIEW<const T_EDGE> edges1,edges2;
    ARRAY_VIEW<const TV> X,X_self_collision_free;
    const T collision_thickness;
    ARRAY_VIEW<const bool> edge1_box_modified,edge2_box_modified;
    const HASHTABLE<VECTOR<int,2*d-2> >& intersecting_edge_edge_pairs;
    MPI_SOLIDS<TV>* mpi_solids;

    TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR(ARRAY<VECTOR<int,2*d-2> >& pairs_internal, ARRAY<VECTOR<int,2*d-2> >& pairs_external,
        const STRUCTURE_INTERACTION_GEOMETRY<TV>& edge_structure1,const STRUCTURE_INTERACTION_GEOMETRY<TV>& edge_structure2,
        const TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry,const T collision_thickness,MPI_SOLIDS<TV>* mpi_solids);

    ~TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR();

    bool Cull_Self(const int edge) const
    {return !edge1_box_modified(edge);}

    bool Cull(const int edge1,const int edge2) const
    {return !edge1_box_modified(edge1) && !edge2_box_modified(edge2);}

    void Store_Helper(const int segment1,const int segment2,VECTOR<T,3>) PHYSBAM_ALWAYS_INLINE
    {const VECTOR<int,2> &segment1_nodes=edges1(segment1),&segment2_nodes=edges2(segment2);
    if(!segment1_nodes.Contains_Any(segment2_nodes))
        Store_Helper_Helper(segment1,segment2);} // split the large part off so that we can always inline the first bit

    void Store(const int edge1,const int edge2) PHYSBAM_ALWAYS_INLINE
    {Store_Helper(edge1,edge2,TV());}

//#####################################################################
    void Store_Helper_Helper(const int segment1,const int segment2);
    void Store_Helper(const int point1,const int point2,VECTOR<T,2>);
//#####################################################################
};
}
#endif

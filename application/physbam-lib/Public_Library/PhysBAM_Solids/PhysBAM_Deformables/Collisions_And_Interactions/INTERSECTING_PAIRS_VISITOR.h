//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERSECTING_PAIRS_VISITOR
//#####################################################################
#ifndef __INTERSECTING_PAIRS_VISITOR__
#define __INTERSECTING_PAIRS_VISITOR__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>

namespace PhysBAM{
template<class TV> class STRUCTURE_INTERACTION_GEOMETRY;
template<class TV>
struct INTERSECTING_PAIRS_VISITOR
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};

    HASHTABLE<VECTOR<int,d+1> >& intersecting_point_face_pairs;
    HASHTABLE<VECTOR<int,2*d-2> >& intersecting_edge_edge_pairs;
    STRUCTURE_INTERACTION_GEOMETRY<TV>& face_structure;
    const ARRAY<VECTOR<int,2> >& segments;
    const ARRAY<VECTOR<int,d> >& faces;
    const ARRAY<TV>& X_self_collision_free;
    const T thickness_over_2;
    int& count;

    INTERSECTING_PAIRS_VISITOR(HASHTABLE<VECTOR<int,d+1> >& intersecting_point_face_pairs,HASHTABLE<VECTOR<int,2*d-2> >& intersecting_edge_edge_pairs,
        const STRUCTURE_INTERACTION_GEOMETRY<TV>& segment_structure,STRUCTURE_INTERACTION_GEOMETRY<TV>& face_structure,
        const ARRAY<TV>& X_self_collision_free,const T thickness_over_2,int& count);

    ~INTERSECTING_PAIRS_VISITOR();

    bool Cull_Self(const int segment) const
    {return false;}

    bool Cull(const int segment,const int face) const
    {return false;}

    void Store(const int segment,const int face)
    {if(segments(segment).Contains_Any(faces(face))) return;
    Store_Helper(segment,face,(TV*)NULL);}

//#####################################################################
    void Store_Helper(const int segment1,const int segment2,const VECTOR<T,2>* dimension);
    void Store_Helper(const int segment,const int triangle,const VECTOR<T,3>* dimension);
//#####################################################################
};
}
#endif

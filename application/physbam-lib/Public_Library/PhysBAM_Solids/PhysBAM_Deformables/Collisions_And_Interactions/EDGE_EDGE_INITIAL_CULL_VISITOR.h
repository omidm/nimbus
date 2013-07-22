//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EDGE_EDGE_INITIAL_CULL_VISITOR
//#####################################################################
#ifndef __EDGE_EDGE_INITIAL_CULL_VISITOR__
#define __EDGE_EDGE_INITIAL_CULL_VISITOR__
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>

namespace PhysBAM{
template<class TV> class SEGMENT_ADHESION;
template<class TV>
struct EDGE_EDGE_INITIAL_CULL_VISITOR
{
    typedef typename TV::SCALAR T;
    SEGMENT_ADHESION<TV>& adhesion;
    const ARRAY<int> &mesh1_indices,&mesh2_indices;

    EDGE_EDGE_INITIAL_CULL_VISITOR(SEGMENT_ADHESION<TV>& adhesion_input,const ARRAY<int>& mesh1_indices_input,const ARRAY<int>& mesh2_indices_input)
        :adhesion(adhesion_input),mesh1_indices(mesh1_indices_input),mesh2_indices(mesh2_indices_input)
    {}

    bool Cull_Self(const int edge) const
    {return false;}

    bool Cull(const int edge1,const int edge2) const
    {return false;}

    void Store(const int segment1_local_index,const int segment2_local_index);
};
}
#endif

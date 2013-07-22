//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_REPULSIONS_POINT_FACE_VISITOR
//#####################################################################
#ifndef __TRIANGLE_REPULSIONS_POINT_FACE_VISITOR__
#define __TRIANGLE_REPULSIONS_POINT_FACE_VISITOR__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLES_COLLISIONS_FORWARD.h>

namespace PhysBAM{
template<class TV> struct POINT_FACE_REPULSION_PAIR;
template<class TV> class STRUCTURE_INTERACTION_GEOMETRY;
template<class TV>
struct TRIANGLE_REPULSIONS_POINT_FACE_VISITOR
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
    typedef typename BASIC_SIMPLEX_POLICY<TV,d-1>::SIMPLEX T_FACE;
    typedef typename MESH_POLICY<d-1>::MESH T_MESH;

    ARRAY<POINT_FACE_REPULSION_PAIR<TV> >& pairs;
    ARRAY_VIEW<const int> particle_active_indices;
    ARRAY_VIEW<const VECTOR<int,d> > faces;
    ARRAY_VIEW<const TV> X_other,X_self_collision_free;
    ARRAY_VIEW<const T> repulsion_thickness;
    const T thickness_multiplier;
    int& pruned;
    const bool perform_attractions;
    const HASHTABLE<VECTOR<int,d+1> > &intersecting_point_face_pairs,&omit_point_face_repulsion_pairs;

    TRIANGLE_REPULSIONS_POINT_FACE_VISITOR(ARRAY<POINT_FACE_REPULSION_PAIR<TV> >& pairs,const STRUCTURE_INTERACTION_GEOMETRY<TV>& particle_structure,
        const STRUCTURE_INTERACTION_GEOMETRY<TV>& face_structure,ARRAY_VIEW<const TV> X_other,const TRIANGLE_REPULSIONS<TV>& repulsions,int& pruned);

    ~TRIANGLE_REPULSIONS_POINT_FACE_VISITOR();

    bool Cull_Self(const int) const
    {PHYSBAM_FATAL_ERROR();}

    bool Cull(const int point,const int face) const
    {return false;}

    void Store(const int point_index,const int face_index);
};
}
#endif

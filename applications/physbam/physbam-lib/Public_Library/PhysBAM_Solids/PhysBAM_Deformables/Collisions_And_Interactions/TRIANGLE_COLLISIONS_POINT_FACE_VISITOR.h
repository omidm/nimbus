//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_COLLISIONS_POINT_FACE_VISITOR
//#####################################################################
#ifndef __TRIANGLE_COLLISIONS_POINT_FACE_VISITOR__
#define __TRIANGLE_COLLISIONS_POINT_FACE_VISITOR__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>

namespace PhysBAM{
template<class TV> class STRUCTURE_INTERACTION_GEOMETRY;
template<class TV> class TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY;
template<class TV> class MPI_SOLIDS;
template<class TV>
struct TRIANGLE_COLLISIONS_POINT_FACE_VISITOR
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};

    ARRAY<VECTOR<int,d+1> > &pairs_internal,&pairs_external;
    ARRAY_VIEW<const int> particle_active_indices;
    ARRAY_VIEW<const VECTOR<int,d> > faces;
    ARRAY_VIEW<const TV> X,X_self_collision_free;
    const T collision_thickness;
    ARRAY_VIEW<const bool> point_box_modified,face_box_modified;
    const HASHTABLE<VECTOR<int,d+1> >& intersecting_point_face_pairs;
    MPI_SOLIDS<TV>* mpi_solids;

    TRIANGLE_COLLISIONS_POINT_FACE_VISITOR(ARRAY<VECTOR<int,d+1> >& pairs_internal,ARRAY<VECTOR<int,d+1> >& pairs_external,
        const STRUCTURE_INTERACTION_GEOMETRY<TV>& particle_structure,const STRUCTURE_INTERACTION_GEOMETRY<TV>& face_structure,
        const TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry,const T collision_thickness,MPI_SOLIDS<TV>* mpi_solids);

    ~TRIANGLE_COLLISIONS_POINT_FACE_VISITOR();

    bool Cull_Self(const int) const
    {PHYSBAM_FATAL_ERROR();}

    bool Cull(const int point,const int face) const
    {return !point_box_modified(point) && !face_box_modified(face);}

//#####################################################################
    void Store(const int point_index,const int face_index);
//#####################################################################
};
}
#endif

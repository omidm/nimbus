//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR
//#####################################################################
#ifndef __RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR__
#define __RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>

namespace PhysBAM{
template<class TV> class RIGID_STRUCTURE_INTERACTION_GEOMETRY;
template<class TV> class RIGID_TRIANGLE_COLLISIONS_GEOMETRY;
template<class TV> class MPI_RIGIDS;
template<class TV>
struct RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};

    ARRAY<VECTOR<int,d+1> > &pairs_internal,&pairs_external;
    ARRAY_VIEW<const int> particle_active_indices;
    ARRAY_VIEW<const VECTOR<int,d> > faces;
    const RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& particle_structure,&face_structure;
    const T collision_thickness;
    MPI_RIGIDS<TV>* mpi_solids;

    RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR(ARRAY<VECTOR<int,d+1> >& pairs_internal,ARRAY<VECTOR<int,d+1> >& pairs_external,
        const RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& particle_structure,const RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& face_structure,
        const RIGID_TRIANGLE_COLLISIONS_GEOMETRY<TV>& geometry,const T collision_thickness,MPI_RIGIDS<TV>* mpi_solids);

    ~RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR();

    bool Cull_Self(const int) const
    {PHYSBAM_FATAL_ERROR();}

    bool Cull(const int point,const int face) const
    {return false;}

//#####################################################################
    void Store(const int point_index,const int face_index);
//#####################################################################
};
}
#endif

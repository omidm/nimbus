//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR
//#####################################################################
#ifndef __RIGID_TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR__
#define __RIGID_TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Parallel_Computation/MPI_RIGIDS.h>

namespace PhysBAM{
template<class TV> class RIGID_STRUCTURE_INTERACTION_GEOMETRY;
template<class TV> class RIGID_TRIANGLE_COLLISIONS_GEOMETRY;
template<class TV>
struct RIGID_TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
    typedef typename IF<d==2,int,VECTOR<int,2> >::TYPE T_EDGE;

    ARRAY<VECTOR<int,2*d-2> > &pairs_internal,&pairs_external;
    ARRAY_VIEW<const T_EDGE> edges1,edges2;
    const RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& edge_structure1,&edge_structure2;
    const T collision_thickness;
    MPI_RIGIDS<TV>* mpi_solids;

    RIGID_TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR(ARRAY<VECTOR<int,2*d-2> >& pairs_internal, ARRAY<VECTOR<int,2*d-2> >& pairs_external,
        const RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& edge_structure1,const RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& edge_structure2,
        const RIGID_TRIANGLE_COLLISIONS_GEOMETRY<TV>& geometry,const T collision_thickness,MPI_RIGIDS<TV>* mpi_solids);

    ~RIGID_TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR();

    bool Cull_Self(const int edge) const
    {return false;}

    bool Cull(const int edge1,const int edge2) const
    {return false;}

    void Store_Helper(const int segment1,const int segment2,VECTOR<T,3>) PHYSBAM_ALWAYS_INLINE
    {Store_Helper_Helper(segment1,segment2);} // split the large part off so that we can always inline the first bit

    void Store(const int edge1,const int edge2) PHYSBAM_ALWAYS_INLINE
    {Store_Helper(edge1,edge2,TV());}

//#####################################################################
    void Store_Helper_Helper(const int segment1,const int segment2);
    void Store_Helper(const int point1,const int point2,VECTOR<T,2>);
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_TRIANGLE_COLLISIONS
//#####################################################################
#ifndef __RIGID_TRIANGLE_COLLISIONS__
#define __RIGID_TRIANGLE_COLLISIONS__

#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_TRIANGLE_COLLISIONS_GEOMETRY.h>
namespace PhysBAM{

template<class TV> class MPI_RIGIDS;

template<class TV>
class RIGID_TRIANGLE_COLLISIONS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
    typedef typename BASIC_SIMPLEX_POLICY<TV,d-1>::SIMPLEX T_FACE;
    typedef typename BASIC_SIMPLEX_POLICY<TV,d-1>::SIMPLEX_FACE T_EDGE;
public:
    RIGID_TRIANGLE_COLLISIONS_GEOMETRY<TV>& geometry;
    T collision_thickness;
    bool compute_point_face_collisions,compute_edge_edge_collisions;
    ARRAY<VECTOR<int,d+1> > point_face_pairs_internal, point_face_pairs_external;
    ARRAY<VECTOR<int,2*d-2> > edge_edge_pairs_internal, edge_edge_pairs_external;
    int swap_index;
    // MPI data
    MPI_RIGIDS<TV>* mpi_solids;

    RIGID_TRIANGLE_COLLISIONS(RIGID_TRIANGLE_COLLISIONS_GEOMETRY<TV>& geometry);
    ~RIGID_TRIANGLE_COLLISIONS();

    void Set_Collision_Thickness(const T thickness=1e-6)
    {collision_thickness=thickness;}

    void Compute_Point_Face_Collisions(const bool compute=true)
    {compute_point_face_collisions=compute;}

    void Compute_Edge_Edge_Collisions(const bool compute=true)
    {compute_edge_edge_collisions=compute;}

//#####################################################################
    void Update_Swept_Hierachies();
    void Update_Local_Hierachies();
    void Compute_Pairs(const T detection_thickness,const FRAME<TV>& gframe1,const FRAME<TV>& gframe2,VECTOR<FRAME<TV>,2>* frame1=0,VECTOR<FRAME<TV>,2>* frame2=0);
    void Compute_Pairs(const T detection_thickness,VECTOR<FRAME<TV>,2>* frame1=0,VECTOR<FRAME<TV>,2>* frame2=0);
private:
    void Get_Moving_Faces_Near_Moving_Points(RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1,RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2,ARRAY<VECTOR<int,d+1> >& pairs_internal,ARRAY<VECTOR<int,d+1> >& pairs_external,const T detection_thickness,const FRAME<TV>& gframe1,const FRAME<TV>& gframe2,VECTOR<FRAME<TV>,2>* frame1=0,VECTOR<FRAME<TV>,2>* frame2=0);
    void Get_Moving_Edges_Near_Moving_Edges(RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1,RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2,ARRAY<VECTOR<int,2*d-2> >& pairs_internal,ARRAY<VECTOR<int,2*d-2> >& pairs_external,const T detection_thickness,const FRAME<TV>& gframe1,const FRAME<TV>& gframe2,VECTOR<FRAME<TV>,2>* frame1=0,VECTOR<FRAME<TV>,2>* frame2=0);
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2011, Bo Zhu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __DELAUNAY_TRIANGULATION_H__
#define __DELAUNAY_TRIANGULATION_H__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGULATED_SURFACE_HELPER.h>

namespace PhysBAM{
template <class T>
class DELAUNAY_TRIANGULATION
{
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> INT2;typedef VECTOR<int,3> INT3;
public:
    DELAUNAY_TRIANGULATION(ARRAY<TV> &points_input,ARRAY<INT3> &triangles_input);
    ~DELAUNAY_TRIANGULATION();

    void Run();
    void Set_Remove_Long_Edge_Triangle(T edge_length_threshold_input)
    {edge_length_threshold=edge_length_threshold_input;flag_remove_long_edge_triangle=true;}

    static void Compute_Triangles(ARRAY<TV> &points_input,ARRAY<INT3> &triangles_input)
    {DELAUNAY_TRIANGULATION<T> dt(points_input,triangles_input);dt.Run();}
    static void Compute_Triangles_With_Edge_Threshold(ARRAY<TV> &points_input,ARRAY<INT3> &triangles_input,T edge_threshold)
    {DELAUNAY_TRIANGULATION<T> dt(points_input,triangles_input);dt.Set_Remove_Long_Edge_Triangle(edge_threshold);dt.Run();}
private:
    ARRAY<TV> &points;   ////input
    ARRAY<INT3> &triangles;  ////result
    VECTOR<TV,3> bounding_triangle;
    bool flag_remove_long_edge_triangle;
    T edge_length_threshold;

    int Point_Triangle_Relation(VECTOR<T,2> pos,ARRAY<VECTOR<int,3> > &tris,ARRAY<bool> &flags,/*result*/int &tid,/*result*/VECTOR<int,2> &evid,T tolerance=0);
    void Legalize_Edge(int pr,int pi,int pj,int pr_tri,TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::EDGE_HELPER &edge_helper,ARRAY<INT3> &tris);
    bool Legal_Triangle_Pair(int pr,int pi,int pj,int pk);
};
}

#endif

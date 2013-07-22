//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Neil Molino, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HYPOTHETICAL_CUT_TRIANGLES
//##################################################################### 
#ifndef __HYPOTHETICAL_CUT_TRIANGLES__
#define __HYPOTHETICAL_CUT_TRIANGLES__

#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/HYPOTHETICAL_CUT.h>
namespace PhysBAM{

template<class TV> class EMBEDDED_TRIANGULATED_OBJECT;
template<class TV,int d> class HYPOTHETICAL_NODE;

template<class TV>
class HYPOTHETICAL_CUT_TRIANGLES:public HYPOTHETICAL_CUT<TV,2>
{
    typedef typename TV::SCALAR T;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::HYPERPLANE T_HYPERPLANE;
    typedef HYPOTHETICAL_CUT<TV,2> BASE;
    using BASE::embedded_object;using BASE::cut_quality_metric;using BASE::Position;
public:
    using BASE::hypothetical_nodes;

private:
    bool cuts_ij,cuts_ik,cuts_jk;
    T interpolation_fraction_ij,interpolation_fraction_ik,interpolation_fraction_jk;
    int triangle;
    T_HYPERPLANE hyperplane;
public:
    
    HYPOTHETICAL_CUT_TRIANGLES(EMBEDDED_TRIANGULATED_OBJECT<TV>& embedded_object_input)
        :HYPOTHETICAL_CUT<TV,2>(embedded_object_input)
    {}

    T Quality_Of_Cut() const
    {return cut_quality_metric;}

    bool Valid_Cut() const
    {return hypothetical_nodes.m==2;}

private:
    void Compute_Cuts(TV& xi,TV& xj,TV& xk)
    {cuts_ij=hyperplane.Segment_Intersection(xi,xj,interpolation_fraction_ij);cuts_ik=hyperplane.Segment_Intersection(xi,xk,interpolation_fraction_ik);
    cuts_jk=hyperplane.Segment_Intersection(xj,xk,interpolation_fraction_jk);}
public:

//##################################################################### 
    HYPOTHETICAL_CUT_TRIANGLES& operator=(const HYPOTHETICAL_CUT_TRIANGLES& old_cut);
    bool Initialize_Hypothetical_Cut(const T_HYPERPLANE& hyperplane_input,const int triangle_input);
private:
    bool Segment_Cut_Already_In_Embedded_Triangulated_Object();
    void Initialize_Quality_Metric();
public:
    int Number_Of_Nodes_Shared_With_Existing_Embedded_Curve();
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FRACTURE_TRIANGULATED_OBJECT
//#####################################################################
#ifndef __FRACTURE_TRIANGULATED_OBJECT__
#define __FRACTURE_TRIANGULATED_OBJECT__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDING_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_OBJECT.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/VIRTUAL_NODES.h>
namespace PhysBAM{

template<class TV> class HYPOTHETICAL_CUT_TRIANGLES;

template<class TV>
class FRACTURE_TRIANGULATED_OBJECT:public FRACTURE_OBJECT<TV,2>
{
    typedef typename TV::SCALAR T;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::HYPERPLANE T_HYPERPLANE;
public:
    typedef FRACTURE_OBJECT<TV,2> BASE;
    using BASE::embedded_object;using BASE::fracture_quality_threshold;using BASE::force_edge_connected_fracture;
    using BASE::Get_Phi;using BASE::Fracture_Phi_Index;using BASE::Phi_In_Simplex;using BASE::Positive_Count;using BASE::Initiation_Point;

    FRACTURE_TRIANGULATED_OBJECT(EMBEDDED_TRIANGULATED_OBJECT<TV>& embedded_triangulated_object_input)
        :FRACTURE_OBJECT<TV,2>(embedded_triangulated_object_input)
    {}

//#####################################################################
private:
    void Add_First_Embedded_Segment(const int triangle,const TV& fracture_normal);
    void Add_First_Cut_Based_On_Phi(const int triangle);
    void Add_Second_Embedded_Segment(const int triangle,const TV& fracture_normal);
    int Add_Intersected_Points_To_Embedded_Triangulated_Object(const int triangle,const TV& fracture_normal,HYPOTHETICAL_CUT_TRIANGLES<TV>& hypothetical_cut);
    void Add_Cut(const int triangle,const TV& fracture_normal);
//#####################################################################
};
}
#endif

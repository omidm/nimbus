//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Neil Molino, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FRACTURE_TETRAHEDRALIZED_VOLUME
//#####################################################################
#ifndef __FRACTURE_TETRAHEDRALIZED_VOLUME__
#define __FRACTURE_TETRAHEDRALIZED_VOLUME__

#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_OBJECT.h>
namespace PhysBAM{

template<class T> class HYPOTHETICAL_CUT_TETRAHEDRONS;

template<class T_input>
class FRACTURE_TETRAHEDRALIZED_VOLUME:public FRACTURE_OBJECT<VECTOR<T_input,3>,3>
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef FRACTURE_OBJECT<TV,3> BASE;
    using BASE::embedded_object;using BASE::fracture_quality_threshold;using BASE::force_edge_connected_fracture;
    using BASE::Fracture_Phi_Index;using BASE::Phi_In_Simplex;using BASE::Initiation_Point;using BASE::Get_Phi;

    FRACTURE_TETRAHEDRALIZED_VOLUME(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& embedded_tetrahedralized_volume_input)
        :FRACTURE_OBJECT<TV,3>(embedded_tetrahedralized_volume_input)
    {}

//#####################################################################
    void Add_Cut_Based_On_Phi(const int tetrahedron,const VECTOR<T,4>& tetrahedron_phi) PHYSBAM_OVERRIDE;
    void Add_Second_Cut(const int tetrahedron,const TV& fracture_normal,const VECTOR<T,4>* tetrahedron_phi=0);
private:
    void Add_First_Cut_Based_On_Phi(const int tetrahedron,const VECTOR<T,4>& tetrahedron_phi);
    void Add_Next_Cut_Based_On_Phi(const int tetrahedron,const VECTOR<T,4>& tetrahedron_phi,const bool second_cut);
    void Add_First_Cut(const int tetrahedron,const TV& fracture_normal);
    int Add_Intersected_Points_To_Embedded_Tetrahedralized_Volume(const int tetrahedron,const TV& fracture_normal,HYPOTHETICAL_CUT_TETRAHEDRONS<T>& hypothetical_cut);
    bool Emb_Node_Added_By_This_Tet(const HYPOTHETICAL_CUT_TETRAHEDRONS<T>& hypothetical_cut,const int emb_node) const;
    bool Quad_Diagonal_Is_Correct(const HYPOTHETICAL_CUT_TETRAHEDRONS<T>& hypothetical_cut,const int diagonal_embedded_particle1,const int diagonal_embedded_particle2) const;
    bool Emb_Node_In_Triangle(const EMBEDDED_TETRAHEDRALIZED_VOLUME<T> &embedded_tetrahedralized_volume,const int emb_node,const int emb_tri);
    bool Add_Best_Embedded_Triangle_With_Quad(const TV& fracture_normal,const int tetrahedron,const VECTOR<T,4>* tetrahedron_phi);
    T Interpolation_Fraction_For_Best_Normal(const TV& fracture_normal,const int tetrahedron,const TV& ik,const TV& il,const int i,const int j);
    bool Emb_Node_In_Quad(const EMBEDDED_TETRAHEDRALIZED_VOLUME<T> &embedded_tetrahedralized_volume,const int emb_node,const int emb_tri1,const int emb_tri2);
    void Add_Third_Cut(const int tetrahedron,const TV& fracture_normal,const VECTOR<T,4>* tetrahedron_phi=0);
    void Add_Best_Embedded_Triangle_With_Quad_And_Triangle(const TV& fracture_normal,const int tetrahedron,const VECTOR<T,4>* tetrahedron_phi);
    int Add_Best_Embedded_Triangle_Or_Quad_With_Two_Triangles(const TV& fracture_normal,const int tetrahedron,HYPOTHETICAL_CUT_TETRAHEDRONS<T>& hypothetical_cut);
    void Add_Cut(const int triangle,const TV& fracture_normal) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

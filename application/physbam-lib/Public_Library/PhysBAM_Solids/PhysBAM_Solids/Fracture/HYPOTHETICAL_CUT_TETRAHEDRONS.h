//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Neil Molino, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HYPOTHETICAL_CUT_TETRAHEDRONS
//##################################################################### 
#ifndef __HYPOTHETICAL_CUT_TETRAHEDRONS__
#define __HYPOTHETICAL_CUT_TETRAHEDRONS__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/HYPOTHETICAL_CUT.h>
namespace PhysBAM{

template<class T_input>
class HYPOTHETICAL_CUT_TETRAHEDRONS:public HYPOTHETICAL_CUT<VECTOR<T_input,3>,3>
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef HYPOTHETICAL_CUT<VECTOR<T,3>,3> BASE;
    using BASE::embedded_object;using BASE::hypothetical_nodes;using BASE::cut_quality_metric;using BASE::Position;

    VECTOR<T,3> fracture_normal;
    int cut_index; // 1--7 1=i 2=j 3=k 4=l 5=ij 6=ik 7=il
    int tetrahedron;
    VECTOR<int,2> quad_diagonal_indices;
    
    HYPOTHETICAL_CUT_TETRAHEDRONS(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& embedded_object_input)
        :HYPOTHETICAL_CUT<TV,3>(embedded_object_input)
    {}

//##################################################################### 
    HYPOTHETICAL_CUT_TETRAHEDRONS& operator=(const HYPOTHETICAL_CUT_TETRAHEDRONS& old_cut);
    void Initialize_Hypothetical_Cut(const VECTOR<T,3>& fracture_normal_input,const int cut_index_input,const int tetrahedron_input);
private:
    void Initialize_Triangle_Cut(const int isolated_node);
    void Initialize_Quad_Cut(const int ri,const int rj);
public:
    void Initialize_Initiation_Point_Cut(const PLANE<T>& plane,const int tetrahedron_input);
private:
    bool Valid_Cut(const bool cuts_ij,const bool cuts_ik,const bool cuts_il,const bool cuts_jk,const bool cuts_jl,const bool cuts_kl) const;
    bool Embedded_Edge_Exists(const int emb_node1,const int emb_node2);
    bool Embedded_Edge_Exists(const int ppa1,const int ppa2,const int ppb1,const int ppb2);
    T Interpolation_Fraction_For_Best_Normal(const VECTOR<T,3>& ik,const VECTOR<T,3>& il,const int i,const int j);
    void Interpolation_Fractions_For_Best_Normal(const VECTOR<T,3>& ik,const VECTOR<T,3>& il,const int j,const int k,const int l,T& interpolation_fraction_jk,T& interpolation_fraction_jl);
    void Initialize_Quality_Metric_And_Quad_Diagonal_Indices_For_Initiation_Point();
    void Initialize_Quality_Metric_And_Quad_Diagonal_Indices();
public:
    bool Cut_Already_Exists();
    bool Triangle_Cut_Already_In_Embedded_Tetrahedralized_Volume();
    bool Quad_Cut_Already_In_Embedded_Tetrahedralized_Volume();
    T Quality_Of_Cut(const T extra_quad_penalty_multiplier=(T)1) const;
    bool Edges_Shared_With_Existing_Embedded_Surface();
    bool Would_Orphan_Half_Oct();
//##################################################################### 
};
}
#endif

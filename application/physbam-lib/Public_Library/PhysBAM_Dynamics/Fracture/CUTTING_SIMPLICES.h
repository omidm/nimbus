//#####################################################################
// Copyright 2006-2008, Kevin Der, Geoffrey Irving, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUTTING_SIMPLICES
//##################################################################### 
#ifndef __CUTTING_SIMPLICES__
#define __CUTTING_SIMPLICES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays/PROJECTED_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <PhysBAM_Geometry/Fracture/CUTTING_SIMPLEX.h>

namespace PhysBAM{

template<class T,int d>
class CUTTING_SIMPLICES
{
public:
    ARRAY<CUTTING_SIMPLEX<T,d> > simplices;
    int index_for_last_old_cutting_simplex;

    CUTTING_SIMPLICES()
        :index_for_last_old_cutting_simplex(0)
    {}

    int Add_Simplex(const VECTOR<int,d>& nodes,const typename CUTTING_SIMPLEX<T,d>::SIMPLEX_TYPE type)
    {return simplices.Append(CUTTING_SIMPLEX<T,d>(nodes,type));}

    int Add_Simplex(const VECTOR<int,d>& nodes,const typename CUTTING_SIMPLEX<T,d>::SIMPLEX_TYPE type,const VECTOR<VECTOR<T,d>,d>& simplex_weights,
        const int parent,const int element_owner)
    {return simplices.Append(CUTTING_SIMPLEX<T,d>(nodes,type,simplex_weights,parent,element_owner));}

    int Add_Simplex(const VECTOR<int,d>& nodes,typename CUTTING_SIMPLEX<T,d>::SIMPLEX_TYPE type,T abs_tol,int parent,int element_owner,
        const VECTOR<VECTOR<T,d>,d+1>& element_original_coordinates,const VECTOR<VECTOR<T,d>,d>& simplex_original_coordinates)
    {return simplices.Append(CUTTING_SIMPLEX<T,d>(nodes,type,abs_tol,parent,element_owner,element_original_coordinates,simplex_original_coordinates));}

    bool Simplex_Is_Fake(const int index)
    {return (simplices(index).nodes[d]==0);}
    
    bool Simplex_Is_Parent(const int index)
    {return simplices(index).parent==0;}

    // Note: this actually only converts to global indices if one of the simplices is a local embedded face
    template<class T_ARRAY> T_ARRAY Convert_Simplex_Indices_To_Global_Indices(const T_ARRAY& simplex_indices,bool& converted) const
    {T_ARRAY converted_indices=simplex_indices;INDIRECT_ARRAY<const ARRAY<CUTTING_SIMPLEX<T,d> >,T_ARRAY&> simplices_subset(simplices,simplex_indices);
    if(simplices_subset.template Project<typename CUTTING_SIMPLEX<T,d>::SIMPLEX_TYPE,&CUTTING_SIMPLEX<T,d>::type>().Contains(CUTTING_SIMPLEX<T,d>::LOCAL_EMBEDDING_FACE)){
        converted=true;
        for(int i=1;i<=simplex_indices.m;i++){int t=simplex_indices(i);
            if(simplices(t).type==CUTTING_SIMPLEX<T,d>::LOCAL_CUT_FACE || simplices(t).type==CUTTING_SIMPLEX<T,d>::LOCAL_EMBEDDING_FACE)
                converted_indices(i)=simplices(t).parent;}}
    else converted=false;
    return converted_indices;}

    template<class T_ARRAY> T_ARRAY Convert_Simplex_Indices_To_Global_Indices(const T_ARRAY& simplex_indices) const
    {bool converted;return Convert_Simplex_Indices_To_Global_Indices(simplex_indices,converted);}

    // This finds all the nodes that are shared by /all/ of the cutting simplices
    template<class T_ARRAY>
    void Shared_Nodes_On_Simplices(const T_ARRAY& cutting_simplices,ARRAY<int>& shared_nodes) const
    {for(int i=1;i<=d;i++){int node=simplices(cutting_simplices(1)).nodes[i];if(node==0) continue;
        for(int j=2;j<=cutting_simplices.m;j++) if(!simplices(cutting_simplices(j)).nodes.Contains(node)) goto Next_Node;
        shared_nodes.Append(node);
        Next_Node:;}}

    VECTOR<T,d> Weight_For_Node_In_Simplex(const int simplex_index,const int node)
    {const CUTTING_SIMPLEX<T,d>& simplex=simplices(simplex_index);return simplex.weights(simplex.nodes.Find(node));}

    template<int d2> VECTOR<VECTOR<T,d>,d2> Weight_For_Nodes_In_Simplex(const int simplex_index,const VECTOR<int,d2>& nodes)
    {VECTOR<VECTOR<T,d>,d2> r;for(int i=1;i<=nodes.m;i++) r[i]=Weight_For_Node_In_Simplex(simplex_index,nodes[i]);return r;}

    void Set_Index_For_Last_Old_Cutting_Simplex()
    {index_for_last_old_cutting_simplex=simplices.m;}

    void Print() const
    {std::stringstream ss;for(int i=1;i<=simplices.m;++i){
        const CUTTING_SIMPLEX<T,d>& simplex=simplices(i);
         ss<<"Simplex "<<i<<": element owner="<<simplex.element_owner<<", vertices="<<simplex.nodes
             <<", type ="<<(simplex.type==CUTTING_SIMPLEX<T,d>::GLOBAL_EMBEDDING_FACE?"global_embedding_face":
                 (simplex.type==CUTTING_SIMPLEX<T,d>::LOCAL_EMBEDDING_FACE?"local_embedding_face":(simplex.type==CUTTING_SIMPLEX<T,d>::GLOBAL_CUT_FACE?"global_cut_face":"local_cut_face")))
                  <<", parent = "<<simplex.parent<<", weights = ";
         for(int j=1;j<=d;j++) ss<<simplex.weights(j)<<"; ";
         ss<<std::endl;}
    LOG::filecout(ss.str());}
//#####################################################################    
};
}
#include <PhysBAM_Dynamics/Read_Write/Fracture/READ_WRITE_CUTTING_SIMPLICES.h>
#endif

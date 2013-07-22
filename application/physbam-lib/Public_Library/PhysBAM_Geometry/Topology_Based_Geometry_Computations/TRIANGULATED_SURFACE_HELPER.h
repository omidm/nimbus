//#####################################################################
// Copyright 2011, Bo Zhu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGULATED_SURFACE_HELPER_H__
#define __TRIANGULATED_SURFACE_HELPER_H__
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
namespace PhysBAM{
namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS{
class EDGE_ELEMENT
{
public:
    VECTOR<int,2> vtx;////assume vtx(1)<vtx(2)
    int tri_num;
    VECTOR<int,2> tri;

    EDGE_ELEMENT(int vertex_1_input=-1,int vertex_2_input=-1)
    :vtx(vertex_1_input,vertex_2_input),tri_num(0){}
    EDGE_ELEMENT(const EDGE_ELEMENT &copy)
    {*this=copy;}
    const EDGE_ELEMENT& operator=(const EDGE_ELEMENT &copy)
    {vtx=copy.vtx;tri_num=copy.tri_num;tri=copy.tri;return *this;}
    bool operator==(const EDGE_ELEMENT &copy)
    {return Equal_Edge(copy.vtx(1),copy.vtx(2));}
    void Add_Triangle(int tri_id)
    {assert(tri_num<2);tri(++tri_num)=tri_id;}
    bool Equal_Edge(int vtx_1,int vtx_2)
    {return (vtx(1)==vtx_1&&vtx(2)==vtx_2)||(vtx(1)==vtx_2&&vtx(2)==vtx_1);}
    void Change_Triangle(int tri_id_from,int tri_id_to)
    {for(int i=1;i<=tri_num;i++){if(tri(i)==tri_id_from){tri(i)=tri_id_to;return;}}}
    int Get_Triangle_Besides(int tri_id)
    {if(tri_num==2){if(tri(1)==tri_id) return tri(2);else return tri(1);} return -1;}
};

class EDGE_HELPER
{
    typedef VECTOR<int,3> TV_INT;
public:
    ARRAY<ARRAY<EDGE_ELEMENT> > elements;
    int size;
    ARRAY<TV_INT> &tris;

    EDGE_HELPER(int size_input,ARRAY<TV_INT> &tris_input):size(size_input),tris(tris_input)
    {elements.Resize(size);}

    void Insert_All_Edges_In_Triangle_Array()
    {
        elements.Resize(size,true,false);
        for(int i=1;i<=tris.m;i++){
            TV_INT vtx=tris(i);
            Insert_Edge(vtx(1),vtx(2),i);Insert_Edge(vtx(2),vtx(3),i);Insert_Edge(vtx(3),vtx(1),i);}
    }
    void Insert_Edge(int v1,int v2,int tri_id)
    {
        sort_less(v1,v2);
        EDGE_ELEMENT *e=Find_Edge(v1,v2);
        if(e!=NULL){if(e->tri_num<2)e->Add_Triangle(tri_id);return;}
        elements(v1).Append(EDGE_ELEMENT(v1,v2));
        elements(v1)(elements(v1).m).Add_Triangle(tri_id);
    }
    int Insert_Edge_With_Triangle_Overlap_Check(int v1,int v2,int tri_id)
    {
        sort_less(v1,v2);
        EDGE_ELEMENT *e=Find_Edge(v1,v2);
        if(e!=NULL){
            for(int i=1;i<=e->tri_num;i++){
                if(Same_Triangle(tris(e->tri(i)),tris(tri_id))){
                    return e->tri(i);}}
            assert(e->tri_num<2);
            Add_Triangle_To_Element(tri_id,e);}
        elements(v1).Append(EDGE_ELEMENT(v1,v2));
        Add_Triangle_To_Element(tri_id,&elements(v1)(elements(v1).m));
        return -1;
    }
    EDGE_ELEMENT* Find_Edge(int v1,int v2)
    {
        sort_less(v1,v2);
        for(int i=1;i<=elements(v1).m;i++){
            if(elements(v1)(i).Equal_Edge(v1,v2)){
                return &elements(v1)(i);}}
        return NULL;
    }
    bool Remove_Edge(int v1,int v2)
    {
        sort_less(v1,v2);
        for(int i=1;i<=elements(v1).m;i++){
            if(elements(v1)(i).Equal_Edge(v1,v2)){
                elements(v1).Remove_Index(i);
                return true;}}
        return false;
    }
private:
    void swap(int &a,int &b){int tmp=a;a=b;b=tmp;}
    void sort(VECTOR<int,3> &a){if(a(1)>a(2))swap(a(1),a(2));if(a(2)>a(3))swap(a(2),a(3));if(a(1)>a(2))swap(a(1),a(2));}
    void sort_less(int &a,int &b){if(a>b)swap(a,b);}
    bool Same_Triangle(TV_INT &a,TV_INT &b)
    {VECTOR<int,3> sa=a,sb=b;sort(sa);sort(sb);return sa==sb;}

    void Add_Triangle_To_Element(int tri_id,EDGE_ELEMENT *e)
    {assert(e->tri_num<2);if(e->tri_num<2)e->tri(++(e->tri_num))=tri_id;}
};

inline int Vertex_In_Triangle_Besides(int i1,int i2,int i3,int ei1,int ei2)/*choose i=(i1,i2,i3)-(ei1,ei2)*/
{if(i1!=ei1&&i1!=ei2) return i1;else if(i2!=ei1&&i2!=ei2) return i2;else return i3;}
inline bool Triangle_Contains_Vertex(VECTOR<int,3> &tri,int index)
{return (tri(1)==index||tri(2)==index||tri(3)==index);}
inline bool Triangle_Replace_Vertex(VECTOR<int,3> &tri,int index_from,int index_to)
{for(int i=1;i<=3;i++){if(tri(i)==index_from){tri(i)=index_to;return true;}}return false;}
}
}
#endif

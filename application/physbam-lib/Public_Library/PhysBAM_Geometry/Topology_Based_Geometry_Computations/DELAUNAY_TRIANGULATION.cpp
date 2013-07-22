//#####################################################################
// Copyright 2011, Bo Zhu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "DELAUNAY_TRIANGULATION.h"
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
using namespace PhysBAM;
using namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS;
namespace
{
template<class T>
VECTOR<VECTOR<T,2>,3> Get_Bounding_Triangle(ARRAY<VECTOR<T,2> > points)
{
    typedef VECTOR<T,2> TV;
    RANGE<TV> box=RANGE<TV>::Bounding_Box(points);
    TV center=(box.min_corner+box.max_corner)*0.5;
    T r=(T)3.0*(box.max_corner-box.min_corner).Magnitude()*(T)0.5;
    TV v1=center+TV(-sqrt((T) 3.0)*r,-r),v2=center+TV(sqrt((T) 3.0)*r,-r),v3=center+TV(0,2*r);
    return VECTOR<TV,3>(v1,v2,v3);
}

void Generate_Random_Sequence(ARRAY<int> &index)
{
    RANDOM_NUMBERS<float> random(0); int swap_num=index.m;
    for(int i=1;i<=swap_num;i++){
        VECTOR<int,2> s(random.Get_Uniform_Integer(1,index.m),random.Get_Uniform_Integer(1,index.m));
        int tmp=index(s(1));index(s(1))=index(s(2));index(s(2))=tmp;}
}
}

template<class T> DELAUNAY_TRIANGULATION<T>::
DELAUNAY_TRIANGULATION(ARRAY<TV> &points_input,ARRAY<INT3> &triangles_input)
:points(points_input),triangles(triangles_input),flag_remove_long_edge_triangle(false),edge_length_threshold(0)
{
}

template<class T> DELAUNAY_TRIANGULATION<T>::
~DELAUNAY_TRIANGULATION()
{
}

////1,inside;-1,outside;0,on the edge;
////Brute force method
template<class T> int DELAUNAY_TRIANGULATION<T>::
Point_Triangle_Relation(VECTOR<T,2> pos,ARRAY<VECTOR<int,3> > &tris,ARRAY<bool> &flags,/*result*/int &tid,/*result*/VECTOR<int,2> &evid,T tolerance)
{
    for(int i=1;i<=tris.m;i++){
        if(!flags(i))continue;
        TRIANGLE_2D<T> t(points(tris(i)(1)),points(tris(i)(2)),points(tris(i)(3)));
        if(t.Inside(pos,tolerance)){
            tid=i;
            SEGMENT_2D<T> seg1(points(tris(i)(1)),points(tris(i)(2)));
            SEGMENT_2D<T> seg2(points(tris(i)(2)),points(tris(i)(3)));
            SEGMENT_2D<T> seg3(points(tris(i)(3)),points(tris(i)(1)));
            if(seg1.Distance_From_Point_To_Segment(pos)<=tolerance){
                evid=VECTOR<int,2>(tris(i)(1),tris(i)(2));return 0;}
            if(seg2.Distance_From_Point_To_Segment(pos)<=tolerance){
                evid=VECTOR<int,2>(tris(i)(2),tris(i)(3));return 0;}
            if(seg3.Distance_From_Point_To_Segment(pos)<=tolerance){
                evid=VECTOR<int,2>(tris(i)(3),tris(i)(1));return 0;}
            return 1;}}
    return -1;
}


template<class T> bool DELAUNAY_TRIANGULATION<T>::
Legal_Triangle_Pair(int pr,int pi,int pj,int pk)
{
    VECTOR<T,2> c=TRIANGLE_2D<T>::Circumcenter(points(pi),points(pj),points(pr));
    T r=(c-points(pr)).Magnitude(); T l=(c-points(pk)).Magnitude();
    return r<l;
}

template<class T> void DELAUNAY_TRIANGULATION<T>::
Legalize_Edge(int pr,int pi,int pj,int pr_tid,TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::EDGE_HELPER &edge_helper,ARRAY<INT3> &tris)
{
    int pk_tid,pk;
    EDGE_ELEMENT * e=edge_helper.Find_Edge(pi,pj);
    pk_tid=e->Get_Triangle_Besides(pr_tid);
    if(pk_tid==-1){return;}////edge with single triangle
    pk=Vertex_In_Triangle_Besides(tris(pk_tid)(1),tris(pk_tid)(2),tris(pk_tid)(3),pi,pj);
    if(!Legal_Triangle_Pair(pr,pi,pj,pk)){
        ////swap edge pi-pj to pr-pk
        Triangle_Replace_Vertex(tris(pr_tid),pj,pk);
        Triangle_Replace_Vertex(tris(pk_tid),pi,pr);
        edge_helper.Remove_Edge(pi,pj);
        edge_helper.Insert_Edge(pr,pk,pr_tid);
        edge_helper.Insert_Edge(pr,pk,pk_tid);
        EDGE_ELEMENT *e_rj=edge_helper.Find_Edge(pr,pj);
        if(e_rj)e_rj->Change_Triangle(pr_tid,pk_tid);
        EDGE_ELEMENT *e_ik=edge_helper.Find_Edge(pi,pk);
        if(e_ik)e_ik->Change_Triangle(pk_tid,pr_tid);
        ////legalize edge ik,jk
        Legalize_Edge(pr,pi,pk,pr_tid,edge_helper,tris);
        Legalize_Edge(pr,pj,pk,pk_tid,edge_helper,tris);}
}

template<class T> void DELAUNAY_TRIANGULATION<T>::
Run()
{
    //static int step=1;
    ARRAY<int> index; for(int i=1;i<=points.m;i++){index.Append(i);}
    Generate_Random_Sequence(index);
    bounding_triangle=Get_Bounding_Triangle<T>(points);
    points.Append(bounding_triangle(1));
    points.Append(bounding_triangle(2));
    points.Append(bounding_triangle(3));
    INT3 bounding_triangle_index(points.m-2,points.m-1,points.m);
    ARRAY<INT3> tris;
    tris.Append(INT3(points.m-2,points.m-1,points.m));
    EDGE_HELPER edge_helper(points.m,tris);
    edge_helper.Insert_All_Edges_In_Triangle_Array();
    ARRAY<bool> flag_valid_tri;flag_valid_tri.Resize(3*points.m);flag_valid_tri.Fill(true);
    T tolerance=0;
    for(int i=1;i<=index.m;i++){
        int pr=index(i);
        int tid=-1;INT2 evid;
        int relation=Point_Triangle_Relation(points(pr),tris,flag_valid_tri,tid,evid,tolerance);
        if(relation==1){////inside
            ////1.add three triangles
            tris.Append(INT3(pr,tris(tid)(1),tris(tid)(2)));////id:m-2
            tris.Append(INT3(pr,tris(tid)(2),tris(tid)(3)));////id:m-1
            tris.Append(INT3(pr,tris(tid)(3),tris(tid)(1)));////id:m
            ////2.change the corresponding triangle of edge a,b,c
            EDGE_ELEMENT *e_a,*e_b,*e_c;
            e_a=edge_helper.Find_Edge(tris(tid)(1),tris(tid)(2));
            if(e_a!=NULL)e_a->Change_Triangle(tid,tris.m-2);
            e_b=edge_helper.Find_Edge(tris(tid)(2),tris(tid)(3));
            if(e_b!=NULL)e_b->Change_Triangle(tid,tris.m-1);
            e_c=edge_helper.Find_Edge(tris(tid)(3),tris(tid)(1));
            if(e_c!=NULL)e_c->Change_Triangle(tid,tris.m);
            ////3.add three edges d,e,f
            edge_helper.Insert_Edge(pr,tris(tid)(1),tris.m-2);
            edge_helper.Insert_Edge(pr,tris(tid)(1),tris.m);
            edge_helper.Insert_Edge(pr,tris(tid)(2),tris.m-1);
            edge_helper.Insert_Edge(pr,tris(tid)(2),tris.m-2);
            edge_helper.Insert_Edge(pr,tris(tid)(3),tris.m);
            edge_helper.Insert_Edge(pr,tris(tid)(3),tris.m-1);
            ////4.mark the big triangle as invalid
            flag_valid_tri(tid)=false;
            ////5.flip
            Legalize_Edge(pr,tris(tid)(1),tris(tid)(2),tris.m-2,edge_helper,tris);
            Legalize_Edge(pr,tris(tid)(2),tris(tid)(3),tris.m-1,edge_helper,tris);
            Legalize_Edge(pr,tris(tid)(3),tris(tid)(1),tris.m,edge_helper,tris);}
        else if(relation==0){////on the edge
            int Ta,Tb,Td,Te,pl,pi,pj,pk;
            pi=evid(1);pj=evid(2);
            EDGE_ELEMENT*e=edge_helper.Find_Edge(pi,pj);
            assert(!(e==NULL||e->tri_num<2));
            Tb=tid;Ta=e->Get_Triangle_Besides(Tb);
            pl=Vertex_In_Triangle_Besides(tris(Ta)(1),tris(Ta)(2),tris(Ta)(3),pi,pj);
            pk=Vertex_In_Triangle_Besides(tris(Tb)(1),tris(Tb)(2),tris(Tb)(3),pi,pj);
            ////1.T_lji->T_lri,T_ijk->T_irk
            Triangle_Replace_Vertex(tris(Ta),pj,pr);
            Triangle_Replace_Vertex(tris(Tb),pj,pr);
            ////2.Add triangle T_rlj,T_rjk
            tris.Append(INT3(pr,pl,pj));Td=tris.m;
            tris.Append(INT3(pr,pj,pk));Te=tris.m;
            ////3.e_lj Ta->Td, ejk Tb->Te
            EDGE_ELEMENT *e_lj=edge_helper.Find_Edge(pl,pj);
            if(e_lj)e_lj->Change_Triangle(Ta,Td);
            EDGE_ELEMENT *e_jk=edge_helper.Find_Edge(pj,pk);
            if(e_jk)e_jk->Change_Triangle(Tb,Te);
            ////4.delete e_ij
            edge_helper.Remove_Edge(pi,pj);
            ////5.add e_ri,e_rl,e_rj,e_rk
            edge_helper.Insert_Edge(pr,pi,Ta);edge_helper.Insert_Edge(pr,pi,Tb);
            edge_helper.Insert_Edge(pr,pl,Ta);edge_helper.Insert_Edge(pr,pl,Td);
            edge_helper.Insert_Edge(pr,pj,Td);edge_helper.Insert_Edge(pr,pj,Te);
            edge_helper.Insert_Edge(pr,pk,Te);edge_helper.Insert_Edge(pr,pk,Tb);
            ////6.flip
            Legalize_Edge(pr,pi,pl,Ta,edge_helper,tris);
            Legalize_Edge(pr,pl,pj,Td,edge_helper,tris);
            Legalize_Edge(pr,pj,pk,Te,edge_helper,tris);
            Legalize_Edge(pr,pk,pi,Tb,edge_helper,tris);}
        else if(relation==-1){/*outside*/PHYSBAM_FATAL_ERROR("Point is outside Triangles.");}}
    //step++;
    ////remove the triangles relating to the bounding edges
    for(int i=1;i<=tris.m;i++){
        if(!flag_valid_tri(i))continue;
        if(Triangle_Contains_Vertex(tris(i),bounding_triangle_index(1))||
           Triangle_Contains_Vertex(tris(i),bounding_triangle_index(2))||
           Triangle_Contains_Vertex(tris(i),bounding_triangle_index(3))){
            flag_valid_tri(i)=false;}
        if(flag_remove_long_edge_triangle){
            if((points(tris(i)(1))-points(tris(i)(2))).Magnitude()>edge_length_threshold||
               (points(tris(i)(2))-points(tris(i)(3))).Magnitude()>edge_length_threshold||
               (points(tris(i)(3))-points(tris(i)(1))).Magnitude()>edge_length_threshold){
                flag_valid_tri(i)=false;}}}
    ////remove points of the bounding triangle
    for(int i=1;i<=3;i++)points.Remove_End();
    ////collecting triangles
    triangles.Remove_All();
    for(int i=1;i<=tris.m;i++){if(flag_valid_tri(i))triangles.Append(tris(i));}
}

template class DELAUNAY_TRIANGULATION<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DELAUNAY_TRIANGULATION<double>;
#endif

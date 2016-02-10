//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGULATED_SURFACE_REFRESH__
#define __TRIANGULATED_SURFACE_REFRESH__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/QUEUE.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGULATED_SURFACE_HELPER.h>
namespace PhysBAM{
template<class T> class TRIANGULATED_SURFACE;
template<class T> class TRIANGLE_3D;
template<class T> class TRIANGULATED_AREA;
template<class T> class TRIANGLE_2D;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Function Initialize_Hierarchy
//#####################################################################
template<class T>
void Initialize_Hierarchy(TRIANGULATED_SURFACE<T>& ts,const bool update_boxes,const int triangles_per_group) // creates and updates the boxes as well
{
    delete ts.hierarchy;
    if(ts.triangle_list) ts.hierarchy=new TRIANGLE_HIERARCHY<T>(ts.mesh,ts.particles,*ts.triangle_list,update_boxes,triangles_per_group);
    else ts.hierarchy=new TRIANGLE_HIERARCHY<T>(ts.mesh,ts.particles,update_boxes,triangles_per_group);
}
//#####################################################################
// Function Initialize_Particle_Hierarchy
//#####################################################################
template<class T>
void Initialize_Particle_Hierarchy(TRIANGULATED_SURFACE<T>& ts,const INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<T,3> > >& array_input,const bool update_boxes,const int particles_per_group) // creates and updates the boxes as well
{
    typedef VECTOR<T,3> TV;
    delete ts.particle_hierarchy;
    ts.particle_hierarchy=new PARTICLE_HIERARCHY<TV,INDIRECT_ARRAY<ARRAY_VIEW<TV> > >(array_input,update_boxes,particles_per_group);
}
//#####################################################################
// Function Update_Triangle_List
//#####################################################################
template<class T>
void Update_Triangle_List(TRIANGULATED_SURFACE<T>& ts,ARRAY_VIEW<const VECTOR<T,3> > X)
{
    if(!ts.triangle_list) ts.triangle_list=new ARRAY<TRIANGLE_3D<T> >;
    ts.triangle_list->Resize(ts.mesh.elements.m);
    for(int t=1;t<=ts.mesh.elements.m;t++)
        (*ts.triangle_list)(t)=TRIANGLE_3D<T>(X.Subset(ts.mesh.elements(t)));
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class T>
void Rescale(TRIANGULATED_SURFACE<T>& ts,const T scaling_x,const T scaling_y,const T scaling_z)
{
    typedef VECTOR<T,3> TV;
    if(scaling_x*scaling_y*scaling_z<=0) PHYSBAM_FATAL_ERROR();
    for(int k=1;k<=ts.particles.array_collection->Size();k++) ts.particles.X(k)*=TV(scaling_x,scaling_y,scaling_z);
    if(ts.triangle_list) ts.Update_Triangle_List();if(ts.hierarchy) ts.hierarchy->Update_Boxes();if(ts.bounding_box) ts.Update_Bounding_Box();
}
//#####################################################################
// Function Initialize_Segment_Lengths
//#####################################################################
template<class T>
void Initialize_Segment_Lengths(TRIANGULATED_SURFACE<T>& ts)
{
    bool segment_mesh_defined=ts.mesh.segment_mesh!=0;if(!segment_mesh_defined) ts.mesh.Initialize_Segment_Mesh();
    delete ts.segment_lengths;ts.segment_lengths=new ARRAY<T>(ts.mesh.segment_mesh->elements.m);
    for(int t=1;t<=ts.mesh.segment_mesh->elements.m;t++)
        (*ts.segment_lengths)(t)=(ts.particles.X(ts.mesh.segment_mesh->elements(t)(1))-ts.particles.X(ts.mesh.segment_mesh->elements(t)(2))).Magnitude();
    if(!segment_mesh_defined){delete ts.mesh.segment_mesh;ts.mesh.segment_mesh=0;}
}
//#####################################################################
// Function Initialize_Torus_Mesh_And_Particles
//#####################################################################
template<class T>
void Initialize_Torus_Mesh_And_Particles(TRIANGULATED_SURFACE<T>& ts,const int m,const int n,const T major_radius,const T minor_radius)
{
    typedef VECTOR<T,3> TV;
    T di=T(2*pi)/m,dj=T(2*pi)/n;
    for(int j=1;j<=n;j++){
        T phi=-dj*j,radius=major_radius+minor_radius*cos(phi),z=minor_radius*sin(phi);
        for(int i=1;i<=m;i++){
            int p=ts.particles.array_collection->Add_Element();T theta=di*(i-(T).5*(j&1));
            ts.particles.X(p)=TV(radius*cos(theta),radius*sin(theta),z);}}
    ts.mesh.Initialize_Torus_Mesh(m,n);
}
//#####################################################################
// Function Initialize_Cylinder_Mesh_And_Particles
//#####################################################################
template<class T>
void Initialize_Cylinder_Mesh_And_Particles(TRIANGULATED_SURFACE<T>& ts,const int m,const int n,const T length,const T radius,const bool create_caps)
{
    typedef VECTOR<T,3> TV;
    ts.particles.array_collection->Delete_All_Elements();T dtheta=(T)two_pi/n;T dlength=length/(m-1);
    for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){
        int p=ts.particles.array_collection->Add_Element();T theta=(j-1)*dtheta;
        ts.particles.X(p)=TV(dlength*(i-1),radius*sin(theta),radius*cos(theta));}
    if(create_caps){int p_1=ts.particles.array_collection->Add_Element();int p_2=ts.particles.array_collection->Add_Element();ts.particles.X(p_1)=TV(0,0,0);ts.particles.X(p_2)=TV(length,0,0);}
    ts.mesh.Initialize_Cylinder_Mesh(m,n,create_caps);
}
//#####################################################################
// Function Initialize_Sphere_Mesh_And_Particles
//#####################################################################
template<class T>
void Initialize_Sphere_Mesh_And_Particles(TRIANGULATED_SURFACE<T>& ts,const int m,const int n,const T radius)
{
    typedef VECTOR<T,3> TV;
    ts.particles.array_collection->Delete_All_Elements();T dalpha=(T)pi/m;T dbeta=(T)two_pi/n;
    for(int i=1;i<m;i++) for(int j=1;j<=n;j++){
        int p=ts.particles.array_collection->Add_Element();T alpha=(i-1)*dalpha,beta=(j-1)*dbeta;
        ts.particles.X(p)=TV(radius*cos(alpha),radius*sin(alpha)*cos(beta),radius*sin(alpha)*sin(beta));}
    int p=ts.particles.array_collection->Add_Element();
    ts.particles.X(p)=TV(radius,0,0);
    p=ts.particles.array_collection->Add_Element();
    ts.particles.X(p)=TV(-radius,0,0);
    ts.mesh.Initialize_Sphere_Mesh(m,n);
}
//#####################################################################
//Helper functions for Delaunay_Triangulation_On_Piecewise_Flat_Surface
//#####################################################################
namespace
{

template<class TV>
typename TV::SCALAR Triangle_Edge_Angle(ARRAY<TV> &pos,VECTOR<int,3> &idx,int edge_i1,int edge_i2)
{
    int i3=Vertex_In_Triangle_Besides(idx(1),idx(2),idx(3),edge_i1,edge_i2);
    TV e1=pos(edge_i1)-pos(i3),e2=pos(edge_i2)-pos(i3);
    return TV::Angle_Between(e1,e2);
}
template<class TV>
typename TV::SCALAR Triangle_Edge_Angle_Cot(ARRAY<TV> &pos,VECTOR<int,3> &idx,int edge_i1,int edge_i2)
{
    int i3=Vertex_In_Triangle_Besides(idx(1),idx(2),idx(3),edge_i1,edge_i2);
    TV e1=pos(edge_i1)-pos(i3),e2=pos(edge_i2)-pos(i3);
    return TV::Dot_Product(e1,e2)/TV::Cross_Product(e1,e2).Magnitude();
}

template<class TV>
TV Obtuse_Triangle_Foot_Of_Perpendicular(TV pa/*obtuse vertex*/,TV pb,TV pc)
{
    typedef typename TV::SCALAR T;
    TV b=pc-pa,c=pb-pa,a=pc-pb;T a2=a.Magnitude_Squared();
    T h2=TV::Cross_Product(b,c).Magnitude_Squared()/a2;
    T alpha=sqrt((c.Magnitude_Squared()-h2)/a2);
    return (1-alpha)*pb+alpha*pc;
}

bool Correct_Vertex_Sequence(VECTOR<int,3> v,int p,int ev1,int ev2)
{
    for(int i=1;i<=3;i++){
        if(v(i)==p){if(ev1==v(i%3+1))return true;else return false;}}
    return false;
}

template<class T>
bool Is_Three_Points_Colinear(VECTOR<T,2> p1,VECTOR<T,2> p2,VECTOR<T,2> p3, T tolerance=0)
{
    VECTOR<T,2> v1=p2-p1,v2=p3-p1;
    T a=VECTOR<T,2>::Angle_Between(v1,v2);
    return (a<=tolerance)||(a>=half_pi-tolerance);
}
////two triangles, ev1,ev2,v3; ev1,v4,ev2,cannot handle the condition of three colinear points
template<class T>
int Illegal_Delaunay_Edge_Number(VECTOR<T,2> ev1,VECTOR<T,2> ev2,VECTOR<T,2> v3,VECTOR<T,2> v4)
{
    VECTOR<T,2> c=TRIANGLE_2D<T>::Circumcenter(ev1,ev2,v3);
    T r=(c-ev1).Magnitude(); T l=(c-v4).Magnitude();
    if(r>l) return 1;else if(r==l) return 2;else return 0;
}
}

//#####################################################################
//Function Delaunay_Triangulation_On_Piecewise_Flat_Surface
//#####################################################################
template<class T>
void Delaunay_Triangulation_On_Piecewise_Flat_Surface(TRIANGULATED_SURFACE<T>& ts,T threshold_short_edge,T threshold_sum_angle=1e-5,T threshold_obtuse_angle=half_pi*1.02)
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
    ARRAY<TV> pos;for(int i=1;i<=ts.particles.X.m;i++){pos.Append(ts.particles.X(i));}
    ARRAY<TV_INT> tris=ts.mesh.elements;
    EDGE_HELPER edge_helper(pos.m,tris);
    edge_helper.Insert_All_Edges_In_Triangle_Array();
    QUEUE<EDGE_ELEMENT> eq(tris.m*3);
    for(int vi=1;vi<=edge_helper.elements.m;vi++){
        for(int i=1;i<=edge_helper.elements(vi).m;i++){
            eq.Enqueue(edge_helper.elements(vi)(i));}}
    while(!eq.Empty()){
        EDGE_ELEMENT eh0=eq.Peek();eq.Dequeue();
        EDGE_ELEMENT *eh=edge_helper.Find_Edge(eh0.vtx(1),eh0.vtx(2));
        if(eh==NULL)continue;
        if(eh->tri_num==2){
        T a=(T)0;VECTOR<int,2> p,tri_id,ev;
        tri_id=eh->tri;ev=eh->vtx;
        for(int ti=1;ti<=2;ti++){
            TV_INT v=tris(tri_id(ti));
            p(ti)=Vertex_In_Triangle_Besides(v(1),v(2),v(3),ev(1),ev(2));
            a+=Triangle_Edge_Angle_Cot<TV>(pos,v,ev(1),ev(2));}
        if(a<-threshold_sum_angle){
            ////switch edges of the two triangle
            Triangle_Replace_Vertex(tris(tri_id(1)),ev(2),p(2));
            Triangle_Replace_Vertex(tris(tri_id(2)),ev(1),p(1));
            edge_helper.Remove_Edge(ev(1),ev(2));
            edge_helper.Insert_Edge(p(1),p(2),tri_id(1));
            edge_helper.Insert_Edge(p(1),p(2),tri_id(2));
            EDGE_ELEMENT *ev1_p2=edge_helper.Find_Edge(ev(1),p(2));
            if(ev1_p2!=NULL)ev1_p2->Change_Triangle(tri_id(2),tri_id(1));
                EDGE_ELEMENT *ev2_p1=edge_helper.Find_Edge(ev(2),p(1));
            if(ev2_p1!=NULL)ev2_p1->Change_Triangle(tri_id(1),tri_id(2));
                ////insert new edge to be proceed in eq
                eq.Enqueue(EDGE_ELEMENT(ev(1),p(2)));
                eq.Enqueue(EDGE_ELEMENT(ev(1),p(1)));
                eq.Enqueue(EDGE_ELEMENT(ev(2),p(2)));
                eq.Enqueue(EDGE_ELEMENT(ev(2),p(1)));
                eq.Enqueue(EDGE_ELEMENT(p(1),p(2)));}}}
    ////handle obtuse triangles
    ARRAY<bool> flags;flags.Resize(tris.m);flags.Fill(true);
    for(int vi=1;vi<=edge_helper.elements.m;vi++){
        for(int i=1;i<=edge_helper.elements(vi).m;i++){
            EDGE_ELEMENT eh=edge_helper.elements(vi)(i);
            if(eh.tri_num==1){
                int tri_id=eh.tri(1);TV_INT v=tris(tri_id);VECTOR<int,2> ev=eh.vtx;
                int p=Vertex_In_Triangle_Besides(v(1),v(2),v(3),ev(1),ev(2));
                if(!Correct_Vertex_Sequence(v,p,ev(1),ev(2))){int tmp=ev(1);ev(1)=ev(2);ev(2)=tmp;}
            T a=Triangle_Edge_Angle(pos,v,ev(1),ev(2));
            if(a>threshold_obtuse_angle){
                TV ob_pos=Obtuse_Triangle_Foot_Of_Perpendicular(pos(p),pos(ev(1)),pos(ev(2)));
                if((ob_pos-pos(ev(1))).Magnitude()<threshold_short_edge||
                (ob_pos-pos(ev(2))).Magnitude()<threshold_short_edge||
                (ob_pos-pos(p)).Magnitude()<threshold_short_edge){continue;}
                pos.Append(ob_pos);int ob=pos.m;flags(tri_id)=false;
                tris.Append(TV_INT(ev(1),ob,p));flags.Append(true);
                tris.Append(TV_INT(ob,ev(2),p));flags.Append(true);}}}}
    ts.mesh.elements.Remove_All();
    for(int i=1;i<=tris.m;i++){
        if(flags(i))ts.mesh.elements.Append(tris(i));}
    ts.mesh.number_nodes=ts.mesh.elements.m;
    ts.particles.array_collection->Resize(pos.m);
    ts.particles.X=pos;
}

template<class T>
void Adjust_To_Delaunay_Triangles(TRIANGULATED_AREA<T> & tri_area)
{
    typedef VECTOR<T,2> TV; typedef VECTOR<int,3> INT3; typedef VECTOR<int,2> INT2;
    EDGE_HELPER edge_helper(tri_area.particles.X.m,tri_area.mesh.elements);
    edge_helper.Insert_All_Edges_In_Triangle_Array();
    QUEUE<EDGE_ELEMENT> eq(tri_area.mesh.elements.m*3);
    for(int vi=1;vi<=edge_helper.elements.m;vi++){
        for(int i=1;i<=edge_helper.elements(vi).m;i++){
            eq.Enqueue(edge_helper.elements(vi)(i));}}
    while(!eq.Empty()){
        EDGE_ELEMENT eh0=eq.Peek();eq.Dequeue();
        EDGE_ELEMENT *eh=edge_helper.Find_Edge(eh0.vtx(1),eh0.vtx(2));
        if(eh!=NULL&&eh->tri_num==2){
            INT2 p,ev,tri_id;ev=eh->vtx;tri_id=eh->tri;
            for(int ti=1;ti<=2;ti++){
                INT3 v=tri_area.mesh.elements(tri_id(ti));
                p(ti)=Vertex_In_Triangle_Besides(v(1),v(2),v(3),ev(1),ev(2));}
            if(Illegal_Delaunay_Edge_Number(tri_area.particles.X(ev(1)),tri_area.particles.X(ev(2)),
                                   tri_area.particles.X(p(1)),tri_area.particles.X(p(2)))==1){
                ////switch edges of the two triangle
                Triangle_Replace_Vertex(tri_area.mesh.elements(tri_id(1)),ev(2),p(2));
                Triangle_Replace_Vertex(tri_area.mesh.elements(tri_id(2)),ev(1),p(1));
                edge_helper.Remove_Edge(ev(1),ev(2));   ////eh does not exist after this
                edge_helper.Insert_Edge(p(1),p(2),tri_id(1));
                edge_helper.Insert_Edge(p(1),p(2),tri_id(2));
                EDGE_ELEMENT *ev1_p2=edge_helper.Find_Edge(ev(1),p(2));
                if(ev1_p2!=NULL)ev1_p2->Change_Triangle(tri_id(2),tri_id(1));
                EDGE_ELEMENT *ev2_p1=edge_helper.Find_Edge(ev(2),p(1));
                if(ev2_p1!=NULL)ev2_p1->Change_Triangle(tri_id(1),tri_id(2));
                ////insert new edge to be proceed in eq
                eq.Enqueue(EDGE_ELEMENT(ev(1),p(2)));
                eq.Enqueue(EDGE_ELEMENT(ev(1),p(1)));
                eq.Enqueue(EDGE_ELEMENT(ev(2),p(2)));
                eq.Enqueue(EDGE_ELEMENT(ev(2),p(1)));}}}
}

}
}
#endif

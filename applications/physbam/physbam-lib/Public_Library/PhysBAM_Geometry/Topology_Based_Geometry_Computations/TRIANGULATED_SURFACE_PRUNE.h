//#####################################################################
// Copyright 2011, Bo Zhu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGULATED_SURFACE_PRUNE__
#define __TRIANGULATED_SURFACE_PRUNE__
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGULATED_SURFACE_HELPER.h>
namespace PhysBAM{
namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS{
template<class T>
void Prune_Short_Triangle_Edge(TRIANGULATED_SURFACE<T>& surface,T threshold)
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
    ARRAY<int> redirect_index,new_index;
    ARRAY<bool> vtx_remove_flag,tri_remove_flag,redirected_flag;
    ARRAY<ARRAY<int> >redirected_points;
    int vtx_num=surface.particles.X.m,tri_num=surface.mesh.elements.m;
    GEOMETRY_PARTICLES<TV>&particles=surface.particles;
    redirect_index.Resize(vtx_num);
    for(int i=1;i<=vtx_num; i++)redirect_index(i)=i;
    new_index.Resize(vtx_num);new_index.Fill(-1);
    vtx_remove_flag.Resize(vtx_num);vtx_remove_flag.Fill(false);
    tri_remove_flag.Resize(tri_num);tri_remove_flag.Fill(false);
    redirected_flag.Resize(vtx_num);redirected_flag.Fill(false);
    redirected_points.Resize(vtx_num);
    for(int i=1;i<=tri_num;i++){
        int v1,v2,v3;int rv1,rv2,rv3;
        ////find the root index,v1-->rv1,v2-->rv2,v3-->rv3
        v1=surface.mesh.elements(i)(1);rv1=redirect_index(v1);
        v2=surface.mesh.elements(i)(2);rv2=redirect_index(v2);
        v3=surface.mesh.elements(i)(3);rv3=redirect_index(v3);
        ////if different vertices point to the same particle,remove the triangle
        if(rv1==rv2||rv2==rv3||rv1==rv3){tri_remove_flag(i)=true;continue;}
    ////remove short edges
    int rvid[6];rvid[0]=rv1;rvid[1]=rv2;rvid[2]=rv2;rvid[3]=rv3;rvid[4]=rv3;rvid[5]=rv1;
    for(int j=0;j<3;j++){
        int va=rvid[j*2],vb=rvid[j*2+1];
        if((particles.X(va)-particles.X(vb)).Magnitude()<threshold){
            TV mid=(particles.X(va)+particles.X(vb))/2;
            particles.X(va)=mid;particles.X(vb)=mid;
        ////if b does not contain redirected points,redirect b to a
        if(!redirected_flag(vb)){
            redirect_index(vb)=va;vtx_remove_flag(vb)=true;
            redirected_flag(va)=true;redirected_points(va).Append(vb);}
        ////if a does not contain redirected points,redirect a to b
        else if(!redirected_flag(va)){
            redirect_index(va)=vb;vtx_remove_flag(va)=true;
            redirected_flag(vb)=true;redirected_points(vb).Append(va);}
        ////if both contain redirected points,redirect a to b(and move all vtx pointing to a to b)
        else{
            redirect_index(va)=vb;vtx_remove_flag(va)=true;
            redirected_flag(vb)=true;redirected_points(vb).Append(va);
            for(int k=1; k<=redirected_points(va).m; k++){
                int pid=redirected_points(va)(k);
                redirect_index(pid)=vb;redirected_points(vb).Append(pid);}
            redirected_flag(va)=false;redirected_points(va).Remove_All();}
        tri_remove_flag(i)=true;break;}}}
    ////set new vertex index to triangles, and remove overlapping triangles
    EDGE_HELPER ot_tester(particles.X.m,surface.mesh.elements);
    for(int i=1;i<=tri_num;i++){
        if(tri_remove_flag(i))continue;
        TV_INT v=surface.mesh.elements(i);
        v(1)=redirect_index(v(1));v(2)=redirect_index(v(2));v(3)=redirect_index(v(3));
        surface.mesh.elements(i)=v;////set back to triangles!
        if(v(1)==v(2)||v(2)==v(3)||v(1)==v(3)){tri_remove_flag(i)=true;continue;}
        int vidx[6]= {1,2,2,3,3,1};
        for(int j=0;j<3;j++){
            int idx1=vidx[j*2],idx2=vidx[j*2+1];
            int v1=v(idx1),v2=v(idx2);
            int ti=ot_tester.Insert_Edge_With_Triangle_Overlap_Check(v1,v2,i);
            if(ti!=-1){
                TV_INT tri_ovlp_1=surface.mesh.elements(i);
                TV_INT tri_ovlp_2=surface.mesh.elements(ti);
                for(int k=1;k<=3;k++){
                    tri_ovlp_1(k)=redirect_index(tri_ovlp_1(k));
                    tri_ovlp_2(k)=redirect_index(tri_ovlp_2(k));}
                tri_remove_flag(i)=true;tri_remove_flag(ti)=true;
                vtx_remove_flag(Vertex_In_Triangle_Besides(v(1),v(2),v(3),v1,v2))=true;break;}}}
    ARRAY<TV> new_X;int c=1;
    for(int i=1;i<=vtx_num;i++){
        if(!vtx_remove_flag(i)){new_X.Append(particles.X(i));new_index(i)=c++;}
        else{new_index(i)=-1;}}
    surface.particles.array_collection->Resize(new_X.m);
    surface.particles.X=new_X;
    ARRAY<TV_INT> new_tri;
    for(int i=1;i<=tri_num;i++){
    if(!tri_remove_flag(i)){
        int v1,v2,v3;int rv1,rv2,rv3;
        ////find the root index,v1-->rv1,v2-->rv2,v3-->rv3
        v1=surface.mesh.elements(i)(1);rv1=redirect_index(v1);
        v2=surface.mesh.elements(i)(2);rv2=redirect_index(v2);
        v3=surface.mesh.elements(i)(3);rv3=redirect_index(v3);
        new_tri.Append(TV_INT(new_index(rv1),new_index(rv2),new_index(rv3)));}}
    surface.mesh.number_nodes=new_X.m;
    surface.mesh.elements=new_tri;
}
}
}
#endif

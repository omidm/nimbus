//#####################################################################
// Copyright 2011, Bo Zhu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SEGMENTED_CURVE_PRUNE__
#define __SEGMENTED_CURVE_PRUNE__
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>

namespace PhysBAM
{
namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
template<class TV>
void Prune_Short_Segments(SEGMENTED_CURVE<TV>& curve,typename TV::SCALAR threshold)
{
    typedef VECTOR<int,2> TV_INT;
    ARRAY<TV> vertices;ARRAY<TV_INT> edges;
    vertices=curve.particles.X;edges=curve.mesh.elements;
    ARRAY<int> new_particle_index,remove_edge_flag,remove_particle_flag,merged_flag;
    new_particle_index.Resize(vertices.m);for(int i=1;i<=vertices.m;i++) new_particle_index(i)=i;
    remove_edge_flag.Resize(edges.m);remove_edge_flag.Fill(0);
    remove_particle_flag.Resize(vertices.m);remove_particle_flag.Fill(0);
    merged_flag.Resize(vertices.m);merged_flag.Fill(0);
    for(int i=1;i<=edges.m;i++){
        int v0=new_particle_index(edges(i)(1)),v1=new_particle_index(edges(i)(2));
        edges(i)(1)=v0;edges(i)(2)=v1;
        if((vertices(v0)-vertices(v1)).Magnitude()<threshold){
            remove_edge_flag(i)=1;
            TV mid_point=(vertices(v0)+vertices(v1))/2; vertices(v0)=mid_point;vertices(v1)=mid_point;
            if(merged_flag(v1)==1){new_particle_index(v0)=v1;remove_particle_flag(v0)=1;}
            else{new_particle_index(v1)=v0;remove_particle_flag(v1)=1;merged_flag(v0)=1;}
            remove_edge_flag(i)=1;}}
    ////construct new vertex list
    ARRAY<TV> new_vertices;
    ARRAY<int> idx2;idx2.Resize(vertices.m);int c=1;
    for(int i=1;i<=vertices.m;i++){
        if(remove_particle_flag(i)==0){new_vertices.Append(vertices(i));idx2(i)=c++;}
        else idx2(i)=-1;}
    curve.particles.array_collection->Resize(new_vertices.m);
    curve.particles.X=new_vertices;
    ////construct new edge list
    ARRAY<TV_INT> new_edges;
    for(int i=1;i<=edges.m;i++){
        if(remove_edge_flag(i)==0){
            new_edges.Append(TV_INT(idx2(new_particle_index(edges(i)(1))),idx2(new_particle_index(edges(i)(2)))));}}
    edges=new_edges;
    curve.mesh.number_nodes=new_vertices.m;
    curve.mesh.elements.Resize(new_edges.m);
    curve.mesh.elements=new_edges;
}
}
}
#endif

//#####################################################################
// Copyright 2011, Bo Zhu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MARCHINGSQUARES
//#####################################################################
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHINGSQUARES.h>

using namespace PhysBAM;

namespace{
template<class T> bool Judge_Edge_Has_Vertex(T phi1,T phi2,T iso_value)
{return (phi1<iso_value&&phi2>=iso_value)||(phi2<iso_value&&phi1>=iso_value);}
};

template<class T> void MARCHINGSQUARES<T>::
Marching_Squares()
{
    TV_INT cell_num=grid.Numbers_Of_Cells(),node_num=grid.Numbers_Of_Nodes();
    ARRAY<int,int> vtx_idx_left;vtx_idx_left.Resize(cell_num.y);vtx_idx_left.Fill(-1);
    ARRAY<int,int> vtx_idx_right;vtx_idx_right.Resize(cell_num.y);vtx_idx_right.Fill(-1);
    ARRAY<int,int> vtx_idx_updown;vtx_idx_updown.Resize(node_num.y);vtx_idx_updown.Fill(-1);
    ////calculate column 1 left
    for(int j=1;j<=cell_num.y;j++){
        T val_up=levelset.phi(1,j+1),val_down=levelset.phi(1,j);
        if(Judge_Edge_Has_Vertex(val_up,val_down,contour_value)){
            T alpha=(val_up-contour_value)/(val_up-val_down);TV pos=(1-alpha)*grid.Node(1,j+1)+alpha*grid.Node(1,j);
            vertices.Append(pos);vtx_idx_right(j)=vertices.m;}}
    ////calculate other columns
    for(int i=1;i<=cell_num.x;i++){
        vtx_idx_left=vtx_idx_right;vtx_idx_right.Fill(-1);vtx_idx_updown.Fill(-1);
        ////calculate up,down vertices
        for(int j=1;j<=node_num.y;j++){
            T val_left=levelset.phi(i,j),val_right=levelset.phi(i+1,j);
            if(Judge_Edge_Has_Vertex(val_left,val_right,contour_value)){
            T alpha=(val_left-contour_value)/(val_left-val_right);TV pos=(1-alpha)*grid.Node(i,j)+alpha*grid.Node(i+1,j);
            vertices.Append(pos);vtx_idx_updown(j)=vertices.m;}}
        ////calculate right vertices
        for(int j=1;j<=cell_num.y;j++){
            T val_up=levelset.phi(i+1,j+1),val_down=levelset.phi(i+1,j);
            if(Judge_Edge_Has_Vertex(val_up,val_down,contour_value)){
                T alpha=(val_up-contour_value)/(val_up-val_down);TV pos=(1-alpha)*grid.Node(i+1,j+1)+alpha*grid.Node(i+1,j);
                vertices.Append(pos);vtx_idx_right(j)=vertices.m;}}
        ////calculate edges
        for(int j=1;j<=cell_num.y;j++){
            T value[4]; unsigned int square_type;
            value[0]=levelset.phi(i,j);value[1]=levelset.phi(i+1,j);
            value[2]=levelset.phi(i+1,j+1);value[3]=levelset.phi(i,j+1);
            square_type=Square_Type(value[0],value[1],value[2],value[3],contour_value);
            switch(square_type){
                case 0:
                case 15:/*do nothing*/
                    break;
                case 1:edges.Append(TV_INT(vtx_idx_updown(j),vtx_idx_left(j)));
                    break;
                case 14:edges.Append(TV_INT(vtx_idx_left(j),vtx_idx_updown(j)));
                    break;
                case 2:edges.Append(TV_INT(vtx_idx_right(j),vtx_idx_updown(j)));
                    break;
                case 13:edges.Append(TV_INT(vtx_idx_updown(j),vtx_idx_right(j)));
                    break;
                case 3:edges.Append(TV_INT(vtx_idx_right(j),vtx_idx_left(j)));
                    break;
                case 12:edges.Append(TV_INT(vtx_idx_left(j),vtx_idx_right(j)));
                    break;
                case 4:edges.Append(TV_INT(vtx_idx_updown(j+1),vtx_idx_right(j)));
                    break;
                case 11:edges.Append(TV_INT(vtx_idx_right(j),vtx_idx_updown(j+1)));
                    break;
                case 5:
                    {T avg_v=(T)0.25*(levelset.phi(i,j)+levelset.phi(i+1,j)+levelset.phi(i+1,j+1)+levelset.phi(i,j+1));
                    if(avg_v*levelset.phi(i,j)>contour_value){
                        edges.Append(TV_INT(vtx_idx_updown(j+1),vtx_idx_right(j)));
                        edges.Append(TV_INT(vtx_idx_updown(j),vtx_idx_left(j)));}
                    else{
                        edges.Append(TV_INT(vtx_idx_updown(j+1),vtx_idx_left(j)));
                        edges.Append(TV_INT(vtx_idx_updown(j),vtx_idx_right(j)));}}
                    break;
                case 10:
                    {T avg_v=(T)0.25*(levelset.phi(i,j)+levelset.phi(i+1,j)+levelset.phi(i+1,j+1)+levelset.phi(i,j+1));
                    if(avg_v*levelset.phi(i,j)>contour_value){
                        edges.Append(TV_INT(vtx_idx_right(j),vtx_idx_updown(j+1)));
                        edges.Append(TV_INT(vtx_idx_left(j),vtx_idx_updown(j)));}
                    else{
                        edges.Append(TV_INT(vtx_idx_left(j),vtx_idx_updown(j+1)));
                        edges.Append(TV_INT(vtx_idx_right(j),vtx_idx_updown(j)));}}
                    break;
                case 6:edges.Append(TV_INT(vtx_idx_updown(j+1),vtx_idx_updown(j)));
                    break;
                case 9:edges.Append(TV_INT(vtx_idx_updown(j),vtx_idx_updown(j+1)));
                    break;
                case 7:edges.Append(TV_INT(vtx_idx_updown(j+1),vtx_idx_left(j)));
                    break;
                case 8:edges.Append(TV_INT(vtx_idx_left(j),vtx_idx_updown(j+1)));
                    break;
                default:PHYSBAM_FATAL_ERROR("Invalid marching square type!");
                    break;}}}
}

template<class T> SEGMENTED_CURVE_2D<T>* MARCHINGSQUARES<T>::
Get_Segmented_Curve()
{
    SEGMENTED_CURVE_2D<T>* curve=SEGMENTED_CURVE_2D<T>::Create();
    Get_Segmented_Curve(*curve);
    return curve;
}

template<class T> void MARCHINGSQUARES<T>::
Get_Segmented_Curve(SEGMENTED_CURVE_2D<T>& curve)
{
    int particle_num=vertices.m;int edge_num=edges.m;
    curve.particles.array_collection->Resize(particle_num);
    for(int i=1;i<=particle_num;i++){
        curve.particles.X(i)=vertices(i);}
    curve.mesh.number_nodes=particle_num;
    curve.mesh.elements.Exact_Resize(edge_num);
    for(int t=1;t<=edge_num;t++){
        curve.mesh.elements(t)=edges(t);}
    curve.Update_Segment_List();
}

template<class T> unsigned int MARCHINGSQUARES<T>::
Square_Type(T v0,T v1,T v2,T v3,T iso_value)
{
    unsigned int type=0;
    if(v0<iso_value) type|=1;if(v1<iso_value) type|=2;if(v2<iso_value) type|=4;if(v3<iso_value) type|=8;
    return type;
}

template class MARCHINGSQUARES<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MARCHINGSQUARES<double>;
#endif

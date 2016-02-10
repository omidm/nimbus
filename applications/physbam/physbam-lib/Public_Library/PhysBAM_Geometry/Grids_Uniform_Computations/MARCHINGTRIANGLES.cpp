//#####################################################################
// Copyright 2011, Bo Zhu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MARCHINGTRIANGLES
//#####################################################################
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHINGTRIANGLES.h>

using namespace PhysBAM;

template<class T> void MARCHINGTRIANGLES<T>::
Initialize_Marching_Triangles()
{
    TV_INT x_counts=TV_INT(grid.counts.x-1,grid.counts.y);
    edge_grid_x.Initialize(x_counts,grid.domain,false);
    edge_grid_array_x.Resize(edge_grid_x.Domain_Indices());edge_grid_array_x.Fill(-1);
    TV_INT y_counts=TV_INT(grid.counts.x,grid.counts.y-1);
    edge_grid_y.Initialize(y_counts,grid.domain,false);
    edge_grid_array_y.Resize(edge_grid_y.Domain_Indices());edge_grid_array_y.Fill(-1);
    TV_INT c_counts=TV_INT(grid.counts.x-1,grid.counts.y-1);
    edge_grid_c.Initialize(c_counts,grid.domain,false);
    edge_grid_array_c.Resize(edge_grid_c.Domain_Indices());edge_grid_array_c.Fill(-1);
}

template<class T> int MARCHINGTRIANGLES<T>::
Add_Particle_To_Edge(TV_INT v0,TV_INT v1)
{
    int particle_index=Get_Particle_Index_On_Edge(v0,v1);
    if(particle_index==-1){
        T alpha=(levelset.phi(v0.x,v0.y)-contour_value)/(levelset.phi(v0.x,v0.y)-levelset.phi(v1.x,v1.y));
        TV pos=(1-alpha)*grid.Node(v0.x,v0.y)+alpha*grid.Node(v1.x,v1.y);
        vertices.Append(pos);Set_Particle_Index_On_Edge(v0,v1,vertices.m);return vertices.m;}
    else{return particle_index;}
}

template<class T> void MARCHINGTRIANGLES<T>::
Marching_Triangles()
{
    TV_INT cell_num=grid.Numbers_Of_Cells();
    for(int i=1;i<=cell_num.x;i++){for(int j=1;j<=cell_num.y;j++){
        Process_Triangle(TV_INT(i,j),TV_INT(i+1,j),TV_INT(i+1,j+1));
        Process_Triangle(TV_INT(i+1,j+1),TV_INT(i,j+1),TV_INT(i,j));}}
}

template<class T> void MARCHINGTRIANGLES<T>::
Process_Triangle(TV_INT v0,TV_INT v1,TV_INT v2)
{
    unsigned int type=Get_Triangle_Type(levelset.phi(v0.x,v0.y),levelset.phi(v1.x,v1.y),levelset.phi(v2.x,v2.y),contour_value);
    switch(type){
        case 0:
        case 7:/*do nothing*/
            break;
        case 1:{int p0=Add_Particle_To_Edge(v0,v1),p1=Add_Particle_To_Edge(v0,v2);edges.Append(TV_INT(p0,p1));}
            break;
        case 6:{int p0=Add_Particle_To_Edge(v0,v1),p1=Add_Particle_To_Edge(v0,v2);edges.Append(TV_INT(p1,p0));}
            break;
        case 2:{int p0=Add_Particle_To_Edge(v0,v1),p1=Add_Particle_To_Edge(v1,v2);edges.Append(TV_INT(p1,p0));}
            break;
        case 5:{int p0=Add_Particle_To_Edge(v0,v1),p1=Add_Particle_To_Edge(v1,v2);edges.Append(TV_INT(p0,p1));}
            break;
        case 3:{int p0=Add_Particle_To_Edge(v1,v2),p1=Add_Particle_To_Edge(v0,v2);edges.Append(TV_INT(p0,p1));}
            break;
        case 4:{int p0=Add_Particle_To_Edge(v1,v2),p1=Add_Particle_To_Edge(v0,v2);edges.Append(TV_INT(p1,p0));}
            break;
        default:PHYSBAM_FATAL_ERROR("invalid triangle type");break;}
}

template<class T> int MARCHINGTRIANGLES<T>::
Get_Particle_Index_On_Edge(TV_INT v0,TV_INT v1)
{
    if(v0.x!=v1.x&&v0.y==v1.y){return edge_grid_array_x(min(v0.x,v1.x),v0.y);}
    else if(v0.y!=v1.y&&v0.x==v1.x){return edge_grid_array_y(v0.x,min(v0.y,v1.y));}
    else if(v0.x!=v1.x&&v0.y!=v1.y){return edge_grid_array_c(min(v0.x,v1.x),min(v0.y,v1.y));}
    else{return -1;}
}

template<class T> void MARCHINGTRIANGLES<T>::
Set_Particle_Index_On_Edge(TV_INT v0,TV_INT v1,int input_value)
{
    if(v0.x!=v1.x&&v0.y==v1.y){edge_grid_array_x(min(v0.x,v1.x),v0.y)=input_value;}
    else if(v0.y!=v1.y&&v0.x==v1.x){edge_grid_array_y(v0.x,min(v0.y,v1.y))=input_value;}
    else if(v0.x!=v1.x&&v0.y!=v1.y){edge_grid_array_c(min(v0.x,v1.x),min(v0.y,v1.y))=input_value;}
}

template<class T> SEGMENTED_CURVE_2D<T>* MARCHINGTRIANGLES<T>::
Get_Segmented_Curve()
{
    SEGMENTED_CURVE_2D<T>* curve=SEGMENTED_CURVE_2D<T>::Create();
    Get_Segmented_Curve(*curve);
    return curve;
}

template<class T> void MARCHINGTRIANGLES<T>::
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

template class MARCHINGTRIANGLES<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MARCHINGTRIANGLES<double>;
#endif

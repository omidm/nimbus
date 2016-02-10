//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Boundaries/BOUNDARY_MAC_GRID_SOLID_WALL_SLIP.h>
using namespace PhysBAM;
template<class T_GRID> BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<T_GRID>::
BOUNDARY_MAC_GRID_SOLID_WALL_SLIP(const TV_SIDES& constant_extrapolation)
    :phi(0)
{
    Set_Constant_Extrapolation(constant_extrapolation);
}
template<class T_GRID> BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<T_GRID>::
~BOUNDARY_MAC_GRID_SOLID_WALL_SLIP()
{
}
//#####################################################################
// Function Fill_Ghost_Cells_Face
//#####################################################################
template<class T_GRID> void BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<T_GRID>::
Fill_Ghost_Cells_Face(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& u,T_FACE_ARRAYS_SCALAR& u_ghost,const T time,const int number_of_ghost_cells)
{
    assert(grid.Is_MAC_Grid());
    T_FACE_ARRAYS_SCALAR::Put(u,u_ghost); // interior
    for(int face_axis=1;face_axis<=T_GRID::dimension;face_axis++){
        T_GRID face_grid=grid.Get_Face_Grid(face_axis);
        T_ARRAYS_BASE& u_ghost_component=u_ghost.Component(face_axis);
        ARRAY<RANGE<TV_INT> > regions;Find_Ghost_Regions(face_grid,regions,number_of_ghost_cells);
        for(int side=1;side<=T_GRID::number_of_faces_per_cell;side++){
            if(Constant_Extrapolation(side)) Fill_Single_Ghost_Region(face_grid,u_ghost_component,side,regions(side));
            else Reflect_Single_Ghost_Region(face_axis,face_grid,u_ghost_component,side,regions(side));}}
}
//#####################################################################
// Function Reflect_Single_Ghost_Region
//#####################################################################
template<class T_GRID> void BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<T_GRID>::
Reflect_Single_Ghost_Region(const int face_axis,const T_GRID& face_grid,T_ARRAYS_BASE& u_ghost_component,const int side,const RANGE<TV_INT>& region)
{
    int axis=(side+1)/2,axis_side=side&1?1:2;
    int boundary=Boundary(side,region),reflection_times_two,flip;
    if(face_axis==axis){reflection_times_two=2*boundary;flip=-1;}
    else{reflection_times_two=2*boundary+(axis_side==1?-1:1);flip=1;}
    for(NODE_ITERATOR iterator(face_grid,region);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
        TV_INT reflected_node=node;reflected_node[axis]=reflection_times_two-node[axis];
        u_ghost_component(node)=flip*u_ghost_component(reflected_node);}
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class T_GRID> void BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<T_GRID>::
Apply_Boundary_Condition_Face(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& u,const T time) 
{
    assert(grid.Is_MAC_Grid());
    for(int side=1;side<=T_GRID::number_of_faces_per_cell;side++)
        if(!Constant_Extrapolation(side)) Zero_Single_Boundary_Side(grid,u,side);
}
//#####################################################################
// Function Zero_Single_Boundary_Side
//#####################################################################
template<class T_GRID> void BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<T_GRID>::
Zero_Single_Boundary_Side(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& u,const int side)
{
    int axis=(side+1)/2,axis_side=side&1?1:2;
    FACE_ITERATOR iterator(grid,0,T_GRID::BOUNDARY_REGION,side);
    if(phi){
        TV_INT interior_cell_offset=axis_side==1?TV_INT():-TV_INT::Axis_Vector(axis);
        for(;iterator.Valid();iterator.Next()){TV_INT face=iterator.Face_Index();
            if((*phi)(face+interior_cell_offset)) u.Component(axis)(face)=0;}}
    else for(;iterator.Valid();iterator.Next())u.Component(axis)(iterator.Face_Index())=0;
}
//#####################################################################
template class BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<GRID<VECTOR<float,1> > >;
template class BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<GRID<VECTOR<float,2> > >;
template class BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<GRID<VECTOR<double,1> > >;
template class BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<GRID<VECTOR<double,2> > >;
template class BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<GRID<VECTOR<double,3> > >;
#endif

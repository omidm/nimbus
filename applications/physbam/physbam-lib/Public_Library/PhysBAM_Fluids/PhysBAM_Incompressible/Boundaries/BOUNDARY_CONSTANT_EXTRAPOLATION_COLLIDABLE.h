//#####################################################################
// Copyright 2003-2007, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE
//#####################################################################
#ifndef __BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE__
#define __BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE__

#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
namespace PhysBAM{

template<class T_GRID,class T2>
class BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE:public BOUNDARY_UNIFORM<T_GRID,T2>
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::VECTOR_T TV;
    typedef typename T_GRID::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;
    typedef typename T_GRID::NODE_ITERATOR T_NODE_ITERATOR;typedef typename T_GRID::RIGID_BODY_LIST T_RIGID_BODY_LIST;
public:
    T_RIGID_BODY_LIST& body_list;

    BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE(T_RIGID_BODY_LIST& body_list_input)
        :body_list(body_list_input)
    {}

//#####################################################################
    void Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_T2& u,T_ARRAYS_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) PHYSBAM_OVERRIDE;
    void Apply_Boundary_Condition(const T_GRID& grid,T_ARRAYS_T2& u,const T time) PHYSBAM_OVERRIDE {} // do nothing
    void Collision_Aware_Extrapolate(const T_GRID& grid,T_ARRAYS_T2& u_ghost,const TV_INT& source_index,const TV_INT& ghost_index,const TV& direction,const T ray_length);
//#####################################################################
};
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE<T_GRID,T2>::
Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_T2& u,T_ARRAYS_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells)
{
    T_ARRAYS_T2::Put(u,u_ghost); // interior
    ARRAY<RANGE<TV_INT> > regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells);
    for(int axis=1;axis<=T_GRID::dimension;axis++)for(int axis_side=1;axis_side<=2;axis_side++){
        int side=2*axis+axis_side-2,outward_sign=axis_side?-1:1;
        TV direction=outward_sign*TV::Axis_Vector(axis);
        T signed_dx=outward_sign*grid.DX()[axis];
        int boundary=Boundary(side,regions(side));
        for(T_NODE_ITERATOR iterator(grid,regions(side));iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
            TV_INT boundary_node=node;boundary_node[axis]=boundary;T ray_length=signed_dx*(node[axis]-boundary);
            Collision_Aware_Extrapolate(grid,u_ghost,boundary_node,node,direction,ray_length);}}
}
//#####################################################################
// Function Collision_Aware_Extrapolate
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE<T_GRID,T2>::
Collision_Aware_Extrapolate(const T_GRID& grid,T_ARRAYS_T2& u_ghost,const TV_INT& source_index,const TV_INT& ghost_index,const TV& direction,const T ray_length)
{
    int body_id;
    RAY<TV> ray=RAY<TV>(grid.X(source_index),direction,true);ray.semi_infinite=false;ray.t_max=ray_length;
    PHYSBAM_NOT_IMPLEMENTED(); // update to use fluid collision body list
    if(body_list.Intersection_With_Any_Simplicial_Object(ray,body_id)) u_ghost(ghost_index)=0; // ZERO MIGHT NOT BE GREAT FOR EVERYTHING!
    else u_ghost(ghost_index)=u_ghost(source_index);
}
//#####################################################################
}
#endif

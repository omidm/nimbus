//#####################################################################
// Copyright 2002-2006, Nipun Kwatra, Frank Losasso, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_UNIFORM
//#####################################################################
#ifndef __BOUNDARY_OPEN_WATER__
#define __BOUNDARY_OPEN_WATER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Boundaries/BOUNDARY_FORWARD.h>
namespace PhysBAM{

template<class T_GRID>
class BOUNDARY_OPEN_WATER:public BOUNDARY_UNIFORM<T_GRID,typename T_GRID::SCALAR>
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
public:
    typedef BOUNDARY_UNIFORM<T_GRID,T> BASE;
    using BASE::Constant_Extrapolation;

    ARRAY<bool> open_boundary;
    T attenuate_inflow;
    BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<T_GRID> boundary_mac_grid_solid_wall_slip;

    BOUNDARY_OPEN_WATER(const T attenuate_inflow_input=T(1),const bool left_open_boundary_input=false,const bool right_open_boundary_input=false,const bool bottom_open_boundary_input=false,const bool top_open_boundary_input=false,
        const bool front_open_boundary_input=false,const bool back_open_boundary_input=false)
        :open_boundary(6,true),boundary_mac_grid_solid_wall_slip()
    {
        attenuate_inflow=attenuate_inflow_input;
        open_boundary(1)=left_open_boundary_input;open_boundary(2)=right_open_boundary_input;open_boundary(3)=bottom_open_boundary_input;open_boundary(4)=top_open_boundary_input;open_boundary(5)=front_open_boundary_input;
        open_boundary(6)=back_open_boundary_input;
    }

    ~BOUNDARY_OPEN_WATER(){}

public:

//#####################################################################
    void Fill_Ghost_Cells_Face(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& u,T_FACE_ARRAYS_SCALAR& u_ghost,const T time,const int number_of_ghost_cells=3) PHYSBAM_OVERRIDE;
//#####################################################################
};
//#####################################################################
// Function Fill_Ghost_Cells_Face
//#####################################################################
template<class T_GRID> void BOUNDARY_OPEN_WATER<T_GRID>::
Fill_Ghost_Cells_Face(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& u,T_FACE_ARRAYS_SCALAR& u_ghost,const T time,const int number_of_ghost_cells)
{
    assert(grid.Is_MAC_Grid());
    T_FACE_ARRAYS_SCALAR::Put(u,u_ghost); // interior
    for(int face_axis=1;face_axis<=T_GRID::dimension;face_axis++){
        T_GRID face_grid=grid.Get_Face_Grid(face_axis);
        T_ARRAYS_BASE& u_ghost_component=u_ghost.Component(face_axis);
        ARRAY<RANGE<TV_INT> > regions;Find_Ghost_Regions(face_grid,regions,number_of_ghost_cells);
        for(int side=1;side<=T_GRID::number_of_faces_per_cell;side++){
            RANGE<TV_INT>& region=regions(side);
            int axis=(side+1)/2,boundary=side&1?region.Maximum_Corner()[axis]+1:region.Minimum_Corner()[axis]-1;
            if(open_boundary(side)) for(NODE_ITERATOR iterator(face_grid,region);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index(),boundary_node=node;boundary_node[axis]=boundary; 
                if(side&1) if(u_ghost_component(boundary_node)>0) u_ghost_component(node)=attenuate_inflow*u_ghost_component(boundary_node);else u_ghost_component(node)=u_ghost_component(boundary_node);
                else if(u_ghost_component(boundary_node)<0) u_ghost_component(node)=attenuate_inflow*u_ghost_component(boundary_node);else u_ghost_component(node)=u_ghost_component(boundary_node);}
            else if(Constant_Extrapolation(side)) Fill_Single_Ghost_Region(face_grid,u_ghost_component,side,regions(side));
            else boundary_mac_grid_solid_wall_slip.Reflect_Single_Ghost_Region(face_axis,face_grid,u_ghost_component,side,regions(side));}}
}
}
#endif

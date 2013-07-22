//#####################################################################
// Copyright 2004-2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_SOLID_WALL_SLIP_OUTFLOW
//#####################################################################
#ifndef __BOUNDARY_SOLID_WALL_SLIP_OUTFLOW__
#define __BOUNDARY_SOLID_WALL_SLIP_OUTFLOW__

#include <PhysBAM_Fluids/PhysBAM_Incompressible/Boundaries/BOUNDARY_MAC_GRID_SOLID_WALL_SLIP.h>
namespace PhysBAM{

template<class T_GRID>
class BOUNDARY_SOLID_WALL_SLIP_OUTFLOW:public BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<T_GRID>
{
    typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_GRID::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
public:
    typedef BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<T_GRID> BASE;
    using BASE::Constant_Extrapolation;using BASE::lower_threshold;using BASE::upper_threshold;using BASE::clamp_below;using BASE::clamp_above;

    BOUNDARY_SOLID_WALL_SLIP_OUTFLOW(const bool left_constant_extrapolation_input=false,const bool right_constant_extrapolation_input=false,
        const bool bottom_constant_extrapolation_input=false,const bool top_constant_extrapolation_input=false,
        const bool front_constant_extrapolation_input=false,const bool back_constant_extrapolation_input=false)
        :BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<T_GRID>(left_constant_extrapolation_input,right_constant_extrapolation_input,bottom_constant_extrapolation_input,top_constant_extrapolation_input,
            front_constant_extrapolation_input,back_constant_extrapolation_input)
    {}

//#####################################################################
    void Fill_Ghost_Cells_Face(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& u,T_FACE_ARRAYS_SCALAR& u_ghost,const T time,const int number_of_ghost_cells=3) PHYSBAM_OVERRIDE;
//#####################################################################
};
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class T_GRID> void BOUNDARY_SOLID_WALL_SLIP_OUTFLOW<T_GRID>::
Fill_Ghost_Cells_Face(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& u,T_FACE_ARRAYS_SCALAR& u_ghost,const T time,const int number_of_ghost_cells)
{
    assert(grid.Is_MAC_Grid());
    T_FACE_ARRAYS_SCALAR::Put(u,u_ghost); // interior
    lower_threshold=upper_threshold=0;
    for(int face_axis=1;face_axis<=T_GRID::dimension;face_axis++){
        T_GRID face_grid=grid.Get_Face_Grid(face_axis);
        T_ARRAYS_SCALAR& u_ghost_component=u_ghost.Component(face_axis);
        ARRAY<RANGE<TV_INT> > regions;Find_Ghost_Regions(face_grid,regions,number_of_ghost_cells);
        for(int side=1;side<=T_GRID::number_of_faces_per_cell;side++){
            if(Constant_Extrapolation(side)){
                clamp_below=side&1;clamp_above=!clamp_below;
                Fill_Single_Ghost_Region(face_grid,u_ghost_component,side,regions(side));}
            else Reflect_Single_Ghost_Region(face_axis,face_grid,u_ghost_component,side,regions(side));}}
}
//#####################################################################
}
#endif

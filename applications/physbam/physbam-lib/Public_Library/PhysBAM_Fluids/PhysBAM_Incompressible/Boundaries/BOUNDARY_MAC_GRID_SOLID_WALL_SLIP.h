//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_MAC_GRID_SOLID_WALL_SLIP
//#####################################################################
#ifndef __BOUNDARY_MAC_GRID_SOLID_WALL_SLIP__
#define __BOUNDARY_MAC_GRID_SOLID_WALL_SLIP__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
namespace PhysBAM{

template<class T_GRID>
class BOUNDARY_MAC_GRID_SOLID_WALL_SLIP:public BOUNDARY_UNIFORM<T_GRID,typename T_GRID::SCALAR>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,T_GRID::dimension> TV_SIDES;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;
public:
    typedef BOUNDARY_UNIFORM<T_GRID,T> BASE;
    using BASE::Set_Constant_Extrapolation;using BASE::Constant_Extrapolation;

    const T_ARRAYS_SCALAR* phi;

    BOUNDARY_MAC_GRID_SOLID_WALL_SLIP(const TV_SIDES& constant_extrapolation=TV_SIDES());
    ~BOUNDARY_MAC_GRID_SOLID_WALL_SLIP();

    void Set_Phi(T_ARRAYS_SCALAR& phi_input)
    {phi=&phi_input;}

//#####################################################################
    void Fill_Ghost_Cells_Face(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& u,T_FACE_ARRAYS_SCALAR& u_ghost,const T time,const int number_of_ghost_cells=3) PHYSBAM_OVERRIDE;
    void Apply_Boundary_Condition_Face(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& u,const T time) PHYSBAM_OVERRIDE;
    void Reflect_Single_Ghost_Region(const int face_axis,const T_GRID& face_grid,T_ARRAYS_BASE& u_ghost_component,const int side,const RANGE<TV_INT>& region);
protected:
    void Zero_Single_Boundary_Side(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& u,const int side);
//#####################################################################
};
}
#endif

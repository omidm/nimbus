//#####################################################################
// Copyright 2004, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_SOLID_WALL_SLIP
//#####################################################################
#ifndef __BOUNDARY_SOLID_WALL_SLIP__
#define __BOUNDARY_SOLID_WALL_SLIP__

#include <PhysBAM_Tools/Boundaries/BOUNDARY.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class TV>
class BOUNDARY_SOLID_WALL_SLIP:public BOUNDARY<TV,TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef BOUNDARY<T,TV> BASE;
    using BASE::Set_Constant_Extrapolation;using BASE::Fill_Left_Ghost_Cells;using BASE::Fill_Right_Ghost_Cells;using BASE::Fill_Bottom_Ghost_Cells;
    using BASE::Fill_Top_Ghost_Cells;using BASE::Fill_Front_Ghost_Cells;using BASE::Fill_Back_Ghost_Cells;
    using BASE::left_constant_extrapolation;using BASE::right_constant_extrapolation;using BASE::bottom_constant_extrapolation;
    using BASE::top_constant_extrapolation;using BASE::front_constant_extrapolation;using BASE::back_constant_extrapolation;

    BOUNDARY_SOLID_WALL_SLIP(const bool left_constant_extrapolation_input=false,const bool right_constant_extrapolation_input=false,
        const bool bottom_constant_extrapolation_input=false,const bool top_constant_extrapolation_input=false,
        const bool front_constant_extrapolation_input=false,const bool back_constant_extrapolation_input=false);
    ~BOUNDARY_SOLID_WALL_SLIP();

//#####################################################################
    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAY<TV,VECTOR<int,2> >& V,ARRAY<TV,VECTOR<int,2> >& V_ghost,const T dt,const T time,const int number_of_ghost_cells=3);
    void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAY<TV,VECTOR<int,2> >& V,const T time);
    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAY<TV,VECTOR<int,3> >& V,ARRAY<TV,VECTOR<int,3> >& V_ghost,const T dt,const T time,const int number_of_ghost_cells=3);
    void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAY<TV,VECTOR<int,3> >& V,const T time);
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    void Apply_Boundary_Condition(const QUADTREE_GRID<T>& grid,ARRAY<TV>& V,const T time);
    void Apply_Boundary_Condition(const OCTREE_GRID<T>& grid,ARRAY<TV>& V,const T time);
#endif
//#####################################################################
};
//#####################################################################
}
#endif

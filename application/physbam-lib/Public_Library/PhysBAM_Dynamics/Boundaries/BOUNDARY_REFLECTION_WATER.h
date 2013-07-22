//#####################################################################
// Copyright 2004-2007, Geoffrey Irving, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_REFLECTION_WATER
//#####################################################################
#ifndef __BOUNDARY_REFLECTION_WATER__
#define __BOUNDARY_REFLECTION_WATER__

#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID,class T2>
class BOUNDARY_REFLECTION_WATER:public BOUNDARY_UNIFORM<T_GRID,T2>
{
    typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::FACE_ARRAYS_SCALAR T_FACE_ARRAYS_SCALAR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
public:
    typedef BOUNDARY_UNIFORM<T_GRID,T2> BASE;
    using BASE::Set_Constant_Extrapolation;using BASE::Constant_Extrapolation;

    bool use_extrapolation_mode;
    T tolerance;
private:
    T_ARRAYS_SCALAR* phi;
    T_FACE_ARRAYS_SCALAR* V;
public:

    BOUNDARY_REFLECTION_WATER(const bool left_constant_extrapolation_input=false,const bool right_constant_extrapolation_input=false,const bool bottom_constant_extrapolation_input=false,
        const bool top_constant_extrapolation_input=false,const bool front_constant_extrapolation_input=false,const bool back_constant_extrapolation_input=false)
        :phi(0),V(0)
    {
        Set_Constant_Extrapolation(left_constant_extrapolation_input,right_constant_extrapolation_input,bottom_constant_extrapolation_input,top_constant_extrapolation_input,
            front_constant_extrapolation_input,back_constant_extrapolation_input);
        Use_Extrapolation_Mode(false);
        Set_Tolerance();
    }

    void Use_Extrapolation_Mode(const bool use=true)
    {use_extrapolation_mode=use;}

    void Set_Phi_And_Velocity_Pointers(T_ARRAYS_SCALAR& phi_input,T_FACE_ARRAYS_SCALAR& V_input)
    {phi=&phi_input;V=&V_input;}

    void Set_Tolerance(const T tolerance_input=(T)9.8/24)  // dt*gravity where dt=1/24 is based on the length of a frame
    {tolerance=tolerance_input;}

//#####################################################################
    void Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_T2& u,T_ARRAYS_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) PHYSBAM_OVERRIDE;
//#####################################################################
};
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_REFLECTION_WATER<T_GRID,T2>::
Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_T2& u,T_ARRAYS_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells)
{
    assert(grid.Is_MAC_Grid() && phi && V);T_ARRAYS_T2::Put(u,u_ghost);
    RANGE<TV_INT> domain_indices=grid.Domain_Indices();
    ARRAY<RANGE<TV_INT> > regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells);
    for(int axis=1;axis<=T_GRID::dimension;axis++)for(int axis_side=1;axis_side<=2;axis_side++){
        int side=2*axis+axis_side-2;
        if(use_extrapolation_mode && Constant_Extrapolation(side)) BOUNDARY_UNIFORM<T_GRID,T>::Fill_Single_Ghost_Region(grid,u_ghost,side,regions(side));
        else{ // either phi=phi_object for a wall, or no wall
            int inward_sign=axis_side==1?1:-1;T dx=grid.DX()[axis],half_dx=(T).5*dx;
            int cell_boundary=Boundary(side,regions(side)),face_boundary=cell_boundary+axis_side-1;
            int reflection_times_two=2*cell_boundary+(axis_side==1?-1:1);
            for(CELL_ITERATOR iterator(grid,regions(side));iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
                TV_INT boundary_cell=cell,boundary_face=cell;boundary_cell[axis]=cell_boundary;boundary_face[axis]=face_boundary;
                if((*phi)(boundary_cell) <= 0 && domain_indices.Lazy_Inside(boundary_cell) && inward_sign*V->Component(axis)(boundary_face) > tolerance){
                    TV_INT reflected_cell=cell;reflected_cell[axis]=reflection_times_two-cell[axis];
                    u_ghost(cell)=u_ghost(reflected_cell);}
                else u_ghost(cell)=u_ghost(boundary_cell);}}}
}  
//#####################################################################
}
#endif

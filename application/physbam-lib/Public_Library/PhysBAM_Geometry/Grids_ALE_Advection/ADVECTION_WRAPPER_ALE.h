//#####################################################################
// Copyright 2011, Linhai Qiu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_WRAPPER_ALE
//#####################################################################
#ifndef __ADVECTION_WRAPPER_ALE__
#define __ADVECTION_WRAPPER_ALE__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GRID.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_NESTED_ADVECTION,class T_NESTED_LOOKUP>
class ADVECTION_WRAPPER_ALE:public ADVECTION<T_GRID,T2>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;typedef typename REBIND<T_BOUNDARY,T2>::TYPE T_BOUNDARY_T2;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
public:

    T_NESTED_ADVECTION& nested_advection;
    RIGID_GRID<T_GRID>& rigid_grid;
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2,T_NESTED_LOOKUP> interpolation;
    int number_of_ghost_cells;

    ADVECTION_WRAPPER_ALE(T_NESTED_ADVECTION& nested_advection_input,RIGID_GRID<T_GRID>& rigid_grid_input,int number_of_ghost_cells_input)
        :nested_advection(nested_advection_input),rigid_grid(rigid_grid_input),number_of_ghost_cells(number_of_ghost_cells_input)
    {}

    void Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
        const T_NESTED_LOOKUP& face_velocities,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max)
    {T_GRID cell_center_grid=T_GRID(grid.numbers_of_cells,RANGE<TV>(grid.domain.min_corner+(T).5*grid.dX,grid.domain.max_corner-(T).5*grid.dX));
    T_ARRAYS_VECTOR V(grid.Domain_Indices());
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell=iterator.Cell_Index();TV current_X=iterator.Location();TV X_new_in_old=rigid_grid.Previous_Object_Space_Point(current_X);
        TV grid_velocity=(X_new_in_old-current_X)/dt;TV velocity_in_old_object_space;
        ///////////////////////////////////////////////////////////////////////////DEBUG NUMBER OF GHOST CELLS
        //if(!cell_center_grid.Ghost_Domain(number_of_ghost_cells).Lazy_Inside(X_new_in_old)){
        //    PHYSBAM_FATAL_ERROR("Need more ghost cells in interpolating velocities!");}
        //////////////////////////////////////////////////////////////////////////////
        velocity_in_old_object_space=interpolation.Clamped_To_Array_Face(grid,face_velocities,X_new_in_old);
        V(cell)=velocity_in_old_object_space-grid_velocity;}
    nested_advection.Update_Advection_Equation_Node(cell_center_grid,Z,Z_ghost,V,boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max);}

    void Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_NESTED_LOOKUP& Z_ghost,
        const T_NESTED_LOOKUP& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
        const T_NESTED_LOOKUP* Z_min_ghost,const T_NESTED_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max)
    {T_GRID face_grid;ARRAY<T_ARRAYS_VECTOR> V_advecting(TV::dimension);ARRAY<T_ARRAYS_SCALAR> V_advected(TV::dimension);ARRAY<T_ARRAYS_SCALAR> Z_ghost_components(TV::dimension);
    ARRAY<TV> min_corner_face_grids(TV::dimension);
    for(int axis=1;axis<=TV::dimension;axis++){Z_ghost_components(axis)=Z_ghost.Raw_Data().Component(axis);min_corner_face_grids(axis)=grid.Get_Face_Grid(axis).Domain().Minimum_Corner();}
    for(int i=1;i<=TV::dimension;i++){
        face_grid=grid.Get_Face_Grid(i);
        TV axis_vector=TV::Axis_Vector(i);axis_vector=rigid_grid.Previous_Object_Space_Vector(axis_vector);
        for(int axis=1;axis<=TV::dimension;axis++){V_advected(axis).Resize(face_grid.Domain_Indices());V_advecting(axis).Resize(face_grid.Domain_Indices());}
        for(NODE_ITERATOR iterator(face_grid);iterator.Valid();iterator.Next()){
            TV_INT node=iterator.Node_Index();TV current_X=iterator.Location();TV X_new_in_old=rigid_grid.Previous_Object_Space_Point(current_X);
            TV grid_velocity=(X_new_in_old-current_X)/dt;TV V_advecting_in_old_object_space;
            V_advecting_in_old_object_space=interpolation.Clamped_To_Array_Face(grid,face_velocities,X_new_in_old);
            for(int axis=1;axis<=TV::dimension;axis++){
                V_advecting(axis)(node)=V_advecting_in_old_object_space-grid_velocity+(min_corner_face_grids(axis)-min_corner_face_grids(i))/dt;}}
        for(int axis=1;axis<=TV::dimension;axis++){ 
            nested_advection.Update_Advection_Equation_Node(face_grid,V_advected(axis),Z_ghost_components(axis),V_advecting(axis),boundary,dt,time,0,0,0,0);}
        for(NODE_ITERATOR iterator(face_grid);iterator.Valid();iterator.Next()){
            TV V_advected_in_old_object_space;for(int axis=1;axis<=TV::dimension;axis++) V_advected_in_old_object_space(axis)=V_advected(axis)(iterator.Node_Index());
            Z.Component(i)(iterator.Node_Index())=rigid_grid.Current_Object_Space_Vector(V_advected_in_old_object_space)(i);}}}
    

//#####################################################################
};
}
#endif

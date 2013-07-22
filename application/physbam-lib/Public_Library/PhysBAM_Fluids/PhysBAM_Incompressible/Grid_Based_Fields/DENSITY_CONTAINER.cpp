//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DENSITY_CONTAINER
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Grids_RLE_Boundaries/BOUNDARY_RLE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Geometry/Advection_Collidable/ADVECTION_WRAPPER_COLLIDABLE_CELL.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_DYADIC.h>
#include <PhysBAM_Geometry/Grids_RLE_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_RLE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Grid_Based_Fields/DENSITY_CONTAINER.h>
using namespace PhysBAM;
template<class T_GRID> DENSITY_CONTAINER<T_GRID>::
DENSITY_CONTAINER(T_GRID& grid_input)
    :GRID_AND_ARRAY_CONTAINER<T_GRID,T>(grid_input),density(array),nested_semi_lagrangian_collidable(0),semi_lagrangian_collidable(0)
{
    Set_Ambient_Density();
    boundary_default.Set_Fixed_Boundary(true,ambient_density);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> DENSITY_CONTAINER<T_GRID>::
~DENSITY_CONTAINER()
{
    delete nested_semi_lagrangian_collidable;delete semi_lagrangian_collidable;
}
//#####################################################################
// Function Use_Semi_Lagrangian_Collidable_Advection
//#####################################################################
template<class T_GRID> void DENSITY_CONTAINER<T_GRID>::
Use_Semi_Lagrangian_Collidable_Advection(const T_GRID_BASED_COLLISION_GEOMETRY& body_list,const T_FACE_ARRAYS_BOOL& face_velocities_valid_mask_input)
{
    assert(!nested_semi_lagrangian_collidable&&!semi_lagrangian_collidable);
    nested_semi_lagrangian_collidable=new T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL(body_list,valid_mask_current,valid_mask_next,ambient_density,false);
    semi_lagrangian_collidable=new ADVECTION_WRAPPER_COLLIDABLE_CELL<T_GRID,T,T_FACE_LOOKUP,T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL,T_FACE_LOOKUP_COLLIDABLE>(*nested_semi_lagrangian_collidable,body_list,
        face_velocities_valid_mask_input);
    Set_Custom_Advection(*semi_lagrangian_collidable);
}
//#####################################################################
// Function Initialize_Array
//#####################################################################
template<class T_GRID> void DENSITY_CONTAINER<T_GRID>::
Initialize_Array(const int ghost_cells,const bool initialize_new_elements,const bool copy_existing_elements)
{
    GRID_AND_ARRAY_CONTAINER<T_GRID,T>::Initialize_Array(ghost_cells,initialize_new_elements,copy_existing_elements);
    if(semi_lagrangian_collidable){valid_mask_current.Resize(grid.Cell_Indices(3),true,true,true);valid_mask_next.Resize(grid.Cell_Indices(3),false);} //TODO: Should these be 3 or ghost_cells?
}
//#####################################################################
// Function Fill_Beta_At_Faces
//#####################################################################
template<class T_GRID> void DENSITY_CONTAINER<T_GRID>::
Fill_Beta_At_Faces(const T dt,const T time,T_FACE_ARRAYS_SCALAR& beta_face) const
{
    T_ARRAYS_SCALAR density_ghost(grid.Cell_Indices(1),false);
    boundary->Fill_Ghost_Cells(grid,density,density_ghost,dt,time,1);
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
        TV_INT face_index=iterator.Face_Index();
        TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
        T density_face=(density_ghost(first_cell_index)+density_ghost(second_cell_index))*(T).5;
        beta_face.Component(axis)(face_index)=(T)1/density_face;}
}
//#####################################################################
// Function Get_Ghost_Density
//#####################################################################
template<class T_GRID> void DENSITY_CONTAINER<T_GRID>::
Get_Ghost_Density(const T dt,const T time,const int number_of_ghost_cells,T_ARRAYS_SCALAR& density_ghost) const
{
    boundary->Fill_Ghost_Cells(grid,density,density_ghost,dt,time,number_of_ghost_cells);
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T_GRID> void DENSITY_CONTAINER<T_GRID>::
Euler_Step(const T dt,const T time,const int number_of_ghost_cells)
{  
    GRID_AND_ARRAY_CONTAINER<T_GRID,T>::Euler_Step(dt,time,number_of_ghost_cells);
    array.Clamp_Below(0); // density needs to be non-negative
}
//#####################################################################
template class DENSITY_CONTAINER<GRID<VECTOR<float,1> > >;
template class DENSITY_CONTAINER<GRID<VECTOR<float,2> > >;
template class DENSITY_CONTAINER<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template class DENSITY_CONTAINER<QUADTREE_GRID<float> >;
template class DENSITY_CONTAINER<OCTREE_GRID<float> >;
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template class DENSITY_CONTAINER<RLE_GRID_2D<float> >;
template class DENSITY_CONTAINER<RLE_GRID_3D<float> >;
#endif
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DENSITY_CONTAINER<GRID<VECTOR<double,1> > >;
template class DENSITY_CONTAINER<GRID<VECTOR<double,2> > >;
template class DENSITY_CONTAINER<GRID<VECTOR<double,3> > >;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template class DENSITY_CONTAINER<QUADTREE_GRID<double> >;
template class DENSITY_CONTAINER<OCTREE_GRID<double> >;
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template class DENSITY_CONTAINER<RLE_GRID_2D<double> >;
template class DENSITY_CONTAINER<RLE_GRID_3D<double> >;
#endif
#endif

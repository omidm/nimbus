//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TEMPERATURE_CONTAINER
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Grids_RLE_Boundaries/BOUNDARY_RLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Geometry/Advection_Collidable/ADVECTION_WRAPPER_COLLIDABLE_CELL.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_DYADIC.h>
#include <PhysBAM_Geometry/Grids_RLE_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_RLE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Grid_Based_Fields/TEMPERATURE_CONTAINER.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> TEMPERATURE_CONTAINER<T_GRID>::
TEMPERATURE_CONTAINER(T_GRID& grid_input)
    :GRID_AND_ARRAY_CONTAINER<T_GRID,T>(grid_input),temperature(array),nested_semi_lagrangian_collidable(0),semi_lagrangian_collidable(0)
{
    Set_Ambient_Temperature();
    boundary_default.Set_Fixed_Boundary(true,ambient_temperature);
    Set_Cooling_Constant();
    Set_Hot_Point();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> TEMPERATURE_CONTAINER<T_GRID>::
~TEMPERATURE_CONTAINER()
{
    delete nested_semi_lagrangian_collidable;delete semi_lagrangian_collidable;
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T_GRID> void TEMPERATURE_CONTAINER<T_GRID>::
Euler_Step(const T dt,const T time,const int number_of_ghost_cells)
{
    GRID_AND_ARRAY_CONTAINER<T_GRID,T>::Euler_Step(dt,time,number_of_ghost_cells);
    array.Clamp_Below(ambient_temperature); // temperature needs to be above ambient temperature
    Apply_Cooling(dt,time);
}
//#####################################################################
// Function Apply_Cooling
//#####################################################################
template<class T_GRID> void TEMPERATURE_CONTAINER<T_GRID>::
Apply_Cooling(const T dt,const T time)
{  
    if(!cooling_constant) return;
    T constant=3*cooling_constant*dt/sqr(sqr(hot_point-ambient_temperature));
    for(CELL_ITERATOR iterator(grid,0);iterator.Valid();iterator.Next()) Apply_Individual_Cooling(temperature(iterator.Cell_Index()),constant);
}
// TODO: the following will go away once we figure out a way to merge all the iterators
template<class T,class T_GRID> static void Apply_Cooling_Helper(TEMPERATURE_CONTAINER<T_GRID>& container,const T dt,const T time)
{
    if(!container.cooling_constant) return;
    T constant=3*container.cooling_constant*dt/sqr(sqr(container.hot_point-container.ambient_temperature));
    for(int i=1;i<=container.grid.number_of_cells;i++)container.Apply_Individual_Cooling(container.temperature(i),constant);
}
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template<> void TEMPERATURE_CONTAINER<OCTREE_GRID<float> >::Apply_Cooling(const float dt,const float time){Apply_Cooling_Helper(*this,dt,time);}
template<> void TEMPERATURE_CONTAINER<QUADTREE_GRID<float> >::Apply_Cooling(const float dt,const float time){Apply_Cooling_Helper(*this,dt,time);}
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template<> void TEMPERATURE_CONTAINER<RLE_GRID_2D<float> >::Apply_Cooling(const float dt,const float time){Apply_Cooling_Helper(*this,dt,time);}
template<> void TEMPERATURE_CONTAINER<RLE_GRID_3D<float> >::Apply_Cooling(const float dt,const float time){Apply_Cooling_Helper(*this,dt,time);}
#endif
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template<> void TEMPERATURE_CONTAINER<OCTREE_GRID<double> >::Apply_Cooling(const double dt,const double time){Apply_Cooling_Helper(*this,dt,time);}
template<> void TEMPERATURE_CONTAINER<QUADTREE_GRID<double> >::Apply_Cooling(const double dt,const double time){Apply_Cooling_Helper(*this,dt,time);}
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template<> void TEMPERATURE_CONTAINER<RLE_GRID_2D<double> >::Apply_Cooling(const double dt,const double time){Apply_Cooling_Helper(*this,dt,time);}
template<> void TEMPERATURE_CONTAINER<RLE_GRID_3D<double> >::Apply_Cooling(const double dt,const double time){Apply_Cooling_Helper(*this,dt,time);}
#endif
#endif
//#####################################################################
// Function Apply_Individual_Cooling
//#####################################################################
template<class T_GRID> void TEMPERATURE_CONTAINER<T_GRID>::
Apply_Individual_Cooling(T& temperature,const T constant)
{
    if(abs(temperature-ambient_temperature) > (T).1) temperature=ambient_temperature+pow(1/(constant+1/cube(temperature-ambient_temperature)),(T)one_third);
}
//#####################################################################
// Function Use_Semi_Lagrangian_Collidable_Advection
//#####################################################################
template<class T_GRID> void TEMPERATURE_CONTAINER<T_GRID>::
Use_Semi_Lagrangian_Collidable_Advection(const T_GRID_BASED_COLLISION_GEOMETRY& body_list,const T_FACE_ARRAYS_BOOL& face_velocities_valid_mask_input)
{
    assert(!nested_semi_lagrangian_collidable&&!semi_lagrangian_collidable);
    nested_semi_lagrangian_collidable=new T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL(body_list,valid_mask_current,valid_mask_next,ambient_temperature,false);
    semi_lagrangian_collidable=new ADVECTION_WRAPPER_COLLIDABLE_CELL<T_GRID,T,T_FACE_LOOKUP,T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL,T_FACE_LOOKUP_COLLIDABLE>(*nested_semi_lagrangian_collidable,body_list,
        face_velocities_valid_mask_input);
    Set_Custom_Advection(*semi_lagrangian_collidable);
}
//#####################################################################
template class TEMPERATURE_CONTAINER<GRID<VECTOR<float,1> > >;
template class TEMPERATURE_CONTAINER<GRID<VECTOR<float,2> > >;
template class TEMPERATURE_CONTAINER<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template class TEMPERATURE_CONTAINER<QUADTREE_GRID<float> >;
template class TEMPERATURE_CONTAINER<OCTREE_GRID<float> >;
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template class TEMPERATURE_CONTAINER<RLE_GRID_2D<float> >;
template class TEMPERATURE_CONTAINER<RLE_GRID_3D<float> >;
#endif
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TEMPERATURE_CONTAINER<GRID<VECTOR<double,1> > >;
template class TEMPERATURE_CONTAINER<GRID<VECTOR<double,2> > >;
template class TEMPERATURE_CONTAINER<GRID<VECTOR<double,3> > >;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template class TEMPERATURE_CONTAINER<QUADTREE_GRID<double> >;
template class TEMPERATURE_CONTAINER<OCTREE_GRID<double> >;
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template class TEMPERATURE_CONTAINER<RLE_GRID_2D<double> >;
template class TEMPERATURE_CONTAINER<RLE_GRID_3D<double> >;
#endif
#endif
}

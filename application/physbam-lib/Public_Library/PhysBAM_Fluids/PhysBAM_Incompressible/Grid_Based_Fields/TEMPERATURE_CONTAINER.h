//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TEMPERATURE_CONTAINER
//#####################################################################
#ifndef __TEMPERATURE_CONTAINER__
#define __TEMPERATURE_CONTAINER__

#include <PhysBAM_Geometry/Advection_Collidable/ADVECTION_COLLIDABLE_FORWARD.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Advection_Collidable/ADVECTION_COLLIDABLE_POLICY_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_DYADIC.h>
#include <PhysBAM_Geometry/Grids_RLE_Advection_Collidable/ADVECTION_COLLIDABLE_POLICY_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_RLE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_COLLIDABLE_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Grid_Based_Fields/GRID_AND_ARRAY_CONTAINER.h>
namespace PhysBAM{

template<class T_GRID> struct GRID_ARRAYS_POLICY;

template<class T_GRID>
class TEMPERATURE_CONTAINER:public GRID_AND_ARRAY_CONTAINER<T_GRID,typename T_GRID::SCALAR>
{
    typedef typename T_GRID::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename REBIND<T_ARRAYS_SCALAR,bool>::TYPE T_ARRAYS_BOOL;
    typedef typename ADVECTION_COLLIDABLE_POLICY<T_GRID>::ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_LOOKUP_COLLIDABLE T_FACE_LOOKUP_COLLIDABLE;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef typename REBIND<typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS,bool>::TYPE T_FACE_ARRAYS_BOOL;
public:
    typedef GRID_AND_ARRAY_CONTAINER<T_GRID,T> BASE;
    using BASE::grid;using BASE::array;using BASE::boundary_default;using BASE::boundary;

    T_ARRAYS_SCALAR& temperature;
    T ambient_temperature;
    T cooling_constant;
    T hot_point;

private:
    T_ARRAYS_BOOL valid_mask_current;
    T_ARRAYS_BOOL valid_mask_next;
public:
    T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL* nested_semi_lagrangian_collidable;
    ADVECTION_WRAPPER_COLLIDABLE_CELL<T_GRID,T,T_FACE_LOOKUP,T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL,T_FACE_LOOKUP_COLLIDABLE>* semi_lagrangian_collidable;

    TEMPERATURE_CONTAINER(T_GRID& grid_input);
    ~TEMPERATURE_CONTAINER();

    void Initialize_Array(const int ghost_cells=0,const bool initialize_new_elements=true,const bool copy_existing_elements=true) PHYSBAM_OVERRIDE
    {GRID_AND_ARRAY_CONTAINER<T_GRID,T>::Initialize_Array(ghost_cells,initialize_new_elements,copy_existing_elements);
    if(semi_lagrangian_collidable){valid_mask_current.Resize(grid.Cell_Indices(3),true,true,true);valid_mask_next.Resize(grid.Cell_Indices(3),false);}}
  
    void Set_Ambient_Temperature(const T temperature_input=298)
    {ambient_temperature=temperature_input;boundary->Set_Fixed_Boundary(true,ambient_temperature);}

    void Set_To_Ambient_Temperature()
    {Set_To_Constant_Value(ambient_temperature);}

    void Set_Cooling_Constant(const T cooling_constant_input=3000)
    {cooling_constant=cooling_constant_input;}

    void Set_Hot_Point(const T hot_point_input=3000)
    {hot_point=hot_point_input;}

//#####################################################################
    void Euler_Step(const T dt,const T time,const int number_of_ghost_cells) PHYSBAM_OVERRIDE;
    void Apply_Cooling(const T dt,const T time);
    void Apply_Individual_Cooling(T& temperature,const T constant);
    void Use_Semi_Lagrangian_Collidable_Advection(const T_GRID_BASED_COLLISION_GEOMETRY& body_list,const T_FACE_ARRAYS_BOOL& face_velocities_valid_mask_input);
//#####################################################################
};
}
#endif

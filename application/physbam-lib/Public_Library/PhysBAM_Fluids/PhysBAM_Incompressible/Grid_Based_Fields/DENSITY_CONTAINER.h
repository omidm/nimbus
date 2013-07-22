//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Eran Guendelman, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DENSITY_CONTAINER
//#####################################################################
#ifndef __DENSITY_CONTAINER__    
#define __DENSITY_CONTAINER__

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
class DENSITY_CONTAINER:public GRID_AND_ARRAY_CONTAINER<T_GRID,typename T_GRID::SCALAR>
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename REBIND<T_ARRAYS_SCALAR,bool>::TYPE T_ARRAYS_BOOL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename ADVECTION_COLLIDABLE_POLICY<T_GRID>::ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL;typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_LOOKUP_COLLIDABLE T_FACE_LOOKUP_COLLIDABLE;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef typename REBIND<typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS,bool>::TYPE T_FACE_ARRAYS_BOOL;
public:
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef GRID_AND_ARRAY_CONTAINER<T_GRID,T> BASE;
    using BASE::array;using BASE::grid;
    using BASE::boundary_default;using BASE::boundary;

    T_ARRAYS_SCALAR& density;
    T ambient_density;

    T_ARRAYS_BOOL valid_mask_current;
    T_ARRAYS_BOOL valid_mask_next;
public:
    T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL* nested_semi_lagrangian_collidable;
private:
    ADVECTION_WRAPPER_COLLIDABLE_CELL<T_GRID,T,T_FACE_LOOKUP,T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL,T_FACE_LOOKUP_COLLIDABLE>* semi_lagrangian_collidable;
public:    

    DENSITY_CONTAINER(T_GRID& grid_input);
    ~DENSITY_CONTAINER();

    void Set_Ambient_Density(const T density_input=0)
    {ambient_density=density_input;}

    void Set_To_Ambient_Density()
    {Set_To_Constant_Value(ambient_density);}

//#####################################################################
    void Euler_Step(const T dt,const T time,const int number_of_ghost_cells) PHYSBAM_OVERRIDE;
    void Initialize_Array(const int ghost_cells=0,const bool initialize_new_elements=true,const bool copy_existing_elements=true) PHYSBAM_OVERRIDE;
    void Use_Semi_Lagrangian_Collidable_Advection(const T_GRID_BASED_COLLISION_GEOMETRY& body_list,const T_FACE_ARRAYS_BOOL& face_velocities_valid_mask_input);
    void Fill_Beta_At_Faces(const T dt,const T time,T_FACE_ARRAYS_SCALAR& beta_face) const;
    void Get_Ghost_Density(const T dt,const T time,const int number_of_ghost_cells,T_ARRAYS_SCALAR& density_ghost) const;
//#####################################################################
};      
}
#endif

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2009, Doug Enright, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso, Neil Molino, Avi Robinson-Mosher, Tamar Shinar, Jerry Talton, Elliot English.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_ADVECTION_DYADIC
//##################################################################### 
#ifndef __LEVELSET_ADVECTION_DYADIC__
#define __LEVELSET_ADVECTION_DYADIC__

#include <PhysBAM_Tools/Grids_Dyadic_Arrays/GRID_ARRAYS_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_POLICY_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Advection_Collidable/ADVECTION_COLLIDABLE_POLICY_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_UNIFORM.h>
#include <PhysBAM_Geometry/Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION.h>

namespace PhysBAM {
    
template<class T_GRID>
class LEVELSET_ADVECTION_DYADIC:public LEVELSET_ADVECTION<T_GRID>
{
    typedef LEVELSET_ADVECTION<T_GRID> BASE;
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename ADVECTION_COLLIDABLE_POLICY<T_GRID>::ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL;
    typedef typename ADVECTION_COLLIDABLE_POLICY<T_GRID>::ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_SLIP_CELL T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_SLIP_CELL;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_LOOKUP_COLLIDABLE T_FACE_LOOKUP_COLLIDABLE;
    typedef typename REBIND<typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS,bool>::TYPE T_FACE_ARRAYS_BOOL;typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY_SCALAR;
    typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_LOOKUP_COLLIDABLE_SLIP T_FACE_LOOKUP_COLLIDABLE_SLIP;
    typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
public:
    using BASE::advection;
    using BASE::levelset;

    LEVELSET_ADVECTION_DYADIC(T_LEVELSET* _levelset):BASE(_levelset){Use_Semi_Lagrangian_Advection();};

    void Euler_Step(const ARRAY<T>& face_velocity,const T dt,const T time=0);
    void Use_Semi_Lagrangian_Advection();
};

}
#endif
#endif

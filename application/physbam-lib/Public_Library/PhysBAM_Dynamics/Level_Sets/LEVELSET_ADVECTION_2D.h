//#####################################################################
// Copyright 2009, Elliot English
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_ADVECTION_2D
//##################################################################### 
#ifndef __LEVELSET_ADVECTION_2D__
#define __LEVELSET_ADVECTION_2D__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_COLLIDABLE_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_UNIFORM.h>

namespace PhysBAM {
    
template<class T>
class LEVELSET_ADVECTION_2D:public LEVELSET_ADVECTION_UNIFORM<GRID<VECTOR<T,2> > >
{
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
    typedef LEVELSET_ADVECTION_UNIFORM<GRID<TV> > BASE;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename ADVECTION_COLLIDABLE_POLICY<GRID<TV> >::ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<TV> >::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef typename INTERPOLATION_POLICY<GRID<TV> >::FACE_LOOKUP T_FACE_LOOKUP;typedef typename INTERPOLATION_COLLIDABLE_POLICY<GRID<TV> >::FACE_LOOKUP_COLLIDABLE T_FACE_LOOKUP_COLLIDABLE;
    typedef typename REBIND<typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS,bool>::TYPE T_FACE_ARRAYS_BOOL;typedef typename BOUNDARY_POLICY<GRID<TV> >::BOUNDARY_SCALAR T_BOUNDARY_SCALAR;
    typedef typename INTERPOLATION_COLLIDABLE_POLICY<GRID<TV> >::FACE_LOOKUP_COLLIDABLE_SLIP T_FACE_LOOKUP_COLLIDABLE_SLIP;
    typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;
public:
    using BASE::levelset;
    using BASE::advection;    
    using BASE::reinitialization_cfl;using BASE::reinitialization_runge_kutta_order;using BASE::reinitialization_spatial_order;

    LEVELSET_ADVECTION_2D(T_LEVELSET* _levelset):BASE(_levelset){};
    LEVELSET_ADVECTION_2D():BASE(0){};

    void Euler_Step(const ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocity,const T dt,const T time,const int number_of_ghost_cells);
    void Reinitialize(const int time_steps=10,const T time=0);
private:
    void Euler_Step_Of_Reinitialization(const ARRAY<T,VECTOR<int,2> >& sign_phi,const T dt,const T time);
};

}
#endif

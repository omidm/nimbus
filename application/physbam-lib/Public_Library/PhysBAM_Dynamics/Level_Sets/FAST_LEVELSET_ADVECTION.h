//#####################################################################
// Copyright 2009, Doug Enright, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso, Neil Molino, Avi Robinson-Mosher, Tamar Shinar, Jerry Talton, Elliot English.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FAST_LEVELSET_ADVECTION
//##################################################################### 
#ifndef __FAST_LEVELSET_ADVECTION__
#define __FAST_LEVELSET_ADVECTION__

#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_1D.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_2D.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_3D.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_POLICY.h>

namespace PhysBAM {

template<class T_GRID>
class FAST_LEVELSET_ADVECTION:public LEVELSET_ADVECTION_POLICY<T_GRID>::LEVELSET_ADVECTION
{
    typedef typename LEVELSET_ADVECTION_POLICY<T_GRID>::LEVELSET_ADVECTION BASE;
    typedef typename LEVELSET_POLICY<T_GRID>::FAST_LEVELSET_T T_FAST_LEVELSET;
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    
public:
    using BASE::levelset;
    using BASE::advection;
    using BASE::reinitialization_cfl;using BASE::reinitialization_runge_kutta_order;using BASE::reinitialization_spatial_order;

    int local_advection_spatial_order;
    bool local_semi_lagrangian_advection;

    FAST_LEVELSET_ADVECTION(T_FAST_LEVELSET* fast_levelset):
        BASE(fast_levelset)
    {
        Use_Level_Set_Advection_Method();
    };

    FAST_LEVELSET_ADVECTION():
        BASE(0) {};
            
    void Use_Local_WENO_For_Advection()
    {local_advection_spatial_order=5;local_semi_lagrangian_advection=false;}

    void Use_Local_ENO_For_Advection(const int order=3)
    {local_advection_spatial_order=order;local_semi_lagrangian_advection=false;assert(order >= 1 && order <= 3);}
            
    void Use_Local_Semi_Lagrangian_Advection()
    {local_semi_lagrangian_advection=true;local_advection_spatial_order=0;}

    void Use_Level_Set_Advection_Method()
    {local_semi_lagrangian_advection=false;local_advection_spatial_order=0;}

    void Euler_Step(const T_ARRAYS_VECTOR& velocity,const T dt,const T time,const int number_of_ghost_cells);
    void Euler_Step(const T_FACE_ARRAYS_SCALAR& velocity,const T dt,const T time,const int number_of_ghost_cells);
    void Euler_Step_Of_Reinitialization(const T_ARRAYS_SCALAR& signed_distance,const T_ARRAYS_SCALAR& sign_phi,const T dt,const T time);
    void Reinitialize(const int time_steps=10,const T time=0);
};

}
#endif

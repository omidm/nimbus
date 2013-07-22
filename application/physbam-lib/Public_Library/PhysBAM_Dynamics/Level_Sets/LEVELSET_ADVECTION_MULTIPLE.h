//#####################################################################
// Copyright 2009, Doug Enright, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso, Neil Molino, Avi Robinson-Mosher, Tamar Shinar, Jerry Talton, Elliot English.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_ADVECTION_MULTIPLE
//##################################################################### 
#ifndef __LEVELSET_ADVECTION_MULTIPLE__
#define __LEVELSET_ADVECTION_MULTIPLE__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_UNIFORM.h>
#include <PhysBAM_Geometry/Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_POLICY.h>

namespace PhysBAM {
    
template<class T_GRID>
class LEVELSET_ADVECTION_MULTIPLE
{
    typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET_MULTIPLE T_LEVELSET_MULTIPLE;
    typedef typename LEVELSET_ADVECTION_POLICY<T_GRID>::FAST_LEVELSET_ADVECTION_T T_FAST_LEVELSET_ADVECTION;
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
public:
    T_LEVELSET_MULTIPLE* levelsets;
    ARRAY<T_FAST_LEVELSET_ADVECTION> levelset_advections;

    LEVELSET_ADVECTION_MULTIPLE(T_LEVELSET_MULTIPLE* _levelsets):levelsets(_levelsets)
    {
        for(int i=1;i<=levelsets->levelsets.m;i++)
            levelset_advections.Append(T_FAST_LEVELSET_ADVECTION(levelsets->levelsets(i)));
    }
    
    void Set_Custom_Advection(ADVECTION<T,T,T_GRID>& advection_input)
    {for(int i=1;i<=levelset_advections.m;i++)levelset_advections(i).Set_Custom_Advection(advection_input);}
    
    void Euler_Step(const T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,const int number_of_ghost_cells)
    {for(int i=1;i<=levelset_advections.m;i++) levelset_advections(i).Euler_Step(face_velocities,dt,time,number_of_ghost_cells);}
    
    void Euler_Step(const T_ARRAYS_VECTOR& velocity,const T dt,const T time,const int number_of_ghost_cells)
    {for(int i=1;i<=levelset_advections.m;i++) levelset_advections(i).Euler_Step(velocity,dt,time,number_of_ghost_cells);}
    
    void Reinitialize()
    {for(int i=1;i<=levelset_advections.m;i++) levelset_advections(i).Reinitialize();}
};

}
#endif

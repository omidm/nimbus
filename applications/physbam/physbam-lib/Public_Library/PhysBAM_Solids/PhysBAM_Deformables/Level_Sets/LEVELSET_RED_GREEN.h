#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_RED_GREEN
//#####################################################################
#ifndef __LEVELSET_RED_GREEN__
#define __LEVELSET_RED_GREEN__

#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_POLICY.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_OCTREE.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_QUADTREE.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_GRID_2D.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_GRID_3D.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_POLICY.h>
namespace PhysBAM{

template<class TV> 
class LEVELSET_RED_GREEN:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename RED_GREEN_POLICY<TV>::GRID_T T_RED_GREEN_GRID;
    typedef typename DYADIC_GRID_POLICY<TV>::DYADIC_GRID T_DYADIC_GRID;
    typedef typename LEVELSET_POLICY<T_DYADIC_GRID>::LEVELSET T_LEVELSET_DYADIC;
public:
    T_RED_GREEN_GRID& grid;
    ARRAY<T>& phi;
    T_LEVELSET_DYADIC* overlay_levelset;
    T_DYADIC_GRID* overlay_grid;
    ARRAY<T>* overlay_phi;
    ARRAY<TV> overlay_velocity;
     
    LEVELSET_RED_GREEN(T_RED_GREEN_GRID& grid_input,ARRAY<T>& phi_input);
    ~LEVELSET_RED_GREEN();

    void Tree_Topology_Changed()
    {delete overlay_levelset;overlay_levelset=0;delete overlay_grid;overlay_grid=0;delete overlay_phi;overlay_phi=0;}

    T Phi(const TV& location) const 
    {return grid.Interpolate_Nodes(phi,location);}

//#####################################################################
    void Lazy_Update_Overlay_Levelset(const ARRAY<TV>* velocity=0);
    void Euler_Step(const ARRAY<TV>& velocity,const T dt,const T time=0);
    void Fast_Marching_Method();
//#####################################################################
};
}
#endif
#endif

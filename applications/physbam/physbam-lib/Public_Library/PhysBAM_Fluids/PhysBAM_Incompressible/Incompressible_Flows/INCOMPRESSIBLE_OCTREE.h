#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_OCTREE
//#####################################################################
//
// Solves V_t+V GRAD(V)+GRAD(p)=0.
//
//#####################################################################
#ifndef __INCOMPRESSIBLE_OCTREE__
#define __INCOMPRESSIBLE_OCTREE__

#include <PhysBAM_Tools/Boundaries/BOUNDARY.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_DYADIC.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_DYADIC.h>
namespace PhysBAM{

template<class T_input>
class INCOMPRESSIBLE_OCTREE:public INCOMPRESSIBLE_DYADIC<OCTREE_GRID<T_input> >
{
    typedef T_input T;
public:
    typedef INCOMPRESSIBLE_DYADIC<OCTREE_GRID<T> > BASE;
    using BASE::grid;using BASE::projection;using BASE::viscosity;using BASE::force;using BASE::max_time_step;using BASE::use_force;using BASE::gravity;
    using BASE::downward_direction;using BASE::Set_Custom_Advection;

    INCOMPRESSIBLE_OCTREE(OCTREE_GRID<T>& grid_input,const bool flame=false) 
        :INCOMPRESSIBLE_DYADIC<OCTREE_GRID<T> >(grid_input,flame)
    {}

    ~INCOMPRESSIBLE_OCTREE()
    {}

//#####################################################################
    void Compute_Vorticity_Confinement_Force(OCTREE_GRID<T>& grid,const ARRAY<VECTOR<T,3> >& V_ghost,ARRAY<VECTOR<T,3> >& F) PHYSBAM_OVERRIDE;
    static T Compute_Aggregate_Face_Velocities(void* incompressible,const T v1,const T v2,const T v3,const T v4);
    static void Interpolate_Face_Velocities_Valid_To_Direct_Children(const OCTREE_CELL<T>* cell,ARRAY<bool>& face_values);
//#####################################################################
};
}
#endif
#endif

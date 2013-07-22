#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2002-2005, Doug Enright, Ronald Fedkiw, Geoffrey Irving, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_CHILDREN.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_CHILDREN.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Advection/ADVECTION_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Advection/ADVECTION_SEMI_LAGRANGIAN_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_OCTREE.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_QUADTREE.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_DYADIC.h>

using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T_GRID> void LEVELSET_ADVECTION_DYADIC<T_GRID>::
Euler_Step(const ARRAY<T>& face_velocity,const T dt,const T time)
{
    T_GRID& grid=levelset->grid;
    T_ARRAYS_SCALAR& phi=levelset->phi;
    
    ARRAY<T> phi_ghost(grid.number_of_cells);levelset->boundary->Fill_Ghost_Cells_Cell(grid,phi,phi_ghost,time); 

    if(levelset->curvature_motion) PHYSBAM_NOT_IMPLEMENTED(); // do curvature first - based on phi^n, not supported yet

    advection->Update_Advection_Equation_Cell(grid,phi,phi_ghost,face_velocity,*levelset->boundary,dt,time);
    levelset->boundary->Apply_Boundary_Condition(grid,phi,time+dt); 
}

template<class T_GRID> void LEVELSET_ADVECTION_DYADIC<T_GRID>::
Use_Semi_Lagrangian_Advection()
{
    static typename ADVECTION_POLICY<T_GRID>::ADVECTION_SEMI_LAGRANGIAN_SCALAR semi_lagrangian_advection;
    Set_Custom_Advection(semi_lagrangian_advection);
}

template class LEVELSET_ADVECTION_DYADIC<OCTREE_GRID<float> >;
template class LEVELSET_ADVECTION_DYADIC<QUADTREE_GRID<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_ADVECTION_DYADIC<OCTREE_GRID<double> >;
template class LEVELSET_ADVECTION_DYADIC<QUADTREE_GRID<double> >;
#endif
#endif

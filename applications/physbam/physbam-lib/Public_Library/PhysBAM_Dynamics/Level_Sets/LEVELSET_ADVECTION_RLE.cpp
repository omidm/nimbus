#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2002-2005, Doug Enright, Ronald Fedkiw, Geoffrey Irving, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_RLE_Advection/ADVECTION_POLICY_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Advection/ADVECTION_SEMI_LAGRANGIAN_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Boundaries/BOUNDARY_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/AVERAGING_RLE.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Geometry/Grids_RLE_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_RLE.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_RLE.h>
using namespace PhysBAM;
//#####################################################################
// Function Use_Semi_Lagrangian_Advection
//#####################################################################
template<class T_GRID> void LEVELSET_ADVECTION_RLE<T_GRID>::
Use_Semi_Lagrangian_Advection()
{
    static typename ADVECTION_POLICY<T_GRID>::ADVECTION_SEMI_LAGRANGIAN_SCALAR semi_lagrangian_advection;
    Set_Custom_Advection(semi_lagrangian_advection);
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T_GRID> void LEVELSET_ADVECTION_RLE<T_GRID>::
Euler_Step(const ARRAY<T>& V,const T dt,const T time)
{
    DEBUG_UTILITIES::Debug_Breakpoint();
    
    T_GRID& grid=levelset->grid;
    T_BOUNDARY_SCALAR* boundary=levelset->boundary;
    T_ARRAYS_SCALAR& phi=levelset->phi;
    
    ARRAY<T> phi_ghost(grid.number_of_cells,false);boundary->Fill_Ghost_Cells_Cell(grid,phi,phi_ghost,time);

    if(levelset->curvature_motion){
        levelset->Curvature_Motion(dt,phi_ghost); // do curvature first - based on phi^n
        boundary->Fill_Ghost_Cells_Cell(grid,phi,phi_ghost,time);}

    advection->Update_Advection_Equation_Cell(grid,phi,phi_ghost,V,*boundary,dt,time);
    boundary->Apply_Boundary_Condition(grid,phi,time+dt);
}

template class LEVELSET_ADVECTION_RLE<RLE_GRID_2D<float> >;
template class LEVELSET_ADVECTION_RLE<RLE_GRID_3D<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_ADVECTION_RLE<RLE_GRID_2D<double> >;
template class LEVELSET_ADVECTION_RLE<RLE_GRID_3D<double> >;
#endif
#endif

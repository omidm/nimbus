//#####################################################################
// Copyright 2010, Jon Gretarsson, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEMI_LAGRANGIAN_CONSERVATION  
//##################################################################### 
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UTILITIES.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/BOUNDARY_OBJECT.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/SEMI_LAGRANGIAN_CONSERVATION.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/BOUNDARY_OBJECT_REFLECTION.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER.h>
using namespace PhysBAM;
//#####################################################################
// Function Update_Conservation_Law
//#####################################################################
template<class T_GRID,int d> void SEMI_LAGRANGIAN_CONSERVATION<T_GRID,d>::
Update_Conservation_Law(T_GRID& grid,T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const T_ARRAYS_BOOL& psi,const T dt,
    VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>& eigensystems,VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>& eigensystems_explicit,const T_FACE_ARRAYS_BOOL& psi_N,
    const T_FACE_ARRAYS_SCALAR& face_velocities,const bool thinshell,const TV_BOOL& outflow_boundaries,VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>* eigensystems_auxiliary,
    T_FACE_ARRAYS_DIMENSION_SCALAR* fluxes_auxiliary)
{
    BOUNDARY_UNIFORM<T_GRID,TV_DIMENSION> local_boundary;
    T_ARRAYS_DIMENSION_SCALAR local_U_ghost(grid.Domain_Indices(7));
    local_boundary.Fill_Ghost_Cells(grid,U,local_U_ghost,dt,0,7);
    advection.Update_Advection_Equation_Cell(grid,U,local_U_ghost,face_velocities,local_boundary,dt,0);
}
template class SEMI_LAGRANGIAN_CONSERVATION<GRID<VECTOR<float,1> >,3>;
template class SEMI_LAGRANGIAN_CONSERVATION<GRID<VECTOR<float,2> >,4>;
template class SEMI_LAGRANGIAN_CONSERVATION<GRID<VECTOR<float,3> >,5>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SEMI_LAGRANGIAN_CONSERVATION<GRID<VECTOR<double,1> >,3>;
template class SEMI_LAGRANGIAN_CONSERVATION<GRID<VECTOR<double,2> >,4>;
template class SEMI_LAGRANGIAN_CONSERVATION<GRID<VECTOR<double,3> >,5>;
#endif

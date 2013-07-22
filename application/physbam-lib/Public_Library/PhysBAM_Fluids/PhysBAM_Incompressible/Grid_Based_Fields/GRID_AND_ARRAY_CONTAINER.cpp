//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_AND_ARRAY_CONTAINER
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Advection/ADVECTION_SEMI_LAGRANGIAN_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_BLOCK_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_BLOCK_3D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_3D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_HORIZONTAL.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_Y.h>
#include <PhysBAM_Tools/Grids_RLE_Advection/ADVECTION_SEMI_LAGRANGIAN_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Boundaries/BOUNDARY_POLICY_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Boundaries/BOUNDARY_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/AVERAGING_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/FACE_LOOKUP_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/LINEAR_INTERPOLATION_RLE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_REFLECTION_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Grid_Based_Fields/GRID_AND_ARRAY_CONTAINER.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T2> GRID_AND_ARRAY_CONTAINER<T_GRID,T2>::
GRID_AND_ARRAY_CONTAINER(T_GRID& grid_input)
    :grid(grid_input),advection_maccormack(0),advection_default(*new T_ADVECTION_SEMI_LAGRANGIAN_SCALAR),boundary_default(*new T_BOUNDARY_REFLECTION),face_velocities(0),cell_velocities(0)
{
    advection=&advection_default;
    boundary=&boundary_default;
    Initialize_Array();
    Initialize_Domain_Boundary_Conditions();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID,class T2> GRID_AND_ARRAY_CONTAINER<T_GRID,T2>::
~GRID_AND_ARRAY_CONTAINER()
{
    delete advection_maccormack;
    delete &advection_default;delete &boundary_default;
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T_GRID,class T2> void GRID_AND_ARRAY_CONTAINER<T_GRID,T2>::
Euler_Step(const T dt,const T time,const int number_of_ghost_cells)
{
    T_ARRAYS_T2 array_ghost(grid.Cell_Indices(number_of_ghost_cells),false);boundary->Fill_Ghost_Cells(grid,array,array_ghost,dt,time,number_of_ghost_cells);
    if(cell_velocities) advection->Update_Advection_Equation_Cell(grid,array,array_ghost,*cell_velocities,*boundary,dt,time);
    else advection->Update_Advection_Equation_Cell(grid,array,array_ghost,*face_velocities,*boundary,dt,time);
}
//#####################################################################
template class GRID_AND_ARRAY_CONTAINER<GRID<VECTOR<float,1> >,float>;
template class GRID_AND_ARRAY_CONTAINER<GRID<VECTOR<float,2> >,float>;
template class GRID_AND_ARRAY_CONTAINER<GRID<VECTOR<float,3> >,float>;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template class GRID_AND_ARRAY_CONTAINER<QUADTREE_GRID<float>,float>;
template class GRID_AND_ARRAY_CONTAINER<OCTREE_GRID<float>,float>;
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template class GRID_AND_ARRAY_CONTAINER<RLE_GRID_2D<float>,float>;
template class GRID_AND_ARRAY_CONTAINER<RLE_GRID_3D<float>,float>;
#endif
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_AND_ARRAY_CONTAINER<GRID<VECTOR<double,1> >,double>;
template class GRID_AND_ARRAY_CONTAINER<GRID<VECTOR<double,2> >,double>;
template class GRID_AND_ARRAY_CONTAINER<GRID<VECTOR<double,3> >,double>;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template class GRID_AND_ARRAY_CONTAINER<QUADTREE_GRID<double>,double>;
template class GRID_AND_ARRAY_CONTAINER<OCTREE_GRID<double>,double>;
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template class GRID_AND_ARRAY_CONTAINER<RLE_GRID_2D<double>,double>;
template class GRID_AND_ARRAY_CONTAINER<RLE_GRID_3D<double>,double>;
#endif
#endif

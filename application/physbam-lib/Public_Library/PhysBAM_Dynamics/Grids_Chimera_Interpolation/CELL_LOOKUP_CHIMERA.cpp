//#####################################################################
// Copyright 2005-2010, Eran Guendelman, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Dynamics/Grids_Chimera_Interpolation/CELL_LOOKUP_CHIMERA.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Dynamics/Grids_Chimera_PDE_Linear/Parallel_Computation/LAPLACE_CHIMERA_GRID_MPI.h>

using namespace PhysBAM;

//#####################################################################
// Function operator()
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR CELL_LOOKUP_CHIMERA<T_GRID>::operator()(const GRID_CELL_INDEX& grid_cell_index) const
{
    if(laplace_grid.Local_Grid(grid_cell_index.x)) return (*u(grid_cell_index.x))(grid_cell_index.y);
    else return u_boundary(grid_cell_index.x)(laplace_grid.boundary_cell_indices_to_linear_index(grid_cell_index.x).Get(grid_cell_index.y));
}
//#####################################################################
// Function Valid
//#####################################################################
/*template<class T_GRID> bool CELL_LOOKUP_CHIMERA<T_GRID>::Valid(const GRID_CELL_INDEX& grid_cell_index) const
{
    if(laplace_grid.Local_Grid(grid_cell_index.x)) return u_valid==0 || (*u_valid)(grid_cell_index.x)(grid_cell_index.y);
    else return u_boundary_valid==0 || (*u_boundary_valid)(grid_cell_index.x).Contains(grid_cell_index.y);
    }*/
//#####################################################################

#define INSTANTIATION_HELPER(T,D)               \
    template class CELL_LOOKUP_CHIMERA<GRID<VECTOR<T,D> > >;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif

//#####################################################################
// Copyright 2005-2006, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __AVERAGING_CHIMERA__
#define __AVERAGING_CHIMERA__

#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>

namespace PhysBAM{

template<class T_GRID> class LAPLACE_CHIMERA_MPI;
template<class T_GRID> class LAPLACE_CHIMERA_MPI_CALLBACKS;
template<class T_GRID> class LAPLACE_CHIMERA_GRID_MPI;

template<class T_GRID>
class AVERAGING_CHIMERA
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    typedef FACE_INDEX<TV::dimension> D_FACE_INDEX;
    typedef PAIR<int,TV_INT> GRID_CELL_INDEX;
public:
    const LAPLACE_CHIMERA_GRID_MPI<T_GRID>& laplace_grid;
    SPARSE_MATRIX_FLAT_MXN<T> face_to_cell_matrix;
    VECTOR_ND<T> face_to_cell_vector;
    //v_c=A*v_f+b //v_f is the vector for solved faces

    AVERAGING_CHIMERA(const LAPLACE_CHIMERA_GRID_MPI<T_GRID>& laplace_grid_input):
        laplace_grid(laplace_grid_input) {}
    ~AVERAGING_CHIMERA() {}
    
    //negative grid index indicates voronoi face
    ARRAY<TRIPLE<int,FACE_INDEX<T_GRID::VECTOR_T::dimension>,TV> > Face_To_Cell_Vector_Weights(const GRID_CELL_INDEX& grid_cell_index) const;
    TV Face_To_Cell_Vector(LAPLACE_CHIMERA_MPI<T_GRID>& laplace,LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time,const VECTOR_ND<T>& u,const GRID_CELL_INDEX& grid_cell_index) const;
    void Build_Face_To_Cell_Matrix(LAPLACE_CHIMERA_MPI<T_GRID>& laplace,LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time);
    void Build_Face_To_Cell_Vector(LAPLACE_CHIMERA_MPI<T_GRID>& laplace,LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time);
//#####################################################################
};
}
#endif

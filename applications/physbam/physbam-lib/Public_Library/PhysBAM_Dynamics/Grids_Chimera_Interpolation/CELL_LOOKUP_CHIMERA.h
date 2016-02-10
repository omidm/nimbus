//#####################################################################
// Copyright 2005-2010, Geoffrey Irving, Michael Lentine, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CELL_LOOKUP_CHIMERA
//#####################################################################
#ifndef __CELL_LOOKUP_CHIMERA__
#define __CELL_LOOKUP_CHIMERA__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>

namespace PhysBAM{

template<class T_GRID> class LAPLACE_CHIMERA_GRID_MPI;

template<class T_GRID>
class CELL_LOOKUP_CHIMERA
{
    typedef LAPLACE_CHIMERA_GRID_MPI<T_GRID> T_LAPLACE_GRID;
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_LAPLACE_GRID::GRID_CELL_INDEX GRID_CELL_INDEX;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;

public:
    typedef T ELEMENT;
    
    T_LAPLACE_GRID& laplace_grid;
    const ARRAY<T_ARRAYS_SCALAR*>& u;
    const ARRAY<ARRAY<T> >& u_boundary;
//    const ARRAY<T_ARRAYS_BOOL>* u_valid;
//    const ARRAY<HASHTABLE<TV_INT> >* u_boundary_valid;

    CELL_LOOKUP_CHIMERA(T_LAPLACE_GRID& laplace_grid_input,const ARRAY<T_ARRAYS_SCALAR*>& u_input,const ARRAY<ARRAY<T> >& u_boundary_input/*,const ARRAY<T_ARRAYS_BOOL>* u_valid_input=0,const ARRAY<HASHTABLE<TV_INT> >* u_boundary_valid_input=0*/)
     :laplace_grid(laplace_grid_input),u(u_input),u_boundary(u_boundary_input)//,u_valid(u_valid_input),u_boundary_valid(u_boundary_valid_input)
    {}

    T operator()(const GRID_CELL_INDEX& grid_cell_index) const;
    //bool Valid(const GRID_CELL_INDEX& grid_cell_index) const;
};
}
#endif

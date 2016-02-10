//#####################################################################
// Copyright 2007, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_PROJECTION
//#####################################################################
#ifndef __EULER_PROJECTION__
#define __EULER_PROJECTION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER.h>
namespace PhysBAM{

template<class T_GRID>
class EULER_PROJECTION:public NONCOPYABLE
{
    typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::VECTOR_T TV;
    typedef typename T_GRID::INDEX INDEX;
    typedef VECTOR<T,T_GRID::dimension+2> TV_U;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_U>::TYPE T_ARRAYS_U;

public:
    EULER_PROJECTION()
    {}

    virtual ~EULER_PROJECTION()
    {}
    //#####################################################################

    static void Compute_One_Over_rho_c_squared(const T_GRID& grid,const T_ARRAYS_U& U_ghost,const EOS<T>* eos,T_ARRAYS_SCALAR& one_over_rho_c_squared)
    {
        for(typename T_GRID::CELL_ITERATOR  iterator(grid,grid.number_of_ghost_cells);iterator.Valid();iterator.Next()){
            INDEX cell_index=iterator.Cell_Index();
            const TV_U& U=U_ghost(cell_index);
            T rho=U(1); T c_squared=eos->c_squared(rho,EULER<T_GRID>::e(U));
            one_over_rho_c_squared(cell_index)=(T)1/(rho*c_squared);}
    }
};
}
#endif

//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_INTERIOR_ITERATOR_CELL
//#####################################################################
#ifndef __UNIFORM_INTERIOR_ITERATOR_CELL__
#define __UNIFORM_INTERIOR_ITERATOR_CELL__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>

namespace PhysBAM{

template<class TV> class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO;
template<class TV>
class UNIFORM_INTERIOR_ITERATOR_CELL:public GRID<TV>::CELL_ITERATOR
{
public:
    enum WORKAROUND {d=TV::dimension};
    typedef typename TV::SCALAR T;typedef VECTOR<int,d> TV_INT;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename GRID<TV>::CELL_ITERATOR BASE;
    using BASE::Valid;using BASE::index;

    const T_ARRAYS_BOOL& outside_fluid;

    explicit UNIFORM_INTERIOR_ITERATOR_CELL(const UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info,const int number_of_ghost_cells_input=0,
        const typename GRID<TV>::REGION& region_type_input=GRID<TV>::WHOLE_REGION,const int side_input=0);

    void Next()
    {BASE::Next();while(Valid() && !outside_fluid(index)) BASE::Next();}
};
}
#endif

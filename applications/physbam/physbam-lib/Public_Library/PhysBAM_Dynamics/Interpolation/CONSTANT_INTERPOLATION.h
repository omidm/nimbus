//#####################################################################
// Copyright 2003-2006, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSTANT_INTERPOLATION 
//#####################################################################
#ifndef __CONSTANT_INTERPOLATION__
#define __CONSTANT_INTERPOLATION__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID,class T2>
class CONSTANT_INTERPOLATION:public INTERPOLATION_UNIFORM<T_GRID,T2>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;
public:

    CONSTANT_INTERPOLATION()
    {}

    T2 Clamped_To_Array(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const PHYSBAM_OVERRIDE
    {return u(INTERPOLATION_UNIFORM<T_GRID,T2>::Clamped_Index(grid,u,X));}
    
    T2 From_Base_Node(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X,const TV_INT& index) const PHYSBAM_OVERRIDE
    {return u(index);}

//#####################################################################
};
}
#endif

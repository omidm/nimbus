//#####################################################################
// Copyright 2004-2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXTERNAL_STRAIN_ADJUSTMENT
//##################################################################### 
#ifndef __EXTERNAL_STRAIN_ADJUSTMENT__
#define __EXTERNAL_STRAIN_ADJUSTMENT__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UNIFORM_FORWARD.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T>
class EXTERNAL_STRAIN_ADJUSTMENT
{    
public:
    EXTERNAL_STRAIN_ADJUSTMENT()
    {}

    virtual ~EXTERNAL_STRAIN_ADJUSTMENT()
    {}

//#####################################################################
    virtual void Adjust_Strain(ARRAY<SYMMETRIC_MATRIX<T,2> ,VECTOR<int,2> >& e_ghost,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Adjust_Strain(ARRAY<SYMMETRIC_MATRIX<T,3> ,VECTOR<int,3> >& e_ghost,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
//#####################################################################
};
}
#endif

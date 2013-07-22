//#####################################################################
// Copyright 2012, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MODEL_REDUCTION
// Base class for other model reduction classes
//#####################################################################
#ifndef _MODEL_REDUCTION_
#define _MODEL_REDUCTION_

#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>

namespace PhysBAM{

template<class T>
class MODEL_REDUCTION
{
protected:
    MATRIX_MXN<T> reduced_basis;
public:
    MODEL_REDUCTION()
    {}

    ~MODEL_REDUCTION()
    {}

    virtual MATRIX_MXN<T> Get_Reduced_Basis(){return reduced_basis;}
    virtual void Set_Reduced_Basis(MATRIX_MXN<T>& m_in){reduced_basis=m_in;}
//#####################################################################
    //virtual void Generate_Reduced_Basis();
//##################################################################### 
};
}
#endif

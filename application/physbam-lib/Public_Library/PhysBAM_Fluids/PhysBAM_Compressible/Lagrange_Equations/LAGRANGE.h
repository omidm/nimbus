//#####################################################################
// Copyright 2002, Ronald Fedkiw
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAGRANGE  
//##################################################################### 
//
// Inherited by LAGRANGE_1D, LAGRANGE_2D, and LAGRANGE_3D.
// Initializes the equation of state information.
//
//#####################################################################
#ifndef __LAGRANGE__
#define __LAGRANGE__    

#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS.h>
namespace PhysBAM{

template<class T>
class LAGRANGE
{
public:
    EOS<T>& eos; // needed for equation of state functions

protected:
    LAGRANGE(EOS<T>& eos_input)
            :eos(eos_input)
    {}

//#####################################################################
};
}    
#endif

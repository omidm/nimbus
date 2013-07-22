//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HAMILTONIAN_1D
//#####################################################################
#ifndef __HAMILTONIAN_1D__
#define __HAMILTONIAN_1D__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_1D.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class T>
class HAMILTONIAN_1D
{
    typedef VECTOR<T,1> TV;
protected:
    GRID<TV>& grid;

    HAMILTONIAN_1D(GRID<TV>& grid_input)
        :grid(grid_input)
    {}

public:
//#####################################################################
    virtual ~HAMILTONIAN_1D(){}
    virtual T H(const T phi_x,const int  i=0,const T t=0){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T Maxabs_H1(const T phi_x_1,const T phi_x_2,const int i=0,const T t=0){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
};
}
#endif

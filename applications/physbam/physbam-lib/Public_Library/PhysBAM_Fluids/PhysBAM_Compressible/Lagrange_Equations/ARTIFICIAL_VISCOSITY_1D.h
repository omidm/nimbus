//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARTIFICIAL_VISCOSITY_1D
//#####################################################################
#ifndef __ARTIFICIAL_VISCOSITY_1D__
#define __ARTIFICIAL_VISCOSITY_1D__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/ARTIFICIAL_VISCOSITY.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/GRID_LAGRANGE_1D.h>
#include <iostream>
namespace PhysBAM{

template<class T>
class ARTIFICIAL_VISCOSITY_1D:public ARTIFICIAL_VISCOSITY
{
protected:
    ARTIFICIAL_VISCOSITY_1D()
    {}

public:
//#####################################################################
    virtual ~ARTIFICIAL_VISCOSITY_1D(){}
    virtual void Get_Artificial_Viscosity(EOS<T>& eos,GRID_LAGRANGE_1D<T>& grid,const ARRAY<T,VECTOR<int,1> >& mass,const ARRAY<T,VECTOR<int,1> >& velocity,const ARRAY<T,VECTOR<int,1> >& energy,
        ARRAY<T,VECTOR<int,1> >& Q){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
};
}
#endif

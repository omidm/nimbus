//#####################################################################
// Copyright 2002, Ronald Fedkiw
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARTIFICIAL_VISCOSITY_2D
//#####################################################################
#ifndef __ARTIFICIAL_VISCOSITY_2D__
#define __ARTIFICIAL_VISCOSITY_2D__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/ARTIFICIAL_VISCOSITY.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/GRID_LAGRANGE_2D.h>
#include <iostream>
namespace PhysBAM{

template<class T>
class ARTIFICIAL_VISCOSITY_2D:public ARTIFICIAL_VISCOSITY
{
protected:
    ARTIFICIAL_VISCOSITY_2D()
    {}

public:
//#####################################################################
    virtual ~ARTIFICIAL_VISCOSITY_2D(){}
    virtual void Get_Artificial_Viscosity(EOS<T>& eos,GRID_LAGRANGE_2D<T>& grid,const ARRAY<T,VECTOR<int,2> >& mass,const ARRAY<T,VECTOR<int,2> >& u,const ARRAY<T,VECTOR<int,2> >& v,const ARRAY<T,VECTOR<int,2> >& energy,ARRAY<T,VECTOR<int,2> >& Q1,
        ARRAY<T,VECTOR<int,2> >& Q2,ARRAY<T,VECTOR<int,2> >& Q3,ARRAY<T,VECTOR<int,2> >& Q4){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
};
}
#endif

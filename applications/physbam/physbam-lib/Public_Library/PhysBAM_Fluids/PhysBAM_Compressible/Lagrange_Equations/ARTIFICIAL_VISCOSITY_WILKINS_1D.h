//#####################################################################
// Copyright 2002, Ronald Fedkiw
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARTIFICIAL_VISCOSITY_WILKINS_1D  
//##################################################################### 
//
// Wilkins artifical viscosity used in the HEMP code. First presented by Kuropatenko. 
//
//#####################################################################
#ifndef __ARTIFICIAL_VISCOSITY_WILKINS_1D__
#define __ARTIFICIAL_VISCOSITY_WILKINS_1D__    

#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/ARTIFICIAL_VISCOSITY_1D.h>
namespace PhysBAM{

template<class T>
class ARTIFICIAL_VISCOSITY_WILKINS_1D:public ARTIFICIAL_VISCOSITY_1D<T>
{
private:
    using ARTIFICIAL_VISCOSITY_1D<T>::limiter;

    T linear_constant;    // coefficient of artificial viscosity - linear term
    T quadratic_constant; // coefficient of artificial viscosity - quadratic term

public:
    ARTIFICIAL_VISCOSITY_WILKINS_1D()
    {
        Set_Linear_Constant();
        Set_Quadratic_Constant();
    }

    void Set_Linear_Constant(const T linear_constant_input=.06)
    {linear_constant=linear_constant_input;}
    
    void Set_Quadratic_Constant(const T quadratic_constant_input=1.5)
    {quadratic_constant=quadratic_constant_input;}

//#####################################################################
    void Get_Artificial_Viscosity(EOS<T>& eos,GRID_LAGRANGE_1D<T>& grid,const ARRAY<T,VECTOR<int,1> >& mass,const ARRAY<T,VECTOR<int,1> >& velocity,const ARRAY<T,VECTOR<int,1> >& energy,ARRAY<T,VECTOR<int,1> >& Q) PHYSBAM_OVERRIDE;
//#####################################################################
};
}    
#endif

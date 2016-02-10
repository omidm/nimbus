//#####################################################################
// Copyright 2002, Ronald Fedkiw
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARTIFICIAL_VISCOSITY_VNR_1D  
//##################################################################### 
//
// Von Neumann - Richtmyer artifical viscosity.
//
//#####################################################################
#ifndef __ARTIFICIAL_VISCOSITY_VNR_1D__
#define __ARTIFICIAL_VISCOSITY_VNR_1D__    

#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/ARTIFICIAL_VISCOSITY_1D.h>
namespace PhysBAM{

template<class T>
class ARTIFICIAL_VISCOSITY_VNR_1D:public ARTIFICIAL_VISCOSITY_1D<T>
{
private:
    using ARTIFICIAL_VISCOSITY_1D<T>::limiter;

    T constant; // coefficient of artificial viscosity

public:
    ARTIFICIAL_VISCOSITY_VNR_1D()
    {
        Set_Constant();
    }

    void Set_Constant(const T constant_input=9)
    {constant=constant_input;}

//#####################################################################
    void Get_Artificial_Viscosity(EOS<T>& eos,GRID_LAGRANGE_1D<T>& grid,const ARRAY<T,VECTOR<int,1> >& mass,const ARRAY<T,VECTOR<int,1> >& velocity,const ARRAY<T,VECTOR<int,1> >& energy,ARRAY<T,VECTOR<int,1> >& Q) PHYSBAM_OVERRIDE;
//#####################################################################
};
}    
#endif


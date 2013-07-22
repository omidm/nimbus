//#####################################################################
// Copyright 2002, Ronald Fedkiw
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARTIFICIAL_VISCOSITY_VNR_2D  
//##################################################################### 
//
// Von Neumann - Richtmyer artifical viscosity.
//
//#####################################################################
#ifndef __ARTIFICIAL_VISCOSITY_VNR_2D__
#define __ARTIFICIAL_VISCOSITY_VNR_2D__    

#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/ARTIFICIAL_VISCOSITY_2D.h>
namespace PhysBAM{

template<class T>
class ARTIFICIAL_VISCOSITY_VNR_2D:public ARTIFICIAL_VISCOSITY_2D<T>
{
private:
    using ARTIFICIAL_VISCOSITY_2D<T>::limiter;

    T constant; // coefficient of artificial viscosity

public:
    ARTIFICIAL_VISCOSITY_VNR_2D()
    {
        Set_Constant();
    }

    void Set_Constant(const T constant_input=9)
    {constant=constant_input;}

//#####################################################################
    void Get_Artificial_Viscosity(EOS<T>& eos,GRID_LAGRANGE_2D<T>& grid,const ARRAY<T,VECTOR<int,2> >& mass,const ARRAY<T,VECTOR<int,2> >& u,const ARRAY<T,VECTOR<int,2> >& v,const ARRAY<T,VECTOR<int,2> >& energy,ARRAY<T,VECTOR<int,2> >& Q1,
        ARRAY<T,VECTOR<int,2> >& Q2,ARRAY<T,VECTOR<int,2> >& Q3,ARRAY<T,VECTOR<int,2> >& Q4);
//#####################################################################
};
}    
#endif

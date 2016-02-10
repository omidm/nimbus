//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HAMILTON_JACOBI
//##################################################################### 
//
// Inherited by HAMILTON_JACOBI_1D, HAMILTON_JACOBI_2D, and HAMILTON_JACOBI_3D.
//
//#####################################################################
#ifndef __HAMILTON_JACOBI__
#define __HAMILTON_JACOBI__ 

namespace PhysBAM{

class HAMILTON_JACOBI
{
protected:
    int spatial_order; // 1,2,3, or 5
    int LF_viscosity,LLF_viscosity,LLLF_viscosity;

protected:
    HAMILTON_JACOBI()
    {
        Set_WENO();  
        Use_LLF_Viscosity();
    }

public:
    void Set_WENO() // 5th order
    {spatial_order=5;}

    void Set_ENO(const int order=3) // can be 1, 2, or 3
    {spatial_order=order;}
    
    void Use_LF_Viscosity()
    {LF_viscosity=1;LLF_viscosity=0;LLLF_viscosity=0;}

    void Use_LLF_Viscosity()
    {LF_viscosity=0;LLF_viscosity=1;LLLF_viscosity=0;}

    void Use_LLLF_Viscosity()
    {LF_viscosity=0;LLF_viscosity=0;LLLF_viscosity=1;}

//#####################################################################
};   
}
#endif


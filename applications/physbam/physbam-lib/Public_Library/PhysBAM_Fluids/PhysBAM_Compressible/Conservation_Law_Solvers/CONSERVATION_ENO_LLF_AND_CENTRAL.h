//#####################################################################
// Copyright 2002-2004, Doug Enright, Ronald Fedkiw, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSERVATION_ENO_LLF_AND_CENTRAL  
//#####################################################################
//
// The central scheme can only be 1st or second order with a global alpha,
// regardless of the set order and alpha of ENO-LLF. Note that the code 
// temporarily changes these parameters when applying the central scheme.
//
//#####################################################################
#ifndef __CONSERVATION_ENO_LLF_AND_CENTRAL__
#define __CONSERVATION_ENO_LLF_AND_CENTRAL__   

#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION.h>
namespace PhysBAM{

template<class T_GRID,int d>
class CONSERVATION_ENO_LLF_AND_CENTRAL:public CONSERVATION<T_GRID,d>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;typedef VECTOR<T,d> TV_DIMENSION;
private:
    typedef CONSERVATION<T_GRID,d> BASE;
    using BASE::order;using BASE::field_by_field_alpha;using BASE::amplification_factor;

    int central_order; // 1 or 2, although order=1,2,or 3 for ENO-LLF
    T central_amplification_factor; // for amplifying alpha to increase dissipation

public:
    CONSERVATION_ENO_LLF_AND_CENTRAL()
    {
        Set_Central_Order();
        Amplify_Central_Alpha();
    }

    void Set_Central_Order(const int central_order_input=2)
    {central_order=central_order_input;}

    void Amplify_Central_Alpha(const T central_amplification_factor_input=1)
    {central_amplification_factor=central_amplification_factor_input;}

//#####################################################################
    void Conservation_Solver(const int m,const T dx,const ARRAY<bool,VECTOR<int,1> >& psi,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& Fx,EIGENSYSTEM<T,TV_DIMENSION>& eigensystem,EIGENSYSTEM<T,TV_DIMENSION>& eigensystem_explicit,
        const VECTOR<bool,2>& outflow_boundaries,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_flux=0) PHYSBAM_OVERRIDE;
private:
    template<int order> void Conservation_Solver_Helper(const int m,const T dx,const ARRAY<bool,VECTOR<int,1> >& psi,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& Fx,EIGENSYSTEM<T,TV_DIMENSION>& eigensystem,
        EIGENSYSTEM<T,TV_DIMENSION>& eigensystem_explicit,const VECTOR<bool,2>& outflow_boundaries);
//#####################################################################
};   
}
#endif

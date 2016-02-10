//#####################################################################
// Copyright 2002-2004, Doug Enright, Ronald Fedkiw, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSERVATION_ENO_RF  
//##################################################################### 
#ifndef __CONSERVATION_ENO_RF__
#define __CONSERVATION_ENO_RF__   

#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION.h>
namespace PhysBAM{

template<class T_GRID,int d>
class CONSERVATION_ENO_RF:public CONSERVATION<T_GRID,d>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;typedef VECTOR<T,d> TV_DIMENSION;
public:
    typedef CONSERVATION<T_GRID,d> BASE;
    using BASE::order;using BASE::save_fluxes;using BASE::flux_temp;

    CONSERVATION_ENO_RF()
    {}

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

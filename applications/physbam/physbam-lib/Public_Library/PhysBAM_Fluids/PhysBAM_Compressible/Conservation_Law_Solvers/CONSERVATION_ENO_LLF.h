//#####################################################################
// Copyright 2002-2009, Doug Enright, Ronald Fedkiw, Jon Gretarsson, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSERVATION_ENO_LLF  
//##################################################################### 
#ifndef __CONSERVATION_ENO_LLF__
#define __CONSERVATION_ENO_LLF__   

#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION.h>
namespace PhysBAM{

template<class T_GRID,int d>
class CONSERVATION_ENO_LLF:public CONSERVATION<T_GRID,d>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;typedef VECTOR<T,d> TV_DIMENSION;
public:
    typedef CONSERVATION<T_GRID,d> BASE;
    using BASE::order;using BASE::save_fluxes;using BASE::flux_temp;

    const bool use_global_llf,use_face_velocity_for_fluxes,use_standard_average,use_explicit_eigensystem_for_alphas;

    CONSERVATION_ENO_LLF(bool use_global_llf_input=false,bool use_face_velocity_for_fluxes_input=false,bool use_standard_average_input=true,
        bool use_explicit_eigensystem_for_alphas_input=false):use_global_llf(use_global_llf_input),
        use_face_velocity_for_fluxes(use_face_velocity_for_fluxes_input),use_standard_average(use_standard_average_input),
        use_explicit_eigensystem_for_alphas(use_explicit_eigensystem_for_alphas_input)
    {}

//#####################################################################
    void Conservation_Solver(const int m,const T dx,const ARRAY<bool,VECTOR<int,1> >& psi,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,
        ARRAY<TV_DIMENSION,VECTOR<int,1> >& Fx,EIGENSYSTEM<T,TV_DIMENSION>& eigensystem,EIGENSYSTEM<T,TV_DIMENSION>& eigensystem_explicit,
        const VECTOR<bool,2>& outflow_boundaries,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_flux=0) PHYSBAM_OVERRIDE;
private:
    template<int order> void Conservation_Solver_Helper(const int m,const T dx,const ARRAY<bool,VECTOR<int,1> >& psi,
        const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& Fx,EIGENSYSTEM<T,TV_DIMENSION>& eigensystem,
        EIGENSYSTEM<T,TV_DIMENSION>& eigensystem_explicit,const VECTOR<bool,2>& outflow_boundaries,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_flux);
    template<int order> void Conservation_Solver_Helper_Experimental(const int m,const T dx,const ARRAY<bool,VECTOR<int,1> >& psi,
        const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& Fx,EIGENSYSTEM<T,TV_DIMENSION>& eigensystem,
        EIGENSYSTEM<T,TV_DIMENSION>& eigensystem_explicit,const VECTOR<bool,2>& outflow_boundaries,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_flux);
//#####################################################################
};   
}
#endif

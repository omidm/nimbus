//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Reactive_Euler_Equations/REACTIVE_EULER.h>
#include <cmath>
using namespace PhysBAM;
//#####################################################################
// Function e
//#####################################################################
// internal energy
template<class T_GRID> typename T_GRID::SCALAR REACTIVE_EULER<T_GRID>::
e(const T rho,const T rho_u,const T E)
{
    return E/rho-sqr(rho_u/rho)/2;
}
//#####################################################################
// Function e
//#####################################################################
// internal energy
template<class T_GRID> typename T_GRID::SCALAR REACTIVE_EULER<T_GRID>::
e(const T rho,const T rho_u,const T rho_v,const T E)
{
    return E/rho-(sqr(rho_u/rho)+sqr(rho_v/rho))/2;
}
//#####################################################################
// Function e
//#####################################################################
// internal energy - 3D
template<class T_GRID> typename T_GRID::SCALAR REACTIVE_EULER<T_GRID>::
e(const T rho,const T rho_u,const T rho_v,const T rho_w,const T E)
{
    return E/rho-(sqr(rho_u/rho)+sqr(rho_v/rho)+sqr(rho_w/rho))/2;
}
//#####################################################################
template class REACTIVE_EULER<GRID<VECTOR<float,1> > >;
template class REACTIVE_EULER<GRID<VECTOR<float,2> > >;
template class REACTIVE_EULER<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class REACTIVE_EULER<GRID<VECTOR<double,1> > >;
template class REACTIVE_EULER<GRID<VECTOR<double,2> > >;
template class REACTIVE_EULER<GRID<VECTOR<double,3> > >;
#endif

//#####################################################################
// Copyright 2002-2009, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Parallel_Computation/THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM_DEFINITION.h>

using namespace PhysBAM;

#define INSTANTIATION_HELPER(T,T_GRID,d) \
    template class THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM<T_GRID,SYMMETRIC_MATRIX<T,d>,AVERAGING_UNIFORM<T_GRID,FACE_LOOKUP_UNIFORM<T_GRID > >,LINEAR_INTERPOLATION_UNIFORM<T_GRID,SYMMETRIC_MATRIX<T,d>,FACE_LOOKUP_UNIFORM<T_GRID > > >; \
    template class THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM<T_GRID,T,AVERAGING_UNIFORM<T_GRID,FACE_LOOKUP_UNIFORM<T_GRID > >,LINEAR_INTERPOLATION_UNIFORM<T_GRID,T,FACE_LOOKUP_UNIFORM<T_GRID > > >;
    
//template class THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM<T_GRID,T,AVERAGING_UNIFORM<T_GRID,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<T_GRID > >,LINEAR_INTERPOLATION_UNIFORM<T_GRID,T,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<T_GRID > > >;

#define P(...) __VA_ARGS__
INSTANTIATION_HELPER(float,P(GRID<VECTOR<float,1> >),2)
INSTANTIATION_HELPER(float,P(GRID<VECTOR<float,2> >),2);
INSTANTIATION_HELPER(float,P(GRID<VECTOR<float,3> >),3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,P(GRID<VECTOR<double,1> >),2);
INSTANTIATION_HELPER(double,P(GRID<VECTOR<double,2> >),2);
INSTANTIATION_HELPER(double,P(GRID<VECTOR<double,3> >),3);
#endif
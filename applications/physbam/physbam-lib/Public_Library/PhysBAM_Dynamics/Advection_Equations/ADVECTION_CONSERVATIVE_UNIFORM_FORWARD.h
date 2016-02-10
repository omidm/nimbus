//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header ADVECTION_CONSERVATIVE_UNIFORM_FORWARD
//#####################################################################
#ifndef __ADVECTION_CONSERVATIVE_UNIFORM_FORWARD__
#define __ADVECTION_CONSERVATIVE_UNIFORM_FORWARD__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_POLICY_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_AVERAGING=AVERAGING_UNIFORM<T_GRID>,class T_INTERPOLATION=LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2> > class ADVECTION_CONSERVATIVE_UNIFORM;

template<class T_ADVECTION,class T_NEW> struct REBIND;
template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION,class T_NEW> struct REBIND<ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>,T_NEW>{
    typedef ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T_NEW,T_AVERAGING,typename REBIND<T_INTERPOLATION,T_NEW>::TYPE> TYPE;};

}
#endif

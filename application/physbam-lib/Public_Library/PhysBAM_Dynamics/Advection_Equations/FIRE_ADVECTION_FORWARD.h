//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Avi Robinson-Mosher, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header ADVECTION_FORWARD
//#####################################################################
#ifndef __FIRE_ADVECTION_FORWARD__
#define __FIRE_ADVECTION_FORWARD__

namespace PhysBAM{

template<class TV> class GRID;

template<class T_GRID,class T2,class T_NESTED_ADVECTION> class ADVECTION_WRAPPER_FIRE_DYADIC;
template<class T_GRID,class T2,class T_NESTED_ADVECTION> class ADVECTION_WRAPPER_FIRE_MULTIPHASE_DYADIC;
template<class T_GRID,class T2,class T_NESTED_ADVECTION> class ADVECTION_WRAPPER_FIRE_MULTIPHASE_RLE;
template<class T_GRID,class T2,class T_NESTED_ADVECTION> class ADVECTION_WRAPPER_FIRE_MULTIPHASE_UNIFORM;
template<class T_GRID,class T2,class T_NESTED_ADVECTION> class ADVECTION_WRAPPER_FIRE_RLE;
}
#endif

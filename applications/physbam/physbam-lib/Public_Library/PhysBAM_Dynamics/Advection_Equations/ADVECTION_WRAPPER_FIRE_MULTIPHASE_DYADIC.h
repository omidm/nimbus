#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2005-2006, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_WRAPPER_FIRE_MULTIPHASE_DYADIC
//#####################################################################
#ifndef __ADVECTION_WRAPPER_FIRE_MULTIPHASE_DYADIC__
#define __ADVECTION_WRAPPER_FIRE_MULTIPHASE_DYADIC__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_POLICY_DYADIC.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_POLICY.h>
namespace PhysBAM{

template<class T_GRID> class PROJECTION_DYADIC;
template<class T_GRID,class T2,class T_NESTED_ADVECTION>
class ADVECTION_WRAPPER_FIRE_MULTIPHASE_DYADIC:public ADVECTION<T_GRID,T2>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename INCOMPRESSIBLE_POLICY<T_GRID>::PROJECTION T_PROJECTION;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET_MULTIPLE T_LEVELSET_MULTIPLE;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;typedef typename REBIND<T_BOUNDARY,T2>::TYPE T_BOUNDARY_T2;
public:

    ADVECTION_WRAPPER_FIRE_MULTIPHASE_DYADIC(T_NESTED_ADVECTION&,const PROJECTION_DYADIC<T_GRID>&,const T_LEVELSET_MULTIPLE&){PHYSBAM_NOT_IMPLEMENTED();}

    void Update_Advection_Equation_Node(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
        const T_ARRAYS_VECTOR& V,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const T_ARRAYS_T2* Z_min_ghost=0,const T_ARRAYS_T2* Z_max_ghost=0,T_ARRAYS_T2* Z_min=0,T_ARRAYS_T2* Z_max=0)
    {PHYSBAM_NOT_IMPLEMENTED();}

    void Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
        const FACE_LOOKUP_DYADIC<T_GRID>& face_velocities,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max)
    {PHYSBAM_NOT_IMPLEMENTED();}
    
    void Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const FACE_LOOKUP_DYADIC<T_GRID>& Z_ghost,
        const FACE_LOOKUP_DYADIC<T_GRID>& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
        const FACE_LOOKUP_DYADIC<T_GRID>* Z_min_ghost,const FACE_LOOKUP_DYADIC<T_GRID>* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max)
    {PHYSBAM_NOT_IMPLEMENTED();}

//#####################################################################
};
}
#endif
#endif

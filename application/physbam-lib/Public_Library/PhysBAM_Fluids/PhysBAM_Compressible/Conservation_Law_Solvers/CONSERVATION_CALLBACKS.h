//#####################################################################
// Copyright 2007, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSERVATION_CALLBACKS
//##################################################################### 
#ifndef __CONSERVATION_CALLBACKS__
#define __CONSERVATION_CALLBACKS__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T_GRID,class TV_DIMENSION>
class CONSERVATION_CALLBACKS
{
    typedef typename T_GRID::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename REBIND<T_ARRAYS_SCALAR,bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_FACE_ARRAYS_DIMENSION_SCALAR;
public:
    CONSERVATION_CALLBACKS()
    {}

    virtual ~CONSERVATION_CALLBACKS()
    {}

    virtual void Get_Neumann_Face_Location(const GRID<VECTOR<T,1> >& grid_1d,const int face_index,T& location) const
    {location=grid_1d.Face(1,VECTOR<int,1>(face_index)).x;}

    virtual void Clamp_Dt_Adaptively(T_GRID& grid,const T_ARRAYS_DIMENSION_SCALAR& rhs,const T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_BOOL& psi,T_ARRAYS_SCALAR& rho_dt,T_ARRAYS_SCALAR& e_dt,const T& dt,T& clamp_rho,T& clamp_e)
    {}

    virtual void Clamp_Fluxes(T_GRID& grid,T_FACE_ARRAYS_DIMENSION_SCALAR& fluxes,T_ARRAYS_DIMENSION_SCALAR& rhs,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const T_ARRAYS_BOOL& psi,const T& dt,T& clamp_rho,T& clamp_e)
    {}
};
}
#endif

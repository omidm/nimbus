//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS
//#####################################################################
#ifndef __INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS__
#define __INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS__
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
namespace PhysBAM{

template<class TV>
class INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS
{
    typedef typename TV::SCALAR T;typedef VECTOR<T,TV::dimension+2> TV_DIMENSION;typedef GRID<TV> T_GRID;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename REBIND<T_ARRAYS_SCALAR,bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;

public:
    INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS()
    {}
    virtual ~INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS()
    {}

//#####################################################################
    virtual void Compute_Compressible_Incompressible_Face_Pressures_From_Cell_Pressures(const T_GRID& face_grid,const T_ARRAYS_DIMENSION_SCALAR& U,
        const T_ARRAYS_BOOL& euler_psi,const T_ARRAYS_SCALAR& p_cell,T_FACE_ARRAYS_SCALAR& p_face) const=0;
    virtual void Fill_Incompressible_Beta_Face(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& beta_face) const=0;
//#####################################################################
};  
}   
#endif


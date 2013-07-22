//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_BOUNDARY_CONDITION
//#####################################################################
#ifndef __IMPLICIT_BOUNDARY_CONDITION__
#define __IMPLICIT_BOUNDARY_CONDITION__
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{
template<class TV> class GRID;
template<class TV>
class IMPLICIT_BOUNDARY_CONDITION
{
    typedef VECTOR<int,TV::dimension> TV_INT;typedef typename TV::SCALAR T;
public:
    IMPLICIT_BOUNDARY_CONDITION();
    virtual ~IMPLICIT_BOUNDARY_CONDITION();

    virtual void Update_Boundary_Conditions(const GRID<TV>& grid,ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N,ARRAY<T,TV_INT>& p,
        ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,const T time)=0;
};
}
#endif

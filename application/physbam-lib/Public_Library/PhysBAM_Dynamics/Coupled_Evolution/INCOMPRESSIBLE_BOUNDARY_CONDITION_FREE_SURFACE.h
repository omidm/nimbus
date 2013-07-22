//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_BOUNDARY_CONDITION_FREE_SURFACE
//#####################################################################
#ifndef __INCOMPRESSIBLE_BOUNDARY_CONDITION_FREE_SURFACE__
#define __INCOMPRESSIBLE_BOUNDARY_CONDITION_FREE_SURFACE__
#include <PhysBAM_Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION.h>
namespace PhysBAM{
template<class TV>
class INCOMPRESSIBLE_BOUNDARY_CONDITION_FREE_SURFACE:public IMPLICIT_BOUNDARY_CONDITION<TV>
{
    typedef VECTOR<int,TV::dimension> TV_INT;typedef typename TV::SCALAR T;
    const ARRAY<T,TV_INT>& phi;
public:

    INCOMPRESSIBLE_BOUNDARY_CONDITION_FREE_SURFACE(const ARRAY<T,TV_INT>& phi_input);
    virtual ~INCOMPRESSIBLE_BOUNDARY_CONDITION_FREE_SURFACE();

    void Update_Boundary_Conditions(const GRID<TV>& grid,ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N,ARRAY<T,TV_INT>& p,
        ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,const T time) PHYSBAM_OVERRIDE;
};
}
#endif

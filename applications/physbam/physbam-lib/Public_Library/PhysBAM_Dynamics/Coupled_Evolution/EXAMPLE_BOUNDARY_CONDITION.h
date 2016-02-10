//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXAMPLE_BOUNDARY_CONDITION
//#####################################################################
#ifndef __EXAMPLE_BOUNDARY_CONDITION__
#define __EXAMPLE_BOUNDARY_CONDITION__
#include <PhysBAM_Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION.h>
namespace PhysBAM{
template<class T_GRID> class BOUNDARY_CONDITIONS_CALLBACKS;

template<class TV>
class EXAMPLE_BOUNDARY_CONDITION:public IMPLICIT_BOUNDARY_CONDITION<TV>
{
    typedef VECTOR<int,TV::dimension> TV_INT;typedef typename TV::SCALAR T;typedef GRID<TV> T_GRID;
public:
    BOUNDARY_CONDITIONS_CALLBACKS<TV>* callback;

    EXAMPLE_BOUNDARY_CONDITION(BOUNDARY_CONDITIONS_CALLBACKS<TV>* callback_input);
    virtual ~EXAMPLE_BOUNDARY_CONDITION();

    void Update_Boundary_Conditions(const GRID<TV>& grid,ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N,ARRAY<T,TV_INT>& p,
        ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,const T time) PHYSBAM_OVERRIDE;
};
}
#endif

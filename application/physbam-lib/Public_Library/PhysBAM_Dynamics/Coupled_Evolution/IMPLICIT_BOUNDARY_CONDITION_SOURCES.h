//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_BOUNDARY_CONDITION_SOURCES
//#####################################################################
#ifndef __IMPLICIT_BOUNDARY_CONDITION_SOURCES__
#define __IMPLICIT_BOUNDARY_CONDITION_SOURCES__
#include <PhysBAM_Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION.h>
namespace PhysBAM{
template<class T_GRID> class FLUIDS_PARAMETERS_CALLBACKS;

template<class TV>
class IMPLICIT_BOUNDARY_CONDITION_SOURCES:public IMPLICIT_BOUNDARY_CONDITION<TV>
{
    typedef VECTOR<int,TV::dimension> TV_INT;typedef typename TV::SCALAR T;typedef GRID<TV> T_GRID;
public:
    FLUIDS_PARAMETERS_CALLBACKS<T_GRID>& callbacks;

    IMPLICIT_BOUNDARY_CONDITION_SOURCES(FLUIDS_PARAMETERS_CALLBACKS<T_GRID>& callbacks_input);
    virtual ~IMPLICIT_BOUNDARY_CONDITION_SOURCES();

    void Update_Boundary_Conditions(const GRID<TV>& grid,ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N,ARRAY<T,TV_INT>& p,
        ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,const T time) PHYSBAM_OVERRIDE;
};
}
#endif

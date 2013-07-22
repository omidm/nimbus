//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_PARAMETERS
//#####################################################################
#ifndef __SOLIDS_FLUIDS_PARAMETERS__
#define __SOLIDS_FLUIDS_PARAMETERS__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>

namespace PhysBAM{

template<class TV> class MPI_SOLID_FLUID;
template<class TV> class MPI_SOLID_FLUID_SLIP;
template<class TV> class SOLIDS_FLUIDS_CALLBACKS;

template <class TV>
class SOLIDS_FLUIDS_PARAMETERS:public NONCOPYABLE
{
public:
    SOLIDS_FLUIDS_CALLBACKS<TV>* callbacks;
    MPI_SOLID_FLUID<TV>* mpi_solid_fluid;
    MPI_SOLID_FLUID_SLIP<TV>* mpi_solid_fluid_slip;
    bool use_leakproof_solve,use_fluid_rigid_fracture;

    SOLIDS_FLUIDS_PARAMETERS(SOLIDS_FLUIDS_CALLBACKS<TV>* callbacks);
    virtual ~SOLIDS_FLUIDS_PARAMETERS();

//#####################################################################
};
}
#endif

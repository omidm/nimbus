//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_PARAMETERS
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID_SLIP.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_PARAMETERS.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLIDS_FLUIDS_PARAMETERS<TV>::
SOLIDS_FLUIDS_PARAMETERS(SOLIDS_FLUIDS_CALLBACKS<TV>* callbacks_input)
    :callbacks(callbacks_input),mpi_solid_fluid(0),mpi_solid_fluid_slip(0),use_leakproof_solve(true),use_fluid_rigid_fracture(false)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLIDS_FLUIDS_PARAMETERS<TV>::
~SOLIDS_FLUIDS_PARAMETERS()
{
    delete mpi_solid_fluid;
    //delete mpi_solid_fluid_slip;
}
//#####################################################################
template class SOLIDS_FLUIDS_PARAMETERS<VECTOR<float,1> >;
template class SOLIDS_FLUIDS_PARAMETERS<VECTOR<float,2> >;
template class SOLIDS_FLUIDS_PARAMETERS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SOLIDS_FLUIDS_PARAMETERS<VECTOR<double,1> >;
template class SOLIDS_FLUIDS_PARAMETERS<VECTOR<double,2> >;
template class SOLIDS_FLUIDS_PARAMETERS<VECTOR<double,3> >;
#endif
}

//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/FAST_PROJECTION_DYNAMICS_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_REFINEMENT_UNIFORM.h>
#include <PhysBAM_Dynamics/Parallel_Computation/REFINEMENT_THREADS.h>

using namespace PhysBAM;

template<class TV> void REFINEMENT_TASK<TV>::Run(const int threadid)
{
    //TODO: Prevent repeated allocation
    GRID<TV> local_mac_grid(TV_INT::All_Ones_Vector()*projection->coarse_scale,RANGE<TV>::Centered_Box(),true);
    ARRAY<T,FACE_INDEX<TV::dimension> > local_face_velocities(local_mac_grid);
    FAST_PROJECTION_DYNAMICS_UNIFORM<GRID<TV> > local_projection(projection->coarse_scale);
    projection->Local_Projection_PCG(*fine_face_velocities,local_mac_grid,local_face_velocities,local_projection,dt,time,cell_index);
}
template class REFINEMENT_TASK<VECTOR<float,1> >;
template class REFINEMENT_TASK<VECTOR<float,2> >;
template class REFINEMENT_TASK<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class REFINEMENT_TASK<VECTOR<double,1> >;
template class REFINEMENT_TASK<VECTOR<double,2> >;
template class REFINEMENT_TASK<VECTOR<double,3> >;
#endif

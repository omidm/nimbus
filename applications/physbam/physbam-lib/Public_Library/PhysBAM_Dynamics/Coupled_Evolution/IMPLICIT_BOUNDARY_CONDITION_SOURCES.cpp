//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION_SOURCES.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_CALLBACKS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> IMPLICIT_BOUNDARY_CONDITION_SOURCES<TV>::
IMPLICIT_BOUNDARY_CONDITION_SOURCES(FLUIDS_PARAMETERS_CALLBACKS<T_GRID>& callbacks_input)
    :callbacks(callbacks_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> IMPLICIT_BOUNDARY_CONDITION_SOURCES<TV>::
~IMPLICIT_BOUNDARY_CONDITION_SOURCES()
{
}
//#####################################################################
// Function Update_Boundary_Conditions
//#####################################################################
template<class TV> void IMPLICIT_BOUNDARY_CONDITION_SOURCES<TV>::
Update_Boundary_Conditions(const GRID<TV>& grid,ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N,ARRAY<T,TV_INT>& p,
    ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,const T time)
{
    callbacks.Get_Source_Velocities(face_velocities,psi_N,time);
}
template class IMPLICIT_BOUNDARY_CONDITION_SOURCES<VECTOR<float,1> >;
template class IMPLICIT_BOUNDARY_CONDITION_SOURCES<VECTOR<float,2> >;
template class IMPLICIT_BOUNDARY_CONDITION_SOURCES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class IMPLICIT_BOUNDARY_CONDITION_SOURCES<VECTOR<double,1> >;
template class IMPLICIT_BOUNDARY_CONDITION_SOURCES<VECTOR<double,2> >;
template class IMPLICIT_BOUNDARY_CONDITION_SOURCES<VECTOR<double,3> >;
#endif

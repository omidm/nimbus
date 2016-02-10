//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
using namespace PhysBAM;
//#####################################################################
// Destructor
//#####################################################################
template<class TV> EXAMPLE_FORCES_AND_VELOCITIES<TV>::
~EXAMPLE_FORCES_AND_VELOCITIES()
{}
//#####################################################################
// Function Advance_One_Time_Step_Begin_Callback
//#####################################################################
template<class TV> void EXAMPLE_FORCES_AND_VELOCITIES<TV>::
Advance_One_Time_Step_Begin_Callback(const T dt,const T time)
{}
//#####################################################################
// Function Advance_One_Time_Step_End_Callback
//#####################################################################
template<class TV> void EXAMPLE_FORCES_AND_VELOCITIES<TV>::
Advance_One_Time_Step_End_Callback(const T dt,const T time)
{}
//#####################################################################
// Function Set_PD_Targets
//#####################################################################
template<class TV> void EXAMPLE_FORCES_AND_VELOCITIES<TV>::
Set_PD_Targets(const T dt,const T time)
{}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
template<class TV> void EXAMPLE_FORCES_AND_VELOCITIES<TV>::
Update_Time_Varying_Material_Properties(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Add_External_Impulses_Before
//#####################################################################
template<class TV> void EXAMPLE_FORCES_AND_VELOCITIES<TV>::
Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Add_External_Impulses
//#####################################################################
template<class TV> void EXAMPLE_FORCES_AND_VELOCITIES<TV>::
Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Add_External_Impulse
//#####################################################################
template<class TV> void EXAMPLE_FORCES_AND_VELOCITIES<TV>::
Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}   
//#####################################################################
template class EXAMPLE_FORCES_AND_VELOCITIES<VECTOR<float,1> >;
template class EXAMPLE_FORCES_AND_VELOCITIES<VECTOR<float,2> >;
template class EXAMPLE_FORCES_AND_VELOCITIES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EXAMPLE_FORCES_AND_VELOCITIES<VECTOR<double,1> >;
template class EXAMPLE_FORCES_AND_VELOCITIES<VECTOR<double,2> >;
template class EXAMPLE_FORCES_AND_VELOCITIES<VECTOR<double,3> >;
#endif

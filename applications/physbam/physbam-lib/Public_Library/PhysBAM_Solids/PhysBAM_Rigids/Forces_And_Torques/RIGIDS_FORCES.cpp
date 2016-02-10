//#####################################################################
// Copyright 2007-2008, Michael Lentine, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGIDS_FORCES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGIDS_FORCES<TV>::
RIGIDS_FORCES(RIGID_BODY_COLLECTION<TV>& rigid_body_collection)
    :rigid_body_collection(rigid_body_collection),cfl_number((T)1),allow_external_cfl_number(true),cfl_initialized(false),use_velocity_independent_forces(true),
    use_velocity_dependent_forces(true),use_force_differential(true),use_implicit_velocity_independent_forces(false),use_position_based_state(true),
    unique_id(Get_Unique_Id()),compute_half_forces(false)
{
    Use_Rest_State_For_Strain_Rate(false);
    Limit_Time_Step_By_Strain_Rate();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGIDS_FORCES<TV>::
~RIGIDS_FORCES()
{}
//#####################################################################
// Function Use_Rest_State_For_Strain_Rate
//#####################################################################
template<class TV> void RIGIDS_FORCES<TV>::
Use_Rest_State_For_Strain_Rate(const bool use_rest_state_for_strain_rate_input)
{
    use_rest_state_for_strain_rate=use_rest_state_for_strain_rate_input;
}
//#####################################################################
// Function Limit_Time_Step_By_Strain_Rate
//#####################################################################
template<class TV> void RIGIDS_FORCES<TV>::
Limit_Time_Step_By_Strain_Rate(const bool limit_time_step_by_strain_rate_input,const T max_strain_per_time_step_input)
{
    limit_time_step_by_strain_rate=limit_time_step_by_strain_rate_input;
    assert(max_strain_per_time_step_input>0);max_strain_per_time_step=max_strain_per_time_step_input;
}
template<class TV> int RIGIDS_FORCES<TV>::
Velocity_Dependent_Forces_Size() const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class TV> void RIGIDS_FORCES<TV>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<T> aggregate,const T time) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class TV> void RIGIDS_FORCES<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class TV> void RIGIDS_FORCES<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class TV> void RIGIDS_FORCES<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class TV> typename TV::SCALAR RIGIDS_FORCES<TV>::
Potential_Energy(const T time) const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return 0;
}
template<class TV> void RIGIDS_FORCES<TV>::
Save_Potential_Energy(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void RIGIDS_FORCES<TV>::
Compute_Energy_Error(ARRAY_VIEW<const TWIST<TV> > rigid_velocity_save,const T time,const T dt)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void RIGIDS_FORCES<TV>::
Add_Energy_Correction_Force(ARRAY_VIEW<const TWIST<TV> > rigid_velocity_save,const T time,const T dt)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void RIGIDS_FORCES<TV>::
Compute_Previously_Applied_Forces()
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void RIGIDS_FORCES<TV>::
Store_Velocities()
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template class RIGIDS_FORCES<VECTOR<T,d> >;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif

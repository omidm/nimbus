//#####################################################################
// Copyright 2007-2008, Michael Lentine, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLES_FORCES<TV>::
DEFORMABLES_FORCES(PARTICLES<TV>& particles)
    :particles(particles),cfl_number((T)1),allow_external_cfl_number(true),cfl_initialized(false),use_velocity_independent_forces(true),
    use_velocity_dependent_forces(true),use_force_differential(true),use_implicit_velocity_independent_forces(false),use_position_based_state(true),
    unique_id(Get_Unique_Id()),compute_half_forces(false)
{
    Use_Rest_State_For_Strain_Rate(false);
    Limit_Time_Step_By_Strain_Rate();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DEFORMABLES_FORCES<TV>::
~DEFORMABLES_FORCES()
{}
//#####################################################################
// Function Use_Rest_State_For_Strain_Rate
//#####################################################################
template<class TV> void DEFORMABLES_FORCES<TV>::
Use_Rest_State_For_Strain_Rate(const bool use_rest_state_for_strain_rate_input)
{
    use_rest_state_for_strain_rate=use_rest_state_for_strain_rate_input;
}
//#####################################################################
// Function Limit_Time_Step_By_Strain_Rate
//#####################################################################
template<class TV> void DEFORMABLES_FORCES<TV>::
Limit_Time_Step_By_Strain_Rate(const bool limit_time_step_by_strain_rate_input,const T max_strain_per_time_step_input)
{
    limit_time_step_by_strain_rate=limit_time_step_by_strain_rate_input;
    assert(max_strain_per_time_step_input>0);max_strain_per_time_step=max_strain_per_time_step_input;
}
template<class TV> int DEFORMABLES_FORCES<TV>::
Velocity_Dependent_Forces_Size() const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class TV> typename TV::SCALAR DEFORMABLES_FORCES<TV>::
Potential_Energy(const T time) const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return 0;
}
template<class TV> typename TV::SCALAR DEFORMABLES_FORCES<TV>::
Residual_Energy(const T time) const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return 0;
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name) const
{
};
template<class TV> void DEFORMABLES_FORCES<TV>::
Save_Potential_Energy(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Compute_Energy_Error(ARRAY_VIEW<const TV> velocity_save,const T time,const T dt)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Compute_Energy_Correction_Force(ARRAY_VIEW<const TV> velocity_save,const int max_particle_degree,const T time,const T dt)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Prepare_Energy_Correction_Force()
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Apply_Energy_Correction(const T time,const T dt)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Add_Connectivity(ARRAY<int>& particle_degree)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Compute_Previously_Applied_Forces()
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Setup_Set_Velocity_From_Positions(const T time,const bool is_position_update,const bool reset_alphas)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> int DEFORMABLES_FORCES<TV>::
Get_Element_Count()
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return 0;
}
template<class TV> FORCE_ELEMENTS* DEFORMABLES_FORCES<TV>::
Get_Force_Elements()
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return 0;
}
template<class TV> ARRAY<int>* DEFORMABLES_FORCES<TV>::
Incident_Nodes(const int force_element)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return 0;
}
template<class TV> TV DEFORMABLES_FORCES<TV>::
Get_Direction(const int force_element)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return TV();
}
template<class TV> typename TV::SCALAR DEFORMABLES_FORCES<TV>::
Get_Combined_One_Over_Mass(const int force_element)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return 0;
}
template<class TV> ARRAY<int>* DEFORMABLES_FORCES<TV>::
Incident_Force_Elements(const int particle)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return 0;
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Choose_Solution(const bool use_orig_force,const int force_element,const T dt,const T alpha1,const T alpha2,ARRAY<T>& v_n_hats)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> TV DEFORMABLES_FORCES<TV>::
Get_Force(const int force_element,const int particle,const bool use_original_force)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return TV();
}
template<class TV> typename TV::SCALAR DEFORMABLES_FORCES<TV>::
Get_Force(const int force_element)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return 0;
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Set_Force(const int force_element,const T force)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Get_Damping_Force(const int particle,TV& damping_force,const T dt,const bool use_coefficient)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Update_Residual_Energy(const int force_element,const T residual_energy,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> typename TV::SCALAR DEFORMABLES_FORCES<TV>::
Get_Residual_Energy(const int force_element)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return 0;
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Compute_Quadratic_Contribution_For_Force(T& A,T& a,T&c,const T dt,const int force_element,const T combined_one_over_mass,const bool ignore_PE_terms)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Compute_Quadratic_Contribution_For_Node(T& B,T& C,T&b,const T dt,const int node,const int force_element,const T combined_one_over_mass,const T v_n_correction,const bool ignore_PE_terms)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Compute_Quadratic_Contribution_For_Residual(T& B,T& C,T& a,T&b,T&c,const T dt,const T time,const int force_element,const bool ignore_PE_terms)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Store_Delta_PE(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> typename TV::SCALAR DEFORMABLES_FORCES<TV>::
Get_Total_Delta_PE()
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return 0;
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Save_And_Reset_Elastic_Coefficient()
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Restore_Elastic_Coefficient()
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
template<class TV> void DEFORMABLES_FORCES<TV>::
Store_Velocities()
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template class DEFORMABLES_FORCES<VECTOR<T,d> >;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif

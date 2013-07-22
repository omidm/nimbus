//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Polynomials/QUADRATIC.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLE_GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void DEFORMABLE_GRAVITY<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    if(!gravity) return;
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();F(p)+=gravity*particles.mass(p)*downward_direction;}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_GRAVITY<TV>::
Potential_Energy(int p,const T time) const
{
    TV acceleration=gravity*downward_direction;
    return -particles.mass(p)*TV::Dot_Product(particles.X(p),acceleration);
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_GRAVITY<TV>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
        potential_energy+=Potential_Energy(p,time);}
    return potential_energy;
}
//#####################################################################
// Function Store_Potential_Energy
//#####################################################################
template<class TV> void DEFORMABLE_GRAVITY<TV>::
Save_Potential_Energy(const T time)
{
    potential_energy_save.Resize(particles.array_collection->Size());
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
        potential_energy_save(p)=Potential_Energy(p,time);}
}
//#####################################################################
// Function Residual_Energy
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_GRAVITY<TV>::
Residual_Energy(const T time) const
{
    T residual_energy=0;
    if(residual_PE.m)
        for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
            residual_energy+=residual_PE(p);}
    return residual_energy;
}
//#####################################################################
// Function Setup_Set_Velocity_From_Positions
//#####################################################################
template<class TV> void DEFORMABLE_GRAVITY<TV>::
Setup_Set_Velocity_From_Positions(const T time,const bool is_position_update,const bool reset_alphas)
{
    if(!is_position_update) Update_Position_Based_State(time,is_position_update);

    force_estimates.Resize(particles.array_collection->Size());
    incident_nodes.Resize(particles.array_collection->Size());
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
        incident_nodes(p).Remove_All();
        incident_nodes(p).Append(p);}
    if(is_position_update){
        if(reset_alphas) for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int s=iterator.Data();force_estimates(s)=0;}
        Save_Potential_Energy(time);}
}
//#####################################################################
// Function Get_Element_Count
//#####################################################################
template<class TV> int DEFORMABLE_GRAVITY<TV>::
Get_Element_Count()
{
    return particles.array_collection->Size();
}
//#####################################################################
// Function Get_Force_Elements
//#####################################################################
template<class TV> FORCE_ELEMENTS* DEFORMABLE_GRAVITY<TV>::
Get_Force_Elements()
{
    return &force_particles;
}
//#####################################################################
// Function Incident_Nodes
//#####################################################################
template<class TV> ARRAY<int>* DEFORMABLE_GRAVITY<TV>::
Incident_Nodes(const int force_element)
{
    return &(incident_nodes(force_element));
}
//#####################################################################
// Function Get_Direction
//#####################################################################
template<class TV> TV DEFORMABLE_GRAVITY<TV>::
Get_Direction(const int force_element)
{
    return downward_direction;
}
//#####################################################################
// Function Get_Combined_One_Over_Mass
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_GRAVITY<TV>::
Get_Combined_One_Over_Mass(const int force_element)
{
    return particles.one_over_mass(force_element);
}
//#####################################################################
// Function Incident_Force_Elements
//#####################################################################
template<class TV> ARRAY<int>* DEFORMABLE_GRAVITY<TV>::
Incident_Force_Elements(const int particle)
{
    return &incident_elements(particle);
}
//#####################################################################
// Function Choose_Solution
//#####################################################################
template<class TV> void DEFORMABLE_GRAVITY<TV>::
Choose_Solution(const bool use_orig_force,const int force_element,const T dt,const T alpha1,const T alpha2,ARRAY<T>& v_n_hats)
{
    T original_force=gravity*particles.mass(force_element);

    if(fabs(alpha1-original_force)<fabs(alpha2-original_force)) Set_Force(force_element,alpha1);
    else Set_Force(force_element,alpha2);
}
//#####################################################################
// Function Get_Force
//#####################################################################
template<class TV> TV DEFORMABLE_GRAVITY<TV>::
Get_Force(const int force_element,const int particle,const bool use_original_force)
{
    if(use_original_force) return gravity*particles.mass(particle)*downward_direction;
    else return force_estimates(force_element)*downward_direction;
}
//#####################################################################
// Function Get_Force
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_GRAVITY<TV>::
Get_Force(const int force_element)
{
    return force_estimates(force_element);
}
//#####################################################################
// Function Set_Force
//#####################################################################
template<class TV> void DEFORMABLE_GRAVITY<TV>::
Set_Force(const int force_element,const T force)
{
    force_estimates(force_element)=force;
}
//#####################################################################
// Function Get_Damping_Force
//#####################################################################
template<class TV> void DEFORMABLE_GRAVITY<TV>::
Get_Damping_Force(const int particle,TV& damping_force,const T dt,const bool use_coefficient)
{
}
//#####################################################################
// Function Update_Residual_Energy
//#####################################################################
template<class TV> void DEFORMABLE_GRAVITY<TV>::
Update_Residual_Energy(const int force_element,const T residual_energy,const T time)
{
    residual_PE(force_element)=residual_energy;
}
//#####################################################################
// Function Get_Residual_Energy
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_GRAVITY<TV>::
Get_Residual_Energy(const int force_element)
{
    return residual_PE(force_element);
}
//#####################################################################
// Function Compute_Quadratic_Contribution_For_Node
//#####################################################################
template<class TV> void DEFORMABLE_GRAVITY<TV>::
Compute_Quadratic_Contribution_For_Force(T& A,T& a,T&c,const T dt,const int force_element,const T combined_one_over_mass,const bool ignore_PE_terms)
{
    A=(T).5*dt*dt*combined_one_over_mass;
}
//#####################################################################
// Function Compute_Quadratic_Contribution_For_Node
//#####################################################################
template<class TV> void DEFORMABLE_GRAVITY<TV>::
Compute_Quadratic_Contribution_For_Node(T& B,T& C,T&b,const T dt,const int node,const int force_element,const T combined_one_over_mass,const T v_n_correction,const bool ignore_PE_terms)
{
    T vu=TV::Dot_Product(particles.V(node),downward_direction);

    B+=dt*(vu+(T).5*v_n_correction);

    if(!ignore_PE_terms) C+=-(T)dt*(1/combined_one_over_mass)*gravity*(vu+(T).5*v_n_correction);
}
//#####################################################################
// Function Compute_Quadratic_Contribution_For_Residual
//#####################################################################
template<class TV> void DEFORMABLE_GRAVITY<TV>::
Compute_Quadratic_Contribution_For_Residual(T& B,T& C,T& a,T&b,T&c,const T dt,const T time,const int force_element,const bool ignore_PE_terms)
{
    if(!ignore_PE_terms)
        B+=-(T).5*dt*dt*gravity;
    else C=delta_PE(force_element);

    C+=residual_PE(force_element);
}
//#####################################################################
// Function Store_Delta_PE
//#####################################################################
template<class TV> void DEFORMABLE_GRAVITY<TV>::
Store_Delta_PE(const T time)
{
    delta_PE.Resize(potential_energy_save.m);
    T total_current_PE=0;
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
        delta_PE(p)=Potential_Energy(p,time)-potential_energy_save(p);
        total_current_PE+=Potential_Energy(p,time);}
    total_delta_PE=total_current_PE-ARRAYS_COMPUTATIONS::Sum(potential_energy_save);
}
//#####################################################################
// Function Get_Total_Delta_PE
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_GRAVITY<TV>::
Get_Total_Delta_PE()
{
    return total_delta_PE;
}
//#####################################################################
// Function Save_And_Reset_Elastic_Coefficient
//#####################################################################
template<class TV> void DEFORMABLE_GRAVITY<TV>::
Save_And_Reset_Elastic_Coefficient()
{
}
//#####################################################################
// Function Restore_Elastic_Coefficient
//#####################################################################
template<class TV> void DEFORMABLE_GRAVITY<TV>::
Restore_Elastic_Coefficient()
{
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void DEFORMABLE_GRAVITY<TV>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    if(!residual_PE.m) residual_PE.Resize(particles.array_collection->Size());

    if(!incident_elements.m){
        incident_elements.Resize(particles.array_collection->Size());
        for(int i=1;i<=force_particles.indices.m;i++)
            incident_elements(force_particles.indices(i)).Append(force_particles.indices(i));}
}
//#####################################################################
template class DEFORMABLE_GRAVITY<VECTOR<float,1> >;
template class DEFORMABLE_GRAVITY<VECTOR<float,2> >;
template class DEFORMABLE_GRAVITY<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DEFORMABLE_GRAVITY<VECTOR<double,1> >;
template class DEFORMABLE_GRAVITY<VECTOR<double,2> >;
template class DEFORMABLE_GRAVITY<VECTOR<double,3> >;
#endif

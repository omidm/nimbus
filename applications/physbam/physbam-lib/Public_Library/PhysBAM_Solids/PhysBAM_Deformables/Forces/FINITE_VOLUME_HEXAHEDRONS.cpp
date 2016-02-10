//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Geometry/Basic_Geometry/HEXAHEDRON.h>
#include <PhysBAM_Geometry/Topology/HEXAHEDRON_MESH.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ANISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/DIAGONALIZED_STRESS_DERIVATIVE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/PLASTICITY_MODEL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/STRAIN_MEASURE_HEXAHEDRONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME_HEXAHEDRONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> FINITE_VOLUME_HEXAHEDRONS<T>::
FINITE_VOLUME_HEXAHEDRONS(STRAIN_MEASURE_HEXAHEDRONS<T>& strain_measure,CONSTITUTIVE_MODEL<T,3>& constitutive_model)
    :DEFORMABLES_FORCES<TV>(strain_measure.particles),strain_measure(strain_measure),constitutive_model(constitutive_model),
    dPi_dFe(0),dP_dFe(0),Be_scales_save(0),V(0),node_stiffness(0),edge_stiffness(0),use_uniform_density(false),destroy_data(false),density_list(0)
{
    Update_Be_Scales();
    isotropic_model=dynamic_cast<ISOTROPIC_CONSTITUTIVE_MODEL<T,3>*>(&constitutive_model);
    anisotropic_model=dynamic_cast<ANISOTROPIC_CONSTITUTIVE_MODEL<T,3>*>(&constitutive_model);
    if(anisotropic_model) Save_V();
    strain_measure.mesh.Initialize_Incident_Elements();
    if(use_uniform_density){
        T total_mass=0;
        ARRAY<int> mesh_particles;strain_measure.mesh_object.mesh.elements.Flattened().Get_Unique(mesh_particles);
        for(int i=1;i<=mesh_particles.m;i++) total_mass+=particles.mass(mesh_particles(i));
        density=total_mass/strain_measure.mesh_object.Total_Volume();
        if(density==0) density=TV::dimension==1?1:TV::dimension==2?100:1000;}
    else{
        density_list=new ARRAY<T>(strain_measure.mesh_object.mesh.elements.m);
        for(int i=1;i<=strain_measure.mesh.elements.m;i++){
            const VECTOR<int,8>& nodes=strain_measure.mesh.elements(i);
            T volume=HEXAHEDRON<T>::Signed_Volume(particles.X(nodes(1)),particles.X(nodes(2)),particles.X(nodes(3)),particles.X(nodes(4)),particles.X(nodes(5)),particles.X(nodes(6)),particles.X(nodes(7)),particles.X(nodes(8)));
            for(int j=1;j<=nodes.m;j++) (*density_list)(i)+=particles.mass(nodes(j))/(*strain_measure.mesh.incident_elements)(nodes(j)).m/volume;
            if((*density_list)(i)==0) (*density_list)(i)=TV::dimension==1?1:TV::dimension==2?100:1000;}}
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> FINITE_VOLUME_HEXAHEDRONS<T>::
~FINITE_VOLUME_HEXAHEDRONS()
{
    delete dPi_dFe;delete dP_dFe;delete Be_scales_save;delete V;delete node_stiffness;delete edge_stiffness;
    if(destroy_data){delete &strain_measure;delete &constitutive_model;}
}
//#####################################################################
// Function Save_Stress_Derivative
//#####################################################################
template<class T> void FINITE_VOLUME_HEXAHEDRONS<T>::
Save_Stress_Derivative()
{
    if(isotropic_model || anisotropic_model->use_isotropic_component_of_stress_derivative_only){
        if(!dPi_dFe) dPi_dFe=new ARRAY<VECTOR<DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>,8> >();delete dP_dFe;dP_dFe=0;}
    else {if(!dP_dFe) dP_dFe=new ARRAY<VECTOR<DIAGONALIZED_STRESS_DERIVATIVE<T,3>,8> >();delete dPi_dFe;dPi_dFe=0;}
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class T> void FINITE_VOLUME_HEXAHEDRONS<T>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    strain_measure.mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class T> void FINITE_VOLUME_HEXAHEDRONS<T>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    force_elements.Update(strain_measure.mesh.elements,particle_is_simulated);
    force_particles.Update(strain_measure.mesh.elements.Flattened(),particle_is_simulated);
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T> void FINITE_VOLUME_HEXAHEDRONS<T>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    int elements=strain_measure.H_DmH_inverse.m;
    U.Resize(elements,false,false);De_inverse_hat.Resize(elements,false,false);Fe_hat.Resize(elements,false,false);
    if(dPi_dFe) dPi_dFe->Resize(elements,false,false);if(dP_dFe) dP_dFe->Resize(elements,false,false);if(V) V->Resize(elements,false,false);
    MATRIX<T,3> V_local;
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int e=iterator.Data();for(int g=1;g<=8;g++){
        int gauss_index=8*(e-1)+g;
        strain_measure.F(g,e).Fast_Singular_Value_Decomposition(U(e)(g),Fe_hat(e)(g),V_local);
        if(anisotropic_model) anisotropic_model->Update_State_Dependent_Auxiliary_Variables(Fe_hat(e)(g),V_local,gauss_index);
        else isotropic_model->Update_State_Dependent_Auxiliary_Variables(Fe_hat(e)(g),gauss_index);
        if(dPi_dFe) constitutive_model.Isotropic_Stress_Derivative(Fe_hat(e)(g),(*dPi_dFe)(e)(g),gauss_index);
        if(dP_dFe) anisotropic_model->Stress_Derivative(Fe_hat(e)(g),V_local,(*dP_dFe)(e)(g),gauss_index);
        De_inverse_hat(e)(g).Resize(8);for(int k=1;k<=8;k++) De_inverse_hat(e)(g)(k)=strain_measure.H_DmH_inverse(e)(g)(k)*V_local;
        if(V) (*V)(e)(g)=V_local;}}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T> void FINITE_VOLUME_HEXAHEDRONS<T>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    if(anisotropic_model)
        for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int e=iterator.Data();for(int g=1;g<=8;g++){
            int gauss_index=8*(e-1)+g;
            MATRIX<T,3> U_P_hat=U(e)(g)*anisotropic_model->P_From_Strain(Fe_hat(e)(g),(*V)(e)(g),Be_scales(e)(g),gauss_index);
            for(int k=1;k<=8;k++) F(strain_measure.mesh.elements(e)(k))+=U_P_hat*De_inverse_hat(e)(g)(k);}}
    else
        for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int e=iterator.Data();for(int g=1;g<=8;g++){
            MATRIX<T,3> U_P_hat=U(e)(g)*isotropic_model->P_From_Strain(Fe_hat(e)(g),Be_scales(e)(g),e);
            for(int k=1;k<=8;k++) F(strain_measure.mesh.elements(e)(k))+=U_P_hat*De_inverse_hat(e)(g)(k);}}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class T> void FINITE_VOLUME_HEXAHEDRONS<T>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int e=iterator.Data();for(int g=1;g<=8;g++){
        MATRIX<T,3> Fe_dot_hat=U(e)(g).Transpose_Times(strain_measure.Gradient(V,De_inverse_hat,g,e));
        MATRIX<T,3> U_P_hat=U(e)(g)*constitutive_model.P_From_Strain_Rate(Fe_hat(e)(g),Fe_dot_hat,Be_scales(e)(g),e);
        for(int k=1;k<=8;k++) F(strain_measure.mesh.elements(e)(k))+=U_P_hat*De_inverse_hat(e)(g)(k);}}
}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class T> void FINITE_VOLUME_HEXAHEDRONS<T>::
Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const
{
    if(!dPi_dFe && !dP_dFe) PHYSBAM_FATAL_ERROR();
    if(anisotropic_model)
        for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int e=iterator.Data();for(int g=1;g<=8;g++){
            int gauss_index=8*(e-1)+g;
            MATRIX<T,3> dF_hat=U(e)(g).Transpose_Times(strain_measure.Gradient(dX,De_inverse_hat,g,e));
            MATRIX<T,3> U_dP_hat;
            if(dP_dFe) U_dP_hat=U(e)(g)*anisotropic_model->dP_From_dF(dF_hat,Fe_hat(e)(g),(*V)(e)(g),(*dP_dFe)(e)(g),Be_scales(e)(g),gauss_index);
            else U_dP_hat=U(e)(g)*anisotropic_model->dP_From_dF(dF_hat,Fe_hat(e)(g),(*V)(e)(g),(*dPi_dFe)(e)(g),Be_scales(e)(g),gauss_index);
            for(int k=1;k<=8;k++) dF(strain_measure.mesh.elements(e)(k))+=U_dP_hat*De_inverse_hat(e)(g)(k);}}
    else
        for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int e=iterator.Data();for(int g=1;g<=8;g++){
            int gauss_index=8*(e-1)+g;
            MATRIX<T,3> dF_hat=U(e)(g).Transpose_Times(strain_measure.Gradient(dX,De_inverse_hat,g,e));
            MATRIX<T,3> U_dP_hat=U(e)(g)*isotropic_model->dP_From_dF(dF_hat,(*dPi_dFe)(e)(g),Be_scales(e)(g),gauss_index);
            for(int k=1;k<=8;k++) dF(strain_measure.mesh.elements(e)(k))+=U_dP_hat*De_inverse_hat(e)(g)(k);}}
}
//#####################################################################
// Function Intialize_CFL
//#####################################################################
template<class T> void FINITE_VOLUME_HEXAHEDRONS<T>::
Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)
{
    // TODO: MPI
    T one_over_cfl_number=1/cfl_number,one_over_cfl_number_squared=sqr(one_over_cfl_number);
    T minimum_altitude_squared=sqr(strain_measure.DmH_minimum_altitude); // TODO: use fragment-specific minimum altitude

    ARRAY<FREQUENCY_DATA> fragment_particle_frequency(frequency.Size());
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int e=iterator.Data();
        const VECTOR<int,8>& nodes=strain_measure.mesh.elements(e);
        T one_over_altitude_squared_and_density=(T)1/(minimum_altitude_squared*use_uniform_density?density:(*density_list)(e));
        T elastic_squared=constitutive_model.Maximum_Elastic_Stiffness(e)*one_over_altitude_squared_and_density*one_over_cfl_number_squared;
        T damping=constitutive_model.Maximum_Damping_Stiffness(e)*one_over_altitude_squared_and_density*one_over_cfl_number;
        for(int j=1;j<=nodes.m;j++){FREQUENCY_DATA& data=fragment_particle_frequency(nodes[j]);
            data.elastic_squared=max(data.elastic_squared,elastic_squared);data.damping=max(data.damping,damping);}}
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
        frequency(p).elastic_squared+=fragment_particle_frequency(p).elastic_squared;
        frequency(p).damping+=fragment_particle_frequency(p).damping;}
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class T> T FINITE_VOLUME_HEXAHEDRONS<T>::
CFL_Strain_Rate() const
{
    T max_strain_rate=0;
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int e=iterator.Data();
        for(int g=1;g<=8;g++){max_strain_rate=max(max_strain_rate,strain_measure.Velocity_Gradient(g,e).Max_Abs());}}
    return Robust_Divide(max_strain_per_time_step,max_strain_rate);
}
//#####################################################################
// Function Create_Finite_Volume
//#####################################################################
template<class T> FINITE_VOLUME_HEXAHEDRONS<T>*
Create_Finite_Volume(HEXAHEDRALIZED_VOLUME<T>& hexahedralized_volume,CONSTITUTIVE_MODEL<T,3>* constitutive_model,const bool limit_time_step_by_strain_rate,
    const T max_strain_per_time_step,const bool use_rest_state_for_strain_rate)
{
    STRAIN_MEASURE_HEXAHEDRONS<T>* strain_measure=new STRAIN_MEASURE_HEXAHEDRONS<T>(hexahedralized_volume);
    FINITE_VOLUME_HEXAHEDRONS<T>* fvm=new FINITE_VOLUME_HEXAHEDRONS<T>(*strain_measure,*constitutive_model);
    fvm->Limit_Time_Step_By_Strain_Rate(limit_time_step_by_strain_rate,max_strain_per_time_step);
    fvm->Use_Rest_State_For_Strain_Rate(use_rest_state_for_strain_rate);
    fvm->destroy_data=true;
    return fvm;
}
//#####################################################################
// Function Create_Quasistatic_Finite_Volume
//#####################################################################
template<class T> FINITE_VOLUME_HEXAHEDRONS<T>*
Create_Quasistatic_Finite_Volume(HEXAHEDRALIZED_VOLUME<T>& hexahedralized_volume,CONSTITUTIVE_MODEL<T,3>* constitutive_model)
{
    FINITE_VOLUME_HEXAHEDRONS<T>* fvm=Create_Finite_Volume(hexahedralized_volume,constitutive_model);
    fvm->Use_Quasistatics();return fvm;
}
//#####################################################################
template class FINITE_VOLUME_HEXAHEDRONS<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FINITE_VOLUME_HEXAHEDRONS<double>;
#endif

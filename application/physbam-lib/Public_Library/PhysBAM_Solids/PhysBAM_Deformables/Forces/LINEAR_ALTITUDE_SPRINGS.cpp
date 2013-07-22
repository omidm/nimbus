//#####################################################################
// Copyright 2004-2007, Ronald Fedkiw, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_ALTITUDE_SPRINGS
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/sign.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TETRAHEDRON_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_2D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_S3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> LINEAR_ALTITUDE_SPRINGS<TV,d>::
LINEAR_ALTITUDE_SPRINGS(PARTICLES<TV>& particles,T_MESH& mesh_input)
    :DEFORMABLES_FORCES<TV>(particles),mesh(mesh_input),parameters(mesh.elements.m),use_plasticity(false),cache_strain(false)
{
    Print_Number_Used();
    Use_Springs_Compressed_Beyond_Threshold_Only(false);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d> LINEAR_ALTITUDE_SPRINGS<TV,d>::
~LINEAR_ALTITUDE_SPRINGS()
{
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV,int d> void LINEAR_ALTITUDE_SPRINGS<TV,d>::
Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)
{
    T one_over_d_squared=(T)1/sqr(d);
    T one_over_cfl_number=1/cfl_number,one_over_cfl_number_squared=sqr(one_over_cfl_number);
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        const VECTOR<int,d+1>& nodes=mesh.elements(t);
        for(int i=1;i<=nodes.m;i++) for(int j=1;j<=nodes.m;j++){
            const SPRING_PARAMETER& parameter=parameters(t)(j);
            T one_over_mass_times_restlength=particles.one_over_effective_mass(nodes[i])/parameter.restlength;
            T c=(i==j)?(T)1:one_over_d_squared;
            frequency(nodes[i]).elastic_squared+=4*c*parameter.youngs_modulus*one_over_mass_times_restlength*one_over_cfl_number_squared;
            frequency(nodes[i]).damping+=2*c*parameter.damping*one_over_mass_times_restlength*one_over_cfl_number;}}
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV,int d> void LINEAR_ALTITUDE_SPRINGS<TV,d>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV,int d> void LINEAR_ALTITUDE_SPRINGS<TV,d>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    force_elements.Update(mesh.elements,particle_is_simulated);
    if(cache_strain) strains_of_spring.Resize(mesh.elements.m,false,false);
}
//#####################################################################
// Function Compute_Plasticity
//#####################################################################
template<class TV,int d> void LINEAR_ALTITUDE_SPRINGS<TV,d>::
Compute_Plasticity(const int h,const int t,const T current_length)
{
    SPRING_PARAMETER& parameter=parameters(t)(h);
    PLASTIC_PARAMETER& plastic_parameter=plastic_parameters(t)(h);
    T strain=(current_length-parameter.visual_restlength)/parameter.restlength;T strain_magnitude=abs(strain);
    if(strain_magnitude>plastic_parameter.yield_strain){
        parameter.visual_restlength=clamp(current_length-sign(strain)*plastic_parameter.yield_strain*parameter.restlength,parameter.visual_restlength/plasticity_clamp_ratio,
            parameter.visual_restlength*plasticity_clamp_ratio);
        plastic_parameter.yield_strain+=plastic_parameter.hardening*(strain_magnitude-plastic_parameter.yield_strain);}
}
//#####################################################################
// Function Clamp_Restlength_With_Fraction_Of_Springs
//#####################################################################
template<class TV,int d> void LINEAR_ALTITUDE_SPRINGS<TV,d>::
Clamp_Restlength_With_Fraction_Of_Springs(const T fraction)
{
    ARRAY<T> length((d+1)*parameters.m,false);for(int s=1;s<=parameters.m;s++) for(int k=1;k<=d+1;k++) length((d+1)*(s-1)+k)=parameters(s)(k).restlength;Sort(length);
    T minimum_restlength=length(min((int)(fraction*length.m)+1,length.m));std::stringstream ss;ss<<"Enlarging the restlength of all altitude springs below "<<minimum_restlength<<std::endl;LOG::filecout(ss.str());
    for(int i=1;i<=parameters.m;i++) for(int k=1;k<=d+1;k++) parameters(i)(k).restlength=max(minimum_restlength,parameters(i)(k).restlength);
}
//#####################################################################
// Function Print_Restlength_Statistics
//#####################################################################
template<class TV,int d> void LINEAR_ALTITUDE_SPRINGS<TV,d>::
Print_Restlength_Statistics() const
{
    std::stringstream ss;ss<<"Altitude Springs - Total Springs = "<<(d+1)*parameters.m<<std::endl;
    ARRAY<T> length((d+1)*parameters.m,false);for(int s=1;s<=parameters.m;s++) for(int k=1;k<=d+1;k++) length((d+1)*s+k-(d+1))=parameters(s)(k).restlength;Sort(length);
    ARRAY<T> visual_length((d+1)*parameters.m,false);for(int s=1;s<=parameters.m;s++) for(int k=1;k<=d+1;k++)visual_length((d+1)*s+k-(d+1))=parameters(s)(k).visual_restlength;
    Sort(visual_length);
    int one_percent=(int)(.01*length.m)+1,ten_percent=(int)(.1*length.m)+1,median=(int)(.5*length.m)+1;
    ss<<"Altitude Springs - Smallest Restlength = "<<length(1)<<", Visual Restlength = "<<visual_length(1)<<std::endl;
    ss<<"Altitude Springs - One Percent Restlength = "<<length(one_percent)<<", Visual Restlength = "<<visual_length(one_percent)<<std::endl;
    ss<<"Altitude Springs - Ten Percent Restlength = "<<length(ten_percent)<<", Visual Restlength = "<<visual_length(ten_percent)<<std::endl;
    ss<<"Altitude Springs - Median Restlength = "<<length(median)<<", Visual Restlength = "<<visual_length(median)<<std::endl;LOG::filecout(ss.str());
}
//#####################################################################
// Function Enable_Plasticity
//#####################################################################
template<class TV,int d> void LINEAR_ALTITUDE_SPRINGS<TV,d>::
Enable_Plasticity(const ARRAY<VECTOR<T,d+1> >& plastic_yield_strain_input,const ARRAY<VECTOR<T,d+1> >& plastic_hardening_input,const T plasticity_clamp_ratio_input)
{
    use_plasticity=true;plasticity_clamp_ratio=plasticity_clamp_ratio_input;
    plastic_parameters.Resize(parameters.m,false,false);
    for(int i=1;i<=parameters.m;i++) for(int k=1;k<=d+1;k++){PLASTIC_PARAMETER& plastic_parameter=plastic_parameters(i)(k);
        plastic_parameter.visual_restlength=parameters(i)(k).visual_restlength;
        plastic_parameter.yield_strain=plastic_yield_strain_input(i)(k);
        plastic_parameter.hardening=plastic_hardening_input(i)(k);}
}
//#####################################################################
// Function Enable_Plasticity
//#####################################################################
template<class TV,int d> void LINEAR_ALTITUDE_SPRINGS<TV,d>::
Enable_Plasticity(const T plastic_yield_strain_input,const T plastic_hardening_input,const T plasticity_clamp_ratio_input)
{
    use_plasticity=true;plasticity_clamp_ratio=plasticity_clamp_ratio_input;
    plastic_parameters.Resize(parameters.m,false,false);
    for(int i=1;i<=parameters.m;i++) for(int k=1;k<=d+1;k++){PLASTIC_PARAMETER& plastic_parameter=plastic_parameters(i)(k);
        plastic_parameter.visual_restlength=parameters(i)(k).visual_restlength;
        plastic_parameter.yield_strain=plastic_yield_strain_input;
        plastic_parameter.hardening=plastic_hardening_input;}
}
//#####################################################################
// Function Set_Stiffness
//#####################################################################
template<class TV,int d> void LINEAR_ALTITUDE_SPRINGS<TV,d>::
Set_Stiffness(const T youngs_modulus_input)
{
    for(int i=1;i<=parameters.m;i++) for(int k=1;k<=d+1;k++) parameters(i)(k).youngs_modulus=youngs_modulus_input;Invalidate_CFL();
}
//#####################################################################
// Function Set_Stiffness
//#####################################################################
template<class TV,int d> void LINEAR_ALTITUDE_SPRINGS<TV,d>::
Set_Stiffness(const ARRAY<VECTOR<T,d+1> >& youngs_modulus_input)
{
    for(int i=1;i<=parameters.m;i++) for(int k=1;k<=d+1;k++) parameters(i)(k).youngs_modulus=youngs_modulus_input(i)(k);Invalidate_CFL();
}
//#####################################################################
// Function Set_Restlength
//#####################################################################
template<class TV,int d> void LINEAR_ALTITUDE_SPRINGS<TV,d>::
Set_Restlength(const ARRAY<VECTOR<T,d+1> >& restlength_input)
{
    for(int i=1;i<=parameters.m;i++) for(int k=1;k<=d+1;k++) parameters(i)(k).restlength=restlength_input(i)(k);
}
//#####################################################################
// Function Clamp_Restlength
//#####################################################################
template<class TV,int d> void LINEAR_ALTITUDE_SPRINGS<TV,d>::
Clamp_Restlength(const T clamped_restlength)
{
    Invalidate_CFL();
    for(int i=1;i<=parameters.m;i++) for(int k=1;k<=d+1;k++) parameters(i)(k).restlength=max(parameters(i)(k).visual_restlength,clamped_restlength);
}
//#####################################################################
// Function Set_Damping
//#####################################################################
template<class TV,int d> void LINEAR_ALTITUDE_SPRINGS<TV,d>::
Set_Damping(const ARRAY<VECTOR<T,d+1> >& damping_input)
{
    for(int i=1;i<=parameters.m;i++) for(int k=1;k<=d+1;k++) parameters(i)(k).damping=damping_input(i)(k);Invalidate_CFL();
}
//#####################################################################
// Function Set_Damping
//#####################################################################
template<class TV,int d> void LINEAR_ALTITUDE_SPRINGS<TV,d>::
Set_Damping(const T damping_input)
{
    for(int i=1;i<=parameters.m;i++) for(int k=1;k<=d+1;k++) parameters(i)(k).damping=damping_input;Invalidate_CFL();
}
//#####################################################################
// Function Create_Altitude_Springs_Base
//#####################################################################
template<class T,class TV,class T_MESH> typename SOLIDS_FORCES_POLICY<TV,T_MESH::dimension>::LINEAR_ALTITUDE_SPRINGS* PhysBAM::
Create_Altitude_Springs_Base(PARTICLES<TV>& particles,T_MESH& mesh,const T stiffness,const T overdamping_fraction,
    const bool use_compressed_by_threshold_only,const T fraction_compression,const bool limit_time_step_by_strain_rate,const T max_strain_per_time_step,const bool use_rest_state_for_strain_rate,
    const T restlength_enlargement_fraction,const bool verbose)
{
    typedef typename SOLIDS_FORCES_POLICY<TV,T_MESH::dimension>::LINEAR_ALTITUDE_SPRINGS T_LINEAR_ALTITUDE_SPRINGS;
    T_LINEAR_ALTITUDE_SPRINGS* las=new T_LINEAR_ALTITUDE_SPRINGS(particles,mesh);
    las->Set_Restlength_From_Particles();
    if(restlength_enlargement_fraction) las->Clamp_Restlength_With_Fraction_Of_Springs(restlength_enlargement_fraction);
    las->Set_Stiffness(stiffness);
    las->Set_Overdamping_Fraction(overdamping_fraction);
    las->Use_Springs_Compressed_Beyond_Threshold_Only(use_compressed_by_threshold_only,fraction_compression);
    las->Limit_Time_Step_By_Strain_Rate(limit_time_step_by_strain_rate,max_strain_per_time_step);
    las->Use_Rest_State_For_Strain_Rate(use_rest_state_for_strain_rate);
    if(verbose) las->Print_Restlength_Statistics();
    return las;
}
//#####################################################################
template class LINEAR_ALTITUDE_SPRINGS<VECTOR<float,2>,2>;
template class LINEAR_ALTITUDE_SPRINGS<VECTOR<float,3>,2>;
template class LINEAR_ALTITUDE_SPRINGS<VECTOR<float,3>,3>;
template SOLIDS_FORCES_POLICY<VECTOR<float,3>,TETRAHEDRON_MESH::dimension>::LINEAR_ALTITUDE_SPRINGS* PhysBAM::Create_Altitude_Springs_Base<float,VECTOR<float,3>,
    TETRAHEDRON_MESH>(PARTICLES<VECTOR<float,3> >&,TETRAHEDRON_MESH&,float,float,bool,float,bool,float,bool,float,bool);
template SOLIDS_FORCES_POLICY<VECTOR<float,3>,TRIANGLE_MESH::dimension>::LINEAR_ALTITUDE_SPRINGS* PhysBAM::Create_Altitude_Springs_Base<float,VECTOR<float,3>,
    TRIANGLE_MESH>(PARTICLES<VECTOR<float,3> >&,TRIANGLE_MESH&,float,float,bool,float,bool,float,bool,float,bool);
template SOLIDS_FORCES_POLICY<VECTOR<float,2>,TRIANGLE_MESH::dimension>::LINEAR_ALTITUDE_SPRINGS* PhysBAM::Create_Altitude_Springs_Base<float,VECTOR<float,2>,
    TRIANGLE_MESH>(PARTICLES<VECTOR<float,2> >&,TRIANGLE_MESH&,float,float,bool,float,bool,float,bool,float,bool);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LINEAR_ALTITUDE_SPRINGS<VECTOR<double,2>,2>;
template class LINEAR_ALTITUDE_SPRINGS<VECTOR<double,3>,2>;
template class LINEAR_ALTITUDE_SPRINGS<VECTOR<double,3>,3>;
template SOLIDS_FORCES_POLICY<VECTOR<double,3>,TETRAHEDRON_MESH::dimension>::LINEAR_ALTITUDE_SPRINGS* PhysBAM::Create_Altitude_Springs_Base<double,VECTOR<double,3>,
    TETRAHEDRON_MESH>(PARTICLES<VECTOR<double,3> >&,TETRAHEDRON_MESH&,double,double,bool,double,bool,double,bool,double,bool);
template SOLIDS_FORCES_POLICY<VECTOR<double,3>,TRIANGLE_MESH::dimension>::LINEAR_ALTITUDE_SPRINGS* PhysBAM::Create_Altitude_Springs_Base<double,VECTOR<double,3>,
    TRIANGLE_MESH>(PARTICLES<VECTOR<double,3> >&,TRIANGLE_MESH&,double,double,bool,double,bool,double,bool,double,bool);
template SOLIDS_FORCES_POLICY<VECTOR<double,2>,TRIANGLE_MESH::dimension>::LINEAR_ALTITUDE_SPRINGS* PhysBAM::Create_Altitude_Springs_Base<double,VECTOR<double,2>,
    TRIANGLE_MESH>(PARTICLES<VECTOR<double,2> >&,TRIANGLE_MESH&,double,double,bool,double,bool,double,bool,double,bool);
#endif

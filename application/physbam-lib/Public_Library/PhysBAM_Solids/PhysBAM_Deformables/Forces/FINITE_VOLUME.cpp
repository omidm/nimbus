//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/UPPER_TRIANGULAR_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/UPPER_TRIANGULAR_MATRIX_3X3.h>
#include <PhysBAM_Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ANISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/DIAGONALIZED_STRESS_DERIVATIVE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/PLASTICITY_MODEL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DIAGONALIZED_SEMI_IMPLICIT_ELEMENT.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> FINITE_VOLUME<TV,d>::
FINITE_VOLUME(const bool use_uniform_density_input,STRAIN_MEASURE<TV,d>& strain_measure,CONSTITUTIVE_MODEL<T,d>& constitutive_model,PLASTICITY_MODEL<T,d>* plasticity_model)
    :DEFORMABLES_FORCES<TV>(dynamic_cast<PARTICLES<TV>&>(strain_measure.particles)),strain_measure(strain_measure),constitutive_model(constitutive_model),
    plasticity_model(plasticity_model),dPi_dFe(0),dP_dFe(0),Be_scales_save(0),V(0),twice_max_strain_per_time_step(0),semi_implicit_data(0),use_uniform_density(use_uniform_density_input),
    density_list(0),implicit_surface(0),node_stiffness(0),edge_stiffness(0),force_segments(0),destroy_data(false),half_force_size(0),total_half_force_size(0)
{
    Update_Be_Scales();
    isotropic_model=dynamic_cast<ISOTROPIC_CONSTITUTIVE_MODEL<T,d>*>(&constitutive_model);
    anisotropic_model=dynamic_cast<ANISOTROPIC_CONSTITUTIVE_MODEL<T,d>*>(&constitutive_model);
    if(anisotropic_model) Save_V();
    strain_measure.mesh.Initialize_Incident_Elements();
    if(use_uniform_density){
        T total_mass=0;
        ARRAY<int> mesh_particles;strain_measure.mesh_object.mesh.elements.Flattened().Get_Unique(mesh_particles);
        for(int i=1;i<=mesh_particles.m;i++) total_mass+=particles.mass(mesh_particles(i));
        density=total_mass/strain_measure.mesh_object.Total_Size();
        if(density==0) density=TV::dimension==1?1:TV::dimension==2?100:1000;}
    else{
        density_list=new ARRAY<T>(strain_measure.mesh_object.mesh.elements.m);
        for(int i=1;i<=strain_measure.mesh.elements.m;i++){
            const VECTOR<int,d+1>& nodes=strain_measure.mesh.elements(i);
            T volume=strain_measure.mesh_object.Signed_Size(i);
            for(int j=1;j<=nodes.m;j++) (*density_list)(i)+=particles.mass(nodes(j))/(*strain_measure.mesh.incident_elements)(nodes(j)).m/volume;
            if((*density_list)(i)==0) (*density_list)(i)=TV::dimension==1?1:TV::dimension==2?100:1000;}}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d> FINITE_VOLUME<TV,d>::
~FINITE_VOLUME()
{
    delete dPi_dFe;delete dP_dFe;delete Be_scales_save;delete V;delete semi_implicit_data;delete node_stiffness;delete edge_stiffness;delete density_list;
    if(destroy_data){delete &strain_measure;delete &constitutive_model;delete plasticity_model;}
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,int d> FINITE_VOLUME<TV,d>* FINITE_VOLUME<TV,d>::
Create(T_OBJECT& object,CONSTITUTIVE_MODEL<T,d>* constitutive_model,PLASTICITY_MODEL<T,d>* plasticity_model,const bool verbose,const bool use_uniform_density)
{
    STRAIN_MEASURE<TV,d>* strain_measure=new STRAIN_MEASURE<TV,d>(object);
    if(verbose) strain_measure->Print_Altitude_Statistics();
    FINITE_VOLUME* fvm=new FINITE_VOLUME(use_uniform_density,*strain_measure,*constitutive_model,plasticity_model);
    fvm->destroy_data=true;return fvm;
}
//#####################################################################
// Function Save_V
//#####################################################################
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Save_V()
{
    if(!V) V=new ARRAY<MATRIX<T,d> >();
}
//#####################################################################
// Function Save_Stress_Derivative
//#####################################################################
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Save_Stress_Derivative()
{
    if(isotropic_model || anisotropic_model->use_isotropic_component_of_stress_derivative_only){
        if(!dPi_dFe) dPi_dFe=new ARRAY<DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d> >();delete dP_dFe;dP_dFe=0;}
    else{if(!dP_dFe) dP_dFe=new ARRAY<DIAGONALIZED_STRESS_DERIVATIVE<T,d> >();delete dPi_dFe;dPi_dFe=0;}
}
//#####################################################################
// Function Use_Quasistatics
//#####################################################################
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Use_Quasistatics()
{
    Save_V();Save_Stress_Derivative();
}
//#####################################################################
// Function Use_Stiffness_Matrix
//#####################################################################
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Use_Stiffness_Matrix()
{
    node_stiffness=new ARRAY<SYMMETRIC_MATRIX<T,TV::m> >();edge_stiffness=new ARRAY<MATRIX<T,TV::m> >();
    if(!strain_measure.mesh.segment_mesh) strain_measure.mesh.Initialize_Segment_Mesh();
    if(!strain_measure.mesh.element_edges) strain_measure.mesh.Initialize_Element_Edges();
    if(force_segments) delete force_segments;force_segments=new FORCE_ELEMENTS();
}
//#####################################################################
// Function Update_Be_Scales
//#####################################################################
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Update_Be_Scales()
{
    Be_scales.Resize(strain_measure.Dm_inverse.m,false,false);
    for(int t=1;t<=Be_scales.m;t++) Be_scales(t)=-(T)1/FACTORIAL<d>::value/strain_measure.Dm_inverse(t).Determinant();
}
//#####################################################################
// Function Initialize_Be_Scales_Save
//#####################################################################
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Initialize_Be_Scales_Save()
{
    Be_scales_save=new ARRAY<T>(Be_scales);
}
//#####################################################################
// Function Copy_Be_Scales_Save_Into_Be_Scales
//#####################################################################
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Copy_Be_Scales_Save_Into_Be_Scales(const ARRAY<int>& map)
{
    Be_scales=Be_scales_save->Subset(map);
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    strain_measure.mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    if(!strain_measure.mesh.segment_mesh) strain_measure.mesh.Initialize_Segment_Mesh();
    force_elements.Update(strain_measure.mesh.elements,particle_is_simulated);
    if(force_segments) force_segments->Update(strain_measure.mesh.segment_mesh->elements,particle_is_simulated);
    force_particles.Update(strain_measure.mesh.elements.Flattened(),particle_is_simulated);
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
namespace{
template<class T,int d> void Set_Inversion_Based_On_Implicit_Surface(FINITE_VOLUME<VECTOR<T,d>,d>&,const IMPLICIT_OBJECT<VECTOR<T,d> >&,const int)
{}
template<class T> void Set_Inversion_Based_On_Implicit_Surface(FINITE_VOLUME<VECTOR<T,3>,2>& fvm,const IMPLICIT_OBJECT<VECTOR<T,3> >& implicit_surface,const int triangle)
{
    VECTOR<T,3> normal=implicit_surface.Normal(fvm.strain_measure.mesh_object.Centroid(triangle));
    if(VECTOR<T,3>::Dot_Product(normal,fvm.U(triangle).Weighted_Normal())<0){
        fvm.Fe_hat(triangle).x22*=-1;fvm.U(triangle).Column(2)*=-1;}
}}
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    if(anisotropic_model && !V) PHYSBAM_FATAL_ERROR();
    int elements=strain_measure.Dm_inverse.m;
    U.Resize(elements,false,false);De_inverse_hat.Resize(elements,false,false);Fe_hat.Resize(elements,false,false);
    if(dPi_dFe) dPi_dFe->Resize(elements,false,false);if(dP_dFe) dP_dFe->Resize(elements,false,false);if(V) V->Resize(elements,false,false);
    if(node_stiffness){
        node_stiffness->Resize(strain_measure.mesh_object.particles.array_collection->Size(),false,false);
        ARRAYS_COMPUTATIONS::Fill(*node_stiffness,SYMMETRIC_MATRIX<T,TV::m>());}
    if(edge_stiffness){
        edge_stiffness->Resize(strain_measure.mesh.segment_mesh->elements.m,false,false);
        ARRAYS_COMPUTATIONS::Fill(*edge_stiffness,MATRIX<T,TV::m>());}
    MATRIX<T,d> V_local;
    if(!plasticity_model){
        MATRIX<MATRIX<T,TV::m>,d+1,d+1> dfdx;
        for(FORCE_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            strain_measure.F(t).Fast_Singular_Value_Decomposition(U(t),Fe_hat(t),V_local);
            if(implicit_surface) Set_Inversion_Based_On_Implicit_Surface(*this,*implicit_surface,t);
            if(anisotropic_model) anisotropic_model->Update_State_Dependent_Auxiliary_Variables(Fe_hat(t),V_local,t);
            else isotropic_model->Update_State_Dependent_Auxiliary_Variables(Fe_hat(t),t);
            if(dPi_dFe) constitutive_model.Isotropic_Stress_Derivative(Fe_hat(t),(*dPi_dFe)(t),t);
            if(dP_dFe) anisotropic_model->Stress_Derivative(Fe_hat(t),V_local,(*dP_dFe)(t),t);
            De_inverse_hat(t)=strain_measure.Dm_inverse(t)*V_local;if(V) (*V)(t)=V_local;
            if(node_stiffness && edge_stiffness){
                for(int k=1;k<=TV::m;k++) for(int l=1;l<=d;l++){
                    T_MATRIX dDs,dG;dDs(k,l)=(T)1;
                    if(!dPi_dFe && !dP_dFe) // precompute damping matrix
                        dG=U(t)*constitutive_model.P_From_Strain_Rate(Fe_hat(t),U(t).Transpose_Times(dDs)*De_inverse_hat(t),Be_scales(t),t).Times_Transpose(De_inverse_hat(t));
                    else if(isotropic_model) // precompute stiffness matrix
                        dG=U(t)*isotropic_model->dP_From_dF(U(t).Transpose_Times(dDs)*De_inverse_hat(t),(*dPi_dFe)(t),Be_scales(t),t).Times_Transpose(De_inverse_hat(t));
                    else if(dP_dFe)
                        dG=U(t)*anisotropic_model->dP_From_dF(U(t).Transpose_Times(dDs)*De_inverse_hat(t),Fe_hat(t),(*V)(t),(*dP_dFe)(t),Be_scales(t),t)
                            .Times_Transpose(De_inverse_hat(t));
                    else
                        dG=U(t)*anisotropic_model->dP_From_dF(U(t).Transpose_Times(dDs)*De_inverse_hat(t),Fe_hat(t),(*V)(t),(*dPi_dFe)(t),Be_scales(t),t)
                            .Times_Transpose(De_inverse_hat(t));
                    for(int i=1;i<=TV::m;i++) for(int j=1;j<=d;j++) dfdx(l+1,j+1)(k,i)=dG(i,j);}
                for(int i=2;i<=d+1;i++){dfdx(i,1)=MATRIX<T,TV::m>();for(int j=2;j<=d+1;j++) dfdx(i,1)-=dfdx(i,j);}
                for(int j=1;j<=d+1;j++){dfdx(1,j)=MATRIX<T,TV::m>();for(int i=2;i<=d+1;i++) dfdx(1,j)-=dfdx(i,j);}
                VECTOR<int,d+1> nodes=strain_measure.mesh.elements(t);
                for(int v=1;v<=nodes.m;v++) (*node_stiffness)(nodes[v])+=dfdx(v,v).Symmetric_Part();
                VECTOR<int,d*(d+1)/2> edges=(*strain_measure.mesh.element_edges)(t);
                for(int e=1;e<=edges.m;e++){
                    int i,j;strain_measure.mesh.segment_mesh->elements(edges[e]).Get(i,j);
                    int m=nodes.Find(i),n=nodes.Find(j);assert(m!=n);
                    (*edge_stiffness)(edges[e])+=dfdx(m,n);}}}}
    else{
        if(implicit_surface) PHYSBAM_NOT_IMPLEMENTED();
        if(anisotropic_model) PHYSBAM_NOT_IMPLEMENTED();
        DIAGONAL_MATRIX<T,d> Fe_project_hat;
        for(FORCE_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            T_MATRIX F=strain_measure.F(t);
            (F*plasticity_model->Fp_inverse(t)).Fast_Singular_Value_Decomposition(U(t),Fe_hat(t),V_local);
            if(implicit_surface) Set_Inversion_Based_On_Implicit_Surface(*this,*implicit_surface,t);
            if(plasticity_model->Project_Fe(Fe_hat(t),Fe_project_hat)){
                plasticity_model->Project_Fp(t,Fe_project_hat.Inverse()*U(t).Transpose_Times(F));
                (F*plasticity_model->Fp_inverse(t)).Fast_Singular_Value_Decomposition(U(t),Fe_hat(t),V_local);}
            De_inverse_hat(t)=strain_measure.Dm_inverse(t)*plasticity_model->Fp_inverse(t)*V_local;
            Be_scales(t)=-(T)1/FACTORIAL<d>::value/De_inverse_hat(t).Determinant();}}

    if(compute_half_forces){
        sqrt_Be_scales.Resize(Be_scales.m);
        for(int i=1;i<=Be_scales.m;i++) sqrt_Be_scales(i)=sqrt(-Be_scales(i));
        int n=0;
        for(FORCE_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()) n++;
        half_force_size=constitutive_model.P_From_Strain_Rate_Forces_Size();
        total_half_force_size=half_force_size*n;}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    if(anisotropic_model)
        for(FORCE_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            T_MATRIX forces=U(t)*anisotropic_model->P_From_Strain(Fe_hat(t),(*V)(t),Be_scales(t),t).Times_Transpose(De_inverse_hat(t));
            strain_measure.Distribute_Force(F,t,forces);}
    else
        for(FORCE_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            T_MATRIX forces=U(t)*isotropic_model->P_From_Strain(Fe_hat(t),Be_scales(t),t).Times_Transpose(De_inverse_hat(t));
            strain_measure.Distribute_Force(F,t,forces);}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    //force_elements.Print();
    if(node_stiffness && edge_stiffness && !dPi_dFe && !dP_dFe){
        if(dP_dFe || dPi_dFe) PHYSBAM_FATAL_ERROR();
        for(FORCE_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
            F(p)+=(*node_stiffness)(p)*V(p);}
        for(FORCE_ITERATOR iterator(*force_segments);iterator.Valid();iterator.Next()){int e=iterator.Data();
            int m,n;strain_measure.mesh.segment_mesh->elements(e).Get(m,n);
            F(m)+=(*edge_stiffness)(e)*V(n);F(n)+=(*edge_stiffness)(e).Transpose_Times(V(m));}}
    else{
        for(FORCE_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            MATRIX<T,d> Fe_dot_hat=U(t).Transpose_Times(strain_measure.Ds(V,t))*De_inverse_hat(t);
            T_MATRIX forces=U(t)*constitutive_model.P_From_Strain_Rate(Fe_hat(t),Fe_dot_hat,Be_scales(t),t).Times_Transpose(De_inverse_hat(t));
            strain_measure.Distribute_Force(F,t,forces);}}
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
    constitutive_model.enforce_definiteness=enforce_definiteness_input;
}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const
{
    if(node_stiffness && edge_stiffness){
        for(FORCE_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
            dF(p)+=(*node_stiffness)(p)*dX(p);}
        for(FORCE_ITERATOR iterator(*force_segments);iterator.Valid();iterator.Next()){int e=iterator.Data();
            int m,n;strain_measure.mesh.segment_mesh->elements(e).Get(m,n);
            dF(m)+=(*edge_stiffness)(e)*dX(n);dF(n)+=(*edge_stiffness)(e).Transpose_Times(dX(m));}}
    else if(anisotropic_model){
        if(!dPi_dFe && !dP_dFe) PHYSBAM_FATAL_ERROR();
        for(FORCE_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            T_MATRIX dDs=strain_measure.Ds(dX,t),dG;
            if(dP_dFe) dG=U(t)*anisotropic_model->dP_From_dF(U(t).Transpose_Times(dDs)*De_inverse_hat(t),Fe_hat(t),(*V)(t),(*dP_dFe)(t),Be_scales(t),t).Times_Transpose(De_inverse_hat(t));
            else dG=U(t)*anisotropic_model->dP_From_dF(U(t).Transpose_Times(dDs)*De_inverse_hat(t),Fe_hat(t),(*V)(t),(*dPi_dFe)(t),Be_scales(t),t).Times_Transpose(De_inverse_hat(t));
            strain_measure.Distribute_Force(dF,t,dG);}}
    else{
        if(!dPi_dFe) PHYSBAM_FATAL_ERROR();
        for(FORCE_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            T_MATRIX dDs=strain_measure.Ds(dX,t);
            T_MATRIX dG=U(t)*isotropic_model->dP_From_dF(U(t).Transpose_Times(dDs)*De_inverse_hat(t),(*dPi_dFe)(t),Be_scales(t),t).Times_Transpose(De_inverse_hat(t));
            strain_measure.Distribute_Force(dF,t,dG);}}
}
//#####################################################################
// Function Intialize_CFL
//#####################################################################
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)
{
    // TODO: MPI
    T one_over_cfl_number=1/cfl_number,one_over_cfl_number_squared=sqr(one_over_cfl_number);
    HASHTABLE<int,FREQUENCY_DATA> particle_frequency;
    for(FORCE_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        T one_over_altitude_squared_and_density=(T)1/(sqr(strain_measure.Rest_Altitude(t))*(use_uniform_density?density:(*density_list)(t)));
        T elastic_squared=constitutive_model.Maximum_Elastic_Stiffness(t)*one_over_altitude_squared_and_density*one_over_cfl_number_squared;
        T damping=constitutive_model.Maximum_Damping_Stiffness(t)*one_over_altitude_squared_and_density*one_over_cfl_number;
        const VECTOR<int,d+1>& nodes=strain_measure.mesh.elements(t);
        for(int j=1;j<=nodes.m;j++){FREQUENCY_DATA& frequency_data=particle_frequency.Get_Or_Insert(nodes[j]);
            frequency_data.elastic_squared=max(frequency_data.elastic_squared,elastic_squared);
            frequency_data.damping=max(frequency_data.damping,damping);}}
    for(typename HASHTABLE<int,FREQUENCY_DATA>::ITERATOR iterator(particle_frequency);iterator.Valid();iterator.Next()){
        frequency(iterator.Key()).elastic_squared+=iterator.Data().elastic_squared;
        frequency(iterator.Key()).damping+=iterator.Data().damping;}
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV,int d> typename TV::SCALAR FINITE_VOLUME<TV,d>::
CFL_Strain_Rate() const
{
    T max_strain_rate=0;
    for(FORCE_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        max_strain_rate=max(max_strain_rate,strain_measure.Velocity_Gradient(t).Max_Abs());}
    return Robust_Divide(max_strain_per_time_step,max_strain_rate);
}
//#####################################################################
// Function Semi_Implicit_Impulse_Precomputation
//#####################################################################
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Semi_Implicit_Impulse_Precomputation(const T time,const T max_dt,ARRAY<T>* time_plus_dt,const bool verbose)
{
    assert(!plasticity_model && isotropic_model);
    int elements=strain_measure.Dm_inverse.m;
    if(!semi_implicit_data) semi_implicit_data=new ARRAY<DIAGONALIZED_SEMI_IMPLICIT_ELEMENT<TV,d> >(elements,false);
    else semi_implicit_data->Resize(elements,false,false);
    twice_max_strain_per_time_step=2*max_strain_per_time_step;
    for(FORCE_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        DIAGONALIZED_SEMI_IMPLICIT_ELEMENT<TV,d>& data=(*semi_implicit_data)(t);
        data.nodes=strain_measure.mesh.elements(t);
        data.Bm_scale=Be_scales(t);
        VECTOR<T,d+1> masses(particles.mass.Subset(data.nodes));
        DIAGONAL_MATRIX<T,d> M_inverse=DIAGONAL_MATRIX<T,d>(masses.Remove_Index(1)).Inverse();
        SYMMETRIC_MATRIX<T,d> W=SYMMETRIC_MATRIX<T,d>::Conjugate_With_Transpose(strain_measure.Dm_inverse(t),
            -data.Bm_scale*M_inverse-SYMMETRIC_MATRIX<T,d>::Unit_Matrix(data.Bm_scale/masses[1])); // note sign of scale factor
        MATRIX<T,d> V;W.Fast_Solve_Eigenproblem(data.W,V);
        data.Dm_inverse=strain_measure.Dm_inverse(t)*V;
        data.dt_cfl=min(max_dt,cfl_number/sqrt(constitutive_model.Maximum_Elastic_Stiffness(t)/((use_uniform_density?density:(*density_list)(t))*sqr(strain_measure.Rest_Altitude(t)))));
        data.time=time;if(time_plus_dt) (*time_plus_dt)(t)=time+min(data.dt_cfl,twice_max_strain_per_time_step/strain_measure.Velocity_Gradient(t).Max_Abs());}
    if(verbose){
        T min_W=FLT_MAX,max_W=0,min_dt_cfl=FLT_MAX,max_dt_cfl=0;
        for(FORCE_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            DIAGONALIZED_SEMI_IMPLICIT_ELEMENT<TV,d>& data=(*semi_implicit_data)(t);
            min_W=min(min_W,data.W.Min());max_W=max(max_W,data.W.Max());
            min_dt_cfl=min(min_dt_cfl,data.dt_cfl);max_dt_cfl=max(max_dt_cfl,data.dt_cfl);}
        std::stringstream ss;ss<<"W range = "<<min_W<<" to "<<max_W<<std::endl;
        ss<<"dt_cfl range = "<<min_dt_cfl<<" to "<<max_dt_cfl<<std::endl;LOG::filecout(ss.str());}
}
//#####################################################################
// Function Semi_Implicit_Recompute_Dt
//#####################################################################
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Semi_Implicit_Recompute_Dt(const int element,T& time_plus_dt)
{
    assert(semi_implicit_data);
    DIAGONALIZED_SEMI_IMPLICIT_ELEMENT<TV,d>& data=(*semi_implicit_data)(element);
    T_MATRIX F_dot=STRAIN_MEASURE<TV,d>::Ds(particles.V,data.nodes)*data.Dm_inverse;
    time_plus_dt=data.time+min(data.dt_cfl,twice_max_strain_per_time_step/F_dot.Max_Abs());
}
//#####################################################################
// Function Add_Semi_Implicit_Impulse
//#####################################################################
namespace{
template<class T> inline SYMMETRIC_MATRIX<T,2> Compute_P1_cap(const T dt_beta,const T dt_alpha,const DIAGONAL_MATRIX<T,2>& W,const SYMMETRIC_MATRIX<T,2>& P0_cap)
{
    VECTOR<T,2> P1_cap_diag=((T)1+((T)2*dt_beta+SYMMETRIC_MATRIX<T,2>::Unit_Matrix(dt_alpha))*W).Solve_Linear_System(VECTOR<T,2>(P0_cap.x11,P0_cap.x22));
    T P1_cap_off=P0_cap.x21/((T)1+dt_beta*(W.x11+W.x22));
    return SYMMETRIC_MATRIX<T,2>(P1_cap_diag.x,P1_cap_off,P1_cap_diag.y);
}
template<class T> inline SYMMETRIC_MATRIX<T,3> Compute_P1_cap(const T dt_beta,const T dt_alpha,const DIAGONAL_MATRIX<T,3>& W,const SYMMETRIC_MATRIX<T,3>& P0_cap)
{
    VECTOR<T,3> P1_cap_diag=((T)1+(2*dt_beta+SYMMETRIC_MATRIX<T,3>::Unit_Matrix(dt_alpha))*W).Solve_Linear_System(VECTOR<T,3>(P0_cap.x11,P0_cap.x22,P0_cap.x33));
    VECTOR<T,3> P1_cap_off=((T)1+dt_beta*DIAGONAL_MATRIX<T,3>(W.x11+W.x22,W.x11+W.x33,W.x22+W.x33)).Solve_Linear_System(VECTOR<T,3>(P0_cap.x21,P0_cap.x31,P0_cap.x32));
    return SYMMETRIC_MATRIX<T,3>(P1_cap_diag.x,P1_cap_off.x,P1_cap_off.y,P1_cap_diag.y,P1_cap_off.z,P1_cap_diag.z);
}
}
//#####################################################################
// Function Add_Semi_Implicit_Impulse
//#####################################################################
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Add_Semi_Implicit_Impulse(const int element,const T dt,T* time_plus_dt)
{
    assert(semi_implicit_data);
    DIAGONALIZED_SEMI_IMPLICIT_ELEMENT<TV,d>& data=(*semi_implicit_data)(element);
    // update time
    data.time+=dt;
    // diagonalize F
    T_MATRIX U;MATRIX<T,d> V;DIAGONAL_MATRIX<T,d> F_hat;
    T_MATRIX F=STRAIN_MEASURE<TV,d>::Ds(particles.X,data.nodes)*data.Dm_inverse;
    F.Fast_Singular_Value_Decomposition(U,F_hat,V);
    T_MATRIX Q=U.Times_Transpose(V);
    // compute -s dt P^n
    T dt_scale=dt*data.Bm_scale;
    SYMMETRIC_MATRIX<T,d> P0_cap=SYMMETRIC_MATRIX<T,d>::Conjugate(V,isotropic_model->P_From_Strain(F_hat,dt_scale,element));
    MATRIX<T,d> F_dot_cap=Q.Transpose_Times(STRAIN_MEASURE<TV,d>::Ds(particles.V,data.nodes))*data.Dm_inverse;
    SYMMETRIC_MATRIX<T,d> F_dot_cap_twice_symmetric_part=F_dot_cap.Twice_Symmetric_Part();
    if(time_plus_dt) *time_plus_dt=data.time+min(data.dt_cfl,twice_max_strain_per_time_step/F_dot_cap_twice_symmetric_part.Max_Abs());
    T alpha=isotropic_model->alpha.m?isotropic_model->alpha(element):isotropic_model->constant_alpha;
    T beta=isotropic_model->beta.m?isotropic_model->beta(element):isotropic_model->constant_beta;
    P0_cap+=dt_scale*beta*F_dot_cap_twice_symmetric_part+dt_scale*alpha*F_dot_cap.Trace();
    // solve for -s dt P^{n+1}
    T dt_beta=dt*beta,dt_alpha=dt*alpha;
    SYMMETRIC_MATRIX<T,d> P1_cap=Compute_P1_cap(dt_beta,dt_alpha,data.W,P0_cap);
    // change velocities
    T_MATRIX G=(Q*P1_cap).Times_Transpose(data.Dm_inverse);
    STRAIN_MEASURE<TV,d>::Distribute_Impulse(particles.V,data.nodes,G,particles.one_over_mass);
}
//#####################################################################
// Function Add_Semi_Implicit_Impulses
//#####################################################################
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Add_Semi_Implicit_Impulses(const T dt,ARRAY<T>* time_plus_dt)
{
    if(time_plus_dt) for(FORCE_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();Add_Semi_Implicit_Impulse(t,dt,&(*time_plus_dt)(t));}
    else for(FORCE_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();Add_Semi_Implicit_Impulse(t,dt,0);}
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
template<class TV,int d> int FINITE_VOLUME<TV,d>::
Velocity_Dependent_Forces_Size() const
{
    return total_half_force_size;
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
    PHYSBAM_ASSERT(!(node_stiffness && edge_stiffness && !dPi_dFe && !dP_dFe) && compute_half_forces);
    int start_force=1;
    for(FORCE_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        MATRIX<T,d> Fe_dot_hat=U(t).Transpose_Times(strain_measure.Ds(V,t))*De_inverse_hat(t);
        ARRAY_VIEW<T,int> intermediate(half_force_size,const_cast<T *>(aggregate.Get_Array_Pointer())+Value(start_force-1));
        constitutive_model.P_From_Strain_Rate_First_Half(Fe_hat(t),intermediate,Fe_dot_hat,sqrt_Be_scales(t),t);
        start_force+=half_force_size;}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class TV,int d> void FINITE_VOLUME<TV,d>::
Add_Velocity_Dependent_Forces_Second_Half(const ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
    int start_force=1;
    PHYSBAM_ASSERT(!(node_stiffness && edge_stiffness && !dPi_dFe && !dP_dFe) && compute_half_forces);
    for(FORCE_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        ARRAY_VIEW<T,int> intermediate(half_force_size,const_cast<T *>(aggregate.Get_Array_Pointer())+Value(start_force-1));
        T_MATRIX forces=U(t)*constitutive_model.P_From_Strain_Rate_Second_Half(Fe_hat(t),intermediate,sqrt_Be_scales(t),t).Times_Transpose(De_inverse_hat(t));
        strain_measure.Distribute_Force(F,t,forces);
        start_force+=half_force_size;}
}
template<class T,int d>
T Determinant_Helper(const MATRIX<T,d>& M)
{
    return M.Determinant();
}
template<class T,int m,int n>
T Determinant_Helper(const MATRIX<T,m,n>& M)
{
    PHYSBAM_NOT_IMPLEMENTED();
    return 0;
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV,int d> typename TV::SCALAR FINITE_VOLUME<TV,d>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    if(anisotropic_model) PHYSBAM_NOT_IMPLEMENTED();
    else
        for(FORCE_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            potential_energy-=Be_scales(t)*isotropic_model->Energy_Density(Fe_hat(t),t);}
    return potential_energy;
}
//#####################################################################
template class FINITE_VOLUME<VECTOR<float,2>,2>;
template class FINITE_VOLUME<VECTOR<float,3>,2>;
template class FINITE_VOLUME<VECTOR<float,3>,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FINITE_VOLUME<VECTOR<double,2>,2>;
template class FINITE_VOLUME<VECTOR<double,3>,2>;
template class FINITE_VOLUME<VECTOR<double,3>,3>;
#endif

//#####################################################################
// Copyright 2006-2008, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Craig Schroeder, Andrew Selle, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_BACKWARD_EULER_SYSTEM
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/INNER_PRODUCT.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_MASS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigids_Evolution/RIGIDS_BACKWARD_EULER_SYSTEM.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigids_Evolution/RIGIDS_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigids_Evolution/RIGIDS_VELOCITY.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGIDS_BACKWARD_EULER_SYSTEM<TV>::
RIGIDS_BACKWARD_EULER_SYSTEM(RIGIDS_EVOLUTION<TV>& rigids_evolution_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const T dt,const T current_velocity_time,
    const T current_position_time,ARTICULATED_RIGID_BODY<TV>* arb_input,const bool velocity_update_input,const bool enforce_poststabilization_in_cg)
    :KRYLOV_SYSTEM_BASE<typename TV::SCALAR>(false,true),rigids_evolution(rigids_evolution_input),rigid_body_collection(rigid_body_collection),
    dt(dt),current_velocity_time(current_velocity_time),current_position_time(current_position_time),arb(arb_input),
    velocity_update(velocity_update_input),project_nullspace_frequency(INT_MAX),projection_data(rigid_body_collection)
{
    if(!enforce_poststabilization_in_cg) arb=0;
    if(arb) arb->Initialize_Poststabilization_Projection();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGIDS_BACKWARD_EULER_SYSTEM<TV>::
~RIGIDS_BACKWARD_EULER_SYSTEM()
{
}
//#####################################################################
// Function Force
//#####################################################################
template<class TV> void RIGIDS_BACKWARD_EULER_SYSTEM<TV>::
Force(const VECTOR_T& V,VECTOR_T& F) const
{
    INDIRECT_ARRAY<ARRAY_VIEW<TWIST<TV> >,ARRAY<int>&> rigid_V_subset=F.rigid_V.array.Subset(rigid_body_collection.simulated_rigid_body_particles);
    ARRAYS_COMPUTATIONS::Fill(rigid_V_subset,TWIST<TV>());
    rigid_body_collection.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Velocities(V.rigid_V.array);
    rigid_body_collection.Implicit_Velocity_Independent_Forces(V.rigid_V.array,F.rigid_V.array,dt,current_velocity_time+dt);
    rigid_body_collection.Add_Velocity_Dependent_Forces(V.rigid_V.array,F.rigid_V.array,current_velocity_time+dt);
    rigid_body_collection.rigid_body_cluster_bindings.Distribute_Force_To_Parents(F.rigid_V.array);
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void RIGIDS_BACKWARD_EULER_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF) const
{
    const VECTOR_T& V=debug_cast<const VECTOR_T&>(BV);VECTOR_T& F=debug_cast<VECTOR_T&>(BF);
    Force(V,F);
    for(int i=1;i<=V.rigid_V.Size();i++) F.rigid_V(i)=V.rigid_V(i)-dt*(projection_data.mass.world_space_rigid_mass_inverse(i)*F.rigid_V(i));
}
//#####################################################################
// Function Set_Global_Boundary_Conditions
//#####################################################################
template<class TV> void RIGIDS_BACKWARD_EULER_SYSTEM<TV>::
Set_Global_Boundary_Conditions(VECTOR_T& V,ARRAY<TV>& rigid_X_save,ARRAY<ROTATION<TV> >& rigid_rotation_save,
    ARRAY<TWIST<TV> >& rigid_velocity_save,ARRAY<typename TV::SPIN>& rigid_angular_momentum_save,bool test_system,bool print_matrix) const
{
    rigids_evolution.kinematic_evolution.Set_External_Velocities(V.rigid_V.array,current_velocity_time+dt,current_position_time);
    // TODO: make Solve_Velocities_for_PD take rigid_V.array and call that instead
    if(arb){
        //PHYSBAM_ASSERT(ARRAY_VIEW<TWIST<TV> >::Same_Array(V.rigid_V.array,rigid_body_collection.rigid_body_particle.twist));
        arb->Apply_Poststabilization(test_system,print_matrix); // Do not project out pd directions here
        if(arb->Has_Actuators() && arb->constrain_pd_directions){
            arb->Compute_Position_Based_State(dt,current_velocity_time);arb->Solve_Velocities_for_PD(current_velocity_time,dt,test_system,print_matrix);}}
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void RIGIDS_BACKWARD_EULER_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& BV) const
{
    VECTOR_T& V=debug_cast<VECTOR_T&>(BV);
    // Applying the projections in this order is equivalent to repeating Zero_Out_Enslaved_Velocity_Nodes after Poststabilization_Projection, which is a (mass) symmetric projection.
    rigids_evolution.Zero_Out_Enslaved_Velocity_Nodes(V.rigid_V.array,current_velocity_time+dt,current_position_time);
    for(int i=1;i<=rigids_evolution.rigids_parameters.implicit_solve_parameters.cg_projection_iterations;i++){
        if(arb){
            arb->Poststabilization_Projection(V.rigid_V.array,true);
            rigids_evolution.Zero_Out_Enslaved_Velocity_Nodes(V.rigid_V.array,current_velocity_time+dt,current_position_time);}}
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double RIGIDS_BACKWARD_EULER_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& BV1,const KRYLOV_VECTOR_BASE<T>& BV2) const
{
    const VECTOR_T& V1=debug_cast<const VECTOR_T&>(BV1),&V2=debug_cast<const VECTOR_T&>(BV2);
    double inner_product=ARRAYS_COMPUTATIONS::Inner_Product_Double_Precision(projection_data.mass.world_space_rigid_mass,V1.rigid_V,V2.rigid_V);
    return inner_product;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR RIGIDS_BACKWARD_EULER_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const
{
    const VECTOR_T& R=debug_cast<const VECTOR_T&>(BR);
    T convergence_norm_squared=0;
    for(int p=1;p<=R.rigid_V.Size();p++){
        const TWIST<TV>& twist=R.rigid_V(p);
        const RIGID_BODY_MASS<TV,true> &rigid_mass=projection_data.mass.world_space_rigid_mass(p),&rigid_mass_inverse=projection_data.mass.world_space_rigid_mass_inverse(p);
        convergence_norm_squared=max(convergence_norm_squared,
            twist.linear.Magnitude_Squared()+rigid_mass_inverse.mass*Dot_Product(twist.angular,rigid_mass.inertia_tensor*twist.angular));}
    T convergence_norm=sqrt(convergence_norm_squared);
    return convergence_norm;
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void RIGIDS_BACKWARD_EULER_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const
{
    if(++projection_data.project_nullspace_counter>=project_nullspace_frequency){Project(V);projection_data.project_nullspace_counter=0;}
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void RIGIDS_BACKWARD_EULER_SYSTEM<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const
{
    R.Copy((T)1,V);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGIDS_MASS<TV>::
RIGIDS_MASS(RIGID_BODY_COLLECTION<TV>& rigid_body_collection)
    :rigid_mass(rigid_body_collection.rigid_body_particle.mass,rigid_body_collection.dynamic_rigid_body_particles),
    rigid_inertia_tensor(rigid_body_collection.rigid_body_particle.inertia_tensor,rigid_body_collection.dynamic_rigid_body_particles)
{
    Initialize_World_Space_Masses(rigid_body_collection);
}
//#####################################################################
// Function Initialize_World_Space_Masses
//#####################################################################
template<class TV> void RIGIDS_MASS<TV>::
Initialize_World_Space_Masses(const RIGID_BODY_COLLECTION<TV>& rigid_body_collection)
{
    world_space_rigid_mass.Resize(rigid_body_collection.dynamic_rigid_body_particles.m,false,false);
    world_space_rigid_mass_inverse.Resize(rigid_body_collection.dynamic_rigid_body_particles.m,false,false);
    for(int i=1;i<=rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=rigid_body_collection.dynamic_rigid_body_particles(i);
        world_space_rigid_mass(i)=rigid_body_collection.State(p).World_Space_Rigid_Mass(RIGID_BODY_MASS<TV>(rigid_mass(i),rigid_inertia_tensor(i)));
        world_space_rigid_mass_inverse(i)=rigid_body_collection.State(p).World_Space_Rigid_Mass_Inverse(RIGID_BODY_MASS<TV>(rigid_mass(i),rigid_inertia_tensor(i)));}
}
//#####################################################################
// Function Inverse_Multiply
//#####################################################################
template<class TV> void RIGIDS_MASS<TV>::
Inverse_Multiply(const RIGIDS_VELOCITY<TV>& V,RIGIDS_VELOCITY<TV>& F) const
{
    F.rigid_V=world_space_rigid_mass_inverse*V.rigid_V;
}
//#####################################################################
template class RIGIDS_MASS<VECTOR<float,1> >;
template class RIGIDS_MASS<VECTOR<float,2> >;
template class RIGIDS_MASS<VECTOR<float,3> >;
template class RIGIDS_BACKWARD_EULER_SYSTEM<VECTOR<float,1> >;
template class RIGIDS_BACKWARD_EULER_SYSTEM<VECTOR<float,2> >;
template class RIGIDS_BACKWARD_EULER_SYSTEM<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGIDS_MASS<VECTOR<double,1> >;
template class RIGIDS_MASS<VECTOR<double,2> >;
template class RIGIDS_MASS<VECTOR<double,3> >;
template class RIGIDS_BACKWARD_EULER_SYSTEM<VECTOR<double,1> >;
template class RIGIDS_BACKWARD_EULER_SYSTEM<VECTOR<double,2> >;
template class RIGIDS_BACKWARD_EULER_SYSTEM<VECTOR<double,3> >;
#endif

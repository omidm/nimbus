//#####################################################################
// Copyright 2006-2008, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Craig Schroeder, Andrew Selle, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BACKWARD_EULER_SYSTEM
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/INNER_PRODUCT.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/ASYNCHRONOUS_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BACKWARD_EULER_SYSTEM<TV>::
BACKWARD_EULER_SYSTEM(SOLIDS_EVOLUTION<TV>& solids_evolution_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection,const T dt,const T current_velocity_time,
    const T current_position_time,ARTICULATED_RIGID_BODY<TV>* arb_input,TRIANGLE_REPULSIONS<TV>* repulsions_input,MPI_SOLIDS<TV>* mpi_solids,const bool velocity_update_input)
    :KRYLOV_SYSTEM_BASE<typename TV::SCALAR>(false,true),solids_evolution(solids_evolution_input),solid_body_collection(solid_body_collection),
    dt(dt),current_velocity_time(current_velocity_time),current_position_time(current_position_time),//mass(deformable_object),
    arb(arb_input),mpi_solids(mpi_solids),repulsions(repulsions_input),velocity_update(velocity_update_input),project_nullspace_frequency(INT_MAX),projection_data(solid_body_collection)
{
    if(!solids_evolution.solids_parameters.enforce_poststabilization_in_cg) arb=0;
    if(repulsions) repulsions->Set_Collision_Pairs(projection_data.point_face_precomputed,projection_data.edge_edge_precomputed,projection_data.point_face_pairs,projection_data.edge_edge_pairs,(T)1);
    if(arb) arb->Initialize_Poststabilization_Projection();
    // TODO: compare enforcing constraints at n+1/2 and n+1
    //if(solids_evolution.solids_parameters.use_rigid_deformable_contact && solid_body_collection.deformable_object.collisions.collisions_on)
    //    solids_evolution.rigid_deformable_collisions->Initialize_All_Contact_Projections();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BACKWARD_EULER_SYSTEM<TV>::
~BACKWARD_EULER_SYSTEM()
{
}
//#####################################################################
// Function Force
//#####################################################################
template<class TV> void BACKWARD_EULER_SYSTEM<TV>::
Force(const VECTOR_T& V,VECTOR_T& F) const
{
    if(mpi_solids) mpi_solids->Exchange_Binding_Boundary_Data(V.V.array);

    INDIRECT_ARRAY<ARRAY_VIEW<TV>,ARRAY<int>&> V_subset=F.V.array.Subset(solid_body_collection.deformable_body_collection.simulated_particles);
    INDIRECT_ARRAY<ARRAY_VIEW<TWIST<TV> >,ARRAY<int>&> rigid_V_subset=F.rigid_V.array.Subset(solid_body_collection.rigid_body_collection.simulated_rigid_body_particles);
    ARRAYS_COMPUTATIONS::Fill(V_subset,TV());ARRAYS_COMPUTATIONS::Fill(rigid_V_subset,TWIST<TV>());
    solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Velocities(V.rigid_V.array);

    solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities(V.V.array,V.rigid_V.array);
    if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data(V.V.array);
    if(!velocity_update || !solids_evolution.solids_parameters.implicit_solve_parameters.use_half_fully_implicit) solid_body_collection.Implicit_Velocity_Independent_Forces(V.V.array,V.rigid_V.array,F.V.array,F.rigid_V.array,dt,current_velocity_time+dt);
    // we don't inherit velocity dependent forces to the drifted particles (contrary to the explicit velocity independent ones) because that would compromise symmetry
    solid_body_collection.Add_Velocity_Dependent_Forces(V.V.array,V.rigid_V.array,F.V.array,F.rigid_V.array,current_velocity_time+dt);
    if(mpi_solids) mpi_solids->Exchange_Binding_Boundary_Data(F.V.array);
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Force_To_Parents(F.V.array,F.rigid_V.array);
    solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Distribute_Force_To_Parents(F.rigid_V.array);
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void BACKWARD_EULER_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF) const
{
    const VECTOR_T& V=debug_cast<const VECTOR_T&>(BV);VECTOR_T& F=debug_cast<VECTOR_T&>(BF);
    Force(V,F);
    for(int i=1;i<=V.V.Size();i++) F.V(i)=V.V(i)-dt*projection_data.mass.one_over_mass(i)*F.V(i);
    for(int i=1;i<=V.rigid_V.Size();i++) F.rigid_V(i)=V.rigid_V(i)-dt*(projection_data.mass.world_space_rigid_mass_inverse(i)*F.rigid_V(i));
}
//#####################################################################
// Function Set_Global_Boundary_Conditions
//#####################################################################
template<class TV> void BACKWARD_EULER_SYSTEM<TV>::
Set_Global_Boundary_Conditions(VECTOR_T& V,ARRAY<TV>& X_save,ARRAY<TV>& rigid_X_save,ARRAY<ROTATION<TV> >& rigid_rotation_save,
    ARRAY<TWIST<TV> >& rigid_velocity_save,ARRAY<typename TV::SPIN>& rigid_angular_momentum_save,ARRAY<TV>& V_save,bool test_system,bool print_matrix) const
{
    SOLIDS_PARAMETERS<TV>& solids_parameters=solids_evolution.solids_parameters;
    solids_evolution.Set_External_Velocities(V.V.array,current_velocity_time+dt,current_position_time);
    solids_evolution.kinematic_evolution.Set_External_Velocities(V.rigid_V.array,current_velocity_time+dt,current_position_time);
    // TODO: make Solve_Velocities_for_PD take rigid_V.array and call that instead
    if(arb){
        //PHYSBAM_ASSERT(ARRAY_VIEW<TWIST<TV> >::Same_Array(V.rigid_V.array,solid_body_collection.rigid_body_collection.rigid_body_particle.twist));
        arb->Apply_Poststabilization(test_system,print_matrix); // Do not project out pd directions here
        if(arb->Has_Actuators() && arb->constrain_pd_directions){
            arb->Compute_Position_Based_State(dt,current_velocity_time);arb->Solve_Velocities_for_PD(current_velocity_time,dt,solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);}}
    if(velocity_update){
        //PHYSBAM_ASSERT(ARRAY_VIEW<TV>::Same_Array(V.V.array,solid_body_collection.deformable_body_collection.particles.V) && ARRAY_VIEW<TWIST<TV> >::Same_Array(V.rigid_V.array,solid_body_collection.rigid_body_collection.rigid_body_particle.twist));
        if(solids_parameters.use_post_cg_constraints){// TODO: may just want to call Apply_Constraints in this case too
            if(solids_evolution.solids_parameters.use_rigid_deformable_contact && solid_body_collection.deformable_body_collection.collisions.collisions_on)
            {
                ARRAY_VIEW<const ROTATION<TV> > rigid_rotation_save_array_view(rigid_rotation_save);
                solids_evolution.rigid_deformable_collisions->Set_Collision_Velocities(V.V.array,V.rigid_V.array,X_save,rigid_X_save,rigid_rotation_save_array_view,rigid_velocity_save,rigid_angular_momentum_save,V_save);
            }
            if(repulsions) repulsions->Adjust_Velocity_For_Self_Repulsion_Using_History(dt,false,false);}}
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void BACKWARD_EULER_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& BV) const
{
    VECTOR_T& V=debug_cast<VECTOR_T&>(BV);
    // Applying the projections in this order is equivalent to repeating Zero_Out_Enslaved_Velocity_Nodes after Poststabilization_Projection, which is a (mass) symmetric projection.
    solids_evolution.Zero_Out_Enslaved_Velocity_Nodes(V.V.array,current_velocity_time+dt,current_position_time);
    solids_evolution.Zero_Out_Enslaved_Velocity_Nodes(V.rigid_V.array,current_velocity_time+dt,current_position_time);
    if(NEWMARK_EVOLUTION<TV>* newmark=dynamic_cast<NEWMARK_EVOLUTION<TV>*>(&solids_evolution)){
        if(newmark->asynchronous_evolution){
            for(int i=1;i<=solids_evolution.solids_parameters.implicit_solve_parameters.cg_projection_iterations;i++) newmark->asynchronous_evolution->Project(V);}}
    if(!velocity_update){if(arb) arb->Poststabilization_Projection(V.rigid_V.array,true);return;}
    for(int i=1;i<=solids_evolution.solids_parameters.implicit_solve_parameters.cg_projection_iterations;i++){
        int middle_projection=1;
        if(solids_evolution.solids_parameters.use_rigid_deformable_contact && solid_body_collection.deformable_body_collection.collisions.collisions_on){
            solids_evolution.rigid_deformable_collisions->Project_Contact_Pairs(V.V.array,V.rigid_V.array);middle_projection=2;}
        // TODO: arb projections should also project out changes in pd controlled directions
        if(arb){arb->Poststabilization_Projection(V.rigid_V.array,true);middle_projection=3;}
        if(projection_data.point_face_precomputed.m || projection_data.edge_edge_precomputed.m){
            TRIANGLE_REPULSIONS<TV>::Project_All_Moving_Constraints(projection_data.point_face_precomputed,projection_data.edge_edge_precomputed,V.V.array);middle_projection=4;}
        if(arb && middle_projection!=3) arb->Poststabilization_Projection(V.rigid_V.array,true);
        if(solids_evolution.solids_parameters.use_rigid_deformable_contact && solid_body_collection.deformable_body_collection.collisions.collisions_on && middle_projection!=2)
            solids_evolution.rigid_deformable_collisions->Project_Contact_Pairs(V.V.array,V.rigid_V.array);
        if(middle_projection!=1){
            solids_evolution.Zero_Out_Enslaved_Velocity_Nodes(V.V.array,current_velocity_time+dt,current_position_time);
            solids_evolution.Zero_Out_Enslaved_Velocity_Nodes(V.rigid_V.array,current_velocity_time+dt,current_position_time);}}
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double BACKWARD_EULER_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& BV1,const KRYLOV_VECTOR_BASE<T>& BV2) const
{
    const VECTOR_T& V1=debug_cast<const VECTOR_T&>(BV1),&V2=debug_cast<const VECTOR_T&>(BV2);
    double inner_product=ARRAYS_COMPUTATIONS::Inner_Product_Double_Precision(projection_data.mass.mass,V1.V,V2.V)+
        ARRAYS_COMPUTATIONS::Inner_Product_Double_Precision(projection_data.mass.world_space_rigid_mass,V1.rigid_V,V2.rigid_V);
    if(mpi_solids) inner_product=mpi_solids->Reduce_Add(inner_product);
    return inner_product;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR BACKWARD_EULER_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const
{
    const VECTOR_T& R=debug_cast<const VECTOR_T&>(BR);
    T convergence_norm_squared=sqr(ARRAYS_COMPUTATIONS::Maximum_Magnitude(R.V));
    for(int p=1;p<=R.rigid_V.Size();p++){
        const TWIST<TV>& twist=R.rigid_V(p);
        const RIGID_BODY_MASS<TV,true> &rigid_mass=projection_data.mass.world_space_rigid_mass(p),&rigid_mass_inverse=projection_data.mass.world_space_rigid_mass_inverse(p);
        convergence_norm_squared=max(convergence_norm_squared,
            twist.linear.Magnitude_Squared()+rigid_mass_inverse.mass*Dot_Product(twist.angular,rigid_mass.inertia_tensor*twist.angular));}
    T convergence_norm=sqrt(convergence_norm_squared);
    if(mpi_solids) convergence_norm=mpi_solids->Reduce_Max(convergence_norm);
    return convergence_norm;
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void BACKWARD_EULER_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const
{
    if(++projection_data.project_nullspace_counter>=project_nullspace_frequency){Project(V);projection_data.project_nullspace_counter=0;}
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void BACKWARD_EULER_SYSTEM<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const
{
    R.Copy((T)1,V);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GENERALIZED_MASS<TV>::
GENERALIZED_MASS(SOLID_BODY_COLLECTION<TV>& solid_body_collection)
    :mass(solid_body_collection.deformable_body_collection.particles.mass,solid_body_collection.deformable_body_collection.dynamic_particles),
    one_over_mass(solid_body_collection.deformable_body_collection.particles.one_over_mass,solid_body_collection.deformable_body_collection.dynamic_particles),
    rigid_mass(solid_body_collection.rigid_body_collection.rigid_body_particle.mass,solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles),
    rigid_inertia_tensor(solid_body_collection.rigid_body_collection.rigid_body_particle.inertia_tensor,solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles),
    world_space_rigid_mass_inverse(world_space_rigid_mass_inverse_full,rigid_inertia_tensor.indices)
{
    Initialize_World_Space_Masses(solid_body_collection);
}
template<class TV> GENERALIZED_MASS<TV>::
~GENERALIZED_MASS()
{
}
//#####################################################################
// Function Initialize_World_Space_Masses
//#####################################################################
template<class TV> void GENERALIZED_MASS<TV>::
Initialize_World_Space_Masses(const SOLID_BODY_COLLECTION<TV>& solid_body_collection)
{
    world_space_rigid_mass.Resize(solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m,false,false);
    world_space_rigid_mass_inverse_full.Resize(rigid_mass.array.m,false,false);

    for(int i=1;i<=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles(i);
        RIGID_BODY_MASS<TV> rm(rigid_mass.array(i),rigid_inertia_tensor.array(i));
        world_space_rigid_mass(i)=solid_body_collection.rigid_body_collection.State(p).World_Space_Rigid_Mass(rm);
        world_space_rigid_mass_inverse(i)=solid_body_collection.rigid_body_collection.State(i).World_Space_Rigid_Mass_Inverse(rm);}
}
//#####################################################################
// Function Inverse_Multiply
//#####################################################################
template<class TV> void GENERALIZED_MASS<TV>::
Inverse_Multiply(const GENERALIZED_VELOCITY<TV>& V,GENERALIZED_VELOCITY<TV>& F,bool include_static) const
{
    F.V=one_over_mass*V.V;
    F.rigid_V=world_space_rigid_mass_inverse*V.rigid_V;
    if(include_static) F.kinematic_and_static_rigid_V*=0;
}
//#####################################################################
template class GENERALIZED_MASS<VECTOR<float,1> >;
template class GENERALIZED_MASS<VECTOR<float,2> >;
template class GENERALIZED_MASS<VECTOR<float,3> >;
template class BACKWARD_EULER_SYSTEM<VECTOR<float,1> >;
template class BACKWARD_EULER_SYSTEM<VECTOR<float,2> >;
template class BACKWARD_EULER_SYSTEM<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GENERALIZED_MASS<VECTOR<double,1> >;
template class GENERALIZED_MASS<VECTOR<double,2> >;
template class GENERALIZED_MASS<VECTOR<double,3> >;
template class BACKWARD_EULER_SYSTEM<VECTOR<double,1> >;
template class BACKWARD_EULER_SYSTEM<VECTOR<double,2> >;
template class BACKWARD_EULER_SYSTEM<VECTOR<double,3> >;
#endif

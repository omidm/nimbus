//#####################################################################
// Copyright 2006-2009, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NEWMARK_EVOLUTION
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Krylov_Solvers/LANCZOS_ITERATION.h>
#include <PhysBAM_Tools/Krylov_Solvers/SYMMQMR.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_TRIPLE.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_MATRIX.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BINDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/INCOMPRESSIBLE_FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_1D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGID_LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGIDS_NEWMARK_COLLISION_CALLBACKS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/ASYNCHRONOUS_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <stdexcept>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> NEWMARK_EVOLUTION<TV>::
NEWMARK_EVOLUTION(SOLIDS_PARAMETERS<TV>& solids_parameters_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input)
    :SOLIDS_EVOLUTION<TV>(solids_parameters_input,solid_body_collection_input),rigids_evolution_callbacks(*new RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>(*this)),repulsions(0),
    use_existing_contact(false),asynchronous_evolution(0),print_matrix(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> NEWMARK_EVOLUTION<TV>::
~NEWMARK_EVOLUTION()
{
    delete &rigids_evolution_callbacks;
}
//#####################################################################
// Function Prepare_Backward_Euler_System
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Prepare_Backward_Euler_System(BACKWARD_EULER_SYSTEM<TV>& system,const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update)
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=rigid_body_collection.rigid_body_particle;
    MPI_SOLIDS<TV>* mpi_solids=solid_body_collection.deformable_body_collection.mpi_solids;
    ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body=solid_body_collection.rigid_body_collection.articulated_rigid_body; // Needn't be a pointer
    rigid_body_collection.Update_Angular_Velocity(); // make sure omega = I^{-1} L

    F_full.Resize(particles.array_collection->Size(),false,false);rigid_F_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
    R_full.Resize(particles.array_collection->Size(),false,false);rigid_R_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
    S_full.Resize(particles.array_collection->Size(),false,false);rigid_S_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
    B_full.Resize(particles.array_collection->Size(),false,false);rigid_B_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
    if(solids_parameters.implicit_solve_parameters.evolution_solver_type!=krylov_solver_cg){
        AR_full.Resize(particles.array_collection->Size(),false,false);rigid_AR_full.Resize(rigid_body_particles.array_collection->Size(),false,false);}

    GENERALIZED_VELOCITY<TV> B_all(B_full,rigid_B_full,solid_body_collection);
    GENERALIZED_VELOCITY<TV> F_all(F_full,rigid_F_full,solid_body_collection);
    ARRAY<TWIST<TV> > twist;twist.Resize(rigid_body_particles.V.Size());
    for(int i=1;i<=twist.Size();i++) twist(i)=TWIST<TV>(rigid_body_particles.V(i),rigid_body_particles.angular_velocity(i));
    GENERALIZED_VELOCITY<TV> V_all(particles.V,twist,solid_body_collection);

    INDIRECT_ARRAY<ARRAY<TV>,ARRAY<int>&> B_subset=B_full.Subset(solid_body_collection.deformable_body_collection.simulated_particles);
    ARRAYS_COMPUTATIONS::Fill(B_subset,TV());ARRAYS_COMPUTATIONS::Fill(rigid_B_full,TWIST<TV>());
    solid_body_collection.example_forces_and_velocities->Add_External_Forces(B_full,current_velocity_time+dt);
    solid_body_collection.example_forces_and_velocities->Add_External_Forces(rigid_B_full,current_velocity_time+dt);
    if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data_Global(particles.V);
    if(!solids_parameters.set_velocity_from_positions) solid_body_collection.Add_Velocity_Independent_Forces(B_full,rigid_B_full,current_velocity_time+dt); // this is a nop for binding forces
    if(mpi_solids) mpi_solids->Exchange_Binding_Boundary_Data_Global(B_full);
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Force_To_Parents(B_full,rigid_B_full);
    solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Distribute_Force_To_Parents(rigid_B_full);

    if(solid_body_collection.deformable_body_collection.soft_bindings.Need_Bindings_Mapped()){
        if(mpi_solids) mpi_solids->Exchange_Binding_Boundary_Data_Global(B_full);
        solid_body_collection.deformable_body_collection.soft_bindings.Map_Forces_From_Parents(B_full,rigid_B_full);
        solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(B_full);
        if(mpi_solids) mpi_solids->Exchange_Binding_Boundary_Data_Global(B_full);
        for(int k=1;k<=solid_body_collection.solids_forces.m;k++) if(dynamic_cast<BINDING_SPRINGS<TV>*>(&*solid_body_collection.solids_forces(k)))
            solid_body_collection.solids_forces(k)->Add_Force_Differential(particles.X,B_full,current_velocity_time+dt);
        if(mpi_solids) mpi_solids->Exchange_Binding_Boundary_Data_Global(B_full);
        solid_body_collection.deformable_body_collection.binding_list.Distribute_Force_To_Parents(B_full,rigid_B_full);}

    Initialize_World_Space_Masses();
    if(articulated_rigid_body.constrain_pd_directions){
        for(int i=1;i<=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles(i);rigid_B_full(p)=world_space_rigid_mass_inverse(p)*rigid_B_full(p);}
        articulated_rigid_body.Poststabilization_Projection(rigid_B_full,true);
        for(int i=1;i<=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles(i);rigid_B_full(p)=world_space_rigid_mass(p)*rigid_B_full(p);}}

    if(asynchronous_evolution) asynchronous_evolution->Project(B_all);
    for(int i=1;i<=solid_body_collection.deformable_body_collection.dynamic_particles.m;i++){int p=solid_body_collection.deformable_body_collection.dynamic_particles(i);
        B_full(p)=particles.V(p)+dt*particles.one_over_mass(p)*B_full(p);}
    for(int i=1;i<=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles(i);
        rigid_B_full(p)=rigid_body_particles.rigid_geometry(p)->Twist()+world_space_rigid_mass_inverse(p)*rigid_B_full(p)*dt;}

    if(solids_parameters.set_velocity_from_positions_physbam){
        if(!velocity_update){
            // Compute deltaPE due to just the elastic forces
            V_postelasticforces.Resize(solid_body_collection.deformable_body_collection.particles.array_collection->Size());
            for(int i=1;i<=solid_body_collection.deformable_body_collection.dynamic_particles.m;i++){
                int p=solid_body_collection.deformable_body_collection.dynamic_particles(i);
                particles.V(p)=B_full(p);
                V_postelasticforces(p)=particles.V(p);}
            Euler_Step_Position(dt*2,time);
            postelasticforces_PE=0;
            for(int i=1;i<=solid_body_collection.deformable_body_collection.deformables_forces.m;i++)
                postelasticforces_PE+=solid_body_collection.deformable_body_collection.deformables_forces(i)->Potential_Energy(time);
            // Store just delta PE from elastic
            for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
                DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
                force.Store_Delta_PE(time);}
            Restore_Position(X_save,rigid_X_save,rigid_rotation_save);
            Restore_Velocity();}
        else{
            // Compute deltaKE due to just the elastic forces
            T deltaKE_elastic=0;
            for(int i=1;i<=solid_body_collection.deformable_body_collection.dynamic_particles.m;i++){
                int p=solid_body_collection.deformable_body_collection.dynamic_particles(i);
                deltaKE_elastic+=((T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*TV::Dot_Product(B_full(p),B_full(p))-
                    (T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*TV::Dot_Product(particles.V(p),particles.V(p)));}
            T deltaPE_elastic=0;
            for(int i=1;i<=solid_body_collection.deformable_body_collection.deformables_forces.m;i++)
                deltaPE_elastic+=solid_body_collection.deformable_body_collection.deformables_forces(i)->Get_Total_Delta_PE();
            energy_damped+=-(deltaKE_elastic+deltaPE_elastic);
            
            T A=0,B=0,C=0;
            ARRAY<TV> forces(solid_body_collection.deformable_body_collection.particles.array_collection->Size());
            if(energy_damped > 0) {
                for(int i=1;i<=solid_body_collection.deformable_body_collection.dynamic_particles.m;i++){
                    int p=solid_body_collection.deformable_body_collection.dynamic_particles(i);
                    forces(p)=TV();
                    for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
                        DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
                        ARRAY<int>& incident_force_elements=*force.Incident_Force_Elements(p);
                        for(int ife=1;ife<=incident_force_elements.m;ife++)
                            forces(p)+=dt*solid_body_collection.deformable_body_collection.particles.one_over_mass(p)*force.Get_Force(incident_force_elements(ife),p,true);}
                    B+=solid_body_collection.deformable_body_collection.particles.mass(p)*TV::Dot_Product(B_full(p),forces(p));
                    A+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*TV::Dot_Product(forces(p),forces(p));}}
            else {
                for(int i=1;i<=solid_body_collection.deformable_body_collection.dynamic_particles.m;i++){
                    int p=solid_body_collection.deformable_body_collection.dynamic_particles(i);
                    forces(p)=TV();
                    for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
                        DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
                        force.Get_Damping_Force(p,forces(p),dt,false);}
                    B+=solid_body_collection.deformable_body_collection.particles.mass(p)*TV::Dot_Product(B_full(p),forces(p));
                    A+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*TV::Dot_Product(forces(p),forces(p));}}
            C=-energy_damped;
            T epsilon;
            T discriminant=B*B-4*A*C;
            if(discriminant<0) discriminant=0;
            if(A==0) epsilon=0;
            else if(B>0) epsilon=(-B+sqrt(discriminant))/(2*A);
            else epsilon=(-B-sqrt(discriminant))/(2*A);
            
            for(int i=1;i<=solid_body_collection.deformable_body_collection.dynamic_particles.m;i++){
                int p=solid_body_collection.deformable_body_collection.dynamic_particles(i);
                TV dv=epsilon*forces(p);
                energy_damped-=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*(2*TV::Dot_Product(B_full(p),dv)+TV::Dot_Product(dv,dv));
                B_full(p)+=solids_parameters.set_velocity_from_positions_percent_energy_recovered*dv;}

            // Save post elastic forces KE
            V_postelasticforces.Resize(solid_body_collection.deformable_body_collection.particles.array_collection->Size());
            postelasticforces_KE=0;
            for(int i=1;i<=solid_body_collection.deformable_body_collection.dynamic_particles.m;i++){
                int p=solid_body_collection.deformable_body_collection.dynamic_particles(i);
                postelasticforces_KE+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*TV::Dot_Product(B_full(p),B_full(p));
                V_postelasticforces(p)=B_full(p);}}}

    V_all=B_all;
    rigid_body_collection.Update_Angular_Momentum();
    Diagnostics(dt,current_position_time,0,0,604,"Before boundary conditions");
    system.Set_Global_Boundary_Conditions(V_all,X_save,rigid_X_save,rigid_rotation_save,rigid_velocity_save,rigid_angular_momentum_save,V_save,
        solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);
    if(solids_parameters.set_velocity_from_positions_physbam && velocity_update && !solids_parameters.use_post_cg_constraints){
        T prerepulsions_KE=0;
        if(solids_parameters.set_velocity_from_positions_physbam){
            for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++)
                prerepulsions_KE+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*
                    TV::Dot_Product(solid_body_collection.deformable_body_collection.particles.V(p),solid_body_collection.deformable_body_collection.particles.V(p));}
        Apply_Constraints(dt,current_velocity_time);
        if(solids_parameters.set_velocity_from_positions_physbam){
            T postrepulsions_KE=0;
            for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++)
                postrepulsions_KE+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*
                    TV::Dot_Product(solid_body_collection.deformable_body_collection.particles.V(p),solid_body_collection.deformable_body_collection.particles.V(p));
            energy_lost_to_bodies+=(postrepulsions_KE-prerepulsions_KE);
            if(solids_parameters.set_velocity_from_positions_conserve_exactly) energy_damped+=-(postrepulsions_KE-prerepulsions_KE);
            V_postelasticforces.Resize(solid_body_collection.deformable_body_collection.particles.array_collection->Size());
            postelasticforces_KE=0;
            for(int i=1;i<=solid_body_collection.deformable_body_collection.dynamic_particles.m;i++){
                int p=solid_body_collection.deformable_body_collection.dynamic_particles(i);
                postelasticforces_KE+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*TV::Dot_Product(particles.V(p),particles.V(p));
                V_postelasticforces(p)=particles.V(p);}}}

    // TODO: This is completely broken for trapezoid and likely for BE as well.
    if(solid_body_collection.rigid_body_collection.articulated_rigid_body.Has_Actuators() && solid_body_collection.rigid_body_collection.articulated_rigid_body.constrain_pd_directions){
        for(int i=1;i<=saved_pd.Size();i++) saved_pd(i)=rigid_body_particles.rigid_geometry(i)->Twist();
        solid_body_collection.rigid_body_collection.articulated_rigid_body.Poststabilization_Projection(saved_pd,true);
        for(int i=1;i<=saved_pd.Size();i++) saved_pd(i)=rigid_body_particles.rigid_geometry(i)->Twist()-saved_pd(i);}
    if(!velocity_update) solid_body_collection.example_forces_and_velocities->Add_External_Impulses_Before(B_full,current_position_time,(T)2*dt); // 2*dt is position dt TODO: what time?
    
    for(int i=1;i<=twist.Size();i++){rigid_body_particles.V(i)=twist(i).linear;rigid_body_particles.angular_velocity(i)=twist(i).angular;}
}
//#####################################################################
// Function Finish_Backward_Euler_Step
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Finish_Backward_Euler_Step(KRYLOV_SYSTEM_BASE<T>& system,const T dt,const T current_position_time,const bool velocity_update)
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=rigid_body_collection.rigid_body_particle;
    GENERALIZED_VELOCITY<TV> F_all(F_full,rigid_F_full,solid_body_collection);
    rigid_deformable_collisions->rigid_body_collisions.Remove_Contact_Joints(); // TODO: Fix me.

    if(velocity_update && solids_parameters.use_post_cg_constraints){ // return rhs + dt Fd V^n+1 for friction processing
        if(solids_parameters.no_contact_friction){
            /*ARRAY<TV> V_no_projection;V_no_projection.Resize(particles.array_collection->Size(),false,false);
            ARRAY<TV> V_projected;V_projected.Resize(particles.array_collection->Size(),false,false);
            ARRAY<TWIST<TV> > rigid_V_no_projection;rigid_V_no_projection.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
            ARRAY<TWIST<TV> > rigid_V_projected;rigid_V_projected.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
            for(int i=1;i<=solid_body_collection.deformable_body_collection.dynamic_particles.m;i++){
                int p=solid_body_collection.deformable_body_collection.dynamic_particles(i);
                V_no_projection(p)=B_full(p)+dt*particles.one_over_mass(p)*F_full(p);}
            for(int i=1;i<=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m;i++){
                int p=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles(i);
                rigid_V_no_projection(p)=rigid_B_full(p)+world_space_rigid_mass_inverse(p)*rigid_F_full(p)*dt;}
            V_projected=V_no_projection;
            GENERALIZED_VELOCITY<TV> V_projected_all(V_projected,rigid_V_projected,solid_body_collection);
            system.Project(V_projected_all);*/}
        else{
            if(asynchronous_evolution) asynchronous_evolution->Project(F_all);
            for(int i=1;i<=solid_body_collection.deformable_body_collection.dynamic_particles.m;i++){int p=solid_body_collection.deformable_body_collection.dynamic_particles(i);
                particles.V(p)=B_full(p)+dt*particles.one_over_mass(p)*F_full(p);}
            for(int i=1;i<=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles(i);
                TWIST<TV> twist=rigid_B_full(p)+world_space_rigid_mass_inverse(p)*rigid_F_full(p)*dt;
                rigid_body_particles.V(p)=twist.linear;rigid_body_particles.angular_velocity(p)=twist.angular;}
            ARRAY<TWIST<TV> > twist;twist.Resize(rigid_body_particles.V.Size());
            for(int i=1;i<=twist.Size();i++){twist(i)=TWIST<TV>(rigid_body_particles.V(i),rigid_body_particles.angular_velocity(i));}
            if(solid_body_collection.rigid_body_collection.articulated_rigid_body.Has_Actuators() && solid_body_collection.rigid_body_collection.articulated_rigid_body.constrain_pd_directions){
                solid_body_collection.rigid_body_collection.articulated_rigid_body.Poststabilization_Projection(twist,true);twist+=saved_pd;}
            for(int i=1;i<=twist.Size();i++){rigid_body_particles.V(i)=twist(i).linear;rigid_body_particles.angular_velocity(i)=twist(i).angular;}
            Diagnostics(dt,current_position_time,0,0,610,"After undo projections");}}
    rigid_body_collection.Update_Angular_Momentum();

    if(!velocity_update) solid_body_collection.example_forces_and_velocities->Add_External_Impulses(particles.V,current_position_time,(T)2*dt); // 2*dt is for

    solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Velocities();
    solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities(); // TODO: MPI safe?

    rigid_body_collection.Update_Angular_Momentum(); // make sure L = I omega
}
//#####################################################################
// Function Backward_Euler_Step_Velocity_Helper
//#####################################################################
// assumes all solids_forces are linear in velocity, with a symmetric positive definite Jacobian.
template<class TV> void NEWMARK_EVOLUTION<TV>::
Backward_Euler_Step_Velocity_Helper(const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update)
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=rigid_body_collection.rigid_body_particle;
    MPI_SOLIDS<TV>* mpi_solids=solid_body_collection.deformable_body_collection.mpi_solids;
    BACKWARD_EULER_SYSTEM<TV> system(*this,solid_body_collection,dt,current_velocity_time,current_position_time,&solid_body_collection.rigid_body_collection.articulated_rigid_body,
        (velocity_update && solids_parameters.enforce_repulsions_in_cg)?repulsions:0,mpi_solids,velocity_update);

    Prepare_Backward_Euler_System(system,dt,current_velocity_time,current_position_time,velocity_update);

    static CONJUGATE_GRADIENT<T> cg;
    static CONJUGATE_RESIDUAL<T> cr;
    static SYMMQMR<T> symmqmr;
    KRYLOV_SOLVER<T>* solver=0;
    const char* solver_name=0;
    if(solids_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_cg){solver=&cg;solver_name="CG";}
    else if(solids_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_cr){solver=&cr;solver_name="CONJUGATE_RESIDUAL";}
    else if(solids_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_symmqmr){solver=&symmqmr;solver_name="SYMMQMR";}
    solver->print_diagnostics=solid_body_collection.print_diagnostics;
    solver->print_residuals=solid_body_collection.print_residuals;
    solver->iterations_used=&solid_body_collection.iterations_used_diagnostic;
    solver->restart_iterations=solids_parameters.implicit_solve_parameters.cg_restart_iterations;
    system.project_nullspace_frequency=solids_parameters.implicit_solve_parameters.project_nullspace_frequency;
    ARRAY<TWIST<TV> > twist;twist.Resize(rigid_body_particles.V.Size());
    for(int i=1;i<=twist.Size();i++){twist(i)=TWIST<TV>(rigid_body_particles.V(i),rigid_body_particles.angular_velocity(i));}
    GENERALIZED_VELOCITY<TV> V(particles.V,twist,solid_body_collection),F(F_full,rigid_F_full,solid_body_collection),
        R(R_full,rigid_R_full,solid_body_collection),S(S_full,rigid_S_full,solid_body_collection),B(B_full,rigid_B_full,solid_body_collection),
        AR(AR_full,rigid_AR_full,solid_body_collection);
    GENERALIZED_MASS<TV> mass(solid_body_collection); // TODO: Doing duplicate computation of mass.

    if(solids_parameters.implicit_solve_parameters.spectral_analysis){V=B;LANCZOS_ITERATION<T>::Print_Spectral_Information(system,V,F,S,(T)1e-2*solids_parameters.implicit_solve_parameters.cg_tolerance,solids_parameters.implicit_solve_parameters.lanczos_iterations);}

    LOG::Time(solver_name);
    Diagnostics(dt,current_position_time,0,0,606,"Before solve");
    static int solve_id=0;solve_id++;
    if(print_matrix){
        {std::stringstream ss;ss<<"solve id "<<solve_id<<std::endl;LOG::filecout(ss.str());}
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("M-%i.txt",solve_id).c_str()).Write("M",system,S,R);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("P-%i.txt",solve_id).c_str()).Write_Projection("P",system,S);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("b-%i.txt",solve_id).c_str()).Write("b",B);}
    if(!solver->Solve(system,V,B,F,S,R,AR,solids_parameters.implicit_solve_parameters.cg_tolerance,1,solids_parameters.implicit_solve_parameters.cg_iterations) && solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure)
        throw std::runtime_error("Backward Euler Failed");
    if(print_matrix) OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("x-%i.txt",solve_id).c_str()).Write("x",V);
    Diagnostics(dt,current_position_time,0,0,607,"After solve");
    LOG::Stop_Time();

    if(velocity_update && solids_parameters.use_post_cg_constraints) system.Force(V,F);

    Finish_Backward_Euler_Step(system,dt,current_position_time,velocity_update);
    
    if(solids_parameters.enforce_energy_conservation)
        solid_body_collection.Compute_Previously_Applied_Forces();
    
    for(int i=1;i<=twist.Size();i++){rigid_body_particles.V(i)=twist(i).linear;rigid_body_particles.angular_velocity(i)=twist(i).angular;}
}
//#####################################################################
// Function Average_And_Exchange_Position
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Average_And_Exchange_Position()
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    assert(X_save.m==particles.array_collection->Size());
    assert(rigid_X_save.m==rigid_body_collection.rigid_body_particle.array_collection->Size() && rigid_rotation_save.m==rigid_body_collection.rigid_body_particle.array_collection->Size());
    const ARRAY<int>& simulated_particles=solid_body_collection.deformable_body_collection.simulated_particles;
    const ARRAY<int>& simulated_rigid_body_particles=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles;
    INDIRECT_ARRAY<ARRAY<TV> > X_save_fragment(X_save,simulated_particles);
    INDIRECT_ARRAY<ARRAY_VIEW<TV> > X_fragment(particles.X,simulated_particles);
    for(int i=1;i<=X_fragment.Size();i++){TV X_average=(T).5*(X_fragment(i)+X_save_fragment(i));X_save_fragment(i)=X_fragment(i);X_fragment(i)=X_average;}
    ARRAY<int> rigid_body_indices(simulated_rigid_body_particles);rigid_body_indices.Append_Elements(solid_body_collection.rigid_body_collection.kinematic_rigid_bodies);
    for(int i=1;i<=rigid_body_indices.Size();i++){int p=rigid_body_indices(i);
        TV tmp_X=TV::Interpolate(rigid_body_collection.rigid_body_particle.X(p),rigid_X_save(p),(T).5);
        rigid_X_save(p)=rigid_body_collection.rigid_body_particle.X(p);
        rigid_body_collection.rigid_body_particle.X(p)=tmp_X;
        ROTATION<TV> tmp_rotation=ROTATION<TV>::Spherical_Linear_Interpolation(rigid_body_collection.rigid_body_particle.rotation(p),rigid_rotation_save(p),(T).5);
        rigid_rotation_save(p)=rigid_body_collection.rigid_body_particle.rotation(p);
        rigid_body_collection.rigid_body_particle.rotation(p)=tmp_rotation;}
    for(int i=1;i<=rigid_body_indices.m;i++) rigid_body_collection.Rigid_Body(rigid_body_indices(i)).Update_Angular_Velocity();
}
//#####################################################################
// Function Trapezoidal_Step_Velocity
//#####################################################################
// assumes X is fixed at time+dt/2
template<class TV> void NEWMARK_EVOLUTION<TV>::
Trapezoidal_Step_Velocity(const T dt,const T time)
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    // save V at time
    Save_Velocity();
    // update V implicitly to time+dt/2
    Backward_Euler_Step_Velocity_Helper(dt/2,time,time+dt/2,true);
    // set up rigid_V_save for extrapolation step
    rigid_V_save.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
    for(int i=1;i<=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles(i);
        rigid_V_save(p).linear=rigid_velocity_save(p).linear;
        rigid_V_save(p).angular=rigid_body_collection.Rigid_Body(p).World_Space_Inertia_Tensor_Inverse_Times(rigid_angular_momentum_save(p));}
    // Use V_n instead of V_save rather than copying it around another time.  Also simplifies state dependencies.
    // extrapolate V to time+dt based on V at time and time+dt/2
    ARRAY<TWIST<TV> > twist;twist.Resize(rigid_body_collection.rigid_body_particle.V.Size());
    for(int i=1;i<=twist.Size();i++) twist(i)=TWIST<TV>(rigid_body_collection.rigid_body_particle.V(i),rigid_body_collection.rigid_body_particle.angular_velocity(i));
    GENERALIZED_VELOCITY<TV> V_n(V_save,rigid_V_save,solid_body_collection),V(particles.V,twist,solid_body_collection);
    V*=(T)2;V-=V_n;
    for(int i=1;i<twist.Size();i++){rigid_body_collection.rigid_body_particle.V(i)=twist(i).linear;rigid_body_collection.rigid_body_particle.angular_velocity(i)=twist(i).angular;}

    // enforce boundary conditions again
    if(solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions) solid_body_collection.deformable_body_collection.collisions.Activate_Collisions(false);
    Set_External_Velocities(particles.V,time+dt,time+dt/2);
    kinematic_evolution.Set_External_Velocities(rigid_body_collection.rigid_body_particle.V,rigid_body_collection.rigid_body_particle.angular_velocity,time+dt,time+dt/2);
    solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Velocities();
    solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities();
    rigid_body_collection.Update_Angular_Momentum(rigid_body_collection.simulated_rigid_body_particles);
    rigid_body_collection.Update_Angular_Momentum(solid_body_collection.rigid_body_collection.kinematic_rigid_bodies);
}
//#####################################################################
// Function Backward_Euler_Step_Velocity
//#####################################################################
// assumes X is fixed at time+dt
template<class TV> void NEWMARK_EVOLUTION<TV>::
Backward_Euler_Step_Velocity(const T dt,const T time)
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    Backward_Euler_Step_Velocity_Helper(dt,time,time+dt,true);
    // enforce boundary conditions again
    if(solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions) solid_body_collection.deformable_body_collection.collisions.Activate_Collisions(false);
    Set_External_Velocities(particles.V,time+dt,time+dt);
    rigid_body_collection.Update_Angular_Velocity(rigid_body_collection.simulated_rigid_body_particles);
    kinematic_evolution.Set_External_Velocities(rigid_body_collection.rigid_body_particle.V,rigid_body_collection.rigid_body_particle.angular_velocity,time+dt,time+dt);
    
    // enforce boundary conditions again
    solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Velocities();
    solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities();
    rigid_body_collection.Update_Angular_Momentum(rigid_body_collection.simulated_rigid_body_particles);
}
//#####################################################################
// Function Make_Incompressible
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Make_Incompressible(const T dt,const bool correct_volume)
{
    for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++)
        if(INCOMPRESSIBLE_FINITE_VOLUME_BASE<TV>* fvm=dynamic_cast<INCOMPRESSIBLE_FINITE_VOLUME_BASE<TV>*>(&*solid_body_collection.deformable_body_collection.deformables_forces(f))){
            fvm->Set_Neumann_Boundary_Conditions(&solid_body_collection.deformable_body_collection.collisions.particle_states,repulsions);
            fvm->Make_Incompressible(dt,correct_volume);}
}
//#####################################################################
// Function Advance_One_Time_Step_Position
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Advance_One_Time_Step_Position(const T dt,const T time, const bool solids)
{
    PHYSBAM_ASSERT(solids_parameters.use_post_cg_constraints || !solids_parameters.rigid_body_collision_parameters.enforce_rigid_rigid_contact_in_cg || solids_parameters.set_velocity_from_positions_physbam);
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Advance_One_Time_Step_Position Start dt=%f time=%f",dt,time),2,2);
    MPI_SOLIDS<TV>* mpi_solids=solid_body_collection.deformable_body_collection.mpi_solids;
    EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities=*solid_body_collection.example_forces_and_velocities;
    ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body=solid_body_collection.rigid_body_collection.articulated_rigid_body; // Needn't be a pointer
    const bool advance_rigid_bodies=true; //solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m!=0;  TODO: Fix this.
    const T v_dt=solids_parameters.use_backward_euler_position_update?dt:dt/2;

    if(solids_parameters.enforce_energy_conservation || solids_parameters.set_velocity_from_positions_physbam)
        solid_body_collection.Save_Potential_Energy(time);

    if((solids_parameters.triangle_collision_parameters.perform_self_collision && (solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.m || solids_parameters.triangle_collision_parameters.initialize_collisions_without_objects))
       && (!repulsions && solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions))
        repulsions=&solid_body_collection.deformable_body_collection.triangle_repulsions;

    Diagnostics(dt,time,0,0,1,"begin integration");
    example_forces_and_velocities.Advance_One_Time_Step_Begin_Callback(dt,time);

    solids_evolution_callbacks->Update_Solids_Parameters(time);
    if(solids){
        rigid_body_collisions->Initialize_Data_Structures();

        // save position and velocity for later trapezoidal rule
        Save_Velocity();
        Save_Position(X_save,rigid_X_save,rigid_rotation_save);

        if(solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions) solid_body_collection.deformable_body_collection.collisions.Activate_Collisions(false);}

    if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data_Global(solid_body_collection.deformable_body_collection.particles.X);
    solid_body_collection.Update_Position_Based_State(time+dt,true);

    // get momentum difference for v^n -> v^{n+1/2} udpate
    if(articulated_rigid_body.Has_Actuators()) example_forces_and_velocities.Set_PD_Targets(dt,time);

    if(solids_parameters.set_velocity_from_positions){
        Set_Velocity_From_Positions_Position_Update(dt,time);
        Compute_Momentum_Differences();}
    else{
        if(asynchronous_evolution && asynchronous_evolution->Take_Full_Backward_Euler_Step_For_Position_Update()){
            Backward_Euler_Step_Velocity_Helper(dt,time,time,false);
            asynchronous_evolution->Position_Velocity_Update(dt,time);}
        else Backward_Euler_Step_Velocity_Helper(v_dt,time,time,false); // update V implicitly to time+dt/2

        if(solids_parameters.set_velocity_from_positions_physbam){
            Euler_Step_Position(dt,time);
            T postelasticanddampingforces_PE=0;
            for(int i=1;i<=solid_body_collection.deformable_body_collection.deformables_forces.m;i++)
                postelasticanddampingforces_PE+=solid_body_collection.deformable_body_collection.deformables_forces(i)->Potential_Energy(time);
            Restore_Position(X_save,rigid_X_save,rigid_rotation_save);
            // Get damping (no projection) force
            ARRAY<TV> damping_force(solid_body_collection.deformable_body_collection.particles.array_collection->Size());
            for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
                for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
                    DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
                    force.Get_Damping_Force(p,damping_force(p),v_dt,true);}}
            ARRAY<TV> temp_velocities(solid_body_collection.deformable_body_collection.particles.array_collection->Size());
            for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
                temp_velocities(p)=solid_body_collection.deformable_body_collection.particles.V(p);
                solid_body_collection.deformable_body_collection.particles.V(p)=V_postelasticforces(p)+damping_force(p);}
            Euler_Step_Position(dt,time);
            T postelasticanddampingnoprojectionforces_PE=0;
            for(int i=1;i<=solid_body_collection.deformable_body_collection.deformables_forces.m;i++)
                postelasticanddampingnoprojectionforces_PE+=solid_body_collection.deformable_body_collection.deformables_forces(i)->Potential_Energy(time);
            Restore_Position(X_save,rigid_X_save,rigid_rotation_save);
            for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++)
                solid_body_collection.deformable_body_collection.particles.V(p)=temp_velocities(p);
            T energy_lost_to_damping=postelasticanddampingnoprojectionforces_PE-postelasticforces_PE;
            T energy_lost_to_damping_and_projection=postelasticanddampingforces_PE-postelasticforces_PE;
            T energy_lost_to_projection=energy_lost_to_damping_and_projection-energy_lost_to_damping;
            energy_lost_to_bodies+=energy_lost_to_projection;
            if(solids_parameters.set_velocity_from_positions_conserve_exactly) energy_damped+=-energy_lost_to_projection;
            // Add the damping amount to energy damped to be added in next time around
            energy_damped+=-energy_lost_to_damping;}

        solid_body_collection.Store_Velocities();
        Diagnostics(dt,time,1,0,5,"backward Euler");
        
        if(solids_parameters.use_projections_in_position_update){
            Restore_Velocity();
            Apply_Projections_In_Position_Update(dt,time);
            Diagnostics(dt,time,1,0,6,"apply projections in position update");}
        
        if(solids_parameters.verbose) Print_Maximum_Velocities(time);
        Diagnostics(dt,time,1,0,7,"position step enforce energy conservation");
        if(!solids) return; // early exit for fluids only in parallel
        
        Make_Incompressible(dt,true); // adjust velocity to fix volume
        solids_evolution_callbacks->Filter_Velocities(dt,time+dt,false); // use time+dt since these velocities are used to step to time+dt
        if(repulsions){
            T prerepulsions_PE=0;
            if(solids_parameters.set_velocity_from_positions_physbam){
                Euler_Step_Position(dt,time);
                for(int i=1;i<=solid_body_collection.deformable_body_collection.deformables_forces.m;i++)
                    prerepulsions_PE+=solid_body_collection.deformable_body_collection.deformables_forces(i)->Potential_Energy(time);
                Restore_Position(X_save,rigid_X_save,rigid_rotation_save);}
            repulsions->Adjust_Velocity_For_Self_Repulsion_Using_History(dt,true,false);
            if(solids_parameters.set_velocity_from_positions_physbam){
                Euler_Step_Position(dt,time);
                T postrepulsions_PE=0;
                for(int i=1;i<=solid_body_collection.deformable_body_collection.deformables_forces.m;i++)
                    postrepulsions_PE+=solid_body_collection.deformable_body_collection.deformables_forces(i)->Potential_Energy(time);
                Restore_Position(X_save,rigid_X_save,rigid_rotation_save);
                energy_lost_to_bodies+=(postrepulsions_PE-prerepulsions_PE);
                if(solids_parameters.set_velocity_from_positions_conserve_exactly) energy_damped+=-(postrepulsions_PE-prerepulsions_PE);}}
        Compute_Momentum_Differences();

        // add collision impulses to time n velocities and save
        Euler_Step_Position(dt,time);
        Diagnostics(dt,time,1,2,10,"Euler step position");}
    Exchange_Velocity();
    Diagnostics(dt,time,0,2,11,"restore velocity");

    T precollisions_KE=0;
    TV precollisions_momentum;
    T precollisions_PE=0;
    if(solids_parameters.set_velocity_from_positions_physbam){
        for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
            precollisions_KE+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*
                TV::Dot_Product(solid_body_collection.deformable_body_collection.particles.V(p),solid_body_collection.deformable_body_collection.particles.V(p));
            precollisions_momentum+=solid_body_collection.deformable_body_collection.particles.mass(p)*solid_body_collection.deformable_body_collection.particles.V(p);}
        
        for(int i=1;i<=solid_body_collection.deformable_body_collection.deformables_forces.m;i++) 
            precollisions_PE+=solid_body_collection.deformable_body_collection.deformables_forces(i)->Potential_Energy(time);}

    Process_Collisions(dt,time,advance_rigid_bodies);
    Diagnostics(dt,time,0,2,12,"add elastic collisions");

    Restore_Position(X_save,rigid_X_save,rigid_rotation_save);
    Diagnostics(dt,time,0,0,13,"restore position");
    if(advance_rigid_bodies){
        articulated_rigid_body.Apply_Poststabilization(solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);
        Diagnostics(dt,time,0,0,14,"poststabilization");}
    Save_Velocity();

    // update positions, apply contact, arb prestabilization and push out
    if(solids_parameters.rigid_body_collision_parameters.perform_contact && advance_rigid_bodies) rigid_body_collisions->Compute_Contact_Graph(dt,time,&articulated_rigid_body);
    Update_Velocity_Using_Stored_Differences(v_dt,time);
    Diagnostics(dt,time,1,0,18,"update velocity using stored differences");
    if(advance_rigid_bodies){
        articulated_rigid_body.Apply_Poststabilization(solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);
        Diagnostics(dt,time,0,0,185,"poststabilization");
        if(articulated_rigid_body.Has_Actuators() && !articulated_rigid_body.constrain_pd_directions){
            articulated_rigid_body.Compute_Position_Based_State(dt,time);
            articulated_rigid_body.Solve_Velocities_for_PD(time,v_dt,solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);}
        Diagnostics(dt,time,1,0,19,"solve velocities for pd");}

    Update_Positions_And_Apply_Contact_Forces(dt,time,false);
    Diagnostics(dt,time,1,2,20,"contact, prestabilization");
    if(solids_parameters.rigid_body_collision_parameters.use_push_out){
        if(solids_parameters.rigid_body_collision_parameters.use_legacy_push_out) rigid_body_collisions->Process_Push_Out_Legacy(); 
        else rigid_deformable_collisions->Process_Push_Out();
        Diagnostics(dt,time,1,2,21,"push out");}

    // This is necessary since the velocity update below needs to use soft bound positions, but modifying positions within the velocity update is prohibited.
    if(repulsions) solid_body_collection.deformable_body_collection.soft_bindings.Clamp_Particles_To_Embedded_Positions(true);

    Restore_Velocity();
    Diagnostics(dt,time,0,2,22,"restore velocity");

    // Compute change in momentum and energy due to collisions
    if(solids_parameters.set_velocity_from_positions_physbam){
        T postcollisions_KE=0;
        TV postcollisions_momentum;
        for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
            postcollisions_KE+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*
                TV::Dot_Product(solid_body_collection.deformable_body_collection.particles.V(p),solid_body_collection.deformable_body_collection.particles.V(p));
            postcollisions_momentum+=solid_body_collection.deformable_body_collection.particles.mass(p)*solid_body_collection.deformable_body_collection.particles.V(p);}
        energy_lost_to_bodies+=postcollisions_KE-precollisions_KE;
        if(solids_parameters.set_velocity_from_positions_conserve_exactly) energy_damped+=-(postcollisions_KE-precollisions_KE);
        momentum_lost_to_bodies+=postcollisions_momentum-precollisions_momentum;
        T postcollisions_PE=0;
        for(int i=1;i<=solid_body_collection.deformable_body_collection.deformables_forces.m;i++)
            postcollisions_PE+=solid_body_collection.deformable_body_collection.deformables_forces(i)->Potential_Energy(time);
        energy_lost_to_bodies+=postcollisions_PE-precollisions_PE;
        if(solids_parameters.set_velocity_from_positions_conserve_exactly) energy_damped+=-(postcollisions_PE-precollisions_PE);
        if((postcollisions_KE-precollisions_KE+postcollisions_PE-precollisions_PE)>0){
            T A=0,B=0,C=0;
            ARRAY<TV> damping_forces(solid_body_collection.deformable_body_collection.particles.array_collection->Size());
            for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
                damping_forces(p)=TV();
                for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
                    DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
                    force.Get_Damping_Force(p,damping_forces(p),dt,false);}
                B+=solid_body_collection.deformable_body_collection.particles.mass(p)*TV::Dot_Product(solid_body_collection.deformable_body_collection.particles.V(p),damping_forces(p));
                A+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*TV::Dot_Product(damping_forces(p),damping_forces(p));}
            C=(postcollisions_KE-precollisions_KE+postcollisions_PE-precollisions_PE);
            T epsilon;
            T discriminant=B*B-4*A*C;
            if(discriminant<0) discriminant=0;
            if(B>0) epsilon=(-B+sqrt(discriminant))/(2*A);
            else epsilon=(-B-sqrt(discriminant))/(2*A);
            
            for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++)
                solid_body_collection.deformable_body_collection.particles.V(p)+=epsilon*damping_forces(p);

/*        T A=0,B=0,C=0;
        ARRAY<TV> elastic_forces(solid_body_collection.deformable_body_collection.particles.array_collection->Size());
        for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
            elastic_forces(p)=TV();
            for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
                DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
                ARRAY<int>& incident_force_elements=*force.Incident_Force_Elements(p);
                for(int ife=1;ife<=incident_force_elements.m;ife++)
                    elastic_forces(p)+=dt*solid_body_collection.deformable_body_collection.particles.one_over_mass(p)*force.Get_Force(incident_force_elements(ife),p,false);}
            B+=solid_body_collection.deformable_body_collection.particles.mass(p)*TV::Dot_Product(solid_body_collection.deformable_body_collection.particles.V(p),elastic_forces(p));
            A+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*TV::Dot_Product(elastic_forces(p),elastic_forces(p));}
        C=(postcollisions_KE-precollisions_KE+postcollisions_PE-precollisions_PE);
        T epsilon;
        T discriminant=B*B-4*A*C;
        if(discriminant<0) discriminant=0;
        if(B>0) epsilon=(-B+sqrt(discriminant))/(2*A);
        else epsilon=(-B-sqrt(discriminant))/(2*A);
        
        for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
        solid_body_collection.deformable_body_collection.particles.V(p)+=epsilon*elastic_forces(p);}*/}
        T postpostcollisions_KE=0;
        TV postpostcollisions_momentum;
        for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
            postpostcollisions_KE+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*
                TV::Dot_Product(solid_body_collection.deformable_body_collection.particles.V(p),solid_body_collection.deformable_body_collection.particles.V(p));
            postpostcollisions_momentum+=solid_body_collection.deformable_body_collection.particles.mass(p)*solid_body_collection.deformable_body_collection.particles.V(p);}
        energy_lost_to_bodies+=postpostcollisions_KE-postcollisions_KE;
        if(solids_parameters.set_velocity_from_positions_conserve_exactly) energy_damped+=-(postpostcollisions_KE-postcollisions_KE);}

    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Advance_One_Time_Step_Position End dt=%f time=%f",dt,time),2,2);
}
//#####################################################################
// Function Set_Velocity_From_Positions_Position_Update
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Set_Velocity_From_Positions_Position_Update(const T dt,const T time)
{
    solid_body_collection.deformable_body_collection.Setup_Set_Velocity_From_Positions(time,true,solids_parameters.set_velocity_from_positions_reset_alphas);
    
    ARRAY<ARRAY<T> > convergence_norms(solid_body_collection.deformable_body_collection.deformables_forces.m);
    int iter;
    for(iter=1;iter<=solids_parameters.set_velocity_from_positions_iterations;iter++){
        T convergence_norm=0;
        for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++)
            convergence_norm=max(convergence_norm,ARRAYS_COMPUTATIONS::Maxabs(convergence_norms(f)));
        if(iter>1 && convergence_norm<solids_parameters.set_velocity_from_positions_tolerance) break;
        
        // Iterate over all forces
        ARRAY<ARRAY<VECTOR<T,3> > > force_quadratic_coefficients(solid_body_collection.deformable_body_collection.deformables_forces.m);
        ARRAY<ARRAY<ARRAY<T> > > v_n_hat(solid_body_collection.deformable_body_collection.deformables_forces.m);
        for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
            DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
            force_quadratic_coefficients(f).Resize(force.Get_Element_Count());
            v_n_hat(f).Resize(force.Get_Element_Count());
            FORCE_ELEMENTS& force_elements=*force.Get_Force_Elements();
            for(typename FORCE_ELEMENTS::ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){
                int s=iterator.Data();
                // Update the alpha for this force
                // Get a list of the nodes involved in this force
                T a=0;T b=0;T c=0;
                T A=0,B=0,C=0;
                ARRAY<int>& nodes=*force.Incident_Nodes(s);
                // Get combined mass
                T combined_one_over_mass=force.Get_Combined_One_Over_Mass(s);
                
                force.Compute_Quadratic_Contribution_For_Force(A,a,c,dt,s,combined_one_over_mass,false);
                // Iterate over the nodes
                v_n_hat(f)(s).Resize(nodes.m);
                for(int n=1;n<=nodes.m;n++){
                    int node=nodes(n);
                    TV direction=force.Get_Direction(s);
                    // Compute modified inital velocity
                    T v_n_correction=0;
                    for(int f_incident=1;f_incident<=solid_body_collection.deformable_body_collection.deformables_forces.m;f_incident++){
                        DEFORMABLES_FORCES<TV>& force_incident=*solid_body_collection.deformable_body_collection.deformables_forces(f_incident);
                        ARRAY<int>& incident_force_elements=*force_incident.Incident_Force_Elements(node);
                        for(int ife=1;ife<=incident_force_elements.m;ife++){
                            if(f_incident==f && incident_force_elements(ife)==s) continue;
                            v_n_correction+=dt*solid_body_collection.deformable_body_collection.particles.one_over_mass(node)*
                                TV::Dot_Product(force_incident.Get_Force(incident_force_elements(ife),node,false),direction);}}
                    force.Compute_Quadratic_Contribution_For_Node(B,C,b,dt,node,s,combined_one_over_mass,v_n_correction,false);
                    v_n_hat(f)(s)(n)=TV::Dot_Product(solid_body_collection.deformable_body_collection.particles.V(node),direction)+(T).5*v_n_correction;}
                force.Compute_Quadratic_Contribution_For_Residual(B,C,a,b,c,dt,time,s,false);
                force_quadratic_coefficients(f)(s)=VECTOR<T,3>(A,B,C);}}
            
        // Solve for each force, and get an alpha
        for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
            DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
            FORCE_ELEMENTS& force_elements=*force.Get_Force_Elements();
            for(typename FORCE_ELEMENTS::ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){
                int s=iterator.Data();
                
                T A=force_quadratic_coefficients(f)(s)(1);
                T B=force_quadratic_coefficients(f)(s)(2);
                T C=force_quadratic_coefficients(f)(s)(3);
                
                if((A==(T)0) && (B==(T)0) && (C==(T)0)){
                    force.Set_Force(s,0);
                    continue;}

                T discriminant=B*B-4*A*C;
                if(discriminant<0){
//                    {std::stringstream ss;ss << "No solutions: negative discriminant position update" << std::endl;LOG::filecout(ss.str());}
//                    {std::stringstream ss;ss << "Negative discriminant in position update: " << discriminant << std::endl;LOG::filecout(ss.str());}
                    force.Set_Force(s,-B/(2*A));}
                else{
                    T alpha1=(-B+sqrt(discriminant))/(2*A);
                    T alpha2=(-B-sqrt(discriminant))/(2*A);
                    
                    force.Choose_Solution(solids_parameters.set_velocity_from_positions_use_orig_force,s,(T).5*dt,alpha1,alpha2,v_n_hat(f)(s));}}}
        
        // Iterate over all forces
        for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
            DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
            FORCE_ELEMENTS& force_elements=*force.Get_Force_Elements();
            convergence_norms(f).Resize(force.Get_Element_Count());
            for(typename FORCE_ELEMENTS::ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){
                int s=iterator.Data();
                // Update the alpha for this force
                // Get a list of the nodes involved in this force
                T a=0;T b=0;T c=0;
                T A=0,B=0,C=0;
                ARRAY<int>& nodes=*force.Incident_Nodes(s);
                // Get combined mass
                T combined_one_over_mass=force.Get_Combined_One_Over_Mass(s);
                
                force.Compute_Quadratic_Contribution_For_Force(A,a,c,dt,s,combined_one_over_mass,false);
                // Iterate over the nodes
                for(int n=1;n<=nodes.m;n++){
                    int node=nodes(n);
                    TV direction=force.Get_Direction(s);
                    // Compute modified inital velocity
                    T v_n_correction=0;
                    for(int f_incident=1;f_incident<=solid_body_collection.deformable_body_collection.deformables_forces.m;f_incident++){
                        DEFORMABLES_FORCES<TV>& force_incident=*solid_body_collection.deformable_body_collection.deformables_forces(f_incident);
                        ARRAY<int>& incident_force_elements=*force_incident.Incident_Force_Elements(node);
                        for(int ife=1;ife<=incident_force_elements.m;ife++){
                            if(f_incident==f && incident_force_elements(ife)==s) continue;
                            v_n_correction+=dt*solid_body_collection.deformable_body_collection.particles.one_over_mass(node)*
                                TV::Dot_Product(force_incident.Get_Force(incident_force_elements(ife),node,false),direction);}}
                    force.Compute_Quadratic_Contribution_For_Node(B,C,b,dt,node,s,combined_one_over_mass,v_n_correction,false);}
                force.Compute_Quadratic_Contribution_For_Residual(B,C,a,b,c,dt,time,s,false);
                T new_alpha=force.Get_Force(s);
                convergence_norms(f)(s)=A*new_alpha*new_alpha+B*new_alpha+C;}}}
    // Update positions based on velocities
    for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
        TV vtemp=solid_body_collection.deformable_body_collection.particles.V(p);
        if(!solid_body_collection.deformable_body_collection.particles.one_over_mass(p)) continue;
        for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
            DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
            ARRAY<int>& incident_force_elements=*force.Incident_Force_Elements(p);
            for(int ife=1;ife<=incident_force_elements.m;ife++){
                vtemp+=(T).5*dt*solid_body_collection.deformable_body_collection.particles.one_over_mass(p)*force.Get_Force(incident_force_elements(ife),p,false);}}
        solid_body_collection.deformable_body_collection.particles.V(p)=vtemp;
        solid_body_collection.deformable_body_collection.particles.X(p)=solid_body_collection.deformable_body_collection.particles.X(p)+vtemp*dt;}
    // Store just delta PE from evolution
    for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
        DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
        force.Store_Delta_PE(time);}
    
    Diagnostics(dt,time,0,2,24,"set velocity from positions: position update");
}
//#####################################################################
// Function Apply_Projections_In_Position_Update
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Apply_Projections_In_Position_Update(const T dt,const T time)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;

    // iterate over rigid/rigid
    if(rigid_body_collection.dynamic_rigid_body_particles.m && solids_parameters.rigid_body_collision_parameters.enforce_rigid_rigid_contact_in_cg){
        for(typename HASHTABLE<TRIPLE<int,int,TV> >::ITERATOR iterator(rigid_deformable_collisions->rigid_body_collisions.rigid_body_particle_intersections);iterator.Valid();iterator.Next()){
            const TRIPLE<int,int,TV>& intersection=iterator.Key();
            const RIGID_BODY<TV> &particle_body=rigid_body_collection.Rigid_Body(intersection.x),
                &levelset_body=rigid_body_collection.Rigid_Body(intersection.y);
            TV relative_velocity=RIGID_GEOMETRY<TV>::Relative_Velocity(particle_body,levelset_body,intersection.z);
            TV normal=levelset_body.implicit_object->Extended_Normal(intersection.z);
            if(TV::Dot_Product(relative_velocity,normal)>0) rigid_deformable_collisions->rigid_body_collisions.rigid_body_particle_intersections.Delete(iterator.Key());}}

    // iterate over rigid/deformable
    for(int i=1;i<=rigid_deformable_collisions->precompute_contact_projections.m;i++){
        typename RIGID_DEFORMABLE_COLLISIONS<TV>::PRECOMPUTE_CONTACT_PROJECTION& precompute=*rigid_deformable_collisions->precompute_contact_projections(i);
        const RIGID_BODY<TV>& body=precompute.rigid_body;
        for(int j=1;j<=precompute.particles.m;j++){const int p=precompute.particles(j);
            TV V_rel=body.Pointwise_Object_Velocity(particles.X(p))-particles.V(p);
            precompute.V_rel_target(j)=TV();
            if(TV::Dot_Product(V_rel,precompute.N(j))<0){
                precompute.particles.Remove_Index_Lazy(j);
                j--;}}}

    Backward_Euler_Step_Velocity_Helper(dt,time,time,true);
}
//#####################################################################
// Function Write_Position_Update_Projection_Data
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Write_Position_Update_Projection_Data(const STREAM_TYPE stream_type,const std::string& prefix)
{
    if(rigid_deformable_collisions->rigid_body_collisions.rigid_body_particle_intersections.Size())
        FILE_UTILITIES::Write_To_File(stream_type,prefix+"projection_data_rigid_rigid",rigid_deformable_collisions->rigid_body_collisions.rigid_body_particle_intersections);
    if(!rigid_deformable_collisions->precompute_contact_projections.m) return;
    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(prefix+"projection_data_rigid_deformable");
    TYPED_OSTREAM typed_output(*output,stream_type);
    Write_Binary(typed_output,rigid_deformable_collisions->precompute_contact_projections.m);
    for(int i=1;i<=rigid_deformable_collisions->precompute_contact_projections.m;i++){
        typename RIGID_DEFORMABLE_COLLISIONS<TV>::PRECOMPUTE_CONTACT_PROJECTION& p=*rigid_deformable_collisions->precompute_contact_projections(i);
        Write_Binary(typed_output,p.rigid_body.particle_index,p.particles,p.V_rel_target,p.N_over_NT_K_N,p.r,p.N,p.rN,p.A,p.A_inverted);}
    delete output;
}
//#####################################################################
// Function Read_Position_Update_Projection_Data
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Read_Position_Update_Projection_Data(const STREAM_TYPE stream_type,const std::string& prefix)
{
    if(FILE_UTILITIES::File_Exists(prefix+"projection_data_rigid_rigid"))
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"projection_data_rigid_rigid",rigid_deformable_collisions->rigid_body_collisions.rigid_body_particle_intersections);
    if(!FILE_UTILITIES::File_Exists(prefix+"projection_data_rigid_deformable")) return;
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(prefix+"projection_data_rigid_deformable");
    TYPED_ISTREAM typed_input(*input,stream_type);
    int precompute_contact_projections_size;
    Read_Binary(typed_input,precompute_contact_projections_size);
    for(int i=1;i<=precompute_contact_projections_size;i++){
        int index;
        Read_Binary(typed_input,index);
        typename RIGID_DEFORMABLE_COLLISIONS<TV>::PRECOMPUTE_CONTACT_PROJECTION* p=
            new typename RIGID_DEFORMABLE_COLLISIONS<TV>::PRECOMPUTE_CONTACT_PROJECTION(solid_body_collection.rigid_body_collection.Rigid_Body(index),false);
        Read_Binary(typed_input,p->particles,p->V_rel_target,p->N_over_NT_K_N,p->r,p->N,p->rN,p->A,p->A_inverted);}
    delete input;
}
//#####################################################################
// Function Process_Collisions
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Process_Collisions(const T dt,const T time,const bool advance_rigid_bodies)
{
    if(solids_parameters.use_rigid_deformable_contact)
        rigid_deformable_collisions->Add_Elastic_Collisions(dt,time,rigid_rotation_save,rigid_angular_momentum_difference,rigid_velocity_difference,rigid_X_save,rigid_velocity_save,
            rigid_angular_momentum_save,X_save,V_save);
    else if(solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions){
        if(solids_parameters.rigid_body_collision_parameters.perform_collisions && advance_rigid_bodies) rigid_body_collisions->Add_Elastic_Collisions(dt,time);
        int interactions=solid_body_collection.deformable_body_collection.collisions.Adjust_Nodes_For_Collision_Body_Collisions(solid_body_collection.deformable_body_collection.binding_list,
            solid_body_collection.deformable_body_collection.soft_bindings,X_save,dt,0);
        if(interactions) LOG::Stat("collision body collisions",interactions);}
}
//#####################################################################
// Function Advance_One_Time_Step_Velocity
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Advance_One_Time_Step_Velocity(const T dt,const T time,const bool solids)
{
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Advance_One_Time_Step_Velocity Start dt=%f time=%f",dt,time),2,2);

    MPI_SOLIDS<TV>* mpi_solids=solid_body_collection.deformable_body_collection.mpi_solids;
    EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities=*solid_body_collection.example_forces_and_velocities;
    ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    const bool advance_rigid_bodies=true; //solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m!=0;  TODO: Fix this.

    if(solids_parameters.set_velocity_from_positions){
        Set_Velocity_From_Positions_Velocity_Update(dt,time);
        return;}
        
    if(solids){
        if(advance_rigid_bodies){
            articulated_rigid_body.Apply_Poststabilization(solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);
            Diagnostics(dt,time,0,2,24,"poststabilization");}
        if(solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions) solid_body_collection.deformable_body_collection.collisions.Activate_Collisions(true);

        // initialize data needed for rigid/deformable contact projection in CG
        if(solids_parameters.use_rigid_deformable_contact) rigid_deformable_collisions->Initialize_All_Contact_Projections();
        PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("after creating joints.  before trapezoidal velocities dt=%f time=%f",dt,time),2,2);}

    if(solids_parameters.use_trapezoidal_rule_for_velocities){
        Average_And_Exchange_Position(); // move to positions at time+dt/2 for trapezoidal step
        Diagnostics(dt,time,0,1,25,"average and exchange positions");
        example_forces_and_velocities.Update_Time_Varying_Material_Properties(time+dt/2);

        if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data_Global(solid_body_collection.deformable_body_collection.particles.X);
        solid_body_collection.Update_Position_Based_State(time+dt/2,(solids_parameters.allow_altitude_spring_change_between_updates?true:false));
        Make_Incompressible(dt,false); // make velocity divergence free
        if(articulated_rigid_body.Has_Actuators()) example_forces_and_velocities.Set_PD_Targets(dt,time);
        Trapezoidal_Step_Velocity(dt,time);
        solids_evolution_callbacks->Filter_Velocities(dt,time+dt,true);
        Diagnostics(dt,time,2,1,29,"trazepoid rule");

        Restore_Position(X_save,rigid_X_save,rigid_rotation_save); // move to final positions at time time+dt
        Diagnostics(dt,time,2,2,30,"restore position");
        if(advance_rigid_bodies){
            articulated_rigid_body.Apply_Poststabilization(solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);
            Diagnostics(dt,time,2,2,31,"poststabilization");}
        example_forces_and_velocities.Update_Time_Varying_Material_Properties(time+dt);}
    else{
        example_forces_and_velocities.Update_Time_Varying_Material_Properties(time+dt);

        if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data_Global(solid_body_collection.deformable_body_collection.particles.X);
        solid_body_collection.Update_Position_Based_State(time+dt,(solids_parameters.allow_altitude_spring_change_between_updates?true:false));
        Make_Incompressible(dt,false); // make velocity divergence free
        if(articulated_rigid_body.Has_Actuators()) example_forces_and_velocities.Set_PD_Targets(dt,time);
        Backward_Euler_Step_Velocity(dt,time); // TODO: Tamar & Craig, do you need post stab?
        Diagnostics(dt,time,2,2,29,"backward Euler");
        solids_evolution_callbacks->Filter_Velocities(dt,time+dt,true);

        if(solids_parameters.set_velocity_from_positions_physbam){
            T postelasticanddamping_KE=0;
            for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
                postelasticanddamping_KE+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*
                    TV::Dot_Product(solid_body_collection.deformable_body_collection.particles.V(p),solid_body_collection.deformable_body_collection.particles.V(p));}
            // Compute damping force
            T postelasticanddampingnoprojection_KE=0;
            ARRAY<TV> damping_force(solid_body_collection.deformable_body_collection.particles.array_collection->Size());
            for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
                for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
                    DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
                    force.Get_Damping_Force(p,damping_force(p),dt,true);}}
            for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
                V_postelasticforces(p)+=damping_force(p);
                postelasticanddampingnoprojection_KE+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*
                    TV::Dot_Product(V_postelasticforces(p),V_postelasticforces(p));}
            // Given the total change in energy, subtract off that due to damping gives you the projection amount
            // Add this to the amount lost to the body (same with momentum)
            T energy_lost_to_damping=postelasticanddampingnoprojection_KE-postelasticforces_KE;
            T energy_lost_to_damping_and_projection=postelasticanddamping_KE-postelasticforces_KE;
            T energy_lost_to_projection=energy_lost_to_damping_and_projection-energy_lost_to_damping;
            energy_lost_to_bodies+=energy_lost_to_projection;
            if(solids_parameters.set_velocity_from_positions_conserve_exactly) energy_damped+=-energy_lost_to_projection;
            // Add the damping amount to energy damped to be added in next time around
            energy_damped+=-energy_lost_to_damping;}}

    if(solids_parameters.enforce_energy_conservation){
        solid_body_collection.Add_Energy_Correction_Force(V_save,rigid_velocity_save,solids_parameters.energy_correction_iterations,time+dt,dt);
        solid_body_collection.Compute_Energy_Error(V_save,rigid_velocity_save,time+dt,dt);}
        
    if(solids){
        PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("after removing joints.  after trapezoidal velocities dt=%f time=%f",dt,time),2,2);

        if(!solids_parameters.no_contact_friction){
            if(solids_parameters.use_post_cg_constraints) Apply_Constraints(dt,time);
            else if(repulsions){
                T prerepulsions_KE=0;
                if(solids_parameters.set_velocity_from_positions_physbam){
                    for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++)
                        prerepulsions_KE+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*
                            TV::Dot_Product(solid_body_collection.deformable_body_collection.particles.V(p),solid_body_collection.deformable_body_collection.particles.V(p));}
                repulsions->Adjust_Velocity_For_Self_Repulsion_Using_History(dt,false,true);
                if(solids_parameters.set_velocity_from_positions_physbam){
                    T postrepulsions_KE=0;
                    for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++)
                        postrepulsions_KE+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*
                            TV::Dot_Product(solid_body_collection.deformable_body_collection.particles.V(p),solid_body_collection.deformable_body_collection.particles.V(p));
                    energy_lost_to_bodies+=(postrepulsions_KE-prerepulsions_KE);
                    if(solids_parameters.set_velocity_from_positions_conserve_exactly) energy_damped+=-(postrepulsions_KE-prerepulsions_KE);}
                Diagnostics(dt,time,2,2,40,"self repulsions");}}
        solid_body_collection.deformable_body_collection.collisions.Reset_Object_Collisions();
        if(advance_rigid_bodies && solids_parameters.rigid_body_evolution_parameters.clamp_rigid_body_velocities){
            Clamp_Velocities(); // TODO: Examples should do this during the Advance_One_Time_Step_End_Callback example callback
            Diagnostics(dt,time,2,2,41,"clamp velocities");}
        if(solids_parameters.verbose) Print_Maximum_Velocities(time);}

    example_forces_and_velocities.Advance_One_Time_Step_End_Callback(dt,time);

    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Advance_One_Time_Step_Velocity End dt=%f time=%f",dt,time),2,2);
}
//#####################################################################
// Function Set_Velocity_From_Positions_Velocity_Update
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Set_Velocity_From_Positions_Velocity_Update(const T dt,const T time)
{
    solid_body_collection.deformable_body_collection.Setup_Set_Velocity_From_Positions(time,false,false);
    
    ARRAY<ARRAY<T> > convergence_norms(solid_body_collection.deformable_body_collection.deformables_forces.m);
    int iter;
    for(iter=1;iter<=solids_parameters.set_velocity_from_positions_iterations;iter++){
        T convergence_norm=0;
        for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++)
            convergence_norm=max(convergence_norm,ARRAYS_COMPUTATIONS::Maxabs(convergence_norms(f)));
        if(iter>1 && convergence_norm<solids_parameters.set_velocity_from_positions_tolerance) break;
        
        // Iterate over all forces
        ARRAY<ARRAY<VECTOR<T,3> > > force_quadratic_coefficients(solid_body_collection.deformable_body_collection.deformables_forces.m);
        ARRAY<ARRAY<ARRAY<T> > > v_n_hat(solid_body_collection.deformable_body_collection.deformables_forces.m);
        for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
            DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
            force_quadratic_coefficients(f).Resize(force.Get_Element_Count());
            v_n_hat(f).Resize(force.Get_Element_Count());
            FORCE_ELEMENTS& force_elements=*force.Get_Force_Elements();
            for(typename FORCE_ELEMENTS::ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){
                int s=iterator.Data();
                // Update the alpha for this force
                // Get a list of the nodes involved in this force
                T a=0;T b=0;T c=0;
                T A=0,B=0,C=0;
                ARRAY<int>& nodes=*force.Incident_Nodes(s);
                // Get combined mass
                T combined_one_over_mass=force.Get_Combined_One_Over_Mass(s);
                
                force.Compute_Quadratic_Contribution_For_Force(A,a,c,dt,s,combined_one_over_mass,true);
                // Iterate over the nodes
                v_n_hat(f)(s).Resize(nodes.m);
                for(int n=1;n<=nodes.m;n++){
                    int node=nodes(n);
                    TV direction=force.Get_Direction(s);
                    // Compute modified inital velocity
                    T v_n_correction=0;
                    for(int f_incident=1;f_incident<=solid_body_collection.deformable_body_collection.deformables_forces.m;f_incident++){
                        DEFORMABLES_FORCES<TV>& force_incident=*solid_body_collection.deformable_body_collection.deformables_forces(f_incident);
                        ARRAY<int>& incident_force_elements=*force_incident.Incident_Force_Elements(node);
                        for(int ife=1;ife<=incident_force_elements.m;ife++){
                            if(f_incident==f && incident_force_elements(ife)==s) continue;
                            v_n_correction+=dt*solid_body_collection.deformable_body_collection.particles.one_over_mass(node)*
                                TV::Dot_Product(force_incident.Get_Force(incident_force_elements(ife),node,false),direction);}}
                    force.Compute_Quadratic_Contribution_For_Node(B,C,b,dt,node,s,combined_one_over_mass,v_n_correction,true);
                    v_n_hat(f)(s)(n)=TV::Dot_Product(solid_body_collection.deformable_body_collection.particles.V(node),direction)+v_n_correction;}
                force.Compute_Quadratic_Contribution_For_Residual(B,C,a,b,c,dt,time,s,true);
                force_quadratic_coefficients(f)(s)=VECTOR<T,3>(A,B,C);}}
        
        // Solve for each force, and get an alpha
        for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
            DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
            FORCE_ELEMENTS& force_elements=*force.Get_Force_Elements();
            for(typename FORCE_ELEMENTS::ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){
                int s=iterator.Data();
                
                T A=force_quadratic_coefficients(f)(s)(1);
                T B=force_quadratic_coefficients(f)(s)(2);
                T C=force_quadratic_coefficients(f)(s)(3);
                
                if((A==(T)0) && (B==(T)0) && (C==(T)0)){
                    force.Set_Force(s,0);
                    continue;}
                
                T discriminant=B*B-4*A*C;
                if(discriminant<0){
//                    {std::stringstream ss;ss << "No solutions: negative discriminant position update" << std::endl;LOG::filecout(ss.str());}
//                    {std::stringstream ss;ss << "Negative discriminant in position update: " << discriminant << std::endl;LOG::filecout(ss.str());}
                    force.Set_Force(s,-B/(2*A));}
                else{
                    T alpha1=(-B+sqrt(discriminant))/(2*A);
                    T alpha2=(-B-sqrt(discriminant))/(2*A);
                    
                    force.Choose_Solution(solids_parameters.set_velocity_from_positions_use_orig_force,s,dt,alpha1,alpha2,v_n_hat(f)(s));}}}
        
        // Iterate over all forces
        for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
            DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
            FORCE_ELEMENTS& force_elements=*force.Get_Force_Elements();
            convergence_norms(f).Resize(force.Get_Element_Count());
            for(typename FORCE_ELEMENTS::ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){
                int s=iterator.Data();
                // Update the alpha for this force
                // Get a list of the nodes involved in this force
                T a=0;T b=0;T c=0;
                T A=0,B=0,C=0;
                ARRAY<int>& nodes=*force.Incident_Nodes(s);
                // Get combined mass
                T combined_one_over_mass=force.Get_Combined_One_Over_Mass(s);
                
                force.Compute_Quadratic_Contribution_For_Force(A,a,c,dt,s,combined_one_over_mass,true);
                // Iterate over the nodes
                for(int n=1;n<=nodes.m;n++){
                    int node=nodes(n);
                    TV direction=force.Get_Direction(s);
                    // Compute modified inital velocity
                    T v_n_correction=0;
                    for(int f_incident=1;f_incident<=solid_body_collection.deformable_body_collection.deformables_forces.m;f_incident++){
                        DEFORMABLES_FORCES<TV>& force_incident=*solid_body_collection.deformable_body_collection.deformables_forces(f_incident);
                        ARRAY<int>& incident_force_elements=*force_incident.Incident_Force_Elements(node);
                        for(int ife=1;ife<=incident_force_elements.m;ife++){
                            if(f_incident==f && incident_force_elements(ife)==s) continue;
                            v_n_correction+=dt*solid_body_collection.deformable_body_collection.particles.one_over_mass(node)*
                                TV::Dot_Product(force_incident.Get_Force(incident_force_elements(ife),node,false),direction);}}
                    force.Compute_Quadratic_Contribution_For_Node(B,C,b,dt,node,s,combined_one_over_mass,v_n_correction,true);}
                force.Compute_Quadratic_Contribution_For_Residual(B,C,a,b,c,dt,time,s,true);
                T new_alpha=force.Get_Force(s);
                convergence_norms(f)(s)=A*new_alpha*new_alpha+B*new_alpha+C;}}}
    // Update velocities
    for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
        TV vtemp=solid_body_collection.deformable_body_collection.particles.V(p);
        if(!solid_body_collection.deformable_body_collection.particles.one_over_mass(p)) continue;
        for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
            DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
            ARRAY<int>& incident_force_elements=*force.Incident_Force_Elements(p);
            for(int ife=1;ife<=incident_force_elements.m;ife++)
                vtemp+=dt*solid_body_collection.deformable_body_collection.particles.one_over_mass(p)*force.Get_Force(incident_force_elements(ife),p,false);}
        solid_body_collection.deformable_body_collection.particles.V(p)=vtemp;}
    
    // Update residuals
    for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
        DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
        FORCE_ELEMENTS& force_elements=*force.Get_Force_Elements();
        for(typename FORCE_ELEMENTS::ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){
            int s=iterator.Data();
            force.Update_Residual_Energy(s,convergence_norms(f)(s),time);}}

    if(solid_body_collection.print_energy){
        TV total_momentum=TV();
        for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++)
            total_momentum+=solid_body_collection.deformable_body_collection.particles.mass(p)*solid_body_collection.deformable_body_collection.particles.V(p);
        {std::stringstream ss;ss << "Total momentum = " << total_momentum << std::endl;LOG::filecout(ss.str());}}
    
    if(solids_parameters.set_velocity_from_positions_move_RE_to_KE_damping){
        T total_residual_energy=0;
        for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
            DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
            total_residual_energy+=force.Residual_Energy(time);}
        T A=0,B=0,C=0;
        ARRAY<TV> damping_forces(solid_body_collection.deformable_body_collection.particles.array_collection->Size());
        for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
            damping_forces(p)=TV();
            for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
                DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
                force.Get_Damping_Force(p,damping_forces(p),dt,false);}
            B+=solid_body_collection.deformable_body_collection.particles.mass(p)*TV::Dot_Product(solid_body_collection.deformable_body_collection.particles.V(p),damping_forces(p));
            A+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*TV::Dot_Product(damping_forces(p),damping_forces(p));}
        C=total_residual_energy-energy_damped;
        T epsilon;
        T discriminant=B*B-4*A*C;
        if(B>0) epsilon=(-B+sqrt(discriminant))/(2*A);
        else epsilon=(-B-sqrt(discriminant))/(2*A);
        
        for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++)
            solid_body_collection.deformable_body_collection.particles.V(p)+=epsilon*damping_forces(p);
        
        for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
            DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
            FORCE_ELEMENTS& force_elements=*force.Get_Force_Elements();
            for(typename FORCE_ELEMENTS::ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){
                int s=iterator.Data();
                force.Update_Residual_Energy(s,0,time);}}}
    else if(solids_parameters.set_velocity_from_positions_move_RE_to_KE_elastic){
        T total_residual_energy=0;
        for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
            DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
            total_residual_energy+=force.Residual_Energy(time);}
        
        T total_kinetic_energy=0;
        for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++)
            total_kinetic_energy+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*
                TV::Dot_Product(solid_body_collection.deformable_body_collection.particles.V(p),solid_body_collection.deformable_body_collection.particles.V(p));
        
        T A=0,B=0,C=0;
        ARRAY<TV> damping_forces(solid_body_collection.deformable_body_collection.particles.array_collection->Size());
        for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
            damping_forces(p)=TV();
            for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
                DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
                ARRAY<int>& incident_force_elements=*force.Incident_Force_Elements(p);
                for(int ife=1;ife<=incident_force_elements.m;ife++)
                    damping_forces(p)+=dt*solid_body_collection.deformable_body_collection.particles.one_over_mass(p)*force.Get_Force(incident_force_elements(ife),p,false);}
            B+=solid_body_collection.deformable_body_collection.particles.mass(p)*TV::Dot_Product(solid_body_collection.deformable_body_collection.particles.V(p),damping_forces(p));
            A+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*TV::Dot_Product(damping_forces(p),damping_forces(p));}
        C=total_residual_energy-energy_damped;
        T epsilon;
        T discriminant=B*B-4*A*C;
        if(discriminant<0) discriminant=0;
        if(B>0) epsilon=(-B+sqrt(discriminant))/(2*A);
        else epsilon=(-B-sqrt(discriminant))/(2*A);
        
        T new_total_kinetic_energy=0;
        for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
            solid_body_collection.deformable_body_collection.particles.V(p)+=epsilon*damping_forces(p);
            new_total_kinetic_energy+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*
                TV::Dot_Product(solid_body_collection.deformable_body_collection.particles.V(p),solid_body_collection.deformable_body_collection.particles.V(p));}
        T residual_to_save=(new_total_kinetic_energy-total_kinetic_energy+total_residual_energy-energy_damped)/(total_residual_energy);
        
        for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
            DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
            FORCE_ELEMENTS& force_elements=*force.Get_Force_Elements();
            for(typename FORCE_ELEMENTS::ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){
                int s=iterator.Data();
                force.Update_Residual_Energy(s,residual_to_save*force.Get_Residual_Energy(s),time);}}}

    if(solids_parameters.set_velocity_from_positions_damping){
        // PhysBAM damping
        T total_kinetic_energy_before_damping=0;
        TV total_momentum_before_damping=TV();
        for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
            total_kinetic_energy_before_damping+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*
                TV::Dot_Product(solid_body_collection.deformable_body_collection.particles.V(p),solid_body_collection.deformable_body_collection.particles.V(p));
            total_momentum_before_damping+=solid_body_collection.deformable_body_collection.particles.mass(p)*solid_body_collection.deformable_body_collection.particles.V(p);}

        ARRAY<TV> pre_damping_velocities(solid_body_collection.deformable_body_collection.particles.array_collection->Size());
        for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
            pre_damping_velocities(p)=solid_body_collection.deformable_body_collection.particles.V(p);}

        for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++)
            solid_body_collection.deformable_body_collection.deformables_forces(f)->Save_And_Reset_Elastic_Coefficient();
            
        // CG solve for damping
        Backward_Euler_Step_Velocity_Helper(dt,time,time,false);

        for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++)
            solid_body_collection.deformable_body_collection.deformables_forces(f)->Restore_Elastic_Coefficient();

        T total_kinetic_energy_after_damping=0;
        TV total_momentum_after_damping=TV();
        for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
            total_kinetic_energy_after_damping+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*
                TV::Dot_Product(solid_body_collection.deformable_body_collection.particles.V(p),solid_body_collection.deformable_body_collection.particles.V(p));
            total_momentum_after_damping+=solid_body_collection.deformable_body_collection.particles.mass(p)*solid_body_collection.deformable_body_collection.particles.V(p);}

        // Compute damping force
        T total_kinetic_energy_no_projection=0;
        ARRAY<TV> damping_force(solid_body_collection.deformable_body_collection.particles.array_collection->Size());
        for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
            for(int f=1;f<=solid_body_collection.deformable_body_collection.deformables_forces.m;f++){
                DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(f);
                force.Get_Damping_Force(p,damping_force(p),dt,true);}}
        for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++){
            pre_damping_velocities(p)+=damping_force(p);
            total_kinetic_energy_no_projection+=(T).5*solid_body_collection.deformable_body_collection.particles.mass(p)*
                TV::Dot_Product(pre_damping_velocities(p),pre_damping_velocities(p));}
        // Given the total change in energy, subtract off that due to damping gives you the projection amount
        // Add this to the amount lost to the body (same with momentum)
        T energy_lost_to_damping=total_kinetic_energy_no_projection-total_kinetic_energy_before_damping;
        T energy_lost_to_damping_and_projection=total_kinetic_energy_after_damping-total_kinetic_energy_before_damping;
        T energy_lost_to_projection=energy_lost_to_damping_and_projection-energy_lost_to_damping;
        energy_lost_to_bodies+=energy_lost_to_projection;
        momentum_lost_to_bodies+=total_momentum_before_damping-total_momentum_after_damping;
        // Add the damping amount to energy damped to be added in next time around
        energy_damped=-energy_lost_to_damping;}
    
    if(solid_body_collection.print_energy){
        TV total_momentum=TV();
        for(int p=1;p<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();p++)
            total_momentum+=solid_body_collection.deformable_body_collection.particles.mass(p)*solid_body_collection.deformable_body_collection.particles.V(p);
        {std::stringstream ss;ss << "Total momentum = " << total_momentum << std::endl;LOG::filecout(ss.str());}}
}
//#####################################################################
// Function Apply_Constraints
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Apply_Constraints(const T dt,const T time)
{
    ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    const bool advance_rigid_bodies=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m!=0;
    if(advance_rigid_bodies){
        articulated_rigid_body.Apply_Poststabilization(solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);
        Diagnostics(dt,time,0,0,345,"poststabilization");
        if(articulated_rigid_body.Has_Actuators() && !articulated_rigid_body.constrain_pd_directions){
            articulated_rigid_body.Compute_Position_Based_State(dt,time);
            articulated_rigid_body.Solve_Velocities_for_PD(time,dt,solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);}
        Diagnostics(dt,time,2,2,35,"solve velocities for pd");}
    Save_Position(X_save_for_constraints,rigid_X_save_for_constraints,rigid_rotation_save_for_constraints);
    if(solids_parameters.rigid_body_collision_parameters.use_persistant_contact){rigid_deformable_collisions->Process_Precomputed_Contact_With_Rigid_Bodies();use_existing_contact=true;}
    Diagnostics(dt,time,2,2,137,"consistent contact");
    Update_Positions_And_Apply_Contact_Forces(dt,time,true);use_existing_contact=false;
    Diagnostics(dt,time,4,2,37,"contact, prestabilization");
    if(rigid_body_collisions->prune_stacks_from_contact) rigid_body_collisions->Construct_Stacks();
    if(rigid_body_collisions->prune_contact_using_velocity) rigid_body_collisions->Compute_Contact_Frequency();
    Restore_Position(X_save_for_constraints,rigid_X_save_for_constraints,rigid_rotation_save_for_constraints);
    Diagnostics(dt,time,2,2,38,"restore position");
    if(advance_rigid_bodies){
        articulated_rigid_body.Apply_Poststabilization(solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);
        Diagnostics(dt,time,2,2,39,"poststabilization");}
    // modify velocity with inelastic and friction repulsions
    if(repulsions) repulsions->Adjust_Velocity_For_Self_Repulsion_Using_History(dt,false,true);
}
//#####################################################################
// Function Print_Maximum_Velocities
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Print_Maximum_Velocities(const T time) const
{
    {std::stringstream ss;ss<<"time = "<<time<<std::endl;LOG::filecout(ss.str());}
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particle;
    int max_index=0;T max_magnitude_squared=-FLT_MAX;const INDIRECT_ARRAY<ARRAY_VIEW<TV>,ARRAY<int>&>& V=particles.V.Subset(solid_body_collection.deformable_body_collection.dynamic_particles);
    for(int i=1;i<=V.Size();i++){T magnitude_squared=V(i).Magnitude_Squared();if(magnitude_squared>max_magnitude_squared){max_magnitude_squared=magnitude_squared;max_index=i;}}
    if(max_index){
        int p=solid_body_collection.deformable_body_collection.dynamic_particles(max_index);T max_magnitude=sqrt(max_magnitude_squared),max_magnitude_global=max_magnitude;
        if(solid_body_collection.deformable_body_collection.mpi_solids) max_magnitude_global=solid_body_collection.deformable_body_collection.mpi_solids->Reduce_Max_Global(max_magnitude_global);
        {std::stringstream ss;ss<<"maximum velocity = "<<max_magnitude_global;LOG::filecout(ss.str());}
        if(solid_body_collection.deformable_body_collection.mpi_solids) {std::stringstream ss;ss<<", local = "<<max_magnitude;LOG::filecout(ss.str());}
        {std::stringstream ss;ss<<" ("<<p<<")"<<std::endl;LOG::filecout(ss.str());}}
    int max_linear_index=0,max_angular_index=0;T max_linear_magnitude_squared=-FLT_MAX,max_angular_magnitude_squared=-FLT_MAX;
    for(int i=1;i<=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m;i++){const int p=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles(i);
        T linear_magnitude_squared=rigid_body_particles.V(p).Magnitude_Squared(),angular_magnitude_squared=rigid_body_particles.angular_velocity(p).Magnitude_Squared();
        if(linear_magnitude_squared>max_linear_magnitude_squared){max_linear_magnitude_squared=linear_magnitude_squared;max_linear_index=p;}
        if(angular_magnitude_squared>max_angular_magnitude_squared){max_angular_magnitude_squared=angular_magnitude_squared;max_angular_index=p;}}
    if(max_linear_index){
        T max_linear_magnitude=sqrt(max_linear_magnitude_squared),max_angular_magnitude=sqrt(max_angular_magnitude_squared);
        T max_linear_magnitude_global=max_linear_magnitude,max_angular_magnitude_global=max_angular_magnitude;
        if(solid_body_collection.deformable_body_collection.mpi_solids){
            max_linear_magnitude_global=solid_body_collection.deformable_body_collection.mpi_solids->Reduce_Max_Global(max_linear_magnitude_global);
            max_angular_magnitude_global=solid_body_collection.deformable_body_collection.mpi_solids->Reduce_Max_Global(max_angular_magnitude_global);}
        {std::stringstream ss;ss<<"maximum rigid linear velocity = "<<max_linear_magnitude_global;LOG::filecout(ss.str());}
        if(solid_body_collection.deformable_body_collection.mpi_solids) {std::stringstream ss;ss<<", local = "<<max_linear_magnitude;LOG::filecout(ss.str());}
            {std::stringstream ss;ss<<" ("<<max_linear_index<<")\n";LOG::filecout(ss.str());}
            {std::stringstream ss;ss<<"maximum rigid angular velocity = "<<max_angular_magnitude_global;LOG::filecout(ss.str());}
        if(solid_body_collection.deformable_body_collection.mpi_solids) {std::stringstream ss;ss<<", local = "<<max_angular_magnitude;LOG::filecout(ss.str());}
        {std::stringstream ss;ss<<" ("<<max_angular_index<<")"<<std::endl;LOG::filecout(ss.str());}}
}
//#####################################################################
// Function Diagnostics
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Diagnostics(const T dt,const T time,const int velocity_time,const int position_time,int step,const char* description)
{
    static const char* time_table[]={"n","(n+1/2)","(n+1)","(n+3/2)","(n+2)"};
    solid_body_collection.Print_Energy(time+position_time*(T).5*time,step);
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Finished step %i (%s).  State: x^%s  v^%s.   dt=%f time=%f",step,description,
        time_table[position_time],time_table[velocity_time],dt,time),2,3);
}
//#####################################################################
// Function Update_Positions_And_Apply_Contact_Forces
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Update_Positions_And_Apply_Contact_Forces(const T dt,const T time,const bool use_saved_pairs)
{
    ARTICULATED_RIGID_BODY<TV>* articulated_rigid_body=&solid_body_collection.rigid_body_collection.articulated_rigid_body;

    Euler_Step_Position(dt,time);

    for(int i=1;i<=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m;i++){
        RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(solid_body_collection.rigid_body_collection.simulated_rigid_body_particles(i));
        if(!rigid_body.is_static) rigid_body.Update_Bounding_Box();}
    for(int i=1;i<=solid_body_collection.rigid_body_collection.kinematic_rigid_bodies.m;i++){
        RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(solid_body_collection.rigid_body_collection.kinematic_rigid_bodies(i));
        rigid_body.Update_Bounding_Box();}

    if(solids_parameters.use_rigid_deformable_contact){
        rigid_deformable_collisions->Process_Contact(dt,time,articulated_rigid_body,use_saved_pairs,use_existing_contact,rigid_X_save,
            rigid_rotation_save,rigid_velocity_difference,rigid_angular_momentum_difference,X_save,solids_parameters.rigid_body_collision_parameters.collision_body_thickness);
        // TODO: rigid/deformable shock propagation step
        if(solids_parameters.rigid_body_collision_parameters.perform_contact && solids_parameters.rigid_body_collision_parameters.use_shock_propagation)
            rigid_body_collisions->Shock_Propagation_Using_Graph(dt,time,articulated_rigid_body,use_saved_pairs);}
    else{
        // TODO: move into rigid/deformable collisions and get interleaved with rigid/deformable contact
        if(solids_parameters.rigid_body_collision_parameters.perform_contact)
            rigid_body_collisions->Process_Contact_Using_Graph(dt,time,articulated_rigid_body,solids_parameters.rigid_body_evolution_parameters.correct_contact_energy,use_saved_pairs);

        // rigid/rigid shock propagation
        if(solids_parameters.rigid_body_collision_parameters.perform_contact && solids_parameters.rigid_body_collision_parameters.use_shock_propagation)
            rigid_body_collisions->Shock_Propagation_Using_Graph(dt,time,articulated_rigid_body,use_saved_pairs);

        if(solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions){
            int interactions=0;
            if(use_existing_contact)
                interactions+=solid_body_collection.deformable_body_collection.collisions.Adjust_Existing_Nodes_For_Collision_Body_Collisions(solid_body_collection.deformable_body_collection.binding_list,solid_body_collection.deformable_body_collection.soft_bindings,X_save,dt,0);
            else 
                interactions+=solid_body_collection.deformable_body_collection.collisions.Adjust_Nodes_For_Collision_Body_Collisions(solid_body_collection.deformable_body_collection.binding_list,solid_body_collection.deformable_body_collection.soft_bindings,X_save,dt,0);
            if(interactions) LOG::Stat("collision body collisions",interactions);}}

}
//#####################################################################
// Function Update_Velocity_Using_Stored_Differences
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Update_Velocity_Using_Stored_Differences(const T dt,const T time,const int p)
{
    RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(p);
    if(rigid_body.is_static) return;
    if(solid_body_collection.rigid_body_collection.rigid_body_particle.kinematic(p)) kinematic_evolution.Set_External_Velocities(rigid_body.V(),rigid_body.Angular_Velocity(),time+dt,p);
    rigid_body.V()+=rigid_velocity_difference(p);
    rigid_body.Angular_Momentum()+=rigid_angular_momentum_difference(p);
    rigid_body.Update_Angular_Velocity();
}
//#####################################################################
// Function Update_Velocity_Using_Stored_Differences
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Update_Velocity_Using_Stored_Differences(const T dt,const T time)
{
    for(int i=1;i<=solid_body_collection.deformable_body_collection.simulated_particles.m;i++){int p=solid_body_collection.deformable_body_collection.simulated_particles(i);
        solid_body_collection.deformable_body_collection.particles.V(p)+=V_difference(p);}
    for(int i=1;i<=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles(i);
        Update_Velocity_Using_Stored_Differences(dt,time,p);}
    for(int i=1;i<=solid_body_collection.rigid_body_collection.kinematic_rigid_bodies.m;i++){int p=solid_body_collection.rigid_body_collection.kinematic_rigid_bodies(i);
        kinematic_evolution.Set_External_Velocities(solid_body_collection.rigid_body_collection.rigid_body_particle.V(p),solid_body_collection.rigid_body_collection.rigid_body_particle.angular_velocity(p),time+dt,p);}
}
//#####################################################################
// Function Compute_Momentum_Differences
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Compute_Momentum_Differences()
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    rigid_velocity_difference.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size());rigid_angular_momentum_difference.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size());
    V_difference.Resize(particles.array_collection->Size());
    for(int i=1;i<=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles(i);
        RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(p);
        rigid_velocity_difference(p)=rigid_body.V()-rigid_velocity_save(p).linear;
        rigid_angular_momentum_difference(p)=rigid_body.Angular_Momentum()-rigid_angular_momentum_save(p);}
    for(int i=1;i<=solid_body_collection.deformable_body_collection.simulated_particles.m;i++){
        int p=solid_body_collection.deformable_body_collection.simulated_particles(i);V_difference(p)=particles.V(p)-V_save(p);}
}
//#####################################################################
// Function Save_Velocity
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Save_Velocity()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    V_save.Resize(particles.array_collection->Size(),false,false);
    rigid_velocity_save.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
    rigid_angular_momentum_save.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
    V_save.Subset(solid_body_collection.deformable_body_collection.simulated_particles)=particles.V.Subset(solid_body_collection.deformable_body_collection.simulated_particles);
    for(int i=1;i<=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles(i);
        rigid_velocity_save(p).linear=rigid_body_collection.rigid_body_particle.V(p);rigid_velocity_save(p).angular=rigid_body_collection.rigid_body_particle.angular_velocity(p);rigid_angular_momentum_save(p)=rigid_body_collection.rigid_body_particle.angular_momentum(p);}
    for(int i=1;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++) if(rigid_body_collection.Is_Active(i)){RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(i);
        if(!body.Is_Simulated()){rigid_velocity_save(i).linear=rigid_body_collection.rigid_body_particle.V(i);rigid_velocity_save(i).angular=rigid_body_collection.rigid_body_particle.angular_velocity(i);rigid_angular_momentum_save(i)=rigid_body_collection.rigid_body_particle.angular_momentum(i);}}
}
//#####################################################################
// Function Restore_Velocity
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Restore_Velocity() const
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    particles.V.Subset(solid_body_collection.deformable_body_collection.simulated_particles)=V_save.Subset(solid_body_collection.deformable_body_collection.simulated_particles);
    for(int i=1;i<=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles(i);
        rigid_body_collection.rigid_body_particle.V(p)=rigid_velocity_save(p).linear;
        rigid_body_collection.rigid_body_particle.angular_momentum(p)=rigid_angular_momentum_save(p);rigid_body_collection.Rigid_Body(p).Update_Angular_Velocity();}
    for(int i=1;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++) if(rigid_body_collection.Is_Active(i)){RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(i);
        if(!body.Is_Simulated()){
            rigid_body_collection.rigid_body_particle.V(i)=rigid_velocity_save(i).linear;rigid_body_collection.rigid_body_particle.angular_momentum(i)=rigid_angular_momentum_save(i);body.Update_Angular_Velocity();}}
}
//#####################################################################
// Function Exchange_Velocity
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Exchange_Velocity()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particle;
    V_save.Resize(particles.array_collection->Size(),false,false);
    rigid_velocity_save.Resize(rigid_body_particles.array_collection->Size(),false,false);
    rigid_angular_momentum_save.Resize(rigid_body_particles.array_collection->Size(),false,false);
    for(int i=1;i<=solid_body_collection.deformable_body_collection.simulated_particles.m;i++){int p=solid_body_collection.deformable_body_collection.simulated_particles(i);
        exchange(V_save(p),deformable_body_collection.particles.V(p));}
    for(int i=1;i<=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles(i);
        exchange(rigid_velocity_save(p).linear,rigid_body_particles.V(p));exchange(rigid_velocity_save(p).angular,rigid_body_particles.angular_velocity(p));
        exchange(rigid_angular_momentum_save(p),rigid_body_particles.angular_momentum(p));}
    for(int i=1;i<=rigid_body_particles.array_collection->Size();i++) if(solid_body_collection.rigid_body_collection.Is_Active(i)){RIGID_BODY<TV>& body=solid_body_collection.rigid_body_collection.Rigid_Body(i);
        if(!body.Is_Simulated()){
            exchange(rigid_velocity_save(i).linear,rigid_body_particles.V(i));exchange(rigid_velocity_save(i).angular,rigid_body_particles.angular_velocity(i));
            exchange(rigid_angular_momentum_save(i),rigid_body_particles.angular_momentum(i));}}
}
//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Initialize_Rigid_Bodies(const T frame_rate, const bool restart)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    // initialize kinematic object positions and velocities
    if(!restart){
        kinematic_evolution.Get_Current_Kinematic_Keyframes(1/frame_rate,time);
        kinematic_evolution.Set_External_Positions(rigid_body_collection.rigid_body_particle.X,rigid_body_collection.rigid_body_particle.rotation,time);
        kinematic_evolution.Set_External_Velocities(rigid_body_collection.rigid_body_particle.V,rigid_body_collection.rigid_body_particle.angular_velocity,time,time);
        rigid_body_collection.Update_Angular_Momentum();
        for(int i(1);i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++) if(rigid_body_collection.Is_Active(i)){rigid_body_collection.rigid_body_particle.rotation(i).Normalize();}}

    RIGID_BODY_COLLISIONS<TV>::Adjust_Bounding_Boxes(rigid_body_collection);
    // rigid body collisions
    rigid_body_collisions=new RIGID_BODY_COLLISIONS<TV>(rigid_body_collection,solids_parameters.rigid_body_collision_parameters,rigids_evolution_callbacks,
        *solid_body_collection.example_forces_and_velocities);
    rigid_body_collisions->spatial_partition->Compute_Voxel_Size(solids_parameters.rigid_body_collision_parameters.rigid_collisions_spatial_partition_voxel_size_heuristic,
        solids_parameters.rigid_body_collision_parameters.rigid_collisions_spatial_partition_number_of_cells,solids_parameters.rigid_body_collision_parameters.rigid_collisions_spatial_partition_voxel_size_scale_factor);
    rigid_body_collisions->verbose=solids_parameters.verbose;

    // partitions and hierarchies
    if(solids_parameters.rigid_body_collision_parameters.rigid_collisions_use_particle_partition){
        rigid_body_collisions->intersections.Use_Particle_Partition(true,solids_parameters.rigid_body_collision_parameters.rigid_collisions_particle_partition_size*VECTOR<int,d>::All_Ones_Vector());
        rigid_body_collisions->intersections.Use_Particle_Partition_Center_Phi_Test(solids_parameters.rigid_body_collision_parameters.rigid_collisions_use_particle_partition_center_phi_test);}
    if(solids_parameters.rigid_body_collision_parameters.rigid_collisions_use_triangle_hierarchy){
        rigid_body_collisions->intersections.Use_Triangle_Hierarchy();
        if(solids_parameters.rigid_body_collision_parameters.rigid_collisions_use_triangle_hierarchy_center_phi_test) rigid_body_collisions->intersections.Use_Triangle_Hierarchy_Center_Phi_Test();
        if(solids_parameters.rigid_body_collision_parameters.rigid_collisions_use_edge_intersection) rigid_body_collisions->intersections.Use_Edge_Intersection();}

    // rigid deformable collisions
    rigid_deformable_collisions=new RIGID_DEFORMABLE_COLLISIONS<TV>(solid_body_collection,*rigid_body_collisions,solids_parameters);

    // dynamics
    solids_parameters.rigid_body_evolution_parameters.rigid_cfl=solids_parameters.cfl;
    solids_parameters.rigid_body_evolution_parameters.rigid_minimum_dt=solids_parameters.rigid_body_evolution_parameters.minimum_rigid_body_time_step_fraction/frame_rate;
    solids_parameters.rigid_body_evolution_parameters.rigid_maximum_dt=solids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction/frame_rate;
    solids_evolution_callbacks->Update_Solids_Parameters(time);
    rigid_body_collisions->Initialize_Data_Structures(); // Must be done before we check interpenetration statistics

    // check for bad initial data (need to give it a chance to set up the collision manager first)
    if(solids_parameters.rigid_body_collision_parameters.rigid_collisions_print_interpenetration_statistics) rigid_body_collisions->Print_Interpenetration_Statistics();
    else if(!rigid_body_collisions->Check_For_Any_Interpenetration()) {std::stringstream ss;ss<<"No initial interpenetration"<<std::endl;LOG::filecout(ss.str());}
}
template<class TV> bool NEWMARK_EVOLUTION<TV>::
Use_CFL() const
{
    return true;
}
//#####################################################################
template class NEWMARK_EVOLUTION<VECTOR<float,1> >;
template class NEWMARK_EVOLUTION<VECTOR<float,2> >;
template class NEWMARK_EVOLUTION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class NEWMARK_EVOLUTION<VECTOR<double,1> >;
template class NEWMARK_EVOLUTION<VECTOR<double,2> >;
template class NEWMARK_EVOLUTION<VECTOR<double,3> >;
#endif

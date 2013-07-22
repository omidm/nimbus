//#####################################################################
// Copyright 2006-2009, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_EVOLUTION
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Krylov_Solvers/LANCZOS_ITERATION.h>
#include <PhysBAM_Tools/Krylov_Solvers/SYMMQMR.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGIDS_ONLY_COLLISION_CALLBACKS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigids_Evolution/RIGIDS_BACKWARD_EULER_SYSTEM.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigids_Evolution/RIGIDS_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigids_Evolution/RIGIDS_EVOLUTION_CALLBACKS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigids_Evolution/RIGIDS_VELOCITY.h>
#include <stdexcept>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGIDS_EVOLUTION<TV>::
RIGIDS_EVOLUTION(RIGIDS_PARAMETERS<TV>& rigids_parameters_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
    :time(0),rigids_collision_callbacks(new RIGIDS_ONLY_COLLISION_CALLBACKS<TV>(*this)),rigid_body_collection(rigid_body_collection_input),
    kinematic_evolution(rigid_body_collection,rigids_parameters_input.rigid_body_evolution_parameters.rigid_geometry_evolution_parameters.use_kinematic_keyframes),
    use_existing_contact(false),rigids_parameters(rigids_parameters_input),mpi_rigids(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGIDS_EVOLUTION<TV>::
~RIGIDS_EVOLUTION()
{
    delete rigid_body_collisions;
    delete rigids_collision_callbacks;
}
//#####################################################################
// Function Prepare_Backward_Euler_System
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Prepare_Backward_Euler_System(RIGIDS_BACKWARD_EULER_SYSTEM<TV>& system,const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update)
{
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=rigid_body_collection.rigid_body_particle;
    ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body=rigid_body_collection.articulated_rigid_body; // Needn't be a pointer
    rigid_body_collection.Update_Angular_Velocity(); // make sure omega = I^{-1} L

    rigid_F_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
    rigid_R_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
    rigid_S_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
    rigid_B_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
    if(rigids_parameters.implicit_solve_parameters.evolution_solver_type!=krylov_solver_cg){
        rigid_AR_full.Resize(rigid_body_particles.array_collection->Size(),false,false);}
    RIGIDS_VELOCITY<TV> B_all(rigid_B_full,rigid_body_collection);
    RIGIDS_VELOCITY<TV> F_all(rigid_F_full,rigid_body_collection);
    ARRAY<TWIST<TV> > twist; twist.Resize(rigid_body_particles.V.Size());
    for(int i=1;i<=twist.Size();i++) twist(i)=TWIST<TV>(rigid_body_particles.V(i),rigid_body_particles.angular_velocity(i));
    RIGIDS_VELOCITY<TV> V_all(twist,rigid_body_collection);

    ARRAYS_COMPUTATIONS::Fill(rigid_B_full,TWIST<TV>());
    rigid_body_collection.rigids_example_forces_and_velocities->Add_External_Forces(rigid_B_full,current_velocity_time+dt);
    rigid_body_collection.Add_Velocity_Independent_Forces(rigid_B_full,current_velocity_time+dt); // this is a nop for binding forces
    rigid_body_collection.rigid_body_cluster_bindings.Distribute_Force_To_Parents(rigid_B_full);

    Initialize_World_Space_Masses();
    if(articulated_rigid_body.constrain_pd_directions){
        for(int i=1;i<=rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=rigid_body_collection.dynamic_rigid_body_particles(i);rigid_B_full(p)=world_space_rigid_mass_inverse(p)*rigid_B_full(p);}
        articulated_rigid_body.Poststabilization_Projection(rigid_B_full,true);
        for(int i=1;i<=rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=rigid_body_collection.dynamic_rigid_body_particles(i);rigid_B_full(p)=world_space_rigid_mass(p)*rigid_B_full(p);}}

    for(int i=1;i<=rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=rigid_body_collection.dynamic_rigid_body_particles(i);
        rigid_B_full(p)=twist(p)+world_space_rigid_mass_inverse(p)*rigid_B_full(p)*dt;}

    V_all=B_all;
    rigid_body_collection.Update_Angular_Momentum();
    Diagnostics(dt,current_position_time,0,0,604,"Before boundary conditions");
    system.Set_Global_Boundary_Conditions(V_all,rigid_X_save,rigid_rotation_save,rigid_velocity_save,rigid_angular_momentum_save,rigids_parameters.implicit_solve_parameters.test_system,
        rigids_parameters.implicit_solve_parameters.print_matrix);
    if(velocity_update && !rigids_parameters.use_post_cg_constraints) Apply_Constraints(dt,current_velocity_time);

    // TODO: This is completely broken for trapezoid and likely for BE as well.
    if(rigid_body_collection.articulated_rigid_body.Has_Actuators() && rigid_body_collection.articulated_rigid_body.constrain_pd_directions){
        for(int i=1;i<=rigid_body_particles.rigid_geometry.Size();i++) saved_pd(i)=rigid_body_particles.rigid_geometry(i)->Twist();
        rigid_body_collection.articulated_rigid_body.Poststabilization_Projection(saved_pd,true);
        for(int i=1;i<=rigid_body_particles.rigid_geometry.Size();i++) saved_pd(i)=rigid_body_particles.rigid_geometry(i)->Twist()-saved_pd(i);}

    for(int i=1;i<=twist.Size();i++){rigid_body_particles.V(i)=twist(i).linear;rigid_body_particles.angular_velocity(i)=twist(i).angular;}
    
    //TODO: add callback
    //if(!velocity_update) example_forces_and_velocities->Add_External_Impulses_Before(B_full,current_position_time,(T)2*dt); // 2*dt is position dt TODO: what time?
}
//#####################################################################
// Function Prepare_Backward_Euler_Step
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Finish_Backward_Euler_Step(const T dt,const T current_position_time,const bool velocity_update)
{
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=rigid_body_collection.rigid_body_particle;
    RIGIDS_VELOCITY<TV> F_all(rigid_F_full,rigid_body_collection);

    rigid_body_collisions->Remove_Contact_Joints();

    if(velocity_update && rigids_parameters.use_post_cg_constraints){ // return rhs + dt Fd V^n+1 for friction processing
        for(int i=1;i<=rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=rigid_body_collection.dynamic_rigid_body_particles(i);
            TWIST<TV> twist=rigid_B_full(p)+world_space_rigid_mass_inverse(p)*rigid_F_full(p)*dt;
            rigid_body_particles.V(p)=twist.linear;rigid_body_particles.angular_velocity(p)=twist.angular;}

        if(rigid_body_collection.articulated_rigid_body.Has_Actuators() && rigid_body_collection.articulated_rigid_body.constrain_pd_directions){
            ARRAY<TWIST<TV> > twist;twist.Resize(rigid_body_particles.V.Size());
            for(int i=1;i<=twist.Size();i++) twist(i)=TWIST<TV>(rigid_body_particles.V(i),rigid_body_particles.angular_velocity(i));
            rigid_body_collection.articulated_rigid_body.Poststabilization_Projection(twist,true);
            for(int i=1;i<=rigid_body_particles.V.Size();i++){rigid_body_particles.V(i)=twist(i).linear+saved_pd(i).linear;rigid_body_particles.angular_velocity(i)=twist(i).angular+saved_pd(i).angular;}}}
    rigid_body_collection.Update_Angular_Momentum();

    Diagnostics(dt,current_position_time,0,0,610,"After undo projections");
    
    rigid_body_collection.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Velocities();
    
    rigid_body_collection.Update_Angular_Momentum(); // make sure L = I omega
}
//#####################################################################
// Function Backward_Euler_Step_Velocity_Helper
//#####################################################################
// assumes all solids_forces are linear in velocity, with a symmetric positive definite Jacobian.
template<class TV> void RIGIDS_EVOLUTION<TV>::
Backward_Euler_Step_Velocity_Helper(const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update)
{
    LOG::SCOPE scope("backward euler step velocity helper",rigids_parameters.threadid);
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=rigid_body_collection.rigid_body_particle;
    RIGIDS_BACKWARD_EULER_SYSTEM<TV> system(*this,rigid_body_collection,dt,current_velocity_time,current_position_time,&rigid_body_collection.articulated_rigid_body,velocity_update,
        rigids_parameters.enforce_poststabilization_in_cg);

    Prepare_Backward_Euler_System(system,dt,current_velocity_time,current_position_time,velocity_update);

    static CONJUGATE_GRADIENT<T> cg;
    static CONJUGATE_RESIDUAL<T> cr;
    static SYMMQMR<T> symmqmr;
    KRYLOV_SOLVER<T>* solver=0;
    const char* solver_name=0;
    if(rigids_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_cg){solver=&cg;solver_name="CG";}
    else if(rigids_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_cr){solver=&cr;solver_name="CONJUGATE_RESIDUAL";}
    else if(rigids_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_symmqmr){solver=&symmqmr;solver_name="SYMMQMR";}
    solver->print_diagnostics=rigid_body_collection.print_diagnostics;
    solver->print_residuals=rigid_body_collection.print_residuals;
    solver->iterations_used=&rigid_body_collection.iterations_used_diagnostic;
    solver->restart_iterations=rigids_parameters.implicit_solve_parameters.cg_restart_iterations;
    system.project_nullspace_frequency=rigids_parameters.implicit_solve_parameters.project_nullspace_frequency;
    ARRAY<TWIST<TV> > twist; twist.Resize(rigid_body_particles.V.Size());
    for(int i=1;i<=twist.Size();i++) twist(i)=TWIST<TV>(rigid_body_particles.V(i),rigid_body_particles.angular_velocity(i));
    RIGIDS_VELOCITY<TV> V(twist,rigid_body_collection),F(rigid_F_full,rigid_body_collection),
        R(rigid_R_full,rigid_body_collection),S(rigid_S_full,rigid_body_collection),B(rigid_B_full,rigid_body_collection),
        AR(rigid_AR_full,rigid_body_collection);
    RIGIDS_MASS<TV> mass(rigid_body_collection); // TODO: Doing duplicate computation of mass.

    if(rigids_parameters.implicit_solve_parameters.spectral_analysis){V=B;LANCZOS_ITERATION<T>::Print_Spectral_Information(system,V,F,S,(T)1e-2*rigids_parameters.implicit_solve_parameters.cg_tolerance,rigids_parameters.implicit_solve_parameters.lanczos_iterations);}

    LOG::Time(solver_name,rigids_parameters.threadid);
    Diagnostics(dt,current_position_time,0,0,606,"Before solve");
    if(!solver->Solve(system,V,B,F,S,R,AR,rigids_parameters.implicit_solve_parameters.cg_tolerance,1,rigids_parameters.implicit_solve_parameters.cg_iterations) && rigids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure)
        throw std::runtime_error("Backward Euler Failed");
    Diagnostics(dt,current_position_time,0,0,607,"After solve");
    LOG::Stop_Time(rigids_parameters.threadid);

    if(velocity_update && rigids_parameters.use_post_cg_constraints)
        system.Force(V,F);
    Finish_Backward_Euler_Step(dt,current_position_time,velocity_update);

    for(int i=1;i<=twist.Size();i++){rigid_body_particles.V(i)=twist(i).linear;rigid_body_particles.angular_velocity(i)=twist(i).angular;}
}
//#####################################################################
// Function Average_And_Exchange_Position
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Average_And_Exchange_Position()
{
    assert(rigid_X_save.m==rigid_body_collection.rigid_body_particle.array_collection->Size() && rigid_rotation_save.m==rigid_body_collection.rigid_body_particle.array_collection->Size());
    const ARRAY<int>& simulated_rigid_body_particles=rigid_body_collection.simulated_rigid_body_particles;
    ARRAY<int> rigid_body_indices(simulated_rigid_body_particles);rigid_body_indices.Append_Elements(rigid_body_collection.kinematic_rigid_bodies);
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
template<class TV> void RIGIDS_EVOLUTION<TV>::
Trapezoidal_Step_Velocity(const T dt,const T time)
{
    LOG::SCOPE scope("trapezoidal step velocity",rigids_parameters.threadid);
    // save V at time
    Save_Velocity();
    // update V implicitly to time+dt/2
    Backward_Euler_Step_Velocity_Helper(dt/2,time,time+dt/2,true);
    // set up rigid_V_save for extrapolation step
    rigid_V_save.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
    for(int i=1;i<=rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=rigid_body_collection.dynamic_rigid_body_particles(i);
        rigid_V_save(p).linear=rigid_velocity_save(p).linear;
        rigid_V_save(p).angular=rigid_body_collection.Rigid_Body(p).World_Space_Inertia_Tensor_Inverse_Times(rigid_angular_momentum_save(p));}
    // Use V_n instead of V_save rather than copying it around another time.  Also simplifies state dependencies.
    // extrapolate V to time+dt based on V at time and time+dt/2
    ARRAY<TWIST<TV> > twist; twist.Resize(rigid_body_collection.rigid_body_particle.V.Size());
    for(int i=1;i<=twist.Size();i++) twist(i)=TWIST<TV>(rigid_body_collection.rigid_body_particle.V(i),rigid_body_collection.rigid_body_particle.angular_velocity(i));
    RIGIDS_VELOCITY<TV> V_n(rigid_V_save,rigid_body_collection),V(twist,rigid_body_collection);
    V*=(T)2;V-=V_n;

    // enforce boundary conditions again
    kinematic_evolution.Set_External_Velocities(rigid_body_collection.rigid_body_particle.V,rigid_body_collection.rigid_body_particle.angular_velocity,time+dt,time+dt/2);
    rigid_body_collection.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Velocities();
    rigid_body_collection.Update_Angular_Momentum(rigid_body_collection.simulated_rigid_body_particles);
    rigid_body_collection.Update_Angular_Momentum(rigid_body_collection.kinematic_rigid_bodies);
    
    for(int i=1;i<=twist.Size();i++){rigid_body_collection.rigid_body_particle.V(i)=twist(i).linear;rigid_body_collection.rigid_body_particle.angular_velocity(i)=twist(i).angular;}
}
//#####################################################################
// Function Backward_Euler_Step_Velocity
//#####################################################################
// assumes X is fixed at time+dt
template<class TV> void RIGIDS_EVOLUTION<TV>::
Backward_Euler_Step_Velocity(const T dt,const T time)
{
    Backward_Euler_Step_Velocity_Helper(dt,time,time+dt,true);
    // enforce boundary conditions again
    rigid_body_collection.Update_Angular_Velocity(rigid_body_collection.simulated_rigid_body_particles);
    kinematic_evolution.Set_External_Velocities(rigid_body_collection.rigid_body_particle.V,rigid_body_collection.rigid_body_particle.angular_velocity,time+dt,time+dt);
    
    // enforce boundary conditions again
    rigid_body_collection.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Velocities();
    rigid_body_collection.Update_Angular_Momentum(rigid_body_collection.simulated_rigid_body_particles);
}
//#####################################################################
// Function Advance_One_Time_Step_Position
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Advance_One_Time_Step_Position(const T dt,const T time,const bool solids)
{
    LOG::SCOPE scope("solids position update",rigids_parameters.threadid);
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Advance_One_Time_Step_Position Start dt=%f time=%f",dt,time),2,2);
    ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body=rigid_body_collection.articulated_rigid_body; // Needn't be a pointer
    const bool advance_rigid_bodies=true; //solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m!=0;  TODO: Fix this.

    Diagnostics(dt,time,0,0,1,"begin integration");
    //TODO: add callback
    //example_forces_and_velocities.Advance_One_Time_Step_Begin_Callback(dt,time);

    rigids_evolution_callbacks->Update_Rigids_Parameters(time);
    if(solids){
        rigid_body_collisions->Initialize_Data_Structures();

        // save position and velocity for later trapezoidal rule
        Save_Velocity();
        Save_Position(rigid_X_save,rigid_rotation_save);}

    rigid_body_collection.Update_Position_Based_State(time+dt);

    // get momentum difference for v^n -> v^{n+1/2} udpate

    //TODO: add callback
    //if(articulated_rigid_body.Has_Actuators()) example_forces_and_velocities.Set_PD_Targets(dt,time);

    Backward_Euler_Step_Velocity_Helper(dt/2,time,time,false); // update V implicitly to time+dt/2

    if(rigids_parameters.verbose) Print_Maximum_Velocities(time);
    Diagnostics(dt,time,1,0,6,"backward Euler");
    if(!solids) return; // early exit for fluids only in parallel

    Compute_Momentum_Differences();

    // add collision impulses to time n velocities and save
    Euler_Step_Position(dt,time);
    Diagnostics(dt,time,1,2,10,"Euler step position");
    Exchange_Velocity();
    Diagnostics(dt,time,0,2,11,"restore velocity");
 
    if(rigids_parameters.rigid_body_collision_parameters.perform_collisions && advance_rigid_bodies) rigid_body_collisions->Add_Elastic_Collisions(dt,time);
 
    Diagnostics(dt,time,0,2,12,"add elastic collisions");

    Restore_Position(rigid_X_save,rigid_rotation_save);
    Diagnostics(dt,time,0,0,13,"restore position");
    if(advance_rigid_bodies){
        articulated_rigid_body.Apply_Poststabilization(rigids_parameters.implicit_solve_parameters.test_system,rigids_parameters.implicit_solve_parameters.print_matrix);
        Diagnostics(dt,time,0,0,14,"poststabilization");}
    Save_Velocity();

    // update positions, apply contact, arb prestabilization and push out
    if(rigids_parameters.rigid_body_collision_parameters.perform_contact && advance_rigid_bodies) rigid_body_collisions->Compute_Contact_Graph(dt,time,&articulated_rigid_body);
 
    Update_Velocity_Using_Stored_Differences(dt/2,time);
    Diagnostics(dt,time,1,0,18,"update velocity using stored differences");
    if(advance_rigid_bodies){
        articulated_rigid_body.Apply_Poststabilization(rigids_parameters.implicit_solve_parameters.test_system,rigids_parameters.implicit_solve_parameters.print_matrix);
        Diagnostics(dt,time,0,0,185,"poststabilization");
        if(articulated_rigid_body.Has_Actuators() && !articulated_rigid_body.constrain_pd_directions){
            articulated_rigid_body.Compute_Position_Based_State(dt,time);
            articulated_rigid_body.Solve_Velocities_for_PD(time,dt/2,rigids_parameters.implicit_solve_parameters.test_system,rigids_parameters.implicit_solve_parameters.print_matrix);}
        Diagnostics(dt,time,1,0,19,"solve velocities for pd");}

    // update positions, apply contact, arb prestabilization and push out
    Update_Positions_And_Apply_Contact_Forces(dt,time,false);
    Diagnostics(dt,time,1,2,20,"contact, prestabilization");
    if(rigids_parameters.rigid_body_collision_parameters.use_push_out){
        if(mpi_rigids || rigids_parameters.rigid_body_collision_parameters.use_legacy_push_out) rigid_body_collisions->Process_Push_Out_Legacy(); 
        else if(rigids_parameters.rigid_body_collision_parameters.use_projected_gauss_seidel_push_out) rigid_body_collisions->Process_Push_Out_Projected_Gauss_Seidel();
        else rigid_body_collisions->Process_Push_Out(true,rigids_parameters.rigid_body_evolution_parameters.residual_push_out_depth);
        Diagnostics(dt,time,1,2,21,"push out");}
 
    Restore_Velocity();
    Diagnostics(dt,time,0,2,22,"restore velocity");

    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Advance_One_Time_Step_Position End dt=%f time=%f",dt,time),2,2);
}
//#####################################################################
// Function Advance_One_Time_Step_Velocity
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Advance_One_Time_Step_Velocity(const T dt,const T time, const bool solids)
{
    LOG::SCOPE scope("solids velocity update",rigids_parameters.threadid);
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Advance_One_Time_Step_Velocity Start dt=%f time=%f",dt,time),2,2);

    ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body=rigid_body_collection.articulated_rigid_body;
    const bool advance_rigid_bodies=true; //solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m!=0;  TODO: Fix this.

    if(solids){
        if(advance_rigid_bodies){
            articulated_rigid_body.Apply_Poststabilization(rigids_parameters.implicit_solve_parameters.test_system,rigids_parameters.implicit_solve_parameters.print_matrix);
            Diagnostics(dt,time,0,2,24,"poststabilization");}
        // initialize data needed for rigid/deformable contact projection in CG
        
        if(rigids_parameters.rigid_body_collision_parameters.perform_contact) rigid_body_collisions->Initialize_All_Contact_Projections(rigids_parameters.rigid_body_collision_parameters.enforce_rigid_rigid_contact_in_cg);
        
        PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("after creating joints.  before trapezoidal velocities dt=%f time=%f",dt,time),2,2);}

    if(rigids_parameters.use_trapezoidal_rule_for_velocities){
        Average_And_Exchange_Position(); // move to positions at time+dt/2 for trapezoidal step
        Diagnostics(dt,time,0,1,25,"average and exchange positions");
        
        rigid_body_collection.Update_Position_Based_State(time+dt/2);
        
        //TODO: add callbacks
        //if(articulated_rigid_body.Has_Actuators()) example_forces_and_velocities.Set_PD_Targets(dt,time);
        
        Trapezoidal_Step_Velocity(dt,time);
        Diagnostics(dt,time,2,1,29,"trazepoid rule");

        Restore_Position(rigid_X_save,rigid_rotation_save); // move to final positions at time time+dt
        Diagnostics(dt,time,2,2,30,"restore position");
        if(advance_rigid_bodies){
            articulated_rigid_body.Apply_Poststabilization(rigids_parameters.implicit_solve_parameters.test_system,rigids_parameters.implicit_solve_parameters.print_matrix);
            Diagnostics(dt,time,2,2,31,"poststabilization");}
    }else{
        rigid_body_collection.Update_Position_Based_State(time+dt);
        
        //TODO: add callbacks
        //if(articulated_rigid_body.Has_Actuators()) example_forces_and_velocities.Set_PD_Targets(dt,time);
        
        Backward_Euler_Step_Velocity(dt,time); // TODO: Tamar & Craig, do you need post stab?
        Diagnostics(dt,time,2,2,29,"backward Euler");}

    if(solids){
        PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("after removing joints.  after trapezoidal velocities dt=%f time=%f",dt,time),2,2);

        if(rigids_parameters.use_post_cg_constraints) Apply_Constraints(dt,time);
        if(advance_rigid_bodies && rigids_parameters.rigid_body_evolution_parameters.clamp_rigid_body_velocities){
            Clamp_Velocities(); // TODO: Examples should do this during the Advance_One_Time_Step_End_Callback example callback
            Diagnostics(dt,time,2,2,41,"clamp velocities");}
        if(rigids_parameters.verbose) Print_Maximum_Velocities(time);}

    //TODO: add callbacks
    //example_forces_and_velocities.Advance_One_Time_Step_End_Callback(dt,time);

    if(mpi_rigids)
        mpi_rigids->Update_Partitions(rigid_body_collection,*rigid_body_collection.rigid_geometry_collection.collision_body_list->spatial_partition);

    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Advance_One_Time_Step_Velocity End dt=%f time=%f",dt,time),2,2);
}
//#####################################################################
// Function Apply_Constraints
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Apply_Constraints(const T dt,const T time)
{
    LOG::SCOPE scope("Apply_Constraints",rigids_parameters.threadid);
    ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body=rigid_body_collection.articulated_rigid_body;
    const bool advance_rigid_bodies=rigid_body_collection.simulated_rigid_body_particles.m!=0;
    if(advance_rigid_bodies){
        articulated_rigid_body.Apply_Poststabilization(rigids_parameters.implicit_solve_parameters.test_system,rigids_parameters.implicit_solve_parameters.print_matrix);
        Diagnostics(dt,time,0,0,345,"poststabilization");
        if(articulated_rigid_body.Has_Actuators() && !articulated_rigid_body.constrain_pd_directions){
            articulated_rigid_body.Compute_Position_Based_State(dt,time);
            articulated_rigid_body.Solve_Velocities_for_PD(time,dt,rigids_parameters.implicit_solve_parameters.test_system,rigids_parameters.implicit_solve_parameters.print_matrix);}
        Diagnostics(dt,time,2,2,35,"solve velocities for pd");}
    Save_Position(rigid_X_save_for_constraints,rigid_rotation_save_for_constraints);

    Diagnostics(dt,time,2,2,137,"consistent contact");
    ARRAY<TWIST<TV> > rigid_velocity_save_mpi,rigid_velocity_mpi;
    ARRAY<T_SPIN> rigid_angular_momentum_save_mpi,rigid_angular_momentum_mpi;
    ARRAY<TV> rigid_X_save_mpi,rigid_X_mpi;
    ARRAY<ROTATION<TV> > rigid_rotation_save_mpi,rigid_rotation_mpi;
    ARRAY<ARRAY<int>,PARTITION_ID> particles_of_partition;
    ARRAY<PARTITION_ID> partition_id_from_particle_index;
    Update_Positions_And_Apply_Contact_Forces(dt,time,true);
    Diagnostics(dt,time,4,2,37,"contact, prestabilization");
    if(rigid_body_collisions->prune_stacks_from_contact) rigid_body_collisions->Construct_Stacks();
    if(rigid_body_collisions->prune_contact_using_velocity) rigid_body_collisions->Compute_Contact_Frequency();
    Restore_Position(rigid_X_save_for_constraints,rigid_rotation_save_for_constraints);
    Diagnostics(dt,time,2,2,38,"restore position");
    if(advance_rigid_bodies){
        articulated_rigid_body.Apply_Poststabilization(rigids_parameters.implicit_solve_parameters.test_system,rigids_parameters.implicit_solve_parameters.print_matrix);
        Diagnostics(dt,time,2,2,39,"poststabilization");}
    // modify velocity with inelastic and friction repulsions
}
//#####################################################################
// Function Print_Maximum_Velocities
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Print_Maximum_Velocities(const T time) const
{
    {std::stringstream ss;ss<<"time = "<<time<<std::endl;LOG::filecout(ss.str(),rigids_parameters.threadid);}
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=rigid_body_collection.rigid_body_particle;
    int max_linear_index=0,max_angular_index=0;T max_linear_magnitude_squared=-FLT_MAX,max_angular_magnitude_squared=-FLT_MAX;
    for(int i=1;i<=rigid_body_collection.dynamic_rigid_body_particles.m;i++){const int p=rigid_body_collection.dynamic_rigid_body_particles(i);
        T linear_magnitude_squared=rigid_body_particles.V(p).Magnitude_Squared(),angular_magnitude_squared=rigid_body_particles.angular_velocity(p).Magnitude_Squared();
        if(linear_magnitude_squared>max_linear_magnitude_squared){max_linear_magnitude_squared=linear_magnitude_squared;max_linear_index=p;}
        if(angular_magnitude_squared>max_angular_magnitude_squared){max_angular_magnitude_squared=angular_magnitude_squared;max_angular_index=p;}}
    if(max_linear_index){
        T max_linear_magnitude=sqrt(max_linear_magnitude_squared),max_angular_magnitude=sqrt(max_angular_magnitude_squared);
        T max_linear_magnitude_global=max_linear_magnitude,max_angular_magnitude_global=max_angular_magnitude;
        {std::stringstream ss;ss<<"maximum rigid linear velocity = "<<max_linear_magnitude_global;LOG::filecout(ss.str(),rigids_parameters.threadid);}
        {std::stringstream ss;ss<<" ("<<max_linear_index<<")\n";LOG::filecout(ss.str(),rigids_parameters.threadid);}
        {std::stringstream ss;ss<<"maximum rigid angular velocity = "<<max_angular_magnitude_global;LOG::filecout(ss.str(),rigids_parameters.threadid);}
        {std::stringstream ss;ss<<" ("<<max_angular_index<<")"<<std::endl;LOG::filecout(ss.str(),rigids_parameters.threadid);}}
}
//#####################################################################
// Function Diagnostics
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Diagnostics(const T dt,const T time,const int velocity_time,const int position_time,int step,const char* description)
{
    static const char* time_table[]={"n","(n+1/2)","(n+1)","(n+3/2)","(n+2)"};
    rigid_body_collection.Print_Energy(time+position_time*(T).5*time,step);
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Finished step %i (%s).  State: x^%s  v^%s.   dt=%f time=%f",step,description,
        time_table[position_time],time_table[velocity_time],dt,time),2,3);
}
//#####################################################################
// Function Update_Positions_And_Apply_Contact_Forces
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Update_Positions_And_Apply_Contact_Forces(const T dt,const T time,const bool use_saved_pairs)
{
    ARTICULATED_RIGID_BODY<TV>* articulated_rigid_body=&rigid_body_collection.articulated_rigid_body;

    if((!use_saved_pairs || !rigids_parameters.rigid_body_collision_parameters.use_projected_gauss_seidel)) Euler_Step_Position(dt,time);

    if(!rigids_parameters.rigid_body_collision_parameters.perform_contact) return;
   
    for(int i=1;i<=rigid_body_collection.simulated_rigid_body_particles.m;i++){
        RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(rigid_body_collection.simulated_rigid_body_particles(i));
        if(!rigid_body.is_static) rigid_body.Update_Bounding_Box();}
    for(int i=1;i<=rigid_body_collection.kinematic_rigid_bodies.m;i++){
        RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(rigid_body_collection.kinematic_rigid_bodies(i));
        rigid_body.Update_Bounding_Box();}

    rigid_body_collisions->Process_Contact_Using_Graph(dt,time,articulated_rigid_body,rigids_parameters.rigid_body_evolution_parameters.correct_contact_energy,use_saved_pairs);
    
    // rigid/rigid shock propagation
    ARRAY<TWIST<TV> > rigid_velocity_save_mpi,rigid_velocity_mpi;
    ARRAY<T_SPIN> rigid_angular_momentum_save_mpi,rigid_angular_momentum_mpi;
    ARRAY<TV> rigid_X_save_mpi,rigid_X_mpi;
    ARRAY<ROTATION<TV> > rigid_rotation_save_mpi,rigid_rotation_mpi;
    ARRAY<ARRAY<int>,PARTITION_ID> particles_of_partition;
    ARRAY<PARTITION_ID> partition_id_from_particle_index;

    if(rigids_parameters.rigid_body_collision_parameters.use_shock_propagation){
        rigid_body_collisions->Shock_Propagation_Using_Graph(dt,time,articulated_rigid_body,use_saved_pairs);}
}
//#####################################################################
// Function Update_Velocity_Using_Stored_Differences
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Update_Velocity_Using_Stored_Differences(const T dt,const T time,const int p)
{
    RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(p);
    if(rigid_body.is_static) return;
    if(rigid_body_collection.rigid_body_particle.kinematic(p)) kinematic_evolution.Set_External_Velocities(rigid_body.V(),rigid_body.Angular_Velocity(),time+dt,p);
    rigid_body.V()+=rigid_velocity_difference(p);
    rigid_body.Angular_Momentum()+=rigid_angular_momentum_difference(p);
    rigid_body.Update_Angular_Velocity();
}
//#####################################################################
// Function Update_Velocity_Using_Stored_Differences
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Update_Velocity_Using_Stored_Differences(const T dt,const T time)
{
    LOG::SCOPE scope("update velocity using stored differences",rigids_parameters.threadid);
    for(int i=1;i<=rigid_body_collection.simulated_rigid_body_particles.m;i++){int p=rigid_body_collection.simulated_rigid_body_particles(i);
        Update_Velocity_Using_Stored_Differences(dt,time,p);}
    for(int i=1;i<=rigid_body_collection.kinematic_rigid_bodies.m;i++){int p=rigid_body_collection.kinematic_rigid_bodies(i);
        kinematic_evolution.Set_External_Velocities(rigid_body_collection.rigid_body_particle.V(p),rigid_body_collection.rigid_body_particle.angular_velocity(p),time+dt,p);}
}
//#####################################################################
// Function Compute_Momentum_Differences
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Compute_Momentum_Differences()
{
    LOG::SCOPE scope("compute momentum difference",rigids_parameters.threadid);
    rigid_velocity_difference.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size());
    rigid_angular_momentum_difference.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size());
    for(int i=1;i<=rigid_body_collection.simulated_rigid_body_particles.m;i++){int p=rigid_body_collection.simulated_rigid_body_particles(i);
        RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(p);
        rigid_velocity_difference(p)=rigid_body.V()-rigid_velocity_save(p).linear;
        rigid_angular_momentum_difference(p)=rigid_body.Angular_Momentum()-rigid_angular_momentum_save(p);}
}
//#####################################################################
// Function Save_Velocity
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Save_Velocity()
{
    Save_Velocity(rigid_velocity_save,rigid_angular_momentum_save);
}
template<class TV> void RIGIDS_EVOLUTION<TV>::
Save_Velocity(ARRAY<TWIST<TV> >& rigid_velocity_save,ARRAY<T_SPIN>& rigid_angular_momentum_save)
{
    rigid_velocity_save.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
    rigid_angular_momentum_save.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
    for(int i=1;i<=rigid_body_collection.simulated_rigid_body_particles.m;i++){int p=rigid_body_collection.simulated_rigid_body_particles(i);
        rigid_velocity_save(p)=rigid_body_collection.Rigid_Body(p).Twist();rigid_angular_momentum_save(p)=rigid_body_collection.rigid_body_particle.angular_momentum(p);}
    for(int i=1;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++) if(rigid_body_collection.Is_Active(i)){RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(i);
        if(!body.Is_Simulated()){rigid_velocity_save(i)=rigid_body_collection.Rigid_Body(i).Twist();rigid_angular_momentum_save(i)=rigid_body_collection.rigid_body_particle.angular_momentum(i);}}
}
//#####################################################################
// Function Restore_Velocity
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Restore_Velocity()
{
    Restore_Velocity(rigid_velocity_save,rigid_angular_momentum_save);
}
template<class TV> void RIGIDS_EVOLUTION<TV>::
Restore_Velocity(ARRAY<TWIST<TV> >& rigid_velocity_save,ARRAY<T_SPIN>& rigid_angular_momentum_save)
{
    for(int i=1;i<=rigid_body_collection.simulated_rigid_body_particles.m;i++){int p=rigid_body_collection.simulated_rigid_body_particles(i);
        rigid_body_collection.rigid_body_particle.V(p)=rigid_velocity_save(p).linear;
        rigid_body_collection.rigid_body_particle.angular_momentum(p)=rigid_angular_momentum_save(p);rigid_body_collection.Rigid_Body(p).Update_Angular_Velocity();}
    for(int i=1;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++) if(rigid_body_collection.Is_Active(i)){RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(i);
        if(!body.Is_Simulated()){
            rigid_body_collection.rigid_body_particle.V(i)=rigid_velocity_save(i).linear;rigid_body_collection.rigid_body_particle.angular_momentum(i)=rigid_angular_momentum_save(i);body.Update_Angular_Velocity();}}
}
//#####################################################################
// Function Exchange_Velocity
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Exchange_Velocity()
{
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=rigid_body_collection.rigid_body_particle;
    rigid_velocity_save.Resize(rigid_body_particles.array_collection->Size(),false,false);
    rigid_angular_momentum_save.Resize(rigid_body_particles.array_collection->Size(),false,false);
    for(int i=1;i<=rigid_body_collection.simulated_rigid_body_particles.m;i++){int p=rigid_body_collection.simulated_rigid_body_particles(i);
        exchange(rigid_velocity_save(p).linear,rigid_body_particles.V(p));
        exchange(rigid_velocity_save(p).angular,rigid_body_particles.angular_velocity(p));
        exchange(rigid_angular_momentum_save(p),rigid_body_particles.angular_momentum(p));}
    for(int i=1;i<=rigid_body_particles.array_collection->Size();i++) if(rigid_body_collection.Is_Active(i)){RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(i);
        if(!body.Is_Simulated()){
            exchange(rigid_velocity_save(i).linear,rigid_body_particles.V(i));
            exchange(rigid_velocity_save(i).angular,rigid_body_particles.angular_velocity(i));
            exchange(rigid_angular_momentum_save(i),rigid_body_particles.angular_momentum(i));}}
}
//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Initialize_Rigid_Bodies(const T frame_rate, const bool restart)
{
    // initialize kinematic object positions and velocities
    if(!restart){
        kinematic_evolution.Get_Current_Kinematic_Keyframes(1/frame_rate,time);
        kinematic_evolution.Set_External_Positions(rigid_body_collection.rigid_body_particle.X,rigid_body_collection.rigid_body_particle.rotation,time);
        kinematic_evolution.Set_External_Velocities(rigid_body_collection.rigid_body_particle.V,rigid_body_collection.rigid_body_particle.angular_velocity,time,time);
        rigid_body_collection.Update_Angular_Momentum();
        for(int i(1);i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++) if(rigid_body_collection.Is_Active(i)){rigid_body_collection.rigid_body_particle.rotation(i).Normalize();}}

    RIGID_BODY_COLLISIONS<TV>::Adjust_Bounding_Boxes(rigid_body_collection);
    // rigid body collisions
    rigid_body_collisions=new RIGID_BODY_COLLISIONS<TV>(rigid_body_collection,rigids_parameters.rigid_body_collision_parameters,*rigids_collision_callbacks,
        *rigid_body_collection.rigids_example_forces_and_velocities);
    rigid_body_collisions->spatial_partition->Compute_Voxel_Size(rigids_parameters.rigid_body_collision_parameters.rigid_collisions_spatial_partition_voxel_size_heuristic,
        rigids_parameters.rigid_body_collision_parameters.rigid_collisions_spatial_partition_number_of_cells,rigids_parameters.rigid_body_collision_parameters.rigid_collisions_spatial_partition_voxel_size_scale_factor);
    rigid_body_collisions->verbose=rigids_parameters.verbose;
    rigid_body_collisions->mpi_rigids=mpi_rigids;
    if(mpi_rigids){
        rigid_body_collisions->spatial_partition->voxel_size=mpi_rigids->Reduce_Min(rigid_body_collisions->spatial_partition->voxel_size);
        rigid_body_collisions->spatial_partition->one_over_voxel_size=1/rigid_body_collisions->spatial_partition->voxel_size;}

    // partitions and hierarchies
    if(rigids_parameters.rigid_body_collision_parameters.rigid_collisions_use_particle_partition){
        rigid_body_collisions->intersections.Use_Particle_Partition(true,rigids_parameters.rigid_body_collision_parameters.rigid_collisions_particle_partition_size*VECTOR<int,d>::All_Ones_Vector());
        rigid_body_collisions->intersections.Use_Particle_Partition_Center_Phi_Test(rigids_parameters.rigid_body_collision_parameters.rigid_collisions_use_particle_partition_center_phi_test);}
    if(rigids_parameters.rigid_body_collision_parameters.rigid_collisions_use_triangle_hierarchy){
        rigid_body_collisions->intersections.Use_Triangle_Hierarchy();
        if(rigids_parameters.rigid_body_collision_parameters.rigid_collisions_use_triangle_hierarchy_center_phi_test) rigid_body_collisions->intersections.Use_Triangle_Hierarchy_Center_Phi_Test();
        if(rigids_parameters.rigid_body_collision_parameters.rigid_collisions_use_edge_intersection) rigid_body_collisions->intersections.Use_Edge_Intersection();}

    // dynamics
    rigids_parameters.rigid_body_evolution_parameters.rigid_cfl=rigids_parameters.cfl;
    rigids_parameters.rigid_body_evolution_parameters.rigid_minimum_dt=rigids_parameters.rigid_body_evolution_parameters.minimum_rigid_body_time_step_fraction/frame_rate;
    rigids_parameters.rigid_body_evolution_parameters.rigid_maximum_dt=rigids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction/frame_rate;
    rigids_evolution_callbacks->Update_Rigids_Parameters(time);
    rigid_body_collisions->Initialize_Data_Structures(); // Must be done before we check interpenetration statistics

    // check for bad initial data (need to give it a chance to set up the collision manager first)
    if(rigids_parameters.rigid_body_collision_parameters.rigid_collisions_print_interpenetration_statistics) rigid_body_collisions->Print_Interpenetration_Statistics();
    else if(!rigid_body_collisions->Check_For_Any_Interpenetration()) {std::stringstream ss;ss<<"No initial interpenetration"<<std::endl;LOG::filecout(ss.str(),rigids_parameters.threadid);}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time)
{
    INDIRECT_ARRAY<ARRAY_VIEW<TWIST<TV> >,ARRAY<int>&> twist_subset=twist.Subset(rigid_body_collection.static_and_kinematic_rigid_bodies);
    ARRAYS_COMPUTATIONS::Fill(twist_subset,TWIST<TV>());
    rigid_body_collection.rigids_example_forces_and_velocities->Zero_Out_Enslaved_Velocity_Nodes(twist,velocity_time,current_position_time);
}
//#####################################################################
// Function Correct_Orientation_For_Kinetic_Energy
//#####################################################################
template<class T>
inline bool Correct_Orientation_For_Kinetic_Energy(RIGID_BODY<VECTOR<T,3> >& rigid_body,const T KE0,int iteration=1)
{
    const int max_orientation_correction_iterations=40;
    if(iteration>max_orientation_correction_iterations) PHYSBAM_FATAL_ERROR("Exceeded maximum number of iterations during solids evolution energy clamping");
    typedef VECTOR<T,3> TV;const TV &w=rigid_body.Twist().angular,&L=rigid_body.Angular_Momentum();
    int status=Correct_Orientation_For_Kinetic_Energy_Using_Direction(rigid_body,KE0,TV::Cross_Product(L,w),false);
    if(status<=0) return true;
    int component;
    if((T).5*TV::Dot_Product(w,L)>KE0) component=rigid_body.Inertia_Tensor().To_Vector().Arg_Max();
    else component=rigid_body.Inertia_Tensor().To_Vector().Arg_Min();
    TV axis=rigid_body.Rotation().Rotate(TV::Axis_Vector(component));
    status=Correct_Orientation_For_Kinetic_Energy_Using_Direction(rigid_body,KE0,TV::Cross_Product(axis,L),true);
    if(status==3) return Correct_Orientation_For_Kinetic_Energy(rigid_body,KE0,iteration+1);
    return status<=0;
}
template<class T>
inline int Correct_Orientation_For_Kinetic_Energy_Using_Direction(RIGID_BODY<VECTOR<T,3> >& rigid_body,const T KE0,const VECTOR<T,3>& direction,const bool use_extrema)
{
    typedef VECTOR<T,3> TV;ROTATION<TV>& R=rigid_body.Rotation();const TV &w=rigid_body.Twist().angular,&L=rigid_body.Angular_Momentum();
    T KE=(T).5*TV::Dot_Product(w,L),k=KE-KE0,error=abs(k);
    if(error<=4*std::numeric_limits<T>::epsilon()*KE0) return 0;
    if(KE0<(T)sqr(std::numeric_limits<T>::epsilon())) return -1;
    TV s=direction.Normalized();T a=TV::Dot_Product(s,TV::Cross_Product(L,w));
    if(a<0){a=-a;s=-s;}
    if(a<16*std::numeric_limits<T>::epsilon()*KE) return 1;
    TV s_cross_L=TV::Cross_Product(s,L);
    T b=(T).5*TV::Dot_Product(rigid_body.World_Space_Inertia_Tensor_Inverse_Times(s_cross_L),s_cross_L);
    T c=KE-b,r=sqrt(sqr(a)+sqr(c)),n=c-2*k,q2=sqr(a)+4*k*(c-k);
    if(q2>=0){
        T q=a*sqrt(q2),p=sqr(a)+q2+2*q;
        T x=(T).5*sqrt(p+sqr(c+n))/r,y=2*fabs(k)/sqrt(p+sqr(c-n));
        if(k-c*sqr(y)<0) y=-y;
        R=ROTATION<TV>::From_Rotation_Vector(atan2(y,x)*s)*R;R.Normalize();
        rigid_body.Update_Angular_Velocity();
        return -2;}
    if(use_extrema){
        if(error<=16*std::numeric_limits<T>::epsilon()*KE0) return 2;
        T sk=sign(k),e=sk*(sk*c>0?atan2(r+sk*c,a):atan2(a,r-sk*c));
        R=ROTATION<TV>::From_Rotation_Vector(e*s)*R;R.Normalize();
        rigid_body.Update_Angular_Velocity();
        return 3;}
    return 4;
}
//#####################################################################
// Function Update_Rotation_Helper
//#####################################################################
template<class T>
inline void Update_Rotation_Helper(const T dt,RIGID_BODY<VECTOR<T,3> >& rigid_body,bool correct_evolution_energy)
{
    typedef VECTOR<T,3> TV;ROTATION<TV>& R=rigid_body.Rotation();const TV &w=rigid_body.Twist().angular,&L=rigid_body.Angular_Momentum();
    T KE0=(T).5*TV::Dot_Product(w,L);
    TV rotate_amount=w-(T).5*dt*rigid_body.World_Space_Inertia_Tensor_Inverse_Times(TV::Cross_Product(w,L));
    R=ROTATION<TV>::From_Rotation_Vector(dt*rotate_amount)*R;R.Normalize();
    rigid_body.Update_Angular_Velocity(); // Note that the value of w changes here.
    if(correct_evolution_energy) Correct_Orientation_For_Kinetic_Energy(rigid_body,KE0);
}
//#####################################################################
// Function Update_Rotation_Helper
//#####################################################################
template<class T>
inline void Update_Rotation_Helper(const T dt,RIGID_BODY<VECTOR<T,2> >& rigid_body,bool correct_evolution_energy)
{
    rigid_body.Rotation()=ROTATION<VECTOR<T,2> >::From_Rotation_Vector(dt*rigid_body.Twist().angular)*rigid_body.Rotation();rigid_body.Rotation().Normalize();
}
//#####################################################################
// Function Update_Rotation_Helper
//#####################################################################
template<class T>
inline void Update_Rotation_Helper(const T dt,RIGID_BODY<VECTOR<T,1> >& rigid_body,bool correct_evolution_energy)
{}
//#####################################################################
// Euler_Step_Position
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Euler_Step_Position(const T dt,const T time,const int p)
{
    RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(p);
    if(rigid_body.is_static) return;
    else if(rigid_body_collection.rigid_body_particle.kinematic(p)){
        kinematic_evolution.Set_External_Positions(rigid_body.X(),rigid_body.Rotation(),time+dt,p);}
    else{
        rigid_body.X()+=dt*rigid_body.Twist().linear;
        Update_Rotation_Helper(dt,rigid_body,rigids_parameters.rigid_body_evolution_parameters.correct_evolution_energy);}
}
//#####################################################################
// Function Euler_Step_Position
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Euler_Step_Position(const T dt,const T time)
{
    for(int i=1;i<=rigid_body_collection.dynamic_rigid_body_particles.m;i++)
        Euler_Step_Position(dt,time,rigid_body_collection.dynamic_rigid_body_particles(i)); // TODO: avoid unnecessary Update_Angular_Velocity
    kinematic_evolution.Set_External_Positions(rigid_body_collection.rigid_body_particle.X,rigid_body_collection.rigid_body_particle.rotation,time+dt);
    rigid_body_collection.Update_Angular_Velocity(); // Note: Possibly remove as we restore velocities right after this function.
    rigid_body_collection.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Positions();
}
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Initialize_World_Space_Masses()
{
    world_space_rigid_mass.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
    world_space_rigid_mass_inverse.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);

    for(int i=1;i<=rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=rigid_body_collection.dynamic_rigid_body_particles(i);
        world_space_rigid_mass(p)=rigid_body_collection.State(p).World_Space_Rigid_Mass(RIGID_BODY_MASS<TV>(rigid_body_collection.rigid_body_particle.mass(p),rigid_body_collection.rigid_body_particle.inertia_tensor(p)));
        world_space_rigid_mass_inverse(p)=rigid_body_collection.State(p).World_Space_Rigid_Mass_Inverse(RIGID_BODY_MASS<TV>(rigid_body_collection.rigid_body_particle.mass(p),rigid_body_collection.rigid_body_particle.inertia_tensor(p)));}
}
//#####################################################################
// Function Clamp_Velocities
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Clamp_Velocities()
{
    const ARRAY<int>& dynamic_rigid_body_particles=rigid_body_collection.dynamic_rigid_body_particles;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=rigid_body_collection.rigid_body_particle;
    T max_linear_velocity_squared=sqr(rigids_parameters.rigid_body_evolution_parameters.max_rigid_body_linear_velocity),max_angular_velocity_squared=sqr(rigids_parameters.rigid_body_evolution_parameters.max_rigid_body_angular_velocity);
    for(int i=1;i<=dynamic_rigid_body_particles.m;i++){int p=dynamic_rigid_body_particles(i);
        T magnitude_squared=rigid_body_particles.V(p).Magnitude_Squared();
        if(magnitude_squared>max_linear_velocity_squared) rigid_body_particles.V(p)*=rigids_parameters.rigid_body_evolution_parameters.max_rigid_body_linear_velocity/sqrt(magnitude_squared);
        magnitude_squared=rigid_body_particles.angular_velocity(p).Magnitude_Squared();
        if(magnitude_squared>max_angular_velocity_squared){
            rigid_body_particles.angular_velocity(p)*=rigids_parameters.rigid_body_evolution_parameters.max_rigid_body_angular_velocity/sqrt(magnitude_squared);
            rigid_body_collection.Rigid_Body(p).Update_Angular_Momentum();}}
}
//#####################################################################
// Function Save_Position
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Save_Position(ARRAY<TV>& rigid_X,ARRAY<ROTATION<TV> >& rigid_rotation)
{
    const ARRAY<int>& simulated_rigid_body_particles=rigid_body_collection.simulated_rigid_body_particles;
    rigid_X.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
    rigid_rotation.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
    rigid_X.Subset(simulated_rigid_body_particles)=rigid_body_collection.rigid_body_particle.X.Subset(simulated_rigid_body_particles);
    rigid_rotation.Subset(simulated_rigid_body_particles)=rigid_body_collection.rigid_body_particle.rotation.Subset(simulated_rigid_body_particles);
    for(int i=1;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++) if(rigid_body_collection.Is_Active(i)){
        if(!rigid_body_collection.Rigid_Body(i).Is_Simulated()){rigid_X(i)=rigid_body_collection.rigid_body_particle.X(i);rigid_rotation(i)=rigid_body_collection.rigid_body_particle.rotation(i);}}
} 
//#####################################################################
// Function Restore_Position
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
Restore_Position(ARRAY_VIEW<const TV> rigid_X,ARRAY_VIEW<const ROTATION<TV> > rigid_rotation)
{
    const ARRAY<int>& simulated_rigid_body_particles=rigid_body_collection.simulated_rigid_body_particles;
    PHYSBAM_ASSERT(rigid_X.Size()==rigid_body_collection.rigid_body_particle.array_collection->Size());
    PHYSBAM_ASSERT(rigid_rotation.Size()==rigid_body_collection.rigid_body_particle.array_collection->Size());
    rigid_body_collection.rigid_body_particle.X.Subset(simulated_rigid_body_particles)=rigid_X.Subset(simulated_rigid_body_particles);
    rigid_body_collection.rigid_body_particle.rotation.Subset(simulated_rigid_body_particles)=rigid_rotation.Subset(simulated_rigid_body_particles);
    for(int i=1;i<=simulated_rigid_body_particles.m;i++) rigid_body_collection.Rigid_Body(simulated_rigid_body_particles(i)).Update_Angular_Velocity();
    for(int i=1;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++) if(rigid_body_collection.Is_Active(i)){RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(i);
        if(!body.Is_Simulated()){rigid_body_collection.rigid_body_particle.X(i)=rigid_X(i);rigid_body_collection.rigid_body_particle.rotation(i)=rigid_rotation(i);body.Update_Angular_Velocity();}}
}
//#####################################################################
// Function CFL
//#####################################################################
template<class TV> void RIGIDS_EVOLUTION<TV>::
CFL(const bool verbose)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
template class RIGIDS_EVOLUTION<VECTOR<float,1> >;
template class RIGIDS_EVOLUTION<VECTOR<float,2> >;
template class RIGIDS_EVOLUTION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGIDS_EVOLUTION<VECTOR<double,1> >;
template class RIGIDS_EVOLUTION<VECTOR<double,2> >;
template class RIGIDS_EVOLUTION<VECTOR<double,3> >;
#endif

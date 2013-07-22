//#####################################################################
// Copyright 2006-2009, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLES_EVOLUTION
//#####################################################################
#include <PhysBAM_Tools/Arrays/CONSTANT_ARRAY.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Krylov_Solvers/LANCZOS_ITERATION.h>
#include <PhysBAM_Tools/Krylov_Solvers/SYMMQMR.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Geometry/Geometry_Particles/RIGID_GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLES_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformables_Evolution/DEFORMABLES_BACKWARD_EULER_SYSTEM.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformables_Evolution/DEFORMABLES_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformables_Evolution/DEFORMABLES_VELOCITY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BINDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <stdexcept>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLES_EVOLUTION<TV>::
DEFORMABLES_EVOLUTION(DEFORMABLES_PARAMETERS<TV>& deformables_parameters_input,DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection_input,RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection_input)
    :time(0),deformable_body_collection(deformable_body_collection_input),rigid_geometry_collection(rigid_geometry_collection_input),
    kinematic_evolution(rigid_geometry_collection,deformables_parameters_input.rigid_geometry_evolution_parameters.use_kinematic_keyframes),repulsions(0),
    deformables_parameters(deformables_parameters_input),use_existing_contact(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DEFORMABLES_EVOLUTION<TV>::
~DEFORMABLES_EVOLUTION()
{
}
//#####################################################################
// Function Prepare_Backward_Euler_System
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Prepare_Backward_Euler_System(DEFORMABLES_BACKWARD_EULER_SYSTEM<TV>& system,const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update)
{
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_GEOMETRY_PARTICLES<TV>& rigid_geometry_particles=rigid_geometry_collection.particles;
    MPI_SOLIDS<TV>* mpi_solids=deformable_body_collection.mpi_solids;

    F_full.Resize(particles.array_collection->Size(),false,false);rigid_F_full.Resize(rigid_geometry_particles.array_collection->Size(),false,false);
    R_full.Resize(particles.array_collection->Size(),false,false);rigid_R_full.Resize(rigid_geometry_particles.array_collection->Size(),false,false);
    S_full.Resize(particles.array_collection->Size(),false,false);rigid_S_full.Resize(rigid_geometry_particles.array_collection->Size(),false,false);
    B_full.Resize(particles.array_collection->Size(),false,false);rigid_B_full.Resize(rigid_geometry_particles.array_collection->Size(),false,false);
    if(deformables_parameters.implicit_solve_parameters.evolution_solver_type!=krylov_solver_cg){
        AR_full.Resize(particles.array_collection->Size(),false,false);rigid_AR_full.Resize(rigid_geometry_particles.array_collection->Size(),false,false);}

    DEFORMABLES_VELOCITY<TV> B_all(B_full,rigid_B_full,deformable_body_collection,rigid_geometry_collection);
    DEFORMABLES_VELOCITY<TV> F_all(F_full,rigid_F_full,deformable_body_collection,rigid_geometry_collection);
    ARRAY<TWIST<TV> > twist;twist.Resize(rigid_geometry_particles.V.Size());
    for(int i=1;i<=twist.Size();i++) twist(i)=TWIST<TV>(rigid_geometry_particles.V(i),rigid_geometry_particles.angular_velocity(i));
    DEFORMABLES_VELOCITY<TV> V_all(particles.V,twist,deformable_body_collection,rigid_geometry_collection);

    INDIRECT_ARRAY<ARRAY<TV>,ARRAY<int>&> B_subset=B_full.Subset(deformable_body_collection.simulated_particles);
    ARRAYS_COMPUTATIONS::Fill(B_subset,TV());
    deformable_body_collection.deformables_example_forces_and_velocities->Add_External_Forces(B_full,current_velocity_time+dt);
    if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data_Global(particles.V);
    deformable_body_collection.Add_Velocity_Independent_Forces(B_full,current_velocity_time+dt); // this is a nop for binding forces
    if(mpi_solids) mpi_solids->Exchange_Binding_Boundary_Data_Global(B_full);
    deformable_body_collection.binding_list.Distribute_Force_To_Parents(B_full,rigid_B_full);

    if(deformable_body_collection.soft_bindings.Need_Bindings_Mapped()){
        if(mpi_solids) mpi_solids->Exchange_Binding_Boundary_Data_Global(B_full);
        deformable_body_collection.soft_bindings.Map_Forces_From_Parents(B_full,rigid_B_full);
        deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(B_full);
        if(mpi_solids) mpi_solids->Exchange_Binding_Boundary_Data_Global(B_full);
        for(int k=1;k<=deformable_body_collection.deformables_forces.m;k++) if(dynamic_cast<BINDING_SPRINGS<TV>*>(&*deformable_body_collection.deformables_forces(k)))
            deformable_body_collection.deformables_forces(k)->Add_Force_Differential(particles.X,B_full,current_velocity_time+dt);
        if(mpi_solids) mpi_solids->Exchange_Binding_Boundary_Data_Global(B_full);
        deformable_body_collection.binding_list.Distribute_Force_To_Parents(B_full,rigid_B_full);}

    for(int i=1;i<=deformable_body_collection.dynamic_particles.m;i++){int p=deformable_body_collection.dynamic_particles(i);
        B_full(p)=particles.V(p)+dt*particles.one_over_mass(p)*B_full(p);}

    V_all=B_all;
    Diagnostics(dt,current_position_time,0,0,604,"Before boundary conditions");
    system.Set_Global_Boundary_Conditions(V_all);
    if(velocity_update && !deformables_parameters.use_post_cg_constraints) Apply_Constraints(dt,current_velocity_time);
    for(int i=1;i<=twist.Size();i++){rigid_geometry_particles.V(i)=twist(i).linear;rigid_geometry_particles.angular_velocity(i)=twist(i).angular;}
}
//#####################################################################
// Function Prepare_Backward_Euler_System
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Finish_Backward_Euler_Step(const T dt,const T current_position_time,const bool velocity_update)
{
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    DEFORMABLES_VELOCITY<TV> F_all(F_full,rigid_F_full,deformable_body_collection,rigid_geometry_collection);

    if(velocity_update && deformables_parameters.use_post_cg_constraints){ // return rhs + dt Fd V^n+1 for friction processing
        for(int i=1;i<=deformable_body_collection.dynamic_particles.m;i++){int p=deformable_body_collection.dynamic_particles(i);
            particles.V(p)=B_full(p)+dt*particles.one_over_mass(p)*F_full(p);}}

    deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities(); // TODO: MPI safe?
}
//#####################################################################
// Function Backward_Euler_Step_Velocity_Helper
//#####################################################################
// assumes all solids_forces are linear in velocity, with a symmetric positive definite Jacobian.
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Backward_Euler_Step_Velocity_Helper(const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update)
{
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_GEOMETRY_PARTICLES<TV>& rigid_geometry_particles=rigid_geometry_collection.particles;
    MPI_SOLIDS<TV>* mpi_solids=deformable_body_collection.mpi_solids;
    DEFORMABLES_BACKWARD_EULER_SYSTEM<TV> system(*this,deformable_body_collection,dt,current_velocity_time,current_position_time,mpi_solids,
        (velocity_update && deformables_parameters.enforce_repulsions_in_cg)?repulsions:0,velocity_update);

    Prepare_Backward_Euler_System(system,dt,current_velocity_time,current_position_time,velocity_update);

    static CONJUGATE_GRADIENT<T> cg;
    static CONJUGATE_RESIDUAL<T> cr;
    static SYMMQMR<T> symmqmr;
    KRYLOV_SOLVER<T>* solver=0;
    const char* solver_name=0;
    if(deformables_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_cg){solver=&cg;solver_name="CG";}
    else if(deformables_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_cr){solver=&cr;solver_name="CONJUGATE_RESIDUAL";}
    else if(deformables_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_symmqmr){solver=&symmqmr;solver_name="SYMMQMR";}
    solver->print_diagnostics=deformable_body_collection.print_diagnostics;
    solver->print_residuals=deformable_body_collection.print_residuals;
    solver->iterations_used=&deformable_body_collection.iterations_used_diagnostic;
    solver->restart_iterations=deformables_parameters.implicit_solve_parameters.cg_restart_iterations;
    system.project_nullspace_frequency=deformables_parameters.implicit_solve_parameters.project_nullspace_frequency;
    ARRAY<TWIST<TV> > twist;twist.Resize(rigid_geometry_particles.V.Size());
    for(int i=1;i<=twist.Size();i++) twist(i)=TWIST<TV>(rigid_geometry_particles.V(i),rigid_geometry_particles.angular_velocity(i));
    DEFORMABLES_VELOCITY<TV> V(particles.V,twist,deformable_body_collection,rigid_geometry_collection),F(F_full,rigid_F_full,deformable_body_collection,rigid_geometry_collection),
        R(R_full,rigid_R_full,deformable_body_collection,rigid_geometry_collection),S(S_full,rigid_S_full,deformable_body_collection,rigid_geometry_collection),
        B(B_full,rigid_B_full,deformable_body_collection,rigid_geometry_collection),AR(AR_full,rigid_AR_full,deformable_body_collection,rigid_geometry_collection);
    DEFORMABLES_MASS<TV> mass(deformable_body_collection); // TODO: Doing duplicate computation of mass.

    if(deformables_parameters.implicit_solve_parameters.spectral_analysis){V=B;LANCZOS_ITERATION<T>::Print_Spectral_Information(system,V,F,S,(T)1e-2*deformables_parameters.implicit_solve_parameters.cg_tolerance,deformables_parameters.implicit_solve_parameters.lanczos_iterations);}

    LOG::Time(solver_name);
    Diagnostics(dt,current_position_time,0,0,606,"Before solve");
    if(!solver->Solve(system,V,B,F,S,R,AR,deformables_parameters.implicit_solve_parameters.cg_tolerance,1,deformables_parameters.implicit_solve_parameters.cg_iterations) && deformables_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure)
        throw std::runtime_error("Backward Euler Failed");
    Diagnostics(dt,current_position_time,0,0,607,"After solve");
    LOG::Stop_Time();

    if(velocity_update && deformables_parameters.use_post_cg_constraints)
        system.Force(V,F);

    Finish_Backward_Euler_Step(dt,current_position_time,velocity_update);
    for(int i=1;i<=twist.Size();i++){rigid_geometry_particles.V(i)=twist(i).linear;rigid_geometry_particles.angular_velocity(i)=twist(i).angular;}
}
//#####################################################################
// Function Average_And_Exchange_Position
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Average_And_Exchange_Position()
{
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    assert(X_save.m==particles.array_collection->Size());
    assert(rigid_X_save.m==rigid_geometry_collection.particles.array_collection->Size() && rigid_rotation_save.m==rigid_geometry_collection.particles.array_collection->Size());
    const ARRAY<int>& simulated_particles=deformable_body_collection.simulated_particles;
    INDIRECT_ARRAY<ARRAY<TV> > X_save_fragment(X_save,simulated_particles);
    INDIRECT_ARRAY<ARRAY_VIEW<TV> > X_fragment(particles.X,simulated_particles);
    for(int i=1;i<=X_fragment.Size();i++){TV X_average=(T).5*(X_fragment(i)+X_save_fragment(i));X_save_fragment(i)=X_fragment(i);X_fragment(i)=X_average;}
    for(int i=1;i<=rigid_geometry_collection.kinematic_rigid_geometry.Size();i++){int p=rigid_geometry_collection.kinematic_rigid_geometry(i); //this needs to be the kinematic bodies
        TV tmp_X=TV::Interpolate(rigid_geometry_collection.particles.X(p),rigid_X_save(p),(T).5);
        rigid_X_save(p)=rigid_geometry_collection.particles.X(p);
        rigid_geometry_collection.particles.X(p)=tmp_X;
        ROTATION<TV> tmp_rotation=ROTATION<TV>::Spherical_Linear_Interpolation(rigid_geometry_collection.particles.rotation(p),rigid_rotation_save(p),(T).5);
        rigid_rotation_save(p)=rigid_geometry_collection.particles.rotation(p);
        rigid_geometry_collection.particles.rotation(p)=tmp_rotation;}
}
//#####################################################################
// Function Trapezoidal_Step_Velocity
//#####################################################################
// assumes X is fixed at time+dt/2
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Trapezoidal_Step_Velocity(const T dt,const T time)
{
    PARTICLES<TV>& particles=deformable_body_collection.particles;

    // save V at time
    Save_Velocity();
    // update V implicitly to time+dt/2
    Backward_Euler_Step_Velocity_Helper(dt/2,time,time+dt/2,true);
    // Use V_n instead of V_save rather than copying it around another time.  Also simplifies state dependencies.
    // extrapolate V to time+dt based on V at time and time+dt/2
    DEFORMABLES_VELOCITY<TV> V_n(V_save,deformable_body_collection),V(particles.V,deformable_body_collection);
    V*=(T)2;V-=V_n;

    // enforce boundary conditions again
    if(deformables_parameters.deformable_object_collision_parameters.perform_collision_body_collisions) deformable_body_collection.collisions.Activate_Collisions(false);
    Set_External_Velocities(particles.V,time+dt,time+dt/2);
    kinematic_evolution.Set_External_Velocities(rigid_geometry_collection.particles.V,rigid_geometry_collection.particles.angular_velocity,time+dt,time+dt/2);
    deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities();
}
//#####################################################################
// Function Backward_Euler_Step_Velocity
//#####################################################################
// assumes X is fixed at time+dt
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Backward_Euler_Step_Velocity(const T dt,const T time)
{
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    
    Backward_Euler_Step_Velocity_Helper(dt,time,time+dt,true);
    // enforce boundary conditions again
    if(deformables_parameters.deformable_object_collision_parameters.perform_collision_body_collisions) deformable_body_collection.collisions.Activate_Collisions(false);
    Set_External_Velocities(particles.V,time+dt,time+dt);
    kinematic_evolution.Set_External_Velocities(rigid_geometry_collection.particles.V,rigid_geometry_collection.particles.angular_velocity,time+dt,time+dt);
    
    // enforce boundary conditions again
    deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities();
}
//#####################################################################
// Function Advance_One_Time_Step_Position
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Advance_One_Time_Step_Position(const T dt,const T time)
{
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Advance_One_Time_Step_Position Start dt=%f time=%f",dt,time),2,2);
    MPI_SOLIDS<TV>* mpi_solids=deformable_body_collection.mpi_solids;

    if((deformables_parameters.triangle_collision_parameters.perform_self_collision && (deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.m || deformables_parameters.triangle_collision_parameters.initialize_collisions_without_objects))
       && (!repulsions && deformables_parameters.triangle_collision_parameters.perform_per_time_step_repulsions))
        repulsions=&deformable_body_collection.triangle_repulsions;

    Diagnostics(dt,time,0,0,1,"begin integration");

    deformables_evolution_callbacks->Update_Deformables_Parameters(time);
    // save position and velocity for later trapezoidal rule
    Save_Velocity();
    Save_Position(X_save,rigid_X_save,rigid_rotation_save);

    if(deformables_parameters.deformable_object_collision_parameters.perform_collision_body_collisions) deformable_body_collection.collisions.Activate_Collisions(false);

    if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data_Global(deformable_body_collection.particles.X);
    deformable_body_collection.Update_Position_Based_State(time+dt,true);

    Backward_Euler_Step_Velocity_Helper(dt/2,time,time,false); // update V implicitly to time+dt/2

    if(deformables_parameters.verbose) Print_Maximum_Velocities(time);
    Diagnostics(dt,time,1,0,6,"backward Euler");

    if(repulsions) repulsions->Adjust_Velocity_For_Self_Repulsion_Using_History(dt,true,false);
    Compute_Momentum_Differences();

    // add collision impulses to time n velocities and save
    Euler_Step_Position(dt,time);
    Diagnostics(dt,time,1,2,10,"Euler step position");
    Exchange_Velocity();
    Diagnostics(dt,time,0,2,11,"restore velocity");
    if(deformables_parameters.deformable_object_collision_parameters.perform_collision_body_collisions){
        int interactions=deformable_body_collection.collisions.Adjust_Nodes_For_Collision_Body_Collisions(deformable_body_collection.binding_list,
            deformable_body_collection.soft_bindings,X_save,dt,0);
        if(interactions) LOG::Stat("collision body collisions",interactions);}

    Restore_Position(X_save,rigid_X_save,rigid_rotation_save);
    Diagnostics(dt,time,0,0,13,"restore position");
    Save_Velocity();

    // update positions, apply contact, arb prestabilization and push out
    Update_Velocity_Using_Stored_Differences(dt/2,time);
    Diagnostics(dt,time,1,0,18,"update velocity using stored differences");

    Update_Positions_And_Apply_Contact_Forces(dt,time,true);
    Diagnostics(dt,time,1,2,20,"contact, prestabilization");

    // This is necessary since the velocity update below needs to use soft bound positions, but modifying positions within the velocity update is prohibited.
    if(repulsions) deformable_body_collection.soft_bindings.Clamp_Particles_To_Embedded_Positions(true);

    Restore_Velocity();
    Diagnostics(dt,time,0,2,22,"restore velocity");
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Advance_One_Time_Step_Position End dt=%f time=%f",dt,time),2,2);
}
//#####################################################################
// Function Advance_One_Time_Step_Velocity
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Advance_One_Time_Step_Velocity(const T dt,const T time)
{
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Advance_One_Time_Step_Velocity Start dt=%f time=%f",dt,time),2,2);

    MPI_SOLIDS<TV>* mpi_solids=deformable_body_collection.mpi_solids;

    if(deformables_parameters.deformable_object_collision_parameters.perform_collision_body_collisions) deformable_body_collection.collisions.Activate_Collisions(true);

    if(deformables_parameters.use_trapezoidal_rule_for_velocities){
        Average_And_Exchange_Position(); // move to positions at time+dt/2 for trapezoidal step
        Diagnostics(dt,time,0,1,25,"average and exchange positions");

        if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data_Global(deformable_body_collection.particles.X);
        deformable_body_collection.Update_Position_Based_State(time+dt/2,false);
        Trapezoidal_Step_Velocity(dt,time);
        Diagnostics(dt,time,2,1,29,"trazepoid rule");

        Restore_Position(X_save,rigid_X_save,rigid_rotation_save); // move to final positions at time time+dt
        Diagnostics(dt,time,2,2,30,"restore position");}
    else{
        if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data_Global(deformable_body_collection.particles.X);
        deformable_body_collection.Update_Position_Based_State(time+dt,false);
        Backward_Euler_Step_Velocity(dt,time); 
        Diagnostics(dt,time,2,2,29,"backward Euler");}

    if(deformables_parameters.use_post_cg_constraints) Apply_Constraints(dt,time);
    else if(repulsions) repulsions->Adjust_Velocity_For_Self_Repulsion_Using_History(dt,false,true);
    if(repulsions) Diagnostics(dt,time,2,2,40,"self repulsions");
    deformable_body_collection.collisions.Reset_Object_Collisions();
    if(deformables_parameters.verbose) Print_Maximum_Velocities(time);

    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Advance_One_Time_Step_Velocity End dt=%f time=%f",dt,time),2,2);
}
//#####################################################################
// Function Apply_Constraints
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Apply_Constraints(const T dt,const T time)
{
    Save_Position(X_save_for_constraints,rigid_X_save_for_constraints,rigid_rotation_save_for_constraints);
    if(deformables_parameters.deformable_object_collision_parameters.use_existing_contact) use_existing_contact=true;
    Update_Positions_And_Apply_Contact_Forces(dt,time,false);use_existing_contact=false;
    Diagnostics(dt,time,4,2,37,"contact, prestabilization");
    Restore_Position(X_save_for_constraints,rigid_X_save_for_constraints,rigid_rotation_save_for_constraints);
    Diagnostics(dt,time,2,2,38,"restore position");
    // modify velocity with inelastic and friction repulsions
    if(repulsions) repulsions->Adjust_Velocity_For_Self_Repulsion_Using_History(dt,false,true);
}
//#####################################################################
// Function Print_Maximum_Velocities
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Print_Maximum_Velocities(const T time) const
{
    std::stringstream ss;
    ss<<"time = "<<time<<std::endl;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    int max_index=0;T max_magnitude_squared=-FLT_MAX;const INDIRECT_ARRAY<ARRAY_VIEW<TV>,ARRAY<int>&>& V=particles.V.Subset(deformable_body_collection.dynamic_particles);
    for(int i=1;i<=V.Size();i++){T magnitude_squared=V(i).Magnitude_Squared();if(magnitude_squared>max_magnitude_squared){max_magnitude_squared=magnitude_squared;max_index=i;}}
    if(max_index){
        int p=deformable_body_collection.dynamic_particles(max_index);T max_magnitude=sqrt(max_magnitude_squared),max_magnitude_global=max_magnitude;
        if(deformable_body_collection.mpi_solids) max_magnitude_global=deformable_body_collection.mpi_solids->Reduce_Max_Global(max_magnitude_global);
        ss<<"maximum velocity = "<<max_magnitude_global;
        if(deformable_body_collection.mpi_solids) ss<<", local = "<<max_magnitude;
        ss<<" ("<<p<<")"<<std::endl;}
    LOG::filecout(ss.str());
}
//#####################################################################
// Function Diagnostics
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Diagnostics(const T dt,const T time,const int velocity_time,const int position_time,int step,const char* description)
{
    static const char* time_table[]={"n","(n+1/2)","(n+1)","(n+3/2)","(n+2)"};
    //deformable_body_collection.Print_Energy(time+position_time*(T).5*time,step);
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Finished step %i (%s).  State: x^%s  v^%s.   dt=%f time=%f",step,description,
        time_table[position_time],time_table[velocity_time],dt,time),2,3);
}
//#####################################################################
// Function Update_Positions_And_Apply_Contact_Forces
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Update_Positions_And_Apply_Contact_Forces(const T dt,const T time,const bool use_saved_pairs)
{
    Euler_Step_Position(dt,time);
    for(int i=1;i<=rigid_geometry_collection.kinematic_rigid_geometry.m;i++){
        RIGID_GEOMETRY<TV>& rigid_geometry=rigid_geometry_collection.Rigid_Geometry(rigid_geometry_collection.kinematic_rigid_geometry(i));
        rigid_geometry.Update_Bounding_Box();}

    if(deformables_parameters.deformable_object_collision_parameters.perform_collision_body_collisions){
        int interactions=0;
        if(use_existing_contact)
            interactions+=deformable_body_collection.collisions.Adjust_Existing_Nodes_For_Collision_Body_Collisions(deformable_body_collection.binding_list,deformable_body_collection.soft_bindings,X_save,dt,0);
        else 
            interactions+=deformable_body_collection.collisions.Adjust_Nodes_For_Collision_Body_Collisions(deformable_body_collection.binding_list,deformable_body_collection.soft_bindings,X_save,dt,0);
        if(interactions) LOG::Stat("collision body collisions",interactions);}
}
//#####################################################################
// Function Update_Velocity_Using_Stored_Differences
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Update_Velocity_Using_Stored_Differences(const T dt,const T time)
{
    for(int i=1;i<=deformable_body_collection.simulated_particles.m;i++){int p=deformable_body_collection.simulated_particles(i);
        deformable_body_collection.particles.V(p)+=V_difference(p);}
    for(int i=1;i<=rigid_geometry_collection.kinematic_rigid_geometry.m;i++){int p=rigid_geometry_collection.kinematic_rigid_geometry(i);
        kinematic_evolution.Set_External_Velocities(rigid_geometry_collection.particles.V(p),rigid_geometry_collection.particles.angular_velocity(p),time+dt,p);}
}
 //#####################################################################
// Function Update_Velocity_Using_Stored_Differences
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Update_Velocity_Using_Stored_Differences(const T dt,const T time,const int p)
{
    PHYSBAM_NOT_IMPLEMENTED(); 
}
//#####################################################################
// Function Compute_Momentum_Differences
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Compute_Momentum_Differences()
{
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    V_difference.Resize(particles.array_collection->Size());
    for(int i=1;i<=deformable_body_collection.simulated_particles.m;i++){
        int p=deformable_body_collection.simulated_particles(i);V_difference(p)=particles.V(p)-V_save(p);}
}
//#####################################################################
// Function Save_Velocity
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Save_Velocity()
{
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    V_save.Resize(particles.array_collection->Size(),false,false);
    rigid_velocity_save.Resize(rigid_geometry_collection.particles.array_collection->Size(),false,false);
    V_save.Subset(deformable_body_collection.simulated_particles)=particles.V.Subset(deformable_body_collection.simulated_particles);
    for(int i=1;i<=rigid_geometry_collection.particles.array_collection->Size();i++) if(rigid_geometry_collection.Is_Active(i)){
        rigid_velocity_save(i)=rigid_geometry_collection.particles.V(i);}
}
//#####################################################################
// Function Restore_Velocity
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Restore_Velocity() const
{

    PARTICLES<TV>& particles=deformable_body_collection.particles;
    particles.V.Subset(deformable_body_collection.simulated_particles)=V_save.Subset(deformable_body_collection.simulated_particles);
    for(int i=1;i<=rigid_geometry_collection.particles.array_collection->Size();i++) if(rigid_geometry_collection.Is_Active(i)){
        rigid_geometry_collection.particles.V(i)=rigid_velocity_save(i);}
}
//#####################################################################
// Function Exchange_Velocity
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Exchange_Velocity()
{

    PARTICLES<TV>& particles=deformable_body_collection.particles;RIGID_GEOMETRY_PARTICLES<TV>& rigid_geometry_particles=rigid_geometry_collection.particles;
    V_save.Resize(particles.array_collection->Size(),false,false);
    rigid_velocity_save.Resize(rigid_geometry_particles.array_collection->Size(),false,false);
    for(int i=1;i<=deformable_body_collection.simulated_particles.m;i++){int p=deformable_body_collection.simulated_particles(i);
        exchange(V_save(p),deformable_body_collection.particles.V(p));}
    for(int i=1;i<=rigid_geometry_particles.array_collection->Size();i++) if(rigid_geometry_collection.Is_Active(i)){
        exchange(rigid_velocity_save(i),rigid_geometry_particles.V(i));}
}
//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Initialize_Rigid_Bodies(const T frame_rate, const bool restart)
{
    // initialize kinematic object positions and velocities
    if(!restart){
        kinematic_evolution.Get_Current_Kinematic_Keyframes(1/frame_rate,time);
        kinematic_evolution.Set_External_Positions(rigid_geometry_collection.particles.X,rigid_geometry_collection.particles.rotation,time);
        kinematic_evolution.Set_External_Velocities(rigid_geometry_collection.particles.V,rigid_geometry_collection.particles.angular_velocity,time,time);
        for(int i(1);i<=rigid_geometry_collection.particles.array_collection->Size();i++) if(rigid_geometry_collection.Is_Active(i)){rigid_geometry_collection.particles.rotation(i).Normalize();}}
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Set_External_Positions(ARRAY_VIEW<TV> X,const T time)
{
    deformable_body_collection.deformables_example_forces_and_velocities->Set_External_Positions(X,time);
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time)
{
    if(deformable_body_collection.collisions.collisions_on) deformable_body_collection.collisions.Set_Collision_Velocities(V);
    deformable_body_collection.deformables_example_forces_and_velocities->Set_External_Velocities(V,velocity_time,current_position_time);
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time)    
{
    if(deformable_body_collection.collisions.collisions_on) deformable_body_collection.collisions.Zero_Out_Collision_Velocities(V);
    deformable_body_collection.deformables_example_forces_and_velocities->Zero_Out_Enslaved_Velocity_Nodes(V,velocity_time,current_position_time);
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time)    
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Euler_Step_Position
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Euler_Step_Position(const T dt,const T time,const int p)
{
    RIGID_GEOMETRY<TV>& rigid_geometry=rigid_geometry_collection.Rigid_Geometry(p);
    if(rigid_geometry.is_static) return;
    kinematic_evolution.Set_External_Positions(rigid_geometry.X(),rigid_geometry.Rotation(),time+dt,p);
}
//#####################################################################
// Function Euler_Step_Position
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Euler_Step_Position(const T dt,const T time)
{
    MPI_SOLIDS<TV>* mpi_solids=deformable_body_collection.mpi_solids;
    deformable_body_collection.particles.Euler_Step_Position(deformable_body_collection.dynamic_particles,dt);
    Set_External_Positions(deformable_body_collection.particles.X,time+dt);
    kinematic_evolution.Set_External_Positions(rigid_geometry_collection.particles.X,rigid_geometry_collection.particles.rotation,time+dt);
    deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Positions();
    if(mpi_solids){
         mpi_solids->Exchange_Force_Boundary_Data_Global(deformable_body_collection.particles.X);
         mpi_solids->Exchange_Binding_Boundary_Data_Global(deformable_body_collection.particles.X);
/*      TODO: exchange rigid body particles data
        mpi_solids->Exchange_Boundary_Data_Global(solid_body_collection.rigid_body_collection.rigid_body_particle.frame);
        TODO: update angular velocity for received boundary data*/
}
}
//#####################################################################
// Function Initialize_World_Space_Masses
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Initialize_World_Space_Masses()
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Clamp_Velocities
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Clamp_Velocities()
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function CFL
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
CFL(const bool verbose)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Save_Position
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Save_Position(ARRAY<TV>& X,ARRAY<TV>& rigid_X,ARRAY<ROTATION<TV> >& rigid_rotation)
{

    PARTICLES<TV>& particles=deformable_body_collection.particles;
    const ARRAY<int>& simulated_particles=deformable_body_collection.simulated_particles;
    X.Resize(particles.array_collection->Size(),false,false);
    X.Subset(simulated_particles)=particles.X.Subset(simulated_particles);
    rigid_X.Resize(rigid_geometry_collection.particles.array_collection->Size(),false,false);
    rigid_rotation.Resize(rigid_geometry_collection.particles.array_collection->Size(),false,false);
    for(int i=1;i<=rigid_geometry_collection.particles.array_collection->Size();i++) if(rigid_geometry_collection.Is_Active(i)){
        rigid_X(i)=rigid_geometry_collection.particles.X(i);rigid_rotation(i)=rigid_geometry_collection.particles.rotation(i);}
}
//#####################################################################
// Function Restore_Position
//#####################################################################
template<class TV> void DEFORMABLES_EVOLUTION<TV>::
Restore_Position(ARRAY_VIEW<const TV> X,ARRAY_VIEW<const TV> rigid_X,ARRAY_VIEW<const ROTATION<TV> > rigid_rotation)
{ 
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    const ARRAY<int>& simulated_particles=deformable_body_collection.simulated_particles;
    PHYSBAM_ASSERT(X.Size()==particles.array_collection->Size());PHYSBAM_ASSERT(rigid_X.Size()==rigid_geometry_collection.particles.array_collection->Size());
    PHYSBAM_ASSERT(rigid_rotation.Size()==rigid_geometry_collection.particles.array_collection->Size());
    particles.X.Subset(simulated_particles)=X.Subset(simulated_particles);
    for(int i=1;i<=rigid_geometry_collection.particles.array_collection->Size();i++) if(rigid_geometry_collection.Is_Active(i)){
        rigid_geometry_collection.particles.X(i)=rigid_X(i);rigid_geometry_collection.particles.rotation(i)=rigid_rotation(i);}
}
//#####################################################################
// Function Adjust_Velocity_For_Self_Repulsion_And_Self_Collisions
//#####################################################################
template<class TV> bool DEFORMABLES_EVOLUTION<TV>::
Adjust_Velocity_For_Self_Repulsion_And_Self_Collisions(const T dt,const T time,int& repulsions_found,int& collisions_found,const bool exit_early)
{
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry=deformable_body_collection.triangle_repulsions_and_collisions_geometry;
    ARRAY<bool>& modified=geometry.modified_full;ARRAY<TV>& X_self_collision_free=geometry.X_self_collision_free;

    repulsions_found=0;collisions_found=0; // important to initialize these in case repulsions/collisions are turned off
    if(!deformables_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions && deformables_parameters.triangle_collision_parameters.turn_off_all_collisions){
        std::stringstream ss;ss<<"all repulsions and collisions are turned off - nothing to do"<<std::endl;LOG::filecout(ss.str());return false;}

    // update soft bound particle positions and velocities
    deformable_body_collection.soft_bindings.Clamp_Particles_To_Embedded_Positions(true);
    deformable_body_collection.soft_bindings.Clamp_Particles_To_Embedded_Velocities(true);

    if(deformable_body_collection.mpi_solids){
        LOG::SCOPE scope("Collisions All Gather","Collisions All Gather");
        deformable_body_collection.mpi_solids->All_Gather_Particles(particles.X,particles.V);}

    modified=CONSTANT_ARRAY<bool>(particles.array_collection->Size(),false);

    // compute average velocity, but save final velocities in case of no repulsion or collisions
    ARRAY<TV> V_save(particles.V);

    ARRAY<TV> V_averaged(1/dt*(particles.X-X_self_collision_free));
    particles.V=V_averaged;

    // check for triangle pairs that are already intersecting, and label them to be ignored
    if(geometry.allow_intersections) geometry.Compute_Intersecting_Segment_Face_Pairs();

    // self repulsion - adjust time n+1/2 velocity and time n+1 position
    if(deformables_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions){
        if(deformable_body_collection.mpi_solids) PHYSBAM_NOT_IMPLEMENTED("per collision step repulsions under mpi");
        repulsions_found=deformable_body_collection.triangle_repulsions.Adjust_Velocity_For_Self_Repulsion(dt,false);}

    // self collisions - adjust time n+1/2 velocity and time n+1 position
    if(!deformables_parameters.triangle_collision_parameters.turn_off_all_collisions)
        collisions_found=deformable_body_collection.triangle_collisions.Adjust_Velocity_For_Self_Collisions(dt,time,exit_early);

    // update velocities
    if(collisions_found<0) // failed to resolve collisions (should only happen when exit_early=true) - otherwise any collisions that were found have been resolved
        PHYSBAM_NOT_IMPLEMENTED("exit_early=true");
    else if(collisions_found>0){ // restore unmodified velocities
        if(deformable_body_collection.triangle_repulsions_and_collisions_geometry.mass_modifier) particles.V=V_save;
        else for(int p=1;p<=particles.array_collection->Size();p++) if(!modified(p)) particles.V(p)=V_save(p);}
    else if(repulsions_found){ // repulsions only, restore velocities for unmodified and apply velocity delta otherwise
        for(int p=1;p<=particles.array_collection->Size();p++) particles.V(p)=modified(p)?V_save(p)+particles.V(p)-V_averaged(p):V_save(p);}
    else{particles.V=V_save;return false;} // restore all the unmodified velocities

    // propagate any changes from soft bound particles to parents
    deformable_body_collection.Adjust_Mesh_For_Self_Collision();

    return true;
}
//#####################################################################
template class DEFORMABLES_EVOLUTION<VECTOR<float,1> >;
template class DEFORMABLES_EVOLUTION<VECTOR<float,2> >;
template class DEFORMABLES_EVOLUTION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DEFORMABLES_EVOLUTION<VECTOR<double,1> >;
template class DEFORMABLES_EVOLUTION<VECTOR<double,2> >;
template class DEFORMABLES_EVOLUTION<VECTOR<double,3> >;
#endif

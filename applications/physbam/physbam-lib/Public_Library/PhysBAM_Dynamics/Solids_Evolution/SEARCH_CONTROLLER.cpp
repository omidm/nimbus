//#####################################################################
// Copyright 2008-2009, Jon Gretarsson, Michael Lentine, Craig Schroeder, Bridget Vuong.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEARCH_CONTROLLER
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/DIRECTED_GRAPH.h>
#include <PhysBAM_Tools/Data_Structures/STACK.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Tools/Math_Tools/sign.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_ELEMENT_ID.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_ROTATION.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_ND.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions_And_Grids/GRID_BASED_COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER_HASH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/ETHER_DRAG.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/WIND_DRAG.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_SYSTEM.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_Evolution/SEARCH_CONTROLLER.h>
#include <cfloat>
using ::std::sqrt;
using namespace PhysBAM;
//#####################################################################
// Function Advance_One_Time_Step_Position
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Advance_One_Time_Step_Position(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    SOLIDS_PARAMETERS<TV>& solids_parameters=driver->example.solids_parameters;
    SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters=driver->example.solids_fluids_parameters;
    ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body=solid_body_collection.rigid_body_collection.articulated_rigid_body; // Needn't be a pointer
    if(NEWMARK_EVOLUTION<TV>* newmark_evolution=dynamic_cast<NEWMARK_EVOLUTION<TV>*>(driver->example.solids_evolution)){
        const bool advance_rigid_bodies=true; //solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m!=0;  TODO: Fix this.
        const bool solids=!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Solid_Node();

        if(articulated_rigid_body.Has_Actuators()) Set_PD_Targets();
        newmark_evolution->Diagnostics(dt,time,0,0,1,"begin integration");

        newmark_evolution->solids_evolution_callbacks->Update_Solids_Parameters(time);
        if(solids){
            newmark_evolution->rigid_body_collisions->Initialize_Data_Structures();
            if(solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions) solid_body_collection.deformable_body_collection.collisions.Activate_Collisions(false);}

        solid_body_collection.Update_Position_Based_State(time+dt,true);
        if(incorporate_fluids) fluids_parameters->incompressible->Advance_One_Time_Step_Forces(face_velocities,dt,time,fluids_parameters->implicit_viscosity,0,fluids_parameters->number_of_ghost_cells);

        newmark_evolution->Backward_Euler_Step_Velocity_Helper(dt,time,time,false); // update V implicitly to time+dt/2
        if(solids_parameters.verbose) newmark_evolution->Print_Maximum_Velocities(time);
        newmark_evolution->Diagnostics(dt,time,1,0,6,"backward Euler");
        if(!solids) return; // early exit for fluids only in parallel

        if(advance_rigid_bodies){
            if(articulated_rigid_body.Has_Actuators()) articulated_rigid_body.Compute_Position_Based_State(dt,time);
            articulated_rigid_body.Solve_Velocities_for_PD(time,dt,solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);
            newmark_evolution->Diagnostics(dt,time,1,0,19,"solve velocities for pd");}}
    else PHYSBAM_FATAL_ERROR("CANNOT FIND NEWMARK!!!");
}
//####################################################################################
// Function Create_Joint_Function
//####################################################################################
template<class T_GRID> JOINT_FUNCTION<typename T_GRID::VECTOR_T>* SEARCH_CONTROLLER<T_GRID>::
Create_Joint_Function(const JOINT_ID joint_id)
{
    joint_mesh(joint_id)->Set_Joint_Function(new JOINT_FUNCTION<TV>(*joint_mesh(joint_id),*Parent(joint_id),*Child(joint_id)));
    return joint_mesh(joint_id)->joint_function;
}
//#####################################################################
// Function Save_Position
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Save_Position(ARRAY<TV>& X,ARRAY<TV>& rigid_X,ARRAY<ROTATION<TV> >& rigid_rotation)
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    const ARRAY<int>& simulated_particles=solid_body_collection.deformable_body_collection.simulated_particles;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particle;
    X.Resize(particles.array_collection->Size(),false,false);
    X.Subset(simulated_particles)=particles.X.Subset(simulated_particles);
    rigid_X.Resize(rigid_body_particles.array_collection->Size(),false,false);
    rigid_rotation.Resize(rigid_body_particles.array_collection->Size(),false,false);
    for(int i=1;i<=rigid_body_particles.array_collection->Size();i++) if(solid_body_collection.rigid_body_collection.Is_Active(i)){rigid_X(i)=rigid_body_particles.X(i);rigid_rotation(i)=rigid_body_particles.rotation(i);}
}
//#####################################################################
// Function Restore_Position
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Restore_Position(ARRAY_VIEW<const TV> X,ARRAY_VIEW<const TV> rigid_X,ARRAY_VIEW<const ROTATION<TV> > rigid_rotation)
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    const ARRAY<int>& simulated_particles=solid_body_collection.deformable_body_collection.simulated_particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    const ARRAY<int>& simulated_rigid_body_particles=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles;
    PHYSBAM_ASSERT(X.Size()==particles.array_collection->Size());PHYSBAM_ASSERT(rigid_X.Size()==rigid_body_collection.rigid_body_particle.array_collection->Size());
    PHYSBAM_ASSERT(rigid_rotation.Size()==rigid_body_collection.rigid_body_particle.array_collection->Size());
    particles.X.Subset(simulated_particles)=X.Subset(simulated_particles);
    for(int i=1;i<=simulated_rigid_body_particles.m;i++) rigid_body_collection.Rigid_Body(simulated_rigid_body_particles(i)).Update_Angular_Velocity();
    for(int i=1;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++) if(rigid_body_collection.Is_Active(i)){
        RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(i);rigid_body_collection.rigid_body_particle.X(i)=rigid_X(i);
        rigid_body_collection.rigid_body_particle.rotation(i)=rigid_rotation(i);body.Update_Angular_Velocity();}
}
//#####################################################################
// Function Save_Velocity
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Save_Velocity(ARRAY<TV>& V,ARRAY<TV>& rigid_velocity,ARRAY<T_SPIN>& rigid_angular_momentum)
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    V.Resize(particles.array_collection->Size(),false,false);
    rigid_velocity.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);rigid_angular_momentum.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
    V.Subset(solid_body_collection.deformable_body_collection.simulated_particles)=particles.V.Subset(solid_body_collection.deformable_body_collection.simulated_particles);
    for(int i=1;i<=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles(i);
        rigid_velocity(p)=rigid_body_collection.rigid_body_particle.V(p);rigid_angular_momentum(p)=rigid_body_collection.rigid_body_particle.angular_momentum(p);}
    for(int i=1;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++) if(rigid_body_collection.Is_Active(i)){RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(i);
        if(!body.Is_Simulated()){rigid_velocity(i)=rigid_body_collection.rigid_body_particle.V(i);rigid_angular_momentum(i)=rigid_body_collection.rigid_body_particle.angular_momentum(i);}}
}
//#####################################################################
// Function Restore_Velocity
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Restore_Velocity(ARRAY<TV>& V,ARRAY<TV>& rigid_velocity,ARRAY<T_SPIN>& rigid_angular_momentum)
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    particles.V.Subset(solid_body_collection.deformable_body_collection.simulated_particles)=V.Subset(solid_body_collection.deformable_body_collection.simulated_particles);
    for(int i=1;i<=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles(i);
        rigid_body_collection.rigid_body_particle.V(p)=rigid_velocity(p);
        rigid_body_collection.rigid_body_particle.angular_momentum(p)=rigid_angular_momentum(p);rigid_body_collection.Rigid_Body(p).Update_Angular_Velocity();}
    for(int i=1;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++) if(rigid_body_collection.Is_Active(i)){RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(i);
        if(!body.Is_Simulated()){
            rigid_body_collection.rigid_body_particle.V(i)=rigid_velocity(i);rigid_body_collection.rigid_body_particle.angular_momentum(i)=rigid_angular_momentum(i);body.Update_Angular_Velocity();}}
}
//#####################################################################
// Function Save_State
//#####################################################################
template<class T_GRID> static typename T_GRID::SCALAR Get_Time_For_Saved_State_Helper(const GRID_BASED_COLLISION_GEOMETRY<T_GRID>& collision_body_list,const int state){
    typename T_GRID::SCALAR old_state_time=0;
    typedef typename T_GRID::VECTOR_T TV;
    if(RIGID_COLLISION_GEOMETRY<TV>* rigid_body_collision_geometry=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(collision_body_list.collision_geometry_collection.bodies(COLLISION_GEOMETRY_ID(1)))){
        old_state_time=rigid_body_collision_geometry->saved_states(state).y;}
    else if(DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>* deformable_body=dynamic_cast<DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>*>(collision_body_list.collision_geometry_collection.bodies(COLLISION_GEOMETRY_ID(1)))){
        old_state_time=deformable_body->saved_states(state).y;}
    else PHYSBAM_FATAL_ERROR();
    return old_state_time;
}
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Save_State(T_FACE_ARRAYS_SCALAR& face_velocities)
{
    //TODO(mlentine): Do we need to save joint frames?
    PHYSBAM_DEBUG_WRITE_SUBSTEP("At Save_State (search_controller)",0,0);
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings;
    rigid_bindings.Clamp_Particles_To_Embedded_Positions();
    rigid_bindings.Clamp_Particles_To_Embedded_Velocities();
    ARRAY<int> active_clusters;
    rigid_bindings.Deactivate_And_Return_Clusters(active_clusters,0);
    Save_Position(X_save,rigid_X_save,rigid_rotation_save);
    Save_Velocity(V_save,rigid_velocity_save,rigid_angular_momentum_save);
    time_save=driver->Time();
    rigid_bindings.Reactivate_Bindings(active_clusters);            
    if(incorporate_fluids && Fluid_Node()) {
        fluid_affects_solid_save=fluids_parameters->fluid_affects_solid;
        project_at_frame_boundaries_save=driver->project_at_frame_boundaries;
        if(fluids_parameters->collision_bodies_affecting_fluid->collision_geometry_collection.bodies.m){
            GRID_BASED_COLLISION_GEOMETRY<T_GRID>& collision_body_list=*fluids_parameters->collision_bodies_affecting_fluid;
            collision_body_list.Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_SAVED_OLD_STATE,time_save);}
        fluid_pressure_save.Resize(fluids_parameters->grid->Domain_Indices(1));fluid_velocity_save.Resize(*fluids_parameters->grid);
        T_ARRAYS_SCALAR::Copy(fluids_parameters->incompressible->projection.p,fluid_pressure_save);
        T_FACE_ARRAYS_SCALAR::Copy(face_velocities,fluid_velocity_save);}
}
//#####################################################################
// Function Restore_State
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Restore_State(T_FACE_ARRAYS_SCALAR& face_velocities)
{
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Before Restore_State (search_controller)",0,0);
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings;
    ARRAY<int> active_clusters;
    rigid_bindings.Deactivate_And_Return_Clusters(active_clusters,0);
    Restore_Position(X_save,rigid_X_save,rigid_rotation_save);
    Restore_Velocity(V_save,rigid_velocity_save,rigid_angular_momentum_save);
    driver->Set_Time(time_save);
    rigid_bindings.Reactivate_Bindings(active_clusters);            
    if(incorporate_fluids && Fluid_Node()){
        fluids_parameters->fluid_affects_solid=fluid_affects_solid_save;
        driver->project_at_frame_boundaries=project_at_frame_boundaries_save;
        if(fluids_parameters->collision_bodies_affecting_fluid->collision_geometry_collection.bodies.m){
            GRID_BASED_COLLISION_GEOMETRY<T_GRID>& collision_body_list=*fluids_parameters->collision_bodies_affecting_fluid;
            collision_body_list.Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_SAVED_OLD_STATE);
            T old_state_time=Get_Time_For_Saved_State_Helper(collision_body_list,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_SAVED_OLD_STATE);
            collision_body_list.Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,old_state_time);}
        T_ARRAYS_SCALAR::Copy(fluid_pressure_save,fluids_parameters->incompressible->projection.p);
        T_FACE_ARRAYS_SCALAR::Copy(fluid_velocity_save,face_velocities);}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After Restore_State (search_controller)",0,0);
}
//#####################################################################
// Function Save_Nested_State
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Save_Nested_State(T_FACE_ARRAYS_SCALAR& face_velocities)
{
    Save_Velocity(V_pd,rigid_velocity_pd,rigid_angular_momentum_pd);
    if(incorporate_fluids && Fluid_Node()){
        fluid_pressure_pd_save.Resize(fluids_parameters->grid->Domain_Indices(1));
        fluid_velocity_pd_save.Resize(*fluids_parameters->grid);
        T_ARRAYS_SCALAR::Copy(fluids_parameters->incompressible->projection.p,fluid_pressure_pd_save);
        T_FACE_ARRAYS_SCALAR::Copy(face_velocities,fluid_velocity_pd_save);}
}   
//#####################################################################
// Function Restore_Nested_State
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Restore_Nested_State(T_FACE_ARRAYS_SCALAR& face_velocities)
{
    Restore_Velocity(V_pd,rigid_velocity_pd,rigid_angular_momentum_pd);
    if(incorporate_fluids && Fluid_Node()){
        T_FACE_ARRAYS_SCALAR::Copy(fluid_velocity_pd_save,face_velocities);
        T_ARRAYS_SCALAR::Copy(fluid_pressure_pd_save,fluids_parameters->incompressible->projection.p);}
}
//#####################################################################
// Function Save_PD_State
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Save_PD_State()
{
    k_p_save.Resize(joint_mesh.Size(),false,false);
    pd_state_save.Resize(joint_mesh.Size(),false,false);
    joint_function_active_save.Resize(joint_mesh.Size(),false,false);
    constrain_pd_save=solid_body_collection.rigid_body_collection.articulated_rigid_body.constrain_pd_directions;
    for(JOINT_ID i(1);i<=joint_mesh.Size();i++){if(!joint_mesh.Is_Active(i)) continue;
        if(joint_mesh(i)->joint_function){
            k_p_save(i)=joint_mesh(i)->joint_function->k_p;
            joint_function_active_save(i)=joint_mesh(i)->joint_function->active;}
        pd_state_save(i)=joint_mesh(i)->global_post_stabilization;}
}
//#####################################################################
// Function Restore_PD_State
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Restore_PD_State()
{
    solid_body_collection.rigid_body_collection.articulated_rigid_body.constrain_pd_directions=constrain_pd_save;
    for(JOINT_ID i(1);i<=joint_mesh.Size();i++){if(!joint_mesh.Is_Active(i)) continue;
        if(joint_mesh(i)->joint_function){
            joint_mesh(i)->joint_function->Set_k_p(k_p_save(i));
            joint_mesh(i)->joint_function->active=joint_function_active_save(i);}
        joint_mesh(i)->global_post_stabilization=joint_mesh(i)->joint_function && joint_mesh(i)->joint_function->active;}
}
//#####################################################################
// Function Propogate_Solid_Helper
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Propogate_Solid_Helper(ARRAY<int>& cluster_bodies,TV& cluster_translation,TWIST<TV>& cluster_twist)
{
    RIGID_BODY_PARTICLES<TV>& rbp=solid_body_collection.rigid_body_collection.rigid_body_particle;
    for(int i=1;i<=cluster_bodies.m;i++){
        int child=cluster_bodies(i);
        rbp.V(child)=cluster_twist.linear+TV::Cross_Product(cluster_twist.angular,rbp.X(child)-cluster_translation);
        //TODO(mlentine): change this to momentum not velocity
        rbp.angular_velocity(child)=cluster_twist.angular;
        solid_body_collection.rigid_body_collection.Rigid_Body(child).Update_Angular_Momentum();}
}
//#####################################################################
// Function Make_Cluster_List
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Make_Cluster_List(int parent,ARRAY<int>& child_bodies)
{
    child_bodies.Append(parent);
    for(int i=1;i<=bone_hierarchy(parent).m;i++) Make_Cluster_List(bone_hierarchy(parent)(i),child_bodies);
}
//#####################################################################
// Function Project_Solid_Velocities
//#####################################################################
template<class T> void Diagonalize_Inertia_Tensor(typename RIGID_BODY_POLICY<VECTOR<T,1> >::WORLD_SPACE_INERTIA_TENSOR& inertia_tensor,typename RIGID_BODY_POLICY<VECTOR<T,1> >::INERTIA_TENSOR& cluster_inertia_tensor,ROTATION<VECTOR<T,1> >& cluster_rotation)
{PHYSBAM_NOT_IMPLEMENTED();}
template<class T> void Diagonalize_Inertia_Tensor(typename RIGID_BODY_POLICY<VECTOR<T,2> >::WORLD_SPACE_INERTIA_TENSOR& inertia_tensor,typename RIGID_BODY_POLICY<VECTOR<T,2> >::INERTIA_TENSOR& cluster_inertia_tensor,ROTATION<VECTOR<T,2> >& cluster_rotation)
{cluster_inertia_tensor=inertia_tensor;}
template<class T> void Diagonalize_Inertia_Tensor(typename RIGID_BODY_POLICY<VECTOR<T,3> >::WORLD_SPACE_INERTIA_TENSOR& inertia_tensor,typename RIGID_BODY_POLICY<VECTOR<T,3> >::INERTIA_TENSOR& cluster_inertia_tensor,ROTATION<VECTOR<T,3> >& cluster_rotation)
{
    MATRIX<T,3> rotation_matrix;inertia_tensor.Solve_Eigenproblem(cluster_inertia_tensor,rotation_matrix);
    cluster_rotation=ROTATION<VECTOR<T,3> >(rotation_matrix);
}
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Project_Solid_Velocities(const T time)
{
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings;
    RIGID_BODY_PARTICLES<TV>& rbp=solid_body_collection.rigid_body_collection.rigid_body_particle;
    TV cluster_translation;
    ROTATION<TV> cluster_rotation;
    TWIST<TV> cluster_twist;
    T cluster_mass=0;
    T_INERTIA_TENSOR cluster_inertia_tensor;
    T_SPIN cluster_momentum;
    T_WORLD_SPACE_INERTIA_TENSOR inertia_tensor;
    ARRAY<int> cluster_list;
    bool has_kinematic=false;
    for(JOINT_ID i(1);i<=joint_mesh.Size();i++){if(!joint_mesh.Is_Active(i)) continue;
        JOINT<TV>& joint=*joint_mesh(i);
        RIGID_BODY<TV> *parent_body=Parent(joint.id_number),*child_body=Child(joint.id_number);
        if((parent_body->Has_Infinite_Inertia() || child_body->Has_Infinite_Inertia()) && !(parent_body->Has_Infinite_Inertia() && child_body->Has_Infinite_Inertia())){
            if(parent_body->Has_Infinite_Inertia()){
                cluster_rotation=parent_body->Rotation();cluster_translation=parent_body->X();cluster_twist=parent_body->Twist();
                cluster_mass=parent_body->Mass();cluster_inertia_tensor=parent_body->Inertia_Tensor();cluster_momentum=parent_body->Angular_Momentum();}
            else{
                cluster_rotation=child_body->Rotation();cluster_translation=child_body->X();cluster_twist=child_body->Twist();
                cluster_mass=child_body->Mass();cluster_inertia_tensor=child_body->Inertia_Tensor();cluster_momentum=child_body->Angular_Momentum();}
            cluster_list.Remove_All();
            if(use_clusters && joint_clusters(joint.id_number).x){
                T_CLUSTER& cluster=*rigid_bindings.reverse_bindings.Get(joint_clusters(joint.id_number).x);
                for(RIGID_CLUSTER_CONSTITUENT_ID j(1);j<=cluster.children.Size();j++) cluster_list.Append(cluster.children(j));}
            else{
                if(rigid_bindings.Is_Parent(parent_body->particle_index)){
                    T_CLUSTER& cluster=*rigid_bindings.reverse_bindings.Get(rigid_bindings.Get_Parent_Index(parent_body->particle_index));
                    for(RIGID_CLUSTER_CONSTITUENT_ID j(1);j<=cluster.children.Size();j++) cluster_list.Append(cluster.children(j));}
                else cluster_list.Append(parent_body->particle_index);}
            if(use_clusters && joint_clusters(joint.id_number).y){
                T_CLUSTER& cluster=*rigid_bindings.reverse_bindings.Get(joint_clusters(joint.id_number).y);
                for(RIGID_CLUSTER_CONSTITUENT_ID j(1);j<=cluster.children.Size();j++) cluster_list.Append(cluster.children(j));}
            else{
                if(rigid_bindings.Is_Parent(child_body->particle_index)){
                    T_CLUSTER& cluster=*rigid_bindings.reverse_bindings.Get(rigid_bindings.Get_Parent_Index(child_body->particle_index));
                    for(RIGID_CLUSTER_CONSTITUENT_ID j(1);j<=cluster.children.Size();j++) cluster_list.Append(cluster.children(j));}
                else cluster_list.Append(child_body->particle_index);}
            ARRAY<int> active_clusters;
            rigid_bindings.Deactivate_And_Return_Clusters(active_clusters,(incorporate_fluids?fluids_parameters->collision_bodies_affecting_fluid:0));
            Propogate_Solid_Helper(cluster_list,cluster_translation,cluster_twist);
            rigid_bindings.Reactivate_Bindings(active_clusters);            
            has_kinematic=true;}}
    if(has_kinematic) return;
    //TODO (mlentine): This needs to be fixed for no kinematic bodies
    Make_Cluster_List(root_particle_index,cluster_list);
    for(int i=1;i<=cluster_list.m;i++){
        int child=cluster_list(i);
        cluster_mass+=rbp.mass(child);
        cluster_translation+=rbp.X(child)*rbp.mass(child);
        cluster_twist.linear+=rbp.V(child)*rbp.mass(child);}
    cluster_translation/=cluster_mass;
    cluster_twist.linear/=cluster_mass;
    for(int i=1;i<=cluster_list.m;i++){
        int child=cluster_list(i);
        RIGID_BODY<TV>& child_body=solid_body_collection.rigid_body_collection.Rigid_Body(child);
        TV s=(child_body.X()-cluster_translation);
        inertia_tensor+=child_body.World_Space_Inertia_Tensor()+MATRIX_POLICY<TV>::CROSS_PRODUCT_MATRIX::Cross_Product_Matrix(child_body.Mass()*s).Times_Cross_Product_Matrix_Transpose_With_Symmetric_Result(s);
        cluster_momentum+=TV::Cross_Product(child_body.X()-cluster_translation,child_body.Mass()*child_body.Twist().linear)+child_body.Angular_Momentum();}
    Diagonalize_Inertia_Tensor(inertia_tensor,cluster_inertia_tensor,cluster_rotation);
    Propogate_Solid_Helper(cluster_list,cluster_translation,cluster_twist);
}
//#####################################################################
// Function Project_Velocities
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Project_Velocities(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Before Project_Velocities (search_controller)",0,1);
    Project_Solid_Velocities(time);
    if(SOLID_FLUID_COUPLED_EVOLUTION<TV>* solid_fluid_coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION<TV>*>(driver->example.solids_evolution)){
        INCOMPRESSIBLE_UNIFORM<T_GRID>* incompressible=fluids_parameters->incompressible;
        if(Fluid_Node()){
            solid_fluid_coupled_evolution->Apply_Solid_Boundary_Conditions(time,false,face_velocities);
            incompressible->projection.Make_Divergence_Free(face_velocities,dt,time);}
        T target_time=time+steady_state_drag_time;
        assert(steady_state_drag_time>(T)0);
        ARRAY<TRIPLE<int,bool,bool> > saved_rigid_body_state(0);

        ARRAY<int> cluster_bodies;Make_Cluster_List(root_particle_index,cluster_bodies);
        for(int i=1;i<=cluster_bodies.m;++i){
            RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(cluster_bodies(i));rigid_body.V()=TV();rigid_body.Angular_Velocity()=T_SPIN();
            saved_rigid_body_state.Append(TRIPLE<int,bool,bool>(rigid_body.particle_index,rigid_body.is_static,!not_affected_by_fluid.Contains(rigid_body.particle_index)));
            rigid_body.is_static=true;not_affected_by_fluid.Set(rigid_body.particle_index);}

        LOG::SCOPE scope(STRING_UTILITIES::string_sprintf("Solving system to steady state time %f",target_time));
        bool fluid_affects_solid=fluids_parameters->fluid_affects_solid;fluids_parameters->fluid_affects_solid=false;
        bool write_output_files=driver->example.write_output_files;driver->example.write_output_files=false;
        driver->Advance_To_Target_Time(target_time);
        driver->example.write_output_files=write_output_files;fluids_parameters->fluid_affects_solid=fluid_affects_solid;
        driver->Set_Time(time);
        for(int i=1;i<=saved_rigid_body_state.m;++i){
            RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(saved_rigid_body_state(i).x);
            rigid_body.is_static=saved_rigid_body_state(i).y;
            if(saved_rigid_body_state(i).z) not_affected_by_fluid.Delete(rigid_body.particle_index);
            else not_affected_by_fluid.Set(rigid_body.particle_index);}}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After Project_Velocities (search_controller)",0,1);
}
//#####################################################################
// Function Set_PD_Targets
//#####################################################################
template<class T> VECTOR<T,1> Set_PD_Targets_Helper(RIGID_BODY<VECTOR<T,1> > *parent,JOINT<VECTOR<T,1> >* joint)
{PHYSBAM_NOT_IMPLEMENTED();}
template<class T> VECTOR<T,1> Set_PD_Targets_Helper(RIGID_BODY<VECTOR<T,2> > *parent,JOINT<VECTOR<T,2> >* joint)
{return joint->joint_function->Angular_Velocity();}
template<class T> VECTOR<T,3> Set_PD_Targets_Helper(RIGID_BODY<VECTOR<T,3> > *parent,JOINT<VECTOR<T,3> >* joint)
{return (parent->Rotation()*joint->F_pj().r).Inverse_Rotate(joint->joint_function->Angular_Velocity());}
template<class TV> void SEARCH_CONTROLLER<TV>::
Set_PD_Targets()
{
    if(solid_body_collection.rigid_body_collection.articulated_rigid_body.Has_Actuators()){
        if(current_joints.m>0) for(int i=1;i<=current_joints.m;i++){assert(joint_mesh.Is_Active(current_joints(i)));
            JOINT<TV>& joint=*joint_mesh(current_joints(i));
            bool control=false;for(int j=1;j<=T_SPIN::dimension;j++) if(joint.control_dof(j)) control=true;
            if(!joint.joint_function || !control) continue;
            joint.joint_function->active=true;joint.global_post_stabilization=true;
            joint.joint_function->Set_Target_Angular_Velocity(Set_PD_Targets_Helper(solid_body_collection.rigid_body_collection.articulated_rigid_body.Parent(joint.id_number),&joint));}
        else for(JOINT_ID i(1);i<=joint_mesh.Size();i++){if(!joint_mesh.Is_Active(i)) continue;
            JOINT<TV>& joint=*joint_mesh(i);
            bool control=false;for(int j=1;j<=T_SPIN::dimension;j++) if(joint.control_dof(j)) control=true;
            if(!joint.joint_function || !control) continue;
            joint.joint_function->active=true;joint.global_post_stabilization=true;
            joint.joint_function->Set_Target_Angular_Velocity(Set_PD_Targets_Helper(solid_body_collection.rigid_body_collection.articulated_rigid_body.Parent(joint.id_number),&joint));}}
}
//#####################################################################
// Function Zero_PD_Targets
//#####################################################################
template<class TV> void SEARCH_CONTROLLER<TV>::
Zero_PD_Targets()
{
    if(solid_body_collection.rigid_body_collection.articulated_rigid_body.Has_Actuators()){
        for(JOINT_ID i(1);i<=joint_mesh.Size();i++){if(!joint_mesh.Is_Active(i) || !joint_mesh(i)->joint_function) continue;
            JOINT<TV>& joint=*joint_mesh(i);
            joint.joint_function->active=true;joint.global_post_stabilization=true;
            joint.joint_function->Set_Target_Angle(joint.joint_function->Angle());}}
}
//#####################################################################
// Function Create_Clusters_From_Joint_List
//#####################################################################
template<class TV> void SEARCH_CONTROLLER<TV>::
Create_Clusters_From_Joint_List(const ARRAY<bool,JOINT_ID>& blocking_joint,ARRAY<ARRAY<int> >& body_lists,ARRAY<VECTOR<int,2>,JOINT_ID>& adjacent_lists)
{
    UNDIRECTED_GRAPH<int,JOINT_ID>& graph=solid_body_collection.rigid_body_collection.articulated_rigid_body.joint_mesh.undirected_graph;
    // 0 = not seen, o.w. = body list index last visited for
    adjacent_lists.Resize(solid_body_collection.rigid_body_collection.articulated_rigid_body.joint_mesh.Size());
    ARRAY<int,int> done(graph.Last_Node());
    STACK<int> stack;
    for(JOINT_ID j(1);j<=blocking_joint.Size();j++) if(blocking_joint(j)){
        VECTOR<int,2> adjacent_list_j(adjacent_lists(j).x,adjacent_lists(j).y);
        for(int k=1;k<=2;k++) if(int b=(k==1?graph.Edges(j).x:graph.Edges(j).y)){
            RIGID_BODY<TV>& rbody=solid_body_collection.rigid_body_collection.Rigid_Body(b);
            if(!rbody.Has_Infinite_Inertia() && done(b)){adjacent_list_j(k)=done(b);continue;}
            stack.Push(b);
            adjacent_list_j(k)=done(b)=body_lists.Append(ARRAY<int>());
            while(!stack.Empty()){
                int id=stack.Pop();
                RIGID_BODY<TV>& body=solid_body_collection.rigid_body_collection.Rigid_Body(id);
                body_lists.Last().Append(body.particle_index);
                if(body.Has_Infinite_Inertia()) continue;
                for(int i=1;i<=graph.Adjacent_Edges(id).m;i++){
                    JOINT_ID jid=graph.Adjacent_Edges(id)(i);
                    if(blocking_joint(jid)) continue;
                    int oid=graph.Edges(jid).y+(Value(graph.Edges(jid).x)-Value(id));
                    RIGID_BODY<TV>& obody=solid_body_collection.rigid_body_collection.Rigid_Body(oid);
                    if((!obody.Has_Infinite_Inertia() && done(oid)) || done(oid)==body_lists.m) continue;
                    done(oid)=body_lists.m;
                    stack.Push(oid);}}}}
}
//#####################################################################
// Function Create_All_Clusters
//#####################################################################
template<class TV> void SEARCH_CONTROLLER<TV>::
Create_All_Clusters(RIGID_BODY_COLLISION_MANAGER_HASH* collision_manager)
{
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings;
    assert(root_particle_index>=0);
    ARRAY<bool,JOINT_ID> blocking_list;ARRAY<int> done;
    blocking_list.Resize(joint_mesh.Size());joint_clusters.Resize(joint_mesh.Size());
    use_clusters=true;rigid_body_particles_number=solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();

    ARRAY<JOINT<TV>*> joints;
    for(JOINT_ID i(1);i<=joint_mesh.Size();i++){JOINT_ID joint_id=i;
        if(!joint_mesh.Is_Active(i)) continue;
        bool controlled=false;
        for(int j=1;j<=T_SPIN::dimension;j++) if(joint_mesh(i)->control_dof(j)) controlled=true;
        if(!controlled) continue;
        int parent_cluster=Parent(joint_id)->particle_index,child_cluster=Child(joint_id)->particle_index;
        ARRAY<bool,JOINT_ID> blocking_list_single;blocking_list_single.Resize(joint_mesh.Size());
        blocking_list_single(joint_id)=true;blocking_list(joint_id)=true;
        ARRAY<ARRAY<int> > body_lists;ARRAY<VECTOR<int,2>,JOINT_ID> adjacent_lists;
        Create_Clusters_From_Joint_List(blocking_list_single,body_lists,adjacent_lists);
        if(body_lists(adjacent_lists(joint_id).x).Size()>1){
            ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children;
            for(int i=1;i<=body_lists(adjacent_lists(joint_id).x).Size();i++) children.Append(body_lists(adjacent_lists(joint_id).x)(i));
            joint_clusters(joint_id).x=rigid_bindings.Add_Binding(children);
            RIGID_BODY<TV>* rigid_body_cluster=&solid_body_collection.rigid_body_collection.Rigid_Body(joint_clusters(joint_id).x);
            parent_cluster=rigid_body_cluster->particle_index;}
        if(body_lists(adjacent_lists(joint_id).y).Size()>1){
            ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children;
            for(int i=1;i<=body_lists(adjacent_lists(joint_id).y).Size();i++) children.Append(body_lists(adjacent_lists(joint_id).y)(i));
            joint_clusters(joint_id).y=rigid_bindings.Add_Binding(children);
            RIGID_BODY<TV>* rigid_body_cluster=&solid_body_collection.rigid_body_collection.Rigid_Body(joint_clusters(joint_id).y);
            child_cluster=rigid_body_cluster->particle_index;}
        //We created a binding so we need to update structures
        if(joint_clusters(joint_id).x || joint_clusters(joint_id).y){
            //update collisions
            if(collision_manager && !collision_manager->hash.Contains(PAIR<int,int>(parent_cluster,child_cluster))){
                collision_manager->hash.Insert(PAIR<int,int>(parent_cluster,child_cluster));
                collision_manager->hash.Insert(PAIR<int,int>(child_cluster,parent_cluster));}
            //deactive new clusters and joints
            if(joint_clusters(joint_id).x) rigid_bindings.Set_Binding_Active(joint_clusters(joint_id).x,false,(incorporate_fluids?fluids_parameters->collision_bodies_affecting_fluid:0));
            if(joint_clusters(joint_id).y) rigid_bindings.Set_Binding_Active(joint_clusters(joint_id).y,false,(incorporate_fluids?fluids_parameters->collision_bodies_affecting_fluid:0));}}
    
    ARRAY<ARRAY<int> > body_lists;ARRAY<VECTOR<int,2>,JOINT_ID> adjacent_lists;
    Create_Clusters_From_Joint_List(blocking_list,body_lists,adjacent_lists);done.Resize(body_lists.Size());
    for(JOINT_ID i(1);i<=joint_mesh.Size();i++){JOINT_ID joint_id=i;
        if(!joint_mesh.Is_Active(i)) continue;
        int parent_cluster=Parent(joint_id)->particle_index,child_cluster=Child(joint_id)->particle_index;
        if(adjacent_lists(joint_id).x>0 && body_lists(adjacent_lists(joint_id).x).Size()>1){
            int cluster_particle=done(adjacent_lists(joint_id).x);
            if(!cluster_particle){
                ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children;
                for(int i=1;i<=body_lists(adjacent_lists(joint_id).x).Size();i++) children.Append(body_lists(adjacent_lists(joint_id).x)(i));
                cluster_particle=rigid_bindings.Add_Binding(children);
                done(adjacent_lists(joint_id).x)=cluster_particle;
                global_clusters.Append(cluster_particle);}
            RIGID_BODY<TV>* rigid_body_cluster=&solid_body_collection.rigid_body_collection.Rigid_Body(cluster_particle);
            parent_cluster=rigid_body_cluster->particle_index;}
        if(adjacent_lists(joint_id).y>0 && body_lists(adjacent_lists(joint_id).y).Size()>1){
            int cluster_particle=done(adjacent_lists(joint_id).y);
            if(!cluster_particle){
                ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> children;
                for(int i=1;i<=body_lists(adjacent_lists(joint_id).y).Size();i++) children.Append(body_lists(adjacent_lists(joint_id).y)(i));
                cluster_particle=rigid_bindings.Add_Binding(children);
                done(adjacent_lists(joint_id).y)=cluster_particle;
                global_clusters.Append(cluster_particle);}
            RIGID_BODY<TV>* rigid_body_cluster=&solid_body_collection.rigid_body_collection.Rigid_Body(cluster_particle);
            child_cluster=rigid_body_cluster->particle_index;}
        //We created a binding so we need to update structures
        if(parent_cluster!=Parent(joint_id)->particle_index || child_cluster!=Child(joint_id)->particle_index){
            if(collision_manager && !collision_manager->hash.Contains(PAIR<int,int>(parent_cluster,child_cluster))){
                collision_manager->hash.Insert(PAIR<int,int>(parent_cluster,child_cluster));
                collision_manager->hash.Insert(PAIR<int,int>(child_cluster,parent_cluster));}}}
    //deactive new clusters
    for(int i=1;i<=global_clusters.m;i++) rigid_bindings.Set_Binding_Active(global_clusters(i),false,(incorporate_fluids?fluids_parameters->collision_bodies_affecting_fluid:0));
}
//#####################################################################
// Function Remove_Rigid
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Remove_Clusters()
{
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings;
    for(int i=1;i<=global_clusters.m;i++) rigid_bindings.Delete_Binding(global_clusters(i));
    global_clusters.Remove_All();
}
//#####################################################################
// Function Add_Solid_Drag
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Add_Solid_Drag(const T time,ARRAY<TV> &F,ARRAY<TWIST<TV> > &rigid_F)
{
    for(int i=1;i<=solid_body_collection.solids_forces.m;i++){
        if(dynamic_cast<WIND_DRAG<TV>*>(&*solid_body_collection.solids_forces(i))||dynamic_cast<ETHER_DRAG<T_GRID>*>(&*solid_body_collection.solids_forces(i))){
            ARRAY<TWIST<TV> > twist;twist.Resize(solid_body_collection.rigid_body_collection.rigid_body_particle.V.Size());
            for(int j=1;j<=twist.Size();j++) twist(j)=TWIST<TV>(solid_body_collection.rigid_body_collection.rigid_body_particle.V(j),solid_body_collection.rigid_body_collection.rigid_body_particle.angular_velocity(j));
            solid_body_collection.solids_forces(i)->Update_Position_Based_State(time);
            solid_body_collection.solids_forces(i)->Add_Velocity_Independent_Forces(F,rigid_F,time);
            solid_body_collection.solids_forces(i)->Add_Velocity_Dependent_Forces(solid_body_collection.deformable_body_collection.particles.V,twist,F,rigid_F,time);
            for(int j=1;j<=twist.Size();j++){solid_body_collection.rigid_body_collection.rigid_body_particle.V(j)=twist(j).linear;solid_body_collection.rigid_body_collection.rigid_body_particle.angular_velocity(j)=twist(j).angular;}}}
}
//#####################################################################
// Function Add_Fluid_Drag
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Add_Fluid_Drag(const T dt,const T time,ARRAY<TV> &F,ARRAY<TWIST<TV> > &rigid_F)
{
    ARRAY<TV> F_neg(F.m);ARRAY<TWIST<TV> >rigid_F_neg(solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size());
    if(SOLID_FLUID_COUPLED_EVOLUTION<TV>* solid_fluid_coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION<TV>*>(driver->example.solids_evolution)){
        ARRAY<INTERVAL<int> > interior_regions;
        { // Set up matrix indices so that we can compute J, J_rigid, and J_rigid_kinematic
            // TODO(jontg): This would be easier if poisson saved these mappings...
            ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array=solid_fluid_coupled_evolution->matrix_index_to_cell_index_array;
            T_ARRAYS_INT cell_index_to_matrix_index;
            int number_of_regions=0;

            if(Fluid_Node()){
                POISSON_COLLIDABLE_UNIFORM<T_GRID>* poisson=dynamic_cast<POISSON_COLLIDABLE_UNIFORM<T_GRID>*>(fluids_parameters->incompressible->projection.elliptic_solver);
                number_of_regions=poisson->number_of_regions;

                ARRAY<int,VECTOR<int,1> > filled_region_cell_count(-1,number_of_regions);
                matrix_index_to_cell_index_array.Resize(number_of_regions);cell_index_to_matrix_index.Resize(fluids_parameters->grid->Domain_Indices(1));
                for(CELL_ITERATOR iterator(*fluids_parameters->grid,1);iterator.Valid();iterator.Next()) filled_region_cell_count(poisson->filled_region_colors(iterator.Cell_Index()))++;
                for(int color=1;color<=number_of_regions;color++) if(poisson->filled_region_touches_dirichlet(color)||poisson->solve_neumann_regions){
                        matrix_index_to_cell_index_array(color).Resize(filled_region_cell_count(color));}
                filled_region_cell_count.Fill(0); // reusing this array in order to make the indirection arrays

                if(!poisson->mpi_grid) poisson->Compute_Matrix_Indices(filled_region_cell_count,matrix_index_to_cell_index_array,cell_index_to_matrix_index);
                else poisson->laplace_mpi->Find_Matrix_Indices(filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);

                int colors=filled_region_cell_count.domain.max_corner.x;
                interior_regions.Resize(colors);
                for(int color=1;color<=colors;color++) if(filled_region_cell_count(color)>0){
                    if(!poisson->laplace_mpi || color>poisson->laplace_mpi->partitions.m) interior_regions(color)=INTERVAL<int>(1,filled_region_cell_count(color));
                    else interior_regions(color)=poisson->laplace_mpi->partitions(color).interior_indices;}}
            solid_fluid_coupled_evolution->Compute_Coupling_Terms_Rigid(cell_index_to_matrix_index,interior_regions,number_of_regions);
            solid_fluid_coupled_evolution->Compute_Coupling_Terms_Deformable(cell_index_to_matrix_index,interior_regions,number_of_regions);
        }

        ARRAY<TWIST<TV> > kinematic_rigid_bodies(solid_fluid_coupled_evolution->rigid_body_count);
        //TODO(mlentine,jontg,avir,cas43): This needs to be changed along with solid-fluid and multiple super-fragments
        GENERALIZED_VELOCITY<TV> V(F_neg,rigid_F_neg,solid_body_collection);
        if(Fluid_Node()){
            T_ARRAYS_SCALAR& p=fluids_parameters->incompressible->projection.p;p*=dt;
            PHYSBAM_DEBUG_WRITE_SUBSTEP("Rescaling pressures to be viewable for Add_Fluid_Drag (search_controller)",0,0);
            VECTOR_ND<T> x_array;
            GENERALIZED_VELOCITY<TV> V_with_kinematic_rigid(F_neg,kinematic_rigid_bodies,solid_body_collection);
            for(int i=1;i<=solid_fluid_coupled_evolution->J_rigid_kinematic.m;++i){
                x_array.Resize(solid_fluid_coupled_evolution->matrix_index_to_cell_index_array(i).Size());
                for(int j=1;j<=solid_fluid_coupled_evolution->matrix_index_to_cell_index_array(i).Size();j++) x_array(j)=p(solid_fluid_coupled_evolution->matrix_index_to_cell_index_array(i)(j));

                VECTOR_ND<T> x_array_i;x_array_i.Set_Subvector_View(x_array,interior_regions(i));
                SOLID_FLUID_SYSTEM<TV,SPARSE_MATRIX_FLAT_NXN<T> >::Add_J_Rigid_Times_Pressure(solid_fluid_coupled_evolution->J_rigid_kinematic(i),x_array_i,V_with_kinematic_rigid);
                if(i <= solid_fluid_coupled_evolution->J_rigid.m) SOLID_FLUID_SYSTEM<TV,SPARSE_MATRIX_FLAT_NXN<T> >::Add_J_Rigid_Times_Pressure(solid_fluid_coupled_evolution->J_rigid(i),x_array_i,V);
                if(i <= solid_fluid_coupled_evolution->J_deformable.m) SOLID_FLUID_SYSTEM<TV,SPARSE_MATRIX_FLAT_NXN<T> >::Add_J_Deformable_Times_Pressure(solid_fluid_coupled_evolution->J_deformable(i),x_array_i,V);}
            for(int i=1;i<=kinematic_rigid_bodies.m;++i) rigid_F_neg(i)+=kinematic_rigid_bodies(i);}

        if(driver->example.solids_fluids_parameters.mpi_solid_fluid) driver->example.solids_fluids_parameters.mpi_solid_fluid->Aggregate_Lists_To_Solid_Node(V);}
    else PHYSBAM_FATAL_ERROR("Cannot find the coupled driver!");
   
    F-=F_neg;rigid_F-=rigid_F_neg;
}
//#####################################################################
// Function Evaluate_Drag_For_Joint
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR SEARCH_CONTROLLER<T_GRID>::
Evaluate_Drag_For_Joint(const T dt,const T time,ARRAY<JOINT_ID>& joints)
{
    return Evaluate_Drag(dt,time);
}
//#####################################################################
// Function Evaluate_Drag
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR SEARCH_CONTROLLER<T_GRID>::
Evaluate_Drag(const T dt,const T time)
{
    T Force=0;ARRAY<TV> F(solid_body_collection.deformable_body_collection.particles.array_collection->Size());ARRAY<TWIST<TV> >rigid_F(solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size());
    if(!incorporate_fluids) Add_Solid_Drag(time,F,rigid_F);
    if(incorporate_fluids) Add_Fluid_Drag(dt,time,F,rigid_F);
    solid_body_collection.rigid_body_collection.articulated_rigid_body.Initialize_Poststabilization_Projection();
    solid_body_collection.rigid_body_collection.articulated_rigid_body.Poststabilization_Projection(rigid_F,true);
    if(use_drag_direction){
        for(int i=1;i<=F.Size();i++) Force+=TV::Dot_Product(F(i),drag_direction);
        for(int i=1;i<=rigid_F.Size();i++) Force+=TV::Dot_Product(rigid_F(i).linear,drag_direction);}
    else{
        for(int i=1;i<=F.Size();i++) Force+=F(i).Magnitude();
        for(int i=1;i<=rigid_F.Size();i++) Force+=rigid_F(i).linear.Magnitude();}
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Evaluate Drag End (search_controller) time=%f force felt=%f",time,Force),0,0);
    return Force;
}
//#####################################################################
// Function Evaluate_Force_To_Stay_For_Joint
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR SEARCH_CONTROLLER<T_GRID>::
Evaluate_Force_To_Stay_For_Joint(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,ARRAY<JOINT_ID>& joints)
{
    assert(dt>0);
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Evaluate Force Start (search_controller) time=%f",time),0,0);
    T F=0;current_joints=joints;Save_PD_State();
    Save_Nested_State(face_velocities);
    for(JOINT_ID i(1);i<=joint_mesh.Size();i++){if(!joint_mesh.Is_Active(i)) continue;
        joint_mesh(i)->global_post_stabilization=false;
        joint_mesh(i)->joint_function->active=false;}
    for(int i=1;i<=joints.m;i++){assert(joint_mesh.Is_Active(joints(i)));
        JOINT<TV>& joint=*joint_mesh(joints(i));
        joint.impulse_accumulator->Reset();
        joint.joint_function->active=true;joint.global_post_stabilization=true;}
    bool fluid_affects_solid=false;
    Set_PD_Targets();pd_step=true;
    solid_body_collection.rigid_body_collection.articulated_rigid_body.constrain_pd_directions=false;
    if(incorporate_fluids){
        fluid_affects_solid=fluids_parameters->fluid_affects_solid;fluids_parameters->fluid_affects_solid=true;
        driver->Advance_To_Target_Time(time+dt);Advance_One_Time_Step_Position(face_velocities,dt,time+dt);
        fluids_parameters->fluid_affects_solid=fluid_affects_solid;}
    else Advance_One_Time_Step_Position(face_velocities,dt,time);
    pd_step=false;
    Restore_Nested_State(face_velocities);
    Restore_PD_State();current_joints.Remove_All();
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Evaluate Force End (search_controller) time=%f",time),0,0);
    for(int i=1;i<=joints.m;i++){assert(joint_mesh.Is_Active(joints(i)));
        F+=joint_mesh(joints(i))->impulse_accumulator->Energy();}
    return F/joints.m;
}
//#####################################################################
// Function Evaluate_Force_To_Stay
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR SEARCH_CONTROLLER<T_GRID>::
Evaluate_Force_To_Stay(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,int hack)
{
    assert(dt>0);
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Evaluate Force Start (search_controller) time=%f",time),0,0);
    T F=0;Save_PD_State();
    Save_Nested_State(face_velocities);
    for(JOINT_ID i(1);i<=joint_mesh.Size();i++){if(!joint_mesh.Is_Active(i)) continue;
        JOINT<TV>& joint=*joint_mesh(i);
        joint.impulse_accumulator->Reset();
        bool control=false;for(int j=1;j<=T_SPIN::dimension;j++) if(joint.control_dof(j)) control=true;
        if(!joint.joint_function || !control){joint.joint_function->active=false;joint.global_post_stabilization=false;continue;} //Turn off joints we don't control
        else{joint.joint_function->active=true;joint.global_post_stabilization=true;}
        if(hack) joint.joint_function->Set_Target_Angle(joint.joint_function->Angle());}
    if(hack){solid_body_collection.deformable_body_collection.particles.V*=(T)0;solid_body_collection.rigid_body_collection.rigid_body_particle.V*=(T)0;solid_body_collection.rigid_body_collection.rigid_body_particle.angular_velocity*=(T)0;solid_body_collection.rigid_body_collection.rigid_body_particle.angular_momentum*=(T)0;}
    else Set_PD_Targets();
    pd_step=true;
    bool fluid_affects_solid=false;
    solid_body_collection.rigid_body_collection.articulated_rigid_body.constrain_pd_directions=false;
    if(incorporate_fluids){
        fluid_affects_solid=fluids_parameters->fluid_affects_solid;fluids_parameters->fluid_affects_solid=true;
        driver->Advance_To_Target_Time(time+dt);Advance_One_Time_Step_Position(face_velocities,dt,time+dt);
        fluids_parameters->fluid_affects_solid=fluid_affects_solid;}
    else Advance_One_Time_Step_Position(face_velocities,dt,time);
    pd_step=false;
    Restore_Nested_State(face_velocities);
    Restore_PD_State();
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Evaluate Force End (search_controller) time=%f",time),0,0);
    for(JOINT_ID i(1);i<=joint_mesh.Size();i++){if(!joint_mesh.Is_Active(i)) continue;
        F+=joint_mesh(i)->impulse_accumulator->Energy();}
    return F;
}
//#####################################################################
// Function Take_Hypothetical_Step
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Take_Hypothetical_Step(T_FACE_ARRAYS_SCALAR& face_velocities,const T time,bool simulate_fluids)
{
    T_FACE_ARRAYS_SCALAR face_velocities_current;
    T_ARRAYS_SCALAR pressure_current;
    MPI_SOLID_FLUID<TV>* mpi_solid_fluid=driver->example.solids_fluids_parameters.mpi_solid_fluid;
    if(incorporate_fluids && Fluid_Node()){
        if(!simulate_fluids) fluids_parameters->simulate=false;
        face_velocities_current.Resize(*fluids_parameters->grid);pressure_current.Resize(fluids_parameters->grid->Domain_Indices(1));
        T_FACE_ARRAYS_SCALAR::Copy(face_velocities,face_velocities_current);
        T_ARRAYS_SCALAR::Copy(fluids_parameters->incompressible->projection.p,pressure_current);}
    if(!simulate_fluids) driver->example.solids_fluids_parameters.mpi_solid_fluid=0;
    driver->Advance_To_Target_Time(time);
    if(!simulate_fluids) driver->example.solids_fluids_parameters.mpi_solid_fluid=mpi_solid_fluid;
    if(incorporate_fluids && Fluid_Node()){
        if(!simulate_fluids) fluids_parameters->simulate=true;
        T_FACE_ARRAYS_SCALAR::Copy(face_velocities_current,face_velocities);
        T_ARRAYS_SCALAR::Copy(pressure_current,fluids_parameters->incompressible->projection.p);}
    if(driver->example.solids_fluids_parameters.mpi_solid_fluid) driver->example.solids_fluids_parameters.mpi_solid_fluid->Exchange_Solid_Positions_And_Velocities(solid_body_collection);
}
//#####################################################################
// Function Take_Hypothetical_Fluids_Step
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Take_Hypothetical_Fluids_Step(const T time)
{
    ARRAY<TRIPLE<int,bool,bool> > saved_rigid_body_state(0);
    assert(fluids_parameters->solid_affects_fluid);
    bool write_frame_boundaries=driver->project_at_frame_boundaries;
    //TODO(jontg): Remove affected_by_fluid flag
    for(int i(1);i<=solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();++i){
        RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(i);rigid_body.V()=TV();rigid_body.Angular_Velocity()=T_SPIN();
        saved_rigid_body_state.Append(TRIPLE<int,bool,bool>(rigid_body.particle_index,rigid_body.is_static,!not_affected_by_fluid.Contains(rigid_body.particle_index)));
        rigid_body.is_static=true;not_affected_by_fluid.Set(rigid_body.particle_index);}
    driver->project_at_frame_boundaries=true;drag_step=true;
    driver->Advance_To_Target_Time(time);
    for(int i=1;i<=saved_rigid_body_state.m;++i){
        RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(saved_rigid_body_state(i).x);
        if(saved_rigid_body_state(i).z) not_affected_by_fluid.Delete(rigid_body.particle_index);
        else not_affected_by_fluid.Set(rigid_body.particle_index);
        rigid_body.is_static=saved_rigid_body_state(i).y;}
    driver->project_at_frame_boundaries=write_frame_boundaries;drag_step=false;
}
//#####################################################################
// Function Evaluate_Force_For_Joint
//#####################################################################
template<class T_GRID> PAIR<typename T_GRID::SCALAR,typename T_GRID::VECTOR_T::SPIN> SEARCH_CONTROLLER<T_GRID>::
Evaluate_Force_For_Joint(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,const JOINT_ID joint_id,const int dimension,const T dx_multiplier)
{
    JOINT<TV>& joint=*joint_mesh(joint_id);
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings;
    T F=0;hypothetical_step=true;
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Hypothetical Step Start (search_controller) old_time=%f time=%f",time_save,time),0,0);
    
    Save_State(face_velocities);
    
    ARRAY<int> active_clusters;
    if(use_clusters){
        rigid_bindings.Deactivate_And_Return_Clusters(active_clusters,(incorporate_fluids?fluids_parameters->collision_bodies_affecting_fluid:0));
        if(joint_clusters(joint.id_number).x) rigid_bindings.Set_Binding_Active(joint_clusters(joint.id_number).x,true);
        if(joint_clusters(joint.id_number).y) rigid_bindings.Set_Binding_Active(joint_clusters(joint.id_number).y,true);
        solid_body_collection.rigid_body_collection.Update_Angular_Momentum();}
    
    Save_PD_State();

    Set_PD_Targets();
    T_SPIN angles;
    angles(dimension)=dx*dx_multiplier;
    ROTATION<TV> target_delta=ROTATION<TV>::From_Euler_Angles(angles);
    ROTATION<TV> target_angle=target_delta*joint.joint_function->Angle();
    T_SPIN new_euler_angles=target_angle.Euler_Angles();joint.Constrain_Angles(new_euler_angles);
    angles=new_euler_angles-joint.joint_function->Angle().Euler_Angles();
    joint.joint_function->Set_Target_Angle(ROTATION<TV>::From_Euler_Angles(new_euler_angles));
    joint.joint_function->Set_k_p(100000); //high constant
    if(incorporate_fluids){driver->project_at_frame_boundaries=false;fluids_parameters->fluid_affects_solid=false;}
    
    Take_Hypothetical_Step(face_velocities,time,use_projection && incorporate_fluids);
    
    Restore_PD_State();
    assert(!joint.joint_function->active);assert(!joint.global_post_stabilization);
    
    if(use_projection) Project_Velocities(face_velocities,dt_hyp,time);
    if(use_clusters){
        if(joint_clusters(joint.id_number).x) rigid_bindings.Set_Binding_Active(joint_clusters(joint.id_number).x,false,(incorporate_fluids?fluids_parameters->collision_bodies_affecting_fluid:0));
        if(joint_clusters(joint.id_number).y) rigid_bindings.Set_Binding_Active(joint_clusters(joint.id_number).y,false,(incorporate_fluids?fluids_parameters->collision_bodies_affecting_fluid:0));
        rigid_bindings.Reactivate_Bindings(active_clusters);
        solid_body_collection.rigid_body_collection.Update_Angular_Momentum();}

    if(!use_projection && incorporate_fluids) Take_Hypothetical_Fluids_Step(time+dt_hyp);

    ARRAY<JOINT_ID> force_list,drag_list;
    for(int i=1;i<=strain_joints.m;i++) force_list.Append(strain_joints(i)); 
    if(objective(joint_id)==FORCE) force_list.Append_Unique(joint_id);
    if(objective(joint_id)==DRAG){drag_list.Append(joint_id);F=Evaluate_Drag_For_Joint(dt_hyp,time,drag_list);}
    if(force_list.m>0) F+=Evaluate_Force_To_Stay_For_Joint(face_velocities,dt,time,force_list);
 
    Restore_State(face_velocities);
    
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Hypothetical Step End (search controller) old_time=%f time=%f force=%f",time_save,time,F),0,0);
    hypothetical_step=false;
    return PAIR<T,T_SPIN>(F,angles);
}
//#####################################################################
// Function Evaluate_Force
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR SEARCH_CONTROLLER<T_GRID>::
Evaluate_Force(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,const T dx,const T magnitude,const VECTOR_ND<T>* grad)
{
    //TODO(mlentine): Rewrite this to be consistant with Force_For_Joint
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings;
    T F=0;hypothetical_step=true;
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Hypothetical Step Start (search controller) old_time=%f time=%f",time_save,time),0,0);
    
    Save_State(face_velocities);
    
    ARRAY<int> active_clusters;
    if(use_clusters){
        rigid_bindings.Deactivate_And_Return_Clusters(active_clusters,(incorporate_fluids?fluids_parameters->collision_bodies_affecting_fluid:0));
        for(int i=1;i<=global_clusters.m;i++) rigid_bindings.Set_Binding_Active(global_clusters(i),true);
        solid_body_collection.rigid_body_collection.Update_Angular_Momentum();}

    Save_PD_State();
    
    for(JOINT_ID i(1);i<=joint_mesh.Size();i++){if(!joint_mesh.Is_Active(i)) continue;
        JOINT<TV>& joint=*joint_mesh(i);
        bool control=false;for(int j=1;j<=T_SPIN::dimension;j++) if(joint.control_dof(j)) control=true;
        if(!joint.joint_function || !control) continue;
        T_SPIN angles;for(int j=1;j<=T_SPIN::dimension;j++) angles(j)=grad?dx*(*grad)((Value(i)-1)*T_SPIN::dimension+j):0;
        ROTATION<TV> target_delta=ROTATION<TV>::From_Euler_Angles(angles);
        ROTATION<TV> target_angle=target_delta*joint.joint_function->Angle();
        T_SPIN new_euler_angles=target_angle.Euler_Angles();joint.Constrain_Angles(new_euler_angles);
        angles=new_euler_angles-joint.joint_function->Angle().Euler_Angles();
        joint.joint_function->Set_Target_Angle(ROTATION<TV>::Spherical_Linear_Interpolation(joint.joint_function->Angle(),ROTATION<TV>::From_Euler_Angles(new_euler_angles),magnitude));
        joint.joint_function->Set_k_p(100000); //high constant
        joint.joint_function->active=true;joint.global_post_stabilization=true;}
    if(incorporate_fluids){driver->project_at_frame_boundaries=false;fluids_parameters->fluid_affects_solid=false;}
    
    Take_Hypothetical_Step(face_velocities,time,use_projection && incorporate_fluids);

    Restore_PD_State();

    if(use_projection) Project_Velocities(face_velocities,dt_hyp,time);
    if(use_clusters){
        rigid_bindings.Reactivate_Bindings(active_clusters);
        solid_body_collection.rigid_body_collection.Update_Angular_Momentum();}

    if(!use_projection && incorporate_fluids) Take_Hypothetical_Fluids_Step(time+dt_hyp);

    ARRAY<JOINT_ID> force_list,drag_list;
    for(int i=1;i<=strain_joints.m;i++) force_list.Append(strain_joints(i)); 
    for(JOINT_ID i(1);i<=joint_mesh.Size();i++){if(!joint_mesh.Is_Active(i)) continue;
        if(objective(i)==FORCE) force_list.Append_Unique(i);
        if(objective(i)==DRAG) drag_list.Append(i);}
    if(force_list.m>0) F=Evaluate_Force_To_Stay_For_Joint(face_velocities,dt,time,force_list);
    if(drag_list.m>0) F+=Evaluate_Drag_For_Joint(dt,time,drag_list);
    
    Restore_State(face_velocities);

    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Hypothetical Step End (search controller) old_time=%f time=%f force=%f",time_save,time,F),0,0);
    hypothetical_step=false;
    return F;
}
//#####################################################################
// Function Normalize_Gradient
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Normalize_Gradient(VECTOR_ND<T>& grad)
{
    if(grad.Magnitude()) for(JOINT_ID i(1);i<=joint_mesh.Size();i++){
        if(!joint_mesh.Is_Active(i)) continue;
        TV joint_gradient;
        for(int j=1;j<=T_SPIN::dimension;j++) joint_gradient(j)=grad((Value(i)-1)*T_SPIN::dimension+j);
        if(joint_gradient.Magnitude()) joint_gradient.Normalize();
        for(int j=1;j<=T_SPIN::dimension;j++) grad((Value(i)-1)*T_SPIN::dimension+j)=joint_gradient(j);}
}
//#####################################################################
// Function Normalized_Gradient
//#####################################################################
template<class T_GRID> VECTOR_ND<typename T_GRID::SCALAR> SEARCH_CONTROLLER<T_GRID>::
Normalized_Gradient(const VECTOR_ND<T>& grad)
{
    VECTOR_ND<T> normalized_gradient;
    if(grad.Magnitude()) for(JOINT_ID i(1);i<=joint_mesh.Size();i++){
        if(!joint_mesh.Is_Active(i)) continue;
        TV joint_gradient;
        for(int j=1;j<=T_SPIN::dimension;j++) joint_gradient(j)=grad((Value(i)-1)*T_SPIN::dimension+j);
        if(joint_gradient.Magnitude()) joint_gradient.Normalize();
        for(int j=1;j<=T_SPIN::dimension;j++) normalized_gradient((Value(i)-1)*T_SPIN::dimension+j)=joint_gradient(j);}
    return normalized_gradient;
}
//#####################################################################
// Function Compute_Gradient
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Compute_Gradient(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    for(JOINT_ID i(1);i<=joint_mesh.Size();i++){if(!joint_mesh.Is_Active(i)) continue;
        JOINT<TV>& joint=*joint_mesh(i);
        if(!joint.joint_function) continue;
        VECTOR<bool,T_SPIN::dimension> constrain=joint.Angular_Constraints();
        bool control=false;for(int j=1;j<=T_SPIN::dimension;j++) if(joint.control_dof(j) && !constrain(j)) control=true;
        if(control){current_angle(i)=joint.joint_function->Angle();joint.joint_function->active=false;joint.global_post_stabilization=false;}}
    for(JOINT_ID i(1);i<=joint_mesh.Size();i++){if(!joint_mesh.Is_Active(i)) continue;
        JOINT<TV>& joint=*joint_mesh(i);
        if(!joint.joint_function) continue;
        VECTOR<bool,T_SPIN::dimension> constrain=joint.Angular_Constraints();
        RANDOM_NUMBERS<T> random;T n=random.Get_Uniform_Number((T)0,(T)1);
        for(int j=1;j<=T_SPIN::dimension;j++){
            if(constrain(j) || !joint.control_dof(j)) continue; //if we can't move the joint dof don't try
            PAIR<T,T_SPIN> force_forward=Evaluate_Force_For_Joint(face_velocities,dt,time,i,j,1);PAIR<T,T_SPIN> force_backward=Evaluate_Force_For_Joint(face_velocities,dt,time,i,j,-1);
            negative_gradient((Value(i)-1)*T_SPIN::dimension+j)=(minimize?(T)-1:(T)1)*(force_forward.x-force_backward.x)/(force_forward.y(j)-force_backward.y(j));
            if(use_random && n<threshold) negative_gradient((Value(i)-1)*T_SPIN::dimension+j)*=-1;
            if(abs(negative_gradient((Value(i)-1)*T_SPIN::dimension+j))<threshold) negative_gradient((Value(i)-1)*T_SPIN::dimension+j)=0;}}
}
//#####################################################################
// Function Compute_Gradient_From_Graph
//#####################################################################
template<class T_GRID> bool SEARCH_CONTROLLER<T_GRID>::
Compute_Gradient_From_Graph(ENVIRONMENTAL_STATE<T_GRID>* current_state)
{
    T threshold=(T)1e-3;
    //TODO: Read these from the example
    //if(FILE_UTILITIES::File_Exists("/disk2/bvuong/PhysBAM/Tools/preprocess_controller_data/test_graph")) FILE_UTILITIES::Read_From_File<float>("/disk2/bvuong/PhysBAM/Tools/preprocess_controller_data/test_graph",states_graph);
    //if(FILE_UTILITIES::File_Exists("/disk2/bvuong/PhysBAM/Tools/preprocess_controller_data/test_hashtable")) FILE_UTILITIES::Read_From_File<float>("/disk2/bvuong/PhysBAM/Tools/preprocess_controller_data/test_hashtable",hashtable);
    for(int i=1;i<=graph_index_to_state_map.Size();i++){assert(graph_index_to_state_map(i)->incorporate_fluids==incorporate_fluids);} //sanity check
    int hash_number=0;for(int i=1;i<=graph_index_to_state_map.Size();i++){if(*graph_index_to_state_map(i)==*current_state) hash_number=i;}
    if(hash_number==0) return false;
    ARRAY<int>& next_possible_states=states_graph->Children(hash_number);
    T minimized_force=0;ENVIRONMENTAL_STATE<T_GRID>* next_state=NULL;
    for(int i=1;i<=next_possible_states.Size();i++){
        next_possible_states.Append_Elements(states_graph->Children(next_possible_states(i)));
        ENVIRONMENTAL_STATE<T_GRID>* possible_next_state=graph_index_to_state_map(next_possible_states(i));
        if(minimized_force==0 || possible_next_state->force_to_stay<minimized_force){
            bool angle_is_possible=true;
            for(JOINT_ID id(1);id<=possible_next_state->angles.Size();id++){
                JOINT<TV>& joint=*joint_mesh(id);
                ROTATION<TV> angle=possible_next_state->angles(id);T_SPIN euler_angle=angle.Euler_Angles();
                T_SPIN angle_constrained=angle.Euler_Angles();joint.Constrain_Angles(angle_constrained);
                T sum=0;for(int k=1;k<=T_SPIN::dimension;k++){sum+=abs(angle_constrained(k)-euler_angle(k));}
                if(sum>threshold){angle_is_possible=false;break;}}
            if(angle_is_possible){minimized_force=possible_next_state->force_to_stay;next_state=possible_next_state;}}}
    if(next_state==NULL) return false;
    target_angles=&next_state->angles;
    return true;
}
//#####################################################################
// Function Golden_Section_Search
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR SEARCH_CONTROLLER<T_GRID>::
Golden_Section_Search(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,T a,T b,const T dx,T* F)
{
    T tau=(T).5*(sqrt((T)5)-1);
    T x_1=a+(1-tau)*(b-a);
    T x_2=a+tau*(b-a);
    T f_1=Evaluate_Force(face_velocities,dt,time,dx,x_1,&negative_gradient);
    T f_2=Evaluate_Force(face_velocities,dt,time,dx,x_2,&negative_gradient);
    while((b-a)>1e2*threshold){
        if((f_1>f_2&&minimize)||(!minimize&&f_2>f_1)){a=x_1;x_1=x_2;f_1=f_2;x_2=a+tau*(b-a);f_2=Evaluate_Force(face_velocities,dt,time,dx,x_2,&negative_gradient);}
        else{b=x_2;x_2=x_1;f_2=f_1;x_1=a+(1-tau)*(b-a);f_1=Evaluate_Force(face_velocities,dt,time,dx,x_1,&negative_gradient);}}
    if(F) *F=(T).5*(f_1+f_2);
    return (T).5*(b+a);
}
//#####################################################################
// Function Steepest_Descent
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR SEARCH_CONTROLLER<T_GRID>::
Steepest_Descent(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    // TODO(jontg): We could use a reduced-order Hessian approximation to speed convergence up here...
    int iterations=0;
    while(iterations<max_iterations){
        PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Steepest descent (search_controller) iterations=%f time=%f",iterations,time),0,0);
        iterations++;Compute_Gradient(face_velocities,dt,time);
        if(negative_gradient.Magnitude()<threshold) return 0;
        Normalize_Gradient(negative_gradient);
        T step=Golden_Section_Search(face_velocities,dt,time,0,1,dx);
        for(JOINT_ID i(1);i<=joint_mesh.Size();i++){if(!joint_mesh.Is_Active(i)) continue;
            JOINT<TV>& joint=*joint_mesh(i);
            bool control=false;for(int j=1;j<=T_SPIN::dimension;j++) if(joint.control_dof(j)) control=true;
            if(!joint.joint_function || !control) continue;
            T_SPIN angles;for(int j=1;j<=T_SPIN::dimension;j++) angles(j)=dx*negative_gradient((Value(i)-1)*T_SPIN::dimension+j);
            ROTATION<TV> target_angle=ROTATION<TV>::From_Euler_Angles(angles)*current_angle(i);
            T_SPIN target_angle_constrained=target_angle.Euler_Angles();joint.Constrain_Angles(target_angle_constrained);
            current_angle(i)=ROTATION<TV>::Spherical_Linear_Interpolation(joint.joint_function->Angle(),ROTATION<TV>::From_Euler_Angles(target_angle_constrained),step);}}
    return 0;
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Update_Position_Based_State(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    assert(!hypothetical_step);
    SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters=driver->example.solids_fluids_parameters;
    bool solids=Simulate_Solids() && (!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Solid_Node());
    bool use_precomputed_angle=true;
    
    ENVIRONMENTAL_STATE<T_GRID>* current_state=0;
    if(collecting_data || states_graph){
        current_state=new ENVIRONMENTAL_STATE<T_GRID>(incorporate_fluids,fluids_parameters?fluids_parameters->grid:0);

        if(incorporate_fluids){
            T_ARRAYS_SCALAR& p=fluids_parameters->incompressible->projection.p;
            current_state->pressures.Resize(p.Domain_Indices());
            for(CELL_ITERATOR iterator(*fluids_parameters->grid,1);iterator.Valid();iterator.Next()) current_state->pressures(iterator.Cell_Index())=p(iterator.Cell_Index());}
        else{
            ARRAY<TV> F(solid_body_collection.deformable_body_collection.particles.array_collection->Size());ARRAY<TWIST<TV> >rigid_F(solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size());
            ARRAY<TWIST<TV> > twist;twist.Resize(solid_body_collection.rigid_body_collection.rigid_body_particle.V.Size());
            for(int i=1;i<=twist.Size();i++) twist(i)=TWIST<TV>(solid_body_collection.rigid_body_collection.rigid_body_particle.V(i),solid_body_collection.rigid_body_collection.rigid_body_particle.angular_velocity(i));
            for(int i=1;i<=solid_body_collection.solids_forces.m;i++){
                solid_body_collection.solids_forces(i)->Add_Velocity_Independent_Forces(F,rigid_F,time);
                solid_body_collection.solids_forces(i)->Add_Velocity_Dependent_Forces(solid_body_collection.deformable_body_collection.particles.V,twist,F,rigid_F,time);}
            for(int j=1;j<=twist.Size();j++){solid_body_collection.rigid_body_collection.rigid_body_particle.V(j)=twist(j).linear;solid_body_collection.rigid_body_collection.rigid_body_particle.angular_velocity(j)=twist(j).angular;}
            for(int i=1;i<=F.Size();i++) current_state->external_force_dir+=F(i);
            for(int i=1;i<=rigid_F.Size();i++) current_state->external_force_dir+=rigid_F(i).linear;
            current_state->external_force_mag=current_state->external_force_dir.Magnitude();
            current_state->external_force_dir.Normalize();}
       
        current_state->angles.Resize(joint_mesh.Size());
        for(JOINT_ID i(1);i<=joint_mesh.Size();i++){if(!joint_mesh.Is_Active(i)) continue;
            JOINT<TV>& joint=*joint_mesh(i);if(!joint.joint_function) continue;
            current_state->angles(i)=joint.joint_function->Angle();}
        
        if(collecting_data){
            ARRAY<JOINT_ID> force_list,drag_list;
            for(int i=1;i<=strain_joints.m;i++) force_list.Append(strain_joints(i)); 
            for(JOINT_ID i(1);i<=joint_mesh.Size();i++){if(!joint_mesh.Is_Active(i)) continue;
                if(objective(i)==FORCE) force_list.Append_Unique(i);
                if(objective(i)==DRAG) drag_list.Append(i);}
            if(force_list.m>0) current_state->force_to_stay=Evaluate_Force_To_Stay_For_Joint(face_velocities,dt,time,force_list);
            if(drag_list.m>0) current_state->force_to_stay+=Evaluate_Drag_For_Joint(dt,time,drag_list);}}
    
    //Update the search direction
    if(time-last_search_time>=dt_per_search_step){
        last_search_time=time;
        if(!(states_graph && Compute_Gradient_From_Graph(current_state))){
            use_precomputed_angle=false;initialized=true;
            if(!dt_hyp) dt_hyp=dt_per_search_step;
            if(!dt_hyp) dt_hyp=dt;
            T target_time=driver->Time()+dt_hyp;
            if(!current_angle.Size()) current_angle.Resize(joint_mesh.Size());
            if(!dF_array_multipliers.Size()){dF_array_multipliers.Resize(T_SPIN::dimension*Value(joint_mesh.Size()));ARRAYS_COMPUTATIONS::Fill(dF_array_multipliers,PAIR<T,T>((T)0,(T)1));}

            negative_gradient.Resize(T_SPIN::dimension*Value(joint_mesh.Size()));
            if(solve_minimization) real_dx=Steepest_Descent(face_velocities,dt,target_time);
            else{
                Compute_Gradient(face_velocities,dt,target_time);
                Normalize_Gradient(negative_gradient);

                if(line_search) real_dx*=Golden_Section_Search(face_velocities,dt,target_time,0,1,real_dx);
            
                for(JOINT_ID i(1);i<=joint_mesh.Size();i++){if(!joint_mesh.Is_Active(i)) continue; for(int j=1;j<=T_SPIN::dimension;j++){
                    int index=(Value(i)-1)*T_SPIN::dimension+j;
                    if(negative_gradient(index)*dF_array_multipliers(index).x<(T)0) dF_array_multipliers(index)=PAIR<T,T>((T)0,max(dF_array_multipliers(index).y*(T).5,min_multiplier));
                    else if(negative_gradient(index)!=(T)0){
                        if(abs(dF_array_multipliers(index).x)>=2) dF_array_multipliers(index)=PAIR<T,T>((T).5*sign(negative_gradient(index)),min((T)1,2*dF_array_multipliers(index).y));
                        else dF_array_multipliers(index).x+=(T)sign(negative_gradient(index));}
                    else dF_array_multipliers(index).x=(T)0;}}}}}

    if(!negative_gradient.Size() || !solids) return;
    for(JOINT_ID i(1);i<=joint_mesh.Size();i++){if(!joint_mesh.Is_Active(i)) continue;
        JOINT<TV>& joint=*joint_mesh(i);
        bool control=false;for(int j=1;j<=T_SPIN::dimension;j++) if(joint.control_dof(j)) control=true;
        if(!joint.joint_function || !control) continue;
        ROTATION<TV> target_angle;
        if(use_precomputed_angle) target_angle=(*target_angles)(i);
        else{
            T_SPIN angles;for(int j=1;j<=T_SPIN::dimension;j++) angles(j)=real_dx*negative_gradient((Value(i)-1)*T_SPIN::dimension+j)*dF_array_multipliers((Value(i)-1)*T_SPIN::dimension+j).y;
            target_angle=ROTATION<TV>::From_Euler_Angles(angles)*current_angle(i);}
        T_SPIN target_angle_constrained=target_angle.Euler_Angles();joint.Constrain_Angles(target_angle_constrained);
        joint.joint_function->Set_Target_Angle(ROTATION<TV>::From_Euler_Angles(target_angle_constrained));
        joint.joint_function->active=true;joint.global_post_stabilization=true;}
}
//#####################################################################
// Function Preprocess_Drag_Substep
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Preprocess_Drag_Substep(const T time)
{
    PHYSBAM_ASSERT(drag_step);
    if(!Fluid_Node()) return;
    fluids_parameters->collision_bodies_affecting_fluid->Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE,time);
    fluids_parameters->collision_bodies_affecting_fluid->Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_SAVED_OLD_STATE);
    fluids_parameters->collision_bodies_affecting_fluid->Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,time);
    fluids_parameters->collision_bodies_affecting_fluid->Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE);
}
//#####################################################################
// Function Read
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Read(const STREAM_TYPE stream_type,const std::string& directory,const int frame) {
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/environmental_state",directory.c_str(),frame),states);
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(STRING_UTILITIES::string_sprintf("%s/%d/search_controller",directory.c_str(),frame));
    TYPED_ISTREAM typed_input=TYPED_ISTREAM(*input,stream_type);
    Read(typed_input);
    delete input;
}
//#####################################################################
// Function Read
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::    
Read(TYPED_ISTREAM& typed_input) {
    Read_Binary(typed_input,dt_per_search_step);
    Read_Binary(typed_input,last_search_time);
    Read_Binary(typed_input,old_time);
    Read_Binary(typed_input,negative_gradient);
    Read_Binary(typed_input,dF_array_multipliers);
    Read_Binary(typed_input,current_angle);
}
//#####################################################################
// Function Read
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Write(const STREAM_TYPE stream_type,const std::string& directory,const int frame) const {
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/environmental_state",directory.c_str(),frame),states);
    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(STRING_UTILITIES::string_sprintf("%s/%d/search_controller",directory.c_str(),frame));
    TYPED_OSTREAM typed_output=TYPED_OSTREAM(*output,stream_type);
    Write(typed_output);
    delete output;
}
//#####################################################################
// Function Write
//#####################################################################
template<class T_GRID> void SEARCH_CONTROLLER<T_GRID>::
Write(TYPED_OSTREAM& typed_output) const {
    Write_Binary(typed_output,dt_per_search_step);
    Write_Binary(typed_output,last_search_time);
    Write_Binary(typed_output,old_time);
    Write_Binary(typed_output,negative_gradient);
    Write_Binary(typed_output,dF_array_multipliers);
    Write_Binary(typed_output,current_angle);
}
//#####################################################################
template class SEARCH_CONTROLLER<GRID<VECTOR<float,2> > >;
template class SEARCH_CONTROLLER<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SEARCH_CONTROLLER<GRID<VECTOR<double,2> > >;
template class SEARCH_CONTROLLER<GRID<VECTOR<double,3> > >;
#endif

//#####################################################################
// Copyright 2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOUNDED_HORIZONTAL_PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/MULTIBODY_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_CONTACT_GRAPH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_SKIP_COLLISION_CHECK.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGIDS_COLLISION_CALLBACKS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/CONTACT_PRECONDITIONER.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/DISCRETE_CONTACT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/LEVELSET_CONTACT_PAIR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/PARTICLES_IN_IMPLICIT_OBJECT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/PARTICLES_IN_PROXIMITY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/PROJECTED_GAUSS_SEIDEL.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/SOLVE_CONTACT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/SPHERE_PLANE_CONTACT_PAIR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/SPHERE_SPHERE_CONTACT_PAIR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/SURFACE_CONTACT_PAIR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/BOX_BOX_CONTACT_PAIR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/BOX_PLANE_CONTACT_PAIR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Parallel_Computation/MPI_RIGIDS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_MASS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>

namespace PhysBAM
{
namespace SOLVE_CONTACT
{
template<class TV>
bool Solve_Projected_Gauss_Seidel(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,RIGID_TRIANGLE_COLLISIONS<TV>* triangle_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,ARRAY<VECTOR<int,2> >& pairs,HASHTABLE<VECTOR<int,2> >& pairs_processed_by_contact,typename TV::SCALAR desired_separation_distance,typename TV::SCALAR contact_proximity,typename TV::SCALAR dt,typename TV::SCALAR tolerance,int iteration_maximum,const bool thin_shells);
//#####################################################################
// Function Solve
//#####################################################################
template<class TV>
void Solve(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,
    RIGID_BODY_COLLISION_PARAMETERS<TV>& parameters,const bool correct_contact_energy,const bool use_saved_pairs,const typename TV::SCALAR dt,const typename TV::SCALAR time,
    MPI_RIGIDS<TV>* mpi_rigids,ARRAY<TWIST<TV> >& mpi_rigid_velocity_save,ARRAY<typename TV::SPIN>& mpi_rigid_angular_momentum_save)
{
    typedef typename TV::SCALAR T;

    LOG::SCOPE scope("rigid body contact kernel",rigid_body_collisions.parameters.threadid);

    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_body_cluster_bindings=rigid_body_collisions.rigid_body_cluster_bindings;
    ARRAY<ARRAY<VECTOR<int,2> > >& contact_pairs_for_level=use_saved_pairs?rigid_body_collisions.saved_contact_pairs_for_level:rigid_body_collisions.precomputed_contact_pairs_for_level;

    PHYSBAM_ASSERT(!parameters.use_projected_gauss_seidel || !mpi_rigids);
    if(!use_saved_pairs || !parameters.use_projected_gauss_seidel){
        LOG::SCOPE scope("guendelman-bridson-fedkiw contact",rigid_body_collisions.parameters.threadid);

        HASHTABLE<VECTOR<std::string,2>,typename ANALYTICS<TV>::UPDATE_ANALYTIC_CONTACT_PAIR_T> analytic_contact_registry;
        if(parameters.use_analytic_collisions){Register_Analytic_Contacts<TV>(analytic_contact_registry);}

        ARTICULATED_RIGID_BODY<TV>* articulated_rigid_body=&rigid_body_collisions.rigid_body_collection.articulated_rigid_body;
        rigid_body_collisions.skip_collision_check.Reset();bool need_another_iteration=true;int iteration=0;T epsilon_scale=1;
        while(need_another_iteration && ++iteration<=parameters.contact_iterations){
            if(mpi_rigids){
                mpi_rigids->Clear_Impulse_Accumulators(rigid_body_collection);
                mpi_rigid_velocity_save.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
                mpi_rigid_angular_momentum_save.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
                for(int p=1;p<=rigid_body_collection.rigid_body_particle.array_collection->Size();p++) {
                    mpi_rigid_velocity_save(p).linear=rigid_body_collection.rigid_body_particle.V(p);
                    mpi_rigid_velocity_save(p).angular=rigid_body_collection.rigid_body_particle.angular_velocity(p);
                    mpi_rigid_angular_momentum_save(p)=rigid_body_collection.rigid_body_particle.angular_momentum(p);}}

            need_another_iteration=false;
            if(!parameters.use_ccd && parameters.use_epsilon_scaling) epsilon_scale=(T)iteration/parameters.contact_iterations;
            for(int level=1;level<=rigid_body_collisions.contact_graph.Number_Of_Levels();level++){
                ARRAY<VECTOR<int,2> >& pairs=contact_pairs_for_level(level);
                bool need_another_level_iteration=true;int level_iteration=0;
                while(need_another_level_iteration && ++level_iteration<=rigid_body_collisions.contact_level_iterations){need_another_level_iteration=false;
                    if(!parameters.use_ccd && parameters.use_epsilon_scaling_for_level) epsilon_scale=(T)iteration*level_iteration/(parameters.contact_iterations*rigid_body_collisions.contact_level_iterations);
                    for(int i=1;i<=pairs.m;i++){int id_1=pairs(i)(1),id_2=pairs(i)(2);
                        if(rigid_body_collisions.skip_collision_check.Skip_Pair(id_1,id_2)) continue;
                        if(rigid_body_collisions.prune_contact_using_velocity){
                            TRIPLE<int,int,int>& pair_scale=rigid_body_collisions.pairs_scale.Get(pairs(i).Sorted());
                            if((rigid_body_collisions.contact_level_iterations-level_iteration)%pair_scale.y || (parameters.contact_iterations-iteration)%pair_scale.x) continue;
                            if(Update_Contact_Pair(rigid_body_collisions,collision_callbacks,analytic_contact_registry,id_1,id_2,correct_contact_energy,
                                    (int)(rigid_body_collisions.contact_pair_iterations/pair_scale.z),epsilon_scale,dt,time,false)){
                                if(!use_saved_pairs) rigid_body_collisions.saved_contact_pairs_for_level(level).Append_Unique(pairs(i)); // TODO: Potentially inefficient
                                need_another_level_iteration=true;need_another_iteration=true;}}
                        else{
                            bool mpi_one_ghost=mpi_rigids && (mpi_rigids->Is_Dynamic_Ghost_Body(rigid_body_collection.Rigid_Body(id_1)) || mpi_rigids->Is_Dynamic_Ghost_Body(rigid_body_collection.Rigid_Body(id_2)));
                            if(Update_Contact_Pair(rigid_body_collisions,collision_callbacks,analytic_contact_registry,id_1,id_2,correct_contact_energy,
                                    rigid_body_collisions.contact_pair_iterations,epsilon_scale,dt,time,mpi_one_ghost)){
                                if(!use_saved_pairs) rigid_body_collisions.saved_contact_pairs_for_level(level).Append_Unique(pairs(i));
                                need_another_level_iteration=true;need_another_iteration=true;}}}
                    if(articulated_rigid_body)
                        for(int i=1;i<=articulated_rigid_body->contact_level_iterations;i++) for(int j=1;j<=articulated_rigid_body->process_list(level).m;j++){
                            rigid_body_collisions.Apply_Prestabilization_To_Joint(dt,time,*articulated_rigid_body,articulated_rigid_body->process_list(level)(j),epsilon_scale);
                            need_another_level_iteration=need_another_iteration=true;}}}

            if(mpi_rigids){
                int need_another_iteration_int=(int)need_another_iteration;
                need_another_iteration=mpi_rigids->Reduce_Max(need_another_iteration_int)>0?true:false;
                mpi_rigids->Exchange_All_Impulses(rigid_body_collection,mpi_rigid_velocity_save,mpi_rigid_angular_momentum_save,rigid_body_collisions,false,dt,time);}}
        if(rigid_body_collisions.prune_stacks_from_contact) rigid_body_collisions.Apply_Stacking_Contact();
    
        if(parameters.use_ccd){
            T old_size=(T)(rigid_body_collection.rigid_body_particle.X.m);int max_iterations=20;
            for(int level=1;level<=rigid_body_collisions.contact_graph.Number_Of_Levels();level++){
                ARRAY<VECTOR<int,2> >& pairs=contact_pairs_for_level(level);
                for(int i=1;i<=pairs.m;i++){int id_1=pairs(i)(1),id_2=pairs(i)(2);int iterations=0;
                    while(!rigid_body_collisions.skip_collision_check.Skip_Pair(id_1,id_2) && Update_Contact_Pair(rigid_body_collisions,collision_callbacks,analytic_contact_registry,id_1,id_2,correct_contact_energy,
                        rigid_body_collisions.contact_pair_iterations,epsilon_scale,dt,time,(mpi_rigids && 
                        (mpi_rigids->Is_Dynamic_Ghost_Body(rigid_body_collection.Rigid_Body(id_1)) || mpi_rigids->Is_Dynamic_Ghost_Body(rigid_body_collection.Rigid_Body(id_2))))) && iterations<max_iterations){iterations++;}
                    if(iterations==max_iterations) Rigidify_Contact_Pair(rigid_body_collisions,collision_callbacks,id_1,id_2,dt,time,false);
                    else Rigidify_Kinematic_Contact_Pair(rigid_body_collisions,collision_callbacks,id_1,id_2,dt,time,false);}}
                    //else Rigidify_Contact_Pair(rigid_body_collisions,collision_callbacks,id_1,id_2,dt,time,false);}}
            ARRAY<int> parents;
            for(typename HASHTABLE<int,typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER*>::ITERATOR i(rigid_body_cluster_bindings.reverse_bindings);i.Valid();i.Next()) parents.Append(i.Key());
            Sort(parents);
            for(int i=parents.m;i>0;i--){rigid_body_cluster_bindings.Delete_Binding(parents(i));rigid_body_collection.rigid_body_particle.array_collection->Delete_Element(parents(i));}
            collision_callbacks.Restore_Size((int)old_size);
            rigid_body_collection.rigid_body_particle.array_collection->Resize((int)old_size);
            rigid_body_collisions.skip_collision_check.Resize(rigid_body_collection.rigid_body_particle.X.m);

            /*for(int level=1;level<=rigid_body_collisions.contact_graph.Number_Of_Levels();level++){
                ARRAY<VECTOR<int,2> >& pairs=contact_pairs_for_level(level);
                for(int i=1;i<=pairs.m;i++){int id_1=pairs(i)(1),id_2=pairs(i)(2);
                    if(Update_Contact_Pair(rigid_body_collisions,collision_callbacks,analytic_contact_registry,id_1,id_2,correct_contact_energy,
                        rigid_body_collisions.contact_pair_iterations,epsilon_scale,dt,time,(mpi_rigids && 
                        (mpi_rigids->Is_Dynamic_Ghost_Body(rigid_body_collection.Rigid_Body(id_1)) || mpi_rigids->Is_Dynamic_Ghost_Body(rigid_body_collection.Rigid_Body(id_2)))))) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Found penetrating pair 1=%d 2=%d",id_1,id_2));}}*/}}
    else{
        LOG::SCOPE scope("projected-gauss-seidel contact",rigid_body_collisions.parameters.threadid);

        if(!use_saved_pairs) rigid_body_collisions.saved_contact_pairs_for_level=contact_pairs_for_level;
        ARRAY<VECTOR<int,2> > pairs;
        for(int level=1;level<=rigid_body_collisions.contact_graph.Number_Of_Levels();level++)
            pairs.Append_Unique_Elements(contact_pairs_for_level(level));
        int iteration_maximum=parameters.contact_iterations*rigid_body_collisions.contact_level_iterations*rigid_body_collisions.contact_pair_iterations;
        
        Solve_Projected_Gauss_Seidel(rigid_body_collection,rigid_body_collisions.triangle_collisions,collision_callbacks,pairs,rigid_body_collisions.pairs_processed_by_contact,rigid_body_collisions.desired_separation_distance,parameters.contact_proximity,dt,parameters.projected_gauss_seidel_tolerance,iteration_maximum,parameters.use_ccd);

        for(int i=1;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++){
            if(rigid_body_collection.Is_Active(i)){
                RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(i);
                if(!body.Has_Infinite_Inertia()){
                    body.Update_Angular_Momentum();
                    collision_callbacks.Euler_Step_Position(i,dt,time);}}}}
}
//#####################################################################
// Function Rigidify_Contact_Pair
//#####################################################################
template<class TV> void
Rigidify_Contact_Pair(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,const int id_1,const int id_2,const typename TV::SCALAR dt,const typename TV::SCALAR time,const bool mpi_one_ghost)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=rigid_body_collisions.rigid_body_collection;
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_body_cluster_bindings=rigid_body_collisions.rigid_body_cluster_bindings;
    int parent_id_1=rigid_body_cluster_bindings.Get_Parent(rigid_body_collection.Rigid_Body(id_1)).particle_index;
    int parent_id_2=rigid_body_cluster_bindings.Get_Parent(rigid_body_collection.Rigid_Body(id_2)).particle_index;
    int parent=0;bool has_kinematic_body=false;
    collision_callbacks.Restore_Position(parent_id_1);collision_callbacks.Restore_Position(parent_id_2);
    rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Positions();
    if(parent_id_1!=id_1){
        typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER& cluster=*rigid_body_cluster_bindings.reverse_bindings.Get(parent_id_1);
        for(RIGID_CLUSTER_CONSTITUENT_ID i(1);i<=cluster.children.Size();i++){int child=cluster.children(i);
            if(cluster.kinematic_child(i)) collision_callbacks.Restore_Position(child);}}
    if(parent_id_2!=id_2){
        typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER& cluster=*rigid_body_cluster_bindings.reverse_bindings.Get(parent_id_2);
        for(RIGID_CLUSTER_CONSTITUENT_ID i(1);i<=cluster.children.Size();i++){int child=cluster.children(i);
            if(cluster.kinematic_child(i)) collision_callbacks.Restore_Position(child);}}
    rigid_body_cluster_bindings.clamp_kinematic_positions=false;
    if(parent_id_1==id_1 && parent_id_2==id_2){ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> ids;ids.Append(id_1);ids.Append(id_2);parent=rigid_body_collisions.rigid_body_cluster_bindings.Add_Binding(ids);}
    else if(parent_id_1==id_1){parent=parent_id_2;rigid_body_collisions.rigid_body_cluster_bindings.Append_To_Binding(parent_id_2,id_1);}
    else if(parent_id_2==id_2){parent=parent_id_1;rigid_body_collisions.rigid_body_cluster_bindings.Append_To_Binding(parent_id_1,id_2);}
    else if(parent_id_1!=parent_id_2) {parent=parent_id_1;rigid_body_collisions.rigid_body_cluster_bindings.Merge_Bindings(parent_id_1,parent_id_2);}
    else{parent=parent_id_1;
        typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER& cluster=*rigid_body_cluster_bindings.reverse_bindings.Get(parent);
        for(RIGID_CLUSTER_CONSTITUENT_ID i(1);i<=cluster.children.Size();i++){int child=cluster.children(i);
            if(child==id_1 || child==id_2){
                RIGID_BODY<TV>& child_body=rigid_body_collection.Rigid_Body(child);cluster.kinematic_child(i)=false;
                cluster.child_to_parent(i)=rigid_body_collection.Rigid_Body(parent).Frame().Inverse()*child_body.Frame();}
            if(cluster.kinematic_child(i)) has_kinematic_body=true;}}
    rigid_body_collisions.skip_collision_check.Resize(rigid_body_collection.rigid_body_particle.X.m);
    collision_callbacks.Save_Position(parent);collision_callbacks.Euler_Step_Position(parent,dt,time);
    if(has_kinematic_body){
        typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER& cluster=*rigid_body_cluster_bindings.reverse_bindings.Get(parent);
        for(RIGID_CLUSTER_CONSTITUENT_ID i(1);i<=cluster.children.Size();i++){
            int child=cluster.children(i);RIGID_BODY<TV>& child_body=rigid_body_collection.Rigid_Body(child);
            if(cluster.kinematic_child(i)){
                collision_callbacks.Euler_Step_Position(child,dt,time);
                cluster.child_to_parent(i)=rigid_body_collection.Rigid_Body(parent).Frame().Inverse()*child_body.Frame();}}}
    rigid_body_cluster_bindings.clamp_kinematic_positions=true;
    rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Positions();
}
//#####################################################################
// Function Rigidify_Contact_Pair
//#####################################################################
template<class TV> void
Rigidify_Kinematic_Contact_Pair(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,const int id_1,const int id_2,const typename TV::SCALAR dt,const typename TV::SCALAR time,const bool mpi_one_ghost)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=rigid_body_collisions.rigid_body_collection;
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_body_cluster_bindings=rigid_body_collisions.rigid_body_cluster_bindings;
    int parent_id_1=rigid_body_cluster_bindings.Get_Parent(rigid_body_collection.Rigid_Body(id_1)).particle_index;
    int parent_id_2=rigid_body_cluster_bindings.Get_Parent(rigid_body_collection.Rigid_Body(id_2)).particle_index;
    int parent=0;bool has_kinematic_body=true;
    collision_callbacks.Restore_Position(parent_id_1);collision_callbacks.Restore_Position(parent_id_2);
    rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Positions();
    if(parent_id_1!=id_1){
        typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER& cluster=*rigid_body_cluster_bindings.reverse_bindings.Get(parent_id_1);
        for(RIGID_CLUSTER_CONSTITUENT_ID i(1);i<=cluster.children.Size();i++){int child=cluster.children(i);
            if(cluster.kinematic_child(i)) collision_callbacks.Restore_Position(child);}}
    if(parent_id_2!=id_2){
        typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER& cluster=*rigid_body_cluster_bindings.reverse_bindings.Get(parent_id_2);
        for(RIGID_CLUSTER_CONSTITUENT_ID i(1);i<=cluster.children.Size();i++){int child=cluster.children(i);
            if(cluster.kinematic_child(i)) collision_callbacks.Restore_Position(child);}}
    rigid_body_cluster_bindings.clamp_kinematic_positions=false;
    if(parent_id_1==id_1 && parent_id_2==id_2){ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> ids;ids.Append(id_1);ids.Append(id_2);parent=rigid_body_collisions.rigid_body_cluster_bindings.Add_Binding(ids,true);}
    else if(parent_id_1==id_1){parent=parent_id_2;rigid_body_collisions.rigid_body_cluster_bindings.Append_To_Binding(parent_id_2,id_1,true);}
    else if(parent_id_2==id_2){parent=parent_id_1;rigid_body_collisions.rigid_body_cluster_bindings.Append_To_Binding(parent_id_1,id_2,true);}
    else if(parent_id_1!=parent_id_2){parent=parent_id_1;rigid_body_collisions.rigid_body_cluster_bindings.Merge_Bindings(parent_id_1,parent_id_2,true);}
    else parent=parent_id_1;
    rigid_body_collisions.skip_collision_check.Resize(rigid_body_collection.rigid_body_particle.X.m);
    collision_callbacks.Save_Position(parent);collision_callbacks.Euler_Step_Position(parent,dt,time);
    if(has_kinematic_body){
        typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER& cluster=*rigid_body_cluster_bindings.reverse_bindings.Get(parent);
        for(RIGID_CLUSTER_CONSTITUENT_ID i(1);i<=cluster.children.Size();i++){
            int child=cluster.children(i);RIGID_BODY<TV>& child_body=rigid_body_collection.Rigid_Body(child);
            if(cluster.kinematic_child(i)){
                collision_callbacks.Euler_Step_Position(child,dt,time);
                cluster.child_to_parent(i)=rigid_body_collection.Rigid_Body(parent).Frame().Inverse()*child_body.Frame();}}}
    rigid_body_cluster_bindings.clamp_kinematic_positions=true;
    rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Positions();
}
//#####################################################################
// Function Update_Analytic_Multibody_Contact
//#####################################################################
template<class TV> bool 
Update_Analytic_Multibody_Contact(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,
    HASHTABLE<VECTOR<std::string,2>,typename ANALYTICS<TV>::UPDATE_ANALYTIC_CONTACT_PAIR_T>& analytic_contact_registry,
    const int id_1,const int id_2,MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>& multibody,IMPLICIT_OBJECT<TV>& levelset,const bool correct_contact_energy,const int max_iterations,
    const typename TV::SCALAR epsilon_scale,const typename TV::SCALAR dt,const typename TV::SCALAR time,const bool mpi_one_ghost)
{
    bool return_val=false;
    IMPLICIT_OBJECT<TV>* levelset_base=&levelset;
    if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* levelset_transformed=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(levelset_base)) 
        levelset_base=levelset_transformed->object_space_implicit_object;
    for(int i=1;i<=multibody.levelsets->m;i++){
        VECTOR<std::string,2> key(typeid(*dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >&>(*(*multibody.levelsets)(i)).object_space_implicit_object).name(),typeid(*levelset_base).name());
        if(typename ANALYTICS<TV>::UPDATE_ANALYTIC_CONTACT_PAIR_T* contact_function=analytic_contact_registry.Get_Pointer(key.Sorted()))
            return_val|=(*contact_function)(rigid_body_collisions,collision_callbacks,id_1,id_2,(*multibody.levelsets)(i),&levelset,correct_contact_energy,max_iterations,epsilon_scale,dt,time,mpi_one_ghost);
        else PHYSBAM_FATAL_ERROR();}
    return return_val;
}
//#####################################################################
// Function Update_Analytic_Multibody_Contact
//#####################################################################
template<class TV> bool 
Update_Analytic_Multibody_Contact(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,
    HASHTABLE<VECTOR<std::string,2>,typename ANALYTICS<TV>::UPDATE_ANALYTIC_CONTACT_PAIR_T>& analytic_contact_registry,
    RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const bool correct_contact_energy,const int max_iterations,
    const typename TV::SCALAR epsilon_scale,const typename TV::SCALAR dt,const typename TV::SCALAR time,const bool mpi_one_ghost)
{
    bool return_val=false;
    MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>* multibody1=dynamic_cast<MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>*>(body1.implicit_object->object_space_implicit_object);
    MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>* multibody2=dynamic_cast<MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>*>(body2.implicit_object->object_space_implicit_object);
    if(multibody1 && multibody2) for(int i=1;i<=multibody2->levelsets->m;i++) 
        return_val|=Update_Analytic_Multibody_Contact(rigid_body_collisions,collision_callbacks,analytic_contact_registry,body1.particle_index,body2.particle_index,*multibody1,*(*multibody2->levelsets)(i),correct_contact_energy,max_iterations,epsilon_scale,dt,time,mpi_one_ghost);
    else if(multibody1) 
        return_val=Update_Analytic_Multibody_Contact(rigid_body_collisions,collision_callbacks,analytic_contact_registry,body1.particle_index,body2.particle_index,*multibody1,*body2.implicit_object->object_space_implicit_object,correct_contact_energy,max_iterations,epsilon_scale,dt,time,mpi_one_ghost);
    else 
        return_val=Update_Analytic_Multibody_Contact(rigid_body_collisions,collision_callbacks,analytic_contact_registry,body1.particle_index,body2.particle_index,*multibody2,*body1.implicit_object->object_space_implicit_object,correct_contact_energy,max_iterations,epsilon_scale,dt,time,mpi_one_ghost);
    return return_val;
}
//#####################################################################
// Function Update_Contact_Pair
//#####################################################################
template<class TV>
bool Update_Contact_Pair(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,HASHTABLE<VECTOR<std::string,2>,
    typename ANALYTICS<TV>::UPDATE_ANALYTIC_CONTACT_PAIR_T>& analytic_contact_registry,const int id_1,const int id_2,const bool correct_contact_energy,const int max_iterations,
    const typename TV::SCALAR epsilon_scale,const typename TV::SCALAR dt,const typename TV::SCALAR time,const bool mpi_one_ghost)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=rigid_body_collisions.rigid_body_collection;
    RIGID_BODY_COLLISION_PARAMETERS<TV>& parameters=rigid_body_collisions.parameters;

    RIGID_BODY<TV>& body1=rigid_body_collection.Rigid_Body(id_1),&body2=rigid_body_collection.Rigid_Body(id_2);
    VECTOR<std::string,2> key(typeid(*body1.implicit_object->object_space_implicit_object).name(),typeid(*body2.implicit_object->object_space_implicit_object).name());
    if(key.x==typeid(MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>).name()||key.y==typeid(MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>).name()){
        if(!body1.simplicial_object && !body2.simplicial_object) 
            return Update_Analytic_Multibody_Contact(rigid_body_collisions,collision_callbacks,analytic_contact_registry,body1,body2,correct_contact_energy,max_iterations,epsilon_scale,dt,time,mpi_one_ghost);}
    if(typename ANALYTICS<TV>::UPDATE_ANALYTIC_CONTACT_PAIR_T* contact_function=analytic_contact_registry.Get_Pointer(key.Sorted()))
        return (*contact_function)(rigid_body_collisions,collision_callbacks,id_1,id_2,body1.implicit_object->object_space_implicit_object,body2.implicit_object->object_space_implicit_object,correct_contact_energy,max_iterations,epsilon_scale,dt,time,mpi_one_ghost);
    if(parameters.use_ccd) return CONTACT_PAIRS::Update_Surface_Contact_Pair(rigid_body_collisions,collision_callbacks,id_1,id_2,correct_contact_energy,max_iterations,epsilon_scale,dt,time,
        parameters.rigid_collisions_use_triangle_hierarchy,parameters.rigid_collisions_use_edge_intersection,parameters.rigid_collisions_use_triangle_hierarchy_center_phi_test,mpi_one_ghost);
    return CONTACT_PAIRS::Update_Levelset_Contact_Pair(rigid_body_collisions,collision_callbacks,id_1,id_2,correct_contact_energy,max_iterations,epsilon_scale,dt,time,
        parameters.rigid_collisions_use_triangle_hierarchy,parameters.rigid_collisions_use_edge_intersection,parameters.rigid_collisions_use_triangle_hierarchy_center_phi_test,mpi_one_ghost);
}
//#####################################################################
// Function Update_Contact_Pair_Helper
//#####################################################################
template<class TV>
void Update_Contact_Pair_Helper(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,const int id_1,const int id_2,
    const typename TV::SCALAR dt,const typename TV::SCALAR time,const typename TV::SCALAR epsilon_scale,const TV& collision_location,const TV& collision_normal,
    const TV& collision_relative_velocity,const bool correct_contact_energy,const bool rolling_friction,const bool mpi_one_ghost)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=rigid_body_collisions.rigid_body_collection;
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_body_cluster_bindings=rigid_body_collisions.rigid_body_cluster_bindings;

    RIGID_BODY<TV>& body1=rigid_body_collection.Rigid_Body(id_1);
    RIGID_BODY<TV>& body2=rigid_body_collection.Rigid_Body(id_2);
    int parent_id_1=rigid_body_cluster_bindings.Get_Parent(body1).particle_index;
    int parent_id_2=rigid_body_cluster_bindings.Get_Parent(body2).particle_index;
    RIGID_BODY<TV>& parent_body_1=rigid_body_collection.Rigid_Body(parent_id_1);
    RIGID_BODY<TV>& parent_body_2=rigid_body_collection.Rigid_Body(parent_id_2);
    rigid_body_collisions.rigid_body_particle_intersections.Set(Tuple(body1.particle_index,body2.particle_index,collision_location));
    TWIST<TV> saved_v1=rigid_body_collection.Rigid_Body(parent_id_1).Twist(),saved_v2=rigid_body_collection.Rigid_Body(parent_id_2).Twist();
    RIGID_BODY<TV>::Apply_Collision_Impulse(parent_body_1,parent_body_2,body1.Rotation(),body2.Rotation(),collision_location,collision_normal,collision_relative_velocity,
        -1+epsilon_scale,RIGID_BODY<TV>::Coefficient_Of_Friction(parent_body_1,parent_body_2),false,rolling_friction,correct_contact_energy,mpi_one_ghost);
    collision_callbacks.Restore_Position(parent_id_1);collision_callbacks.Restore_Position(parent_id_2); // fix saved values & re-evolve bodies
    Euler_Step_Position(rigid_body_collisions,collision_callbacks,parent_id_1,dt,time);Euler_Step_Position(rigid_body_collisions,collision_callbacks,parent_id_2,dt,time);
    if(parent_id_2!=id_2 || parent_id_1!=id_1){
        rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Positions();
        rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Velocities();}
    if(parent_id_1!=id_1){
        typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER& cluster=*rigid_body_cluster_bindings.reverse_bindings.Get(parent_id_1);
        rigid_body_collection.Rigid_Body(parent_id_1).V()-=saved_v1.linear;
        for(RIGID_CLUSTER_CONSTITUENT_ID i(1);i<=cluster.children.Size();i++){int child=cluster.children(i);
            if(cluster.kinematic_child(i)){
                rigid_body_collection.rigid_body_particle.angular_velocity(child)+=rigid_body_collection.Rigid_Body(parent_id_1).Angular_Velocity()-saved_v1.angular;
                rigid_body_collection.rigid_body_particle.V(child)+=rigid_body_collection.Rigid_Body(parent_id_1).Pointwise_Object_Velocity(rigid_body_collection.rigid_body_particle.X(child));}}
        rigid_body_collection.Rigid_Body(parent_id_1).V()+=saved_v1.linear;}
    if(parent_id_2!=id_2){
        typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER& cluster=*rigid_body_cluster_bindings.reverse_bindings.Get(parent_id_2);
        rigid_body_collection.Rigid_Body(parent_id_2).V()-=saved_v2.linear;
        for(RIGID_CLUSTER_CONSTITUENT_ID i(1);i<=cluster.children.Size();i++){int child=cluster.children(i);
            if(cluster.kinematic_child(i)){
                rigid_body_collection.rigid_body_particle.angular_velocity(child)+=rigid_body_collection.Rigid_Body(parent_id_2).Angular_Velocity()-saved_v2.angular;
                rigid_body_collection.rigid_body_particle.V(child)+=rigid_body_collection.Rigid_Body(parent_id_2).Pointwise_Object_Velocity(rigid_body_collection.rigid_body_particle.X(child));}}
       rigid_body_collection.Rigid_Body(parent_id_2).V()+=saved_v2.linear;}
}
//#####################################################################
// Function Euler_Step_Position
//#####################################################################
template<class TV>
void Euler_Step_Position(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,const int id,const typename TV::SCALAR dt,const typename TV::SCALAR time)
{
    collision_callbacks.Euler_Step_Position(id,dt,time);
    rigid_body_collisions.rigid_body_collection.Rigid_Body(id).Update_Bounding_Box();rigid_body_collisions.skip_collision_check.Set_Last_Moved(id);
}
//#####################################################################
// Function Register_Analytic_Collisions
//#####################################################################
template<class TV>
void Register_Analytic_Contacts(HASHTABLE<VECTOR<std::string,2>,typename ANALYTICS<TV>::UPDATE_ANALYTIC_CONTACT_PAIR_T>& analytic_contact_registry)
{
    const char* sphere=typeid(ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >).name();
    const char* plane=typeid(ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<TV> >).name();
    const char* box=typeid(ANALYTIC_IMPLICIT_OBJECT<BOX<TV> >).name();
    analytic_contact_registry.Set(VECTOR<std::string,2>(sphere,sphere),CONTACT_PAIRS::Update_Sphere_Sphere_Contact_Pair<TV>);
    analytic_contact_registry.Set(VECTOR<std::string,2>(box,box),CONTACT_PAIRS::Update_Box_Box_Contact_Pair<TV>);
    analytic_contact_registry.Set(VECTOR<std::string,2>(sphere,plane).Sorted(),CONTACT_PAIRS::Update_Sphere_Plane_Contact_Pair<TV>);
    analytic_contact_registry.Set(VECTOR<std::string,2>(box,plane).Sorted(),CONTACT_PAIRS::Update_Box_Plane_Contact_Pair<TV>);
}
//#####################################################################
// Function Solve_Projected_Gauss_Seidel
//#####################################################################
template<class TV>
bool Solve_Projected_Gauss_Seidel(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,RIGID_TRIANGLE_COLLISIONS<TV>* triangle_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,ARRAY<VECTOR<int,2> >& pairs,HASHTABLE<VECTOR<int,2> >& pairs_processed_by_contact,typename TV::SCALAR desired_separation_distance,typename TV::SCALAR contact_proximity,typename TV::SCALAR dt,typename TV::SCALAR tolerance,int iteration_maximum,const bool thin_shells)
{
    typedef typename TV::SCALAR T;

    int n_bodies=rigid_body_collection.rigid_body_particle.array_collection->Size();
    ARRAY<CONTACT<TV> > contacts;
    if(thin_shells) DISCRETE_CONTACT::Compute_Contacts(rigid_body_collection,triangle_collisions,contacts,(T)2e-2,(T)1e-2,dt);   
    else Get_Contact_Points(rigid_body_collection,collision_callbacks,pairs,contacts,contact_proximity,dt,false,false);
    ARRAY<TWIST<TV> > velocities(n_bodies);
    ARRAY<bool> has_infinite_inertia(n_bodies);

    for(int i=1;i<=n_bodies;i++) if(rigid_body_collection.Is_Active(i)){
        velocities(i)=rigid_body_collection.Rigid_Body(i).Twist();
        has_infinite_inertia(i)=rigid_body_collection.Rigid_Body(i).Has_Infinite_Inertia();}

    int n_contacts=contacts.m;
    ARRAY<T> lambda_normal(n_contacts);
    ARRAY<VECTOR<T,TV::dimension-1> > lambda_tangent(n_contacts);
    PROJECTED_GAUSS_SEIDEL::Solve(velocities,has_infinite_inertia,contacts,lambda_normal,lambda_tangent,tolerance,iteration_maximum,true);
    for(int i=1;i<=n_bodies;i++) if(rigid_body_collection.Is_Active(i) && !has_infinite_inertia(i)){
        RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(i);
        body.V()=velocities(i).linear;
        body.Angular_Velocity()=velocities(i).angular;
        body.Update_Angular_Momentum();}
    return true;
}
//#####################################################################
// Function Get_Contact_Points
//#####################################################################
template<class TV>
void Get_Contact_Points(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,ARRAY<VECTOR<int,2> >& pairs,ARRAY<CONTACT<TV> >& contacts,typename TV::SCALAR contact_proximity,typename TV::SCALAR dt,const bool stagger_points,const bool use_old_states)
{
    if(use_old_states){
        LOG::SCOPE scope_tn("restoring tn states");
        for(int i=1;i<=rigid_body_collection.simulated_rigid_body_particles.m;i++)
            collision_callbacks.Swap_State(rigid_body_collection.simulated_rigid_body_particles(i));
        scope_tn.Pop();}
    
    LOG::SCOPE scope_contacts("creating contacts");
    for(int i=1;i<=pairs.m;i++){
        int id_1=pairs(i)(1),id_2=pairs(i)(2);
        ARRAY<TV> locations,normals;
        ARRAY<typename TV::SCALAR> distances;
        RIGID_BODY<TV>& body_1=rigid_body_collection.Rigid_Body(id_1);
        RIGID_BODY<TV>& body_2=rigid_body_collection.Rigid_Body(id_2);
        PARTICLES_IN_PROXIMITY::All_Particles_In_Proximity(body_1,body_2,locations,normals,distances,contact_proximity,stagger_points);
        for(int j=1;j<=locations.m;j++) contacts.Append(CONTACT<TV>(body_1,body_2,locations(j),normals(j),distances(j),dt));}
    scope_contacts.Pop();

    if(use_old_states){
        LOG::SCOPE scope_tn1("restoring tn+1 states");
        for(int i=1;i<=rigid_body_collection.simulated_rigid_body_particles.m;i++)
            collision_callbacks.Swap_State(rigid_body_collection.simulated_rigid_body_particles(i));
        scope_tn1.Pop();}
}
//#####################################################################
// Function Push_Out
//#####################################################################
template<class TV>
void Push_Out(RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,RIGID_BODY_COLLISION_PARAMETERS<TV>& parameters,HASHTABLE<VECTOR<int,2> >& pairs_processed_by_contact)
{
    LOG::SCOPE scope("Projected Gauss Seidel Push_Out");

    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    typedef typename RIGID_BODY_POLICY<TV>::WORLD_SPACE_INERTIA_TENSOR T_WORLD_SPACE_INERTIA_TENSOR;

    ARRAY<VECTOR<int,2> > pairs;
    pairs_processed_by_contact.Get_Keys(pairs);

    ARRAY<CONTACT<TV> > contacts;
    Get_Contact_Points(rigid_body_collection,collision_callbacks,pairs,contacts,0,1,false,false);
    for(int i=1;i<contacts.m;i++)
        contacts(i).coefficient_of_friction=0;

    //store linear velocity/update rotational momentum
    //set velocity to 0
    ARRAY<TV> linear_velocities(rigid_body_collection.simulated_rigid_body_particles.m);
    for(int i=1;i<=rigid_body_collection.simulated_rigid_body_particles.m;i++)
    {
        RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(rigid_body_collection.simulated_rigid_body_particles(i));
        linear_velocities(i)=body.V();
        body.Update_Angular_Momentum();
        body.V()=TV();
        body.Angular_Velocity()=T_SPIN();
    }

    T tolerance=(T)1e-3;
    int iteration_maximum=0;

    PROJECTED_GAUSS_SEIDEL::Solve(rigid_body_collection,contacts,tolerance,iteration_maximum);

    for(int i=1;i<=rigid_body_collection.simulated_rigid_body_particles.m;i++)
        if(!rigid_body_collection.rigid_body_particle.kinematic(rigid_body_collection.simulated_rigid_body_particles(i)))
            collision_callbacks.Euler_Step_Position(i,1,0); //should pass the actual time.

    //restore linear velocity/update rotational velocity
    for(int i=1;i<=rigid_body_collection.simulated_rigid_body_particles.m;i++)
    {
        RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(rigid_body_collection.simulated_rigid_body_particles(i));
        body.V()=linear_velocities(i);
        body.Update_Angular_Velocity();
    }
}
//#####################################################################

#define INSTANTIATION_HELPER(T,d) \
    template void Solve<VECTOR<T,d> >(RIGID_BODY_COLLISIONS<VECTOR<T,d> >& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks, \
        RIGID_BODY_COLLECTION<VECTOR<T,d> >& rigid_body_collection,RIGID_BODY_COLLISION_PARAMETERS<VECTOR<T,d> >& parameters,const bool correct_contact_energy,const bool use_saved_pairs, \
        const T dt,const T time,MPI_RIGIDS<VECTOR<T,d> >* mpi_rigids,ARRAY<TWIST<VECTOR<T,d> > >& mpi_rigid_velocity_save,ARRAY<VECTOR<T,d>::SPIN>& mpi_rigid_angular_momentum_save); \
    template bool Update_Analytic_Multibody_Contact<VECTOR<T,d> >(RIGID_BODY_COLLISIONS<VECTOR<T,d> >& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks, \
        HASHTABLE<VECTOR<std::string,2>,ANALYTICS<VECTOR<T,d> >::UPDATE_ANALYTIC_CONTACT_PAIR_T>& analytic_contact_registry, \
        const int id_1,const int id_2,MULTIBODY_LEVELSET_IMPLICIT_OBJECT<VECTOR<T,d> >& multibody,IMPLICIT_OBJECT<VECTOR<T,d> >& levelset,const bool correct_contact_energy, \
        const int max_iterations,const T epsilon_scale,const T dt,const T time,const bool mpi_one_ghost); \
    template bool Update_Analytic_Multibody_Contact<VECTOR<T,d> >(RIGID_BODY_COLLISIONS<VECTOR<T,d> >& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks, \
        HASHTABLE<VECTOR<std::string,2>,ANALYTICS<VECTOR<T,d> >::UPDATE_ANALYTIC_CONTACT_PAIR_T>& analytic_contact_registry, \
        RIGID_BODY<VECTOR<T,d> >& body1,RIGID_BODY<VECTOR<T,d> >& body2,const bool correct_contact_energy,const int max_iterations, \
        const T epsilon_scale,const T dt,const T time,const bool mpi_one_ghost); \
    template bool Update_Contact_Pair<VECTOR<T,d> >(RIGID_BODY_COLLISIONS<VECTOR<T,d> >& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks, \
        HASHTABLE<VECTOR<std::string,2>,ANALYTICS<VECTOR<T,d> >::UPDATE_ANALYTIC_CONTACT_PAIR_T>& analytic_contact_registry,const int id_1,const int id_2,const bool correct_contact_energy, \
        const int max_iterations,const T epsilon_scale,const T dt,const T time,const bool mpi_one_ghost); \
    template void Update_Contact_Pair_Helper<VECTOR<T,d> >(RIGID_BODY_COLLISIONS<VECTOR<T,d> >& rigid_body_collection,RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks, \
        const int id_1,const int id_2,const T dt,const T time,const T epsilon_scale,const VECTOR<T,d> & collision_location,const VECTOR<T,d> & collision_normal, \
        const VECTOR<T,d> & collision_relative_velocity,const bool correct_contact_energy,const bool rolling_friction,const bool mpi_one_ghost); \
    template void Euler_Step_Position<VECTOR<T,d> >(RIGID_BODY_COLLISIONS<VECTOR<T,d> >& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks, \
        const int id,const T dt,const T time); \
    template void Register_Analytic_Contacts<VECTOR<T,d> >(HASHTABLE<VECTOR<std::string,2>,ANALYTICS<VECTOR<T,d> >::UPDATE_ANALYTIC_CONTACT_PAIR_T>& analytic_contact_registry); \
    template void Push_Out(RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks,RIGID_BODY_COLLECTION<VECTOR<T,d> >& rigid_body_collection,RIGID_BODY_COLLISION_PARAMETERS<VECTOR<T,d> >& parameters,HASHTABLE<VECTOR<int,2> >& pairs_processed_by_contact); \
    template bool Solve_Projected_Gauss_Seidel<VECTOR<T,d> >(RIGID_BODY_COLLECTION<VECTOR<T,d> >& rigid_body_collection,RIGID_TRIANGLE_COLLISIONS<VECTOR<T,d> >* triangle_collisions,RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks, \
        ARRAY<VECTOR<int,2> >& pairs,HASHTABLE<VECTOR<int,2> >& pairs_processed_by_contact,T desired_separation_distance,T contact_proximity,T dt,T tolerance,int iteration_maximum,const bool thin_shells);

INSTANTIATION_HELPER(float,1)
INSTANTIATION_HELPER(float,2)
INSTANTIATION_HELPER(float,3)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1)
INSTANTIATION_HELPER(double,2)
INSTANTIATION_HELPER(double,3)
#endif
    }
}

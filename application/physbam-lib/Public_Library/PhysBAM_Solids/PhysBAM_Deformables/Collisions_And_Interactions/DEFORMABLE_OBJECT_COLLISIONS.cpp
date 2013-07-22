//#####################################################################
// Copyright 2006-2008, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_OBJECT_COLLISIONS
//#####################################################################
#include <PhysBAM_Tools/Arrays/CONSTANT_ARRAY.h>
#include <PhysBAM_Tools/Arrays/PROJECTED_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Utilities/Find_Type.h>
#include <PhysBAM_Geometry/Collision_Detection/COLLISION_GEOMETRY_SPATIAL_PARTITION.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_PARTICLE_STATE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COLLISION_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TETRAHEDRON_COLLISION_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_MATERIAL_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLE_OBJECT_COLLISIONS<TV>::
DEFORMABLE_OBJECT_COLLISIONS(PARTICLES<TV>& particles,DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection,ARRAY<STRUCTURE<TV>*>& deformable_object_structures,
    COLLISION_GEOMETRY_COLLECTION<TV>& collision_body_list)
    :particles(particles),deformable_body_collection(deformable_body_collection),deformable_object_structures(deformable_object_structures),collision_body_list(collision_body_list),
    collision_tolerance((T)1e-6),use_spatial_partition(true),disable_multiple_levelset_collisions(true),use_protectors(false),maximum_levelset_collision_projection_velocity(FLT_MAX),
    protection_thickness((T)1e-3),ignore_priorities(false),collisions_on(false),collision_tolerances(0),thickness_table(0),friction_table(0),use_structure_collide_collision_body(false)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DEFORMABLE_OBJECT_COLLISIONS<TV>::
~DEFORMABLE_OBJECT_COLLISIONS()
{}
//#####################################################################
// Function Initialize_Object_Collisions
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Initialize_Object_Collisions(const bool collide_with_interior,const T collision_tolerance_input,
    const bool use_spatial_partition_for_levelset_collisions,const bool disable_multiple_levelset_collisions_input,const T maximum_levelset_collision_projection_velocity_input)
{
    collision_tolerance=collision_tolerance_input;use_spatial_partition=use_spatial_partition_for_levelset_collisions;
    disable_multiple_levelset_collisions=disable_multiple_levelset_collisions_input;maximum_levelset_collision_projection_velocity=maximum_levelset_collision_projection_velocity_input;
    check_collision=CONSTANT_ARRAY<bool>(particles.array_collection->Size(),false);particle_states.Resize(particles.array_collection->Size());particle_to_collision_body_id.Resize(particles.array_collection->Size());
    Reset_Object_Collisions(); // in case collisions already exist
    for(int c=1;c<=collision_structures.m;c++){
        if(TRIANGULATED_AREA<T>* triangulated_area=dynamic_cast<TRIANGULATED_AREA<T>*>(collision_structures(c)))
            Add_Collision_Mesh(triangulated_area->mesh,collide_with_interior);
        else if(TRIANGULATED_SURFACE<T>* triangulated_surface=dynamic_cast<TRIANGULATED_SURFACE<T>*>(collision_structures(c)))
            Add_Collision_Mesh(triangulated_surface->mesh,true);
        else if(SEGMENTED_CURVE<TV>* segmented_curve=dynamic_cast<SEGMENTED_CURVE<TV>*>(collision_structures(c))){
            SEGMENT_MESH& mesh=segmented_curve->mesh;for(int s=1;s<=mesh.elements.m;s++){
                INDIRECT_ARRAY<ARRAY<bool>,SEGMENT_MESH::ELEMENT_TYPE&> subset=check_collision.Subset(mesh.elements(s));
                ARRAYS_COMPUTATIONS::Fill(subset,true);}}
        else if(TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(collision_structures(c)))
            Add_Collision_Mesh(tetrahedralized_volume->mesh,collide_with_interior);
        else if(HEXAHEDRALIZED_VOLUME<T>* hexahedralized_volume=dynamic_cast<HEXAHEDRALIZED_VOLUME<T>*>(collision_structures(c)))
            Add_Collision_Mesh(hexahedralized_volume->mesh,collide_with_interior);
        else if(EMBEDDING<TV>* embedding=dynamic_cast<EMBEDDING<TV>*>(collision_structures(c)))
            if(collide_with_interior && !dynamic_cast<EMBEDDED_MATERIAL_SURFACE<VECTOR<T,3>,2>*>(embedding)) PHYSBAM_FATAL_ERROR();
            else Add_Collision_Mesh(embedding->material_surface_mesh,true);
        else if(FREE_PARTICLES<TV>* free_particles=dynamic_cast<FREE_PARTICLES<TV>*>(collision_structures(c))){
            INDIRECT_ARRAY<ARRAY<bool>,ARRAY<int>&> subset=check_collision.Subset(free_particles->nodes);
            ARRAYS_COMPUTATIONS::Fill(subset,true);}
        else PHYSBAM_NOT_IMPLEMENTED(STRING_UTILITIES::string_sprintf("Collisions with %s",typeid(*collision_structures(c)).name()));}
    INDIRECT_ARRAY<ARRAY<bool>,ARRAY<int>&> ignored_subset=check_collision.Subset(ignored_nodes);
    ARRAYS_COMPUTATIONS::Fill(ignored_subset,false);
}
//#####################################################################
// Function Use_Structure_Skip_Collision_Body
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Use_Structure_Collide_Collision_Body(const bool value)
{
    use_structure_collide_collision_body=value;
    structure_collide_collision_body.Resize(deformable_object_structures.m);
}
//#####################################################################
// Function Compute_Candidate_Nodes_For_Collision_Body_Collisions
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Compute_Candidate_Nodes_For_Collision_Body_Collisions(const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>& bodies)
{
    collision_body_candidate_nodes.Resize(collision_body_list.Size(),false);
    for(COLLISION_GEOMETRY_ID id(1);id<=collision_body_list.Size();id++) collision_body_candidate_nodes(id).Remove_All();
    if(!bodies.m) return;
    if(use_spatial_partition){
        ARRAY<COLLISION_GEOMETRY_ID> collision_body_indices;collision_body_indices.Preallocate(100);
        for(int i=1;i<=deformable_body_collection.simulated_particles.m;i++){int p=deformable_body_collection.simulated_particles(i);
            if(check_collision(p)){
                if(thickness_table && thickness_table->Contains(p)) collision_body_list.spatial_partition->collision_body_thickness=thickness_table->Get(p);
                collision_body_list.spatial_partition->Get_Potential_Collisions(particles.X(p),collision_body_indices,get_potential_collisions_already_added);
                for(int j=1;j<=collision_body_indices.m;j++) collision_body_candidate_nodes(collision_body_list.bodies(collision_body_indices(j))->collision_geometry_id).Append(p);}}}
    else{
        ARRAY<int> general_collision_body_candidate_nodes;
        for(int i=1;i<=deformable_body_collection.simulated_particles.m;i++){int p=deformable_body_collection.simulated_particles(i);
            if(check_collision(p)) general_collision_body_candidate_nodes.Append(p);}
        for(COLLISION_GEOMETRY_ID i(1);i<=bodies.m;i++) if(bodies(i)) collision_body_candidate_nodes(i)=general_collision_body_candidate_nodes;}
    if(use_structure_collide_collision_body) for(COLLISION_GEOMETRY_ID body_id(1);body_id<=bodies.m;body_id++) if(bodies(body_id)){
        for(int i=collision_body_candidate_nodes(body_id).m;i>=1;i--){
            int p=collision_body_candidate_nodes(body_id)(i),structure=particle_to_structure(p);
            if(structure && !structure_collide_collision_body(structure).Contains(body_id)) collision_body_candidate_nodes(body_id).Remove_Index_Lazy(i);}}
}
//#####################################################################
// Function Adjust_Nodes_For_Collision_Body_Collisions
//#####################################################################
// TODO: should be accelerated by using the hierarchy to eliminate tests on parts of the surface that aren't close to collision bodies.
template<class TV> int DEFORMABLE_OBJECT_COLLISIONS<TV>::
Adjust_Nodes_For_Collision_Body_Collisions(BINDING_LIST<TV>& binding_list,SOFT_BINDINGS<TV>& soft_bindings,ARRAY_VIEW<const TV> X_old,const T dt,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>* bodies)
{
    soft_bindings.Clamp_Particles_To_Embedded_Positions(true);soft_bindings.Clamp_Particles_To_Embedded_Velocities(true); // TODO: move this elsewhere
    if(!bodies) bodies=&collision_body_list.bodies;
    if(!bodies->m) return 0;
    int interactions=0;
    Compute_Candidate_Nodes_For_Collision_Body_Collisions(*bodies);

    if(use_protectors){
        ARRAY<bool> is_protected(check_collision.m);ARRAY<TV> X_save(particles.X),V_old(particles.V);
        for(COLLISION_GEOMETRY_ID body_id(1);body_id<=bodies->m;body_id++) if((*bodies)(body_id)){
            COLLISION_GEOMETRY<TV>& collision_body=*(*bodies)(body_id);
            for(int j=collision_body_candidate_nodes(body_id).m;j>=1;j--){int node=collision_body_candidate_nodes(body_id)(j);
                if(protecting_bodies_of_nodes(node).Contains(body_id)){
                    if(!is_protected(node) && collision_body.Implicit_Geometry_Lazy_Inside_Extended_Levelset(X_save(node),(T)protection_thickness)){
                        is_protected(node)=true;particles.X(node)=X_save(node);particles.V(node)=V_old(node);particle_states(node).enforce=false;}}
                else if(is_protected(node)) collision_body_candidate_nodes(body_id).Remove_Index_Lazy(j);}
            interactions+=COLLISION_BODY<TV>::Adjust_Nodes_For_Collisions(collision_body,X_old,particles,soft_bindings,collision_body_candidate_nodes(body_id),check_collision,
                collision_tolerance,particle_states,particle_to_collision_body_id,maximum_levelset_collision_projection_velocity,dt,friction_table,thickness_table);}
        ARRAY<int> collision_count(check_collision.m);
        for(COLLISION_GEOMETRY_ID body_id(1);body_id<=bodies->m;body_id++) if((*bodies)(body_id)){
            COLLISION_GEOMETRY<TV>& collision_body=*(*bodies)(body_id);
            for(int j=1;j<=collision_body_candidate_nodes(body_id).m;j++){int node=collision_body_candidate_nodes(body_id)(j);
                if(particle_states(node).enforce && (!is_protected(node) || protecting_bodies_of_nodes(node).Contains(body_id))
                    && collision_body.Implicit_Geometry_Lazy_Inside(particles.X(node),(thickness_table?thickness_table->Get_Default(node,0):0)-(T)1e-5)){
                    collision_count(node)++;if(collision_count(node)>1){particle_states(node).enforce=false;particles.X(node)=X_save(node);particles.V(node)=V_old(node);}}}}}
    else if(disable_multiple_levelset_collisions){
        ARRAY<TV> X_save(particles.X),V_old(particles.V);
        for(COLLISION_GEOMETRY_ID body_id(1);body_id<=bodies->m;body_id++) if((*bodies)(body_id))
            interactions+=COLLISION_BODY<TV>::Adjust_Nodes_For_Collisions(*(*bodies)(body_id),X_old,particles,soft_bindings,collision_body_candidate_nodes(body_id),
                check_collision,collision_tolerance,particle_states,particle_to_collision_body_id,maximum_levelset_collision_projection_velocity,dt,friction_table,thickness_table);
        ARRAY<int> collision_count(check_collision.m);
        for(COLLISION_GEOMETRY_ID body_id(1);body_id<=bodies->m;body_id++) if((*bodies)(body_id))
            for(int j=1;j<=collision_body_candidate_nodes(body_id).m;j++){int node=collision_body_candidate_nodes(body_id)(j);
                if(particle_states(node).enforce && (*bodies)(body_id)->Implicit_Geometry_Lazy_Inside(particles.X(node),(thickness_table?thickness_table->Get_Default(node,0):0)-(T)1e-5)){
                    collision_count(node)++;if(collision_count(node)>1){particle_states(node).enforce=false;particles.X(node)=X_save(node);particles.V(node)=V_old(node);}}}}
    else for(COLLISION_GEOMETRY_ID body_id=bodies->m;body_id>=COLLISION_GEOMETRY_ID(1);body_id--)
        interactions+=COLLISION_BODY<TV>::Adjust_Nodes_For_Collisions(*(*bodies)(body_id),X_old,particles,soft_bindings,collision_body_candidate_nodes(body_id),
            check_collision,collision_tolerance,particle_states,particle_to_collision_body_id,maximum_levelset_collision_projection_velocity,dt,friction_table,thickness_table);

    binding_list.Clamp_Particles_To_Embedded_Positions();binding_list.Clamp_Particles_To_Embedded_Velocities();
    return interactions;
}
//#####################################################################
// Function Adjust_Existing_Nodes_For_Collision_Body_Collisions
//#####################################################################
template<class TV> int DEFORMABLE_OBJECT_COLLISIONS<TV>::
Adjust_Existing_Nodes_For_Collision_Body_Collisions(BINDING_LIST<TV>& binding_list,SOFT_BINDINGS<TV>& soft_bindings,ARRAY_VIEW<const TV> X_old,const T dt,
    const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>* bodies)
{
    int interactions=0;T depth=0,one_over_dt=1/dt;ARRAY_VIEW<TV> &X=particles.X,&V=particles.V;
    for(int i=1;i<=deformable_body_collection.simulated_particles.m;i++){int p=deformable_body_collection.simulated_particles(i);
        COLLISION_PARTICLE_STATE<TV>& collision=particle_states(p);
        if(!collision.enforce) continue;interactions++;
        COLLISION_BODY_HELPER<TV>::Adjust_Point_For_Collision(collision_body_list(particle_to_collision_body_id(p)),X_old(p),X(p),V(p),particles.mass(p),depth,dt,one_over_dt,maximum_levelset_collision_projection_velocity,collision,collision.friction);}
    return interactions;
}
//#####################################################################
// Function Set_Collision_Velocities
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Set_Collision_Velocities(ARRAY_VIEW<TV> V) // for external forces and velocities
{
    enforced_particles.Remove_All();
    for(int i=1;i<=deformable_body_collection.simulated_particles.m;i++){int p=deformable_body_collection.simulated_particles(i);COLLISION_PARTICLE_STATE<TV>& collision=particle_states(p);
        if(collision.enforce){
            T VN=TV::Dot_Product(V(p),collision.normal);
            V(p)+=(collision.VN-VN)*collision.normal;
            collision.delta_VN=collision.VN-VN;}}
// TODO: put this back for efficiency and test
//            enforced_particles.Append(PAIR<int,COLLISION_PARTICLE_STATE<TV> >(p,collision));
}
//#####################################################################
// Function Zero_Out_Collision_Velocities
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Zero_Out_Collision_Velocities(ARRAY_VIEW<TV> V) // for external forces and velocities
{
    for(int p=1;p<=particle_states.m;p++){COLLISION_PARTICLE_STATE<TV>& collision=particle_states(p);
        if(collision.enforce){
            V(p)-=TV::Dot_Product(V(p),collision.normal)*collision.normal;}}
// TODO: put this back for efficiency and test
//    for(int i=1;i<=enforced_particles.m;i++){
//        const int p=enforced_particles(i).x;const TV& collision_normal=enforced_particles(i).y.normal;
//        V(p)-=TV::Dot_Product(V(p),collision_normal)*collision_normal;}
}
//#####################################################################
// Function Reset_Object_Collisions
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Reset_Object_Collisions() // The unusual signature is not great
{
    collisions_on=false;
    if(particle_states.m){
        typedef INDIRECT_ARRAY<ARRAY<COLLISION_PARTICLE_STATE<TV> >,ARRAY<int>&> ARRAY_TYPE;
        typedef FIELD_PROJECTOR<COLLISION_PARTICLE_STATE<TV>,bool,&COLLISION_PARTICLE_STATE<TV>::enforce> T_FIELD_PROJECTOR;
        INDIRECT_ARRAY<ARRAY<COLLISION_PARTICLE_STATE<TV> >,ARRAY<int>&> subset=particle_states.Subset(deformable_body_collection.simulated_particles);
        PROJECTED_ARRAY<ARRAY_TYPE,T_FIELD_PROJECTOR> projected_array=subset.template Project<bool,&COLLISION_PARTICLE_STATE<TV>::enforce>();
        ARRAYS_COMPUTATIONS::Fill(projected_array,false);}
}
//#####################################################################
// Function Add_Collision_Mesh
//#####################################################################
template<class TV> template<class T_MESH> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Add_Collision_Mesh(T_MESH& mesh,const bool collide_with_interior)
{
    if(collide_with_interior) for(int t=1;t<=mesh.elements.m;t++){
        INDIRECT_ARRAY<ARRAY<bool>,typename T_MESH::ELEMENT_TYPE&> subset=check_collision.Subset(mesh.elements(t));
        ARRAYS_COMPUTATIONS::Fill(subset,true);}
    else{bool boundary_nodes_defined=mesh.boundary_nodes!=0;if(!boundary_nodes_defined) mesh.Initialize_Boundary_Nodes();
        const ARRAY<int>* boundary_nodes=mesh.boundary_nodes;
        for(int p=1;p<=boundary_nodes->m;++p) check_collision((*boundary_nodes)(p))=true;
        if(!boundary_nodes_defined){delete mesh.boundary_nodes;mesh.boundary_nodes=0;}}
}
//#####################################################################
// Function Update_Simulated_Particles
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISIONS<TV>::
Update_Simulated_Particles()
{
    particle_to_structure.Resize(particles.array_collection->Size(),false,false);ARRAYS_COMPUTATIONS::Fill(particle_to_structure,0);
    for(int s=1;s<=deformable_object_structures.m;s++) deformable_object_structures(s)->Mark_Nodes_Referenced(particle_to_structure,s);
}
//#####################################################################
template class DEFORMABLE_OBJECT_COLLISIONS<VECTOR<float,1> >;
template class DEFORMABLE_OBJECT_COLLISIONS<VECTOR<float,2> >;
template class DEFORMABLE_OBJECT_COLLISIONS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DEFORMABLE_OBJECT_COLLISIONS<VECTOR<double,1> >;
template class DEFORMABLE_OBJECT_COLLISIONS<VECTOR<double,2> >;
template class DEFORMABLE_OBJECT_COLLISIONS<VECTOR<double,3> >;
#endif

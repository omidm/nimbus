//#####################################################################
// Copyright 2007-2008, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_DEFORMABLE_COLLISIONS
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Math_Tools/clamp.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Vectors/Dot_Product.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Collision_Detection/COLLISION_GEOMETRY_SPATIAL_PARTITION.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TETRAHEDRON_COLLISION_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_CONTACT_GRAPH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_PARTICLE_INTERSECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_SKIP_COLLISION_CHECK.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGIDS_COLLISION_CALLBACKS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/SOLVE_CONTACT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/NORMAL_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_DEFORMABLE_COLLISIONS<TV>::
RIGID_DEFORMABLE_COLLISIONS(SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions_input,SOLIDS_PARAMETERS<TV>& solids_parameters_input)
    :solid_body_collection(solid_body_collection_input),rigid_body_collisions(rigid_body_collisions_input),solids_parameters(solids_parameters_input),fix_position_gap(false),
    normal_relative_velocity_threshold((T)1e-7)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_DEFORMABLE_COLLISIONS<TV>::
~RIGID_DEFORMABLE_COLLISIONS()
{
    precompute_contact_projections.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Get_Rigid_Bodies_Intersecting_Rigid_Body
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_COLLISIONS<TV>::
Get_Rigid_Bodies_Intersecting_Rigid_Body(const int particle_index,ARRAY<int>& rigid_bodies,ARRAY<TV>& collision_locations,ARRAY<TV>& body_distances,
    ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections) const
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    const RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(particle_index);
    const int rigid_body_id=rigid_body.particle_index;
    const int rigid_body_level=rigid_body_collisions.contact_graph.directed_graph.Level_Of_Node(rigid_body_id);

    const ARRAY<VECTOR<int,2> >& level_contact_pairs=rigid_body_collisions.precomputed_contact_pairs_for_level(rigid_body_level);
    for(int i=1;i<=level_contact_pairs.m;i++) if(level_contact_pairs(i)(1)==rigid_body_id || level_contact_pairs(i)(2)==rigid_body_id){
        const int other_body_id=level_contact_pairs(i)(1)==rigid_body_id?level_contact_pairs(i)(2):level_contact_pairs(i)(1);
        const RIGID_BODY<TV>& other_rigid_body=rigid_body_collection.Rigid_Body(other_body_id);
        if(rigid_body_collisions.skip_collision_check.Skip_Pair(rigid_body_id,other_body_id)
            || (rigid_body_collisions.collision_manager && !rigid_body_collisions.collision_manager->Either_Body_Collides_With_The_Other(other_rigid_body.particle_index,rigid_body.particle_index))) continue;
        T smallest_value;int smallest_index;TV collision_location,collision_normal,collision_relative_velocity;bool ignored_separating;
        if(rigid_body_collisions.Get_Deepest_Intersection_Point(rigid_body_id,other_body_id,particle_intersections,
                smallest_value,smallest_index,collision_location,collision_normal,collision_relative_velocity,false,
                rigid_body_collisions.desired_separation_distance-solids_parameters.rigid_body_evolution_parameters.residual_push_out_depth,ignored_separating)){
            // reverse the direction of separation distance to be relative to rigid_body_id
            if(particle_intersections(smallest_index).levelset_body!=rigid_body_id) smallest_value=-smallest_value;
            rigid_bodies.Append(other_rigid_body.particle_index);collision_locations.Append(collision_location);
            smallest_value-=solids_parameters.rigid_body_evolution_parameters.residual_push_out_depth;
            body_distances.Append(smallest_value*collision_normal);}
            else rigid_body_collisions.skip_collision_check.Set_Last_Checked(rigid_body_id,other_body_id);} // set last checked even if ignored_separating==true
}
//#####################################################################
// Function Get_Particles_Intersecting_Rigid_Body
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_COLLISIONS<TV>::
Get_Particles_Intersecting_Rigid_Body(const RIGID_BODY<TV>& rigid_body,ARRAY<int>& particles,ARRAY<TV>& particle_distances) const
{

    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    particles.Remove_All();particle_distances.Remove_All();
    if(!solids_parameters.use_rigid_deformable_contact) return;
    const HASHTABLE<int>* candidate_particles=particles_collided_with_rigid_body.Get_Pointer(rigid_body.particle_index);
    if(!candidate_particles) return;

    T depth=0;
    for(HASHTABLE<int>::ITERATOR i(*candidate_particles);i.Valid();i.Next()){int p=i.Key();
        if(!solid_body_collection.deformable_body_collection.binding_list.Binding_Index_From_Particle_Index(p)
          && rigid_body.Implicit_Geometry_Lazy_Inside_And_Value(deformable_body_collection.particles.X(p),depth,solids_parameters.rigid_body_collision_parameters.collision_body_thickness)){
            depth-=solids_parameters.rigid_body_collision_parameters.collision_body_thickness;
            particles.Append(p);particle_distances.Append(depth*rigid_body.Implicit_Geometry_Normal(deformable_body_collection.particles.X(p)));}}
}
//#####################################################################
// Function Get_Point_Surface_Element_Pairs_Helper
//#####################################################################
template<class T> void
Get_Point_Surface_Element_Pairs_Helper(const RIGID_DEFORMABLE_COLLISIONS<VECTOR<T,1> >& rigid_deformable_collisions,const int particle_index,ARRAY<VECTOR<int,1> >& elements,
    ARRAY<VECTOR<T,1> >& weights,ARRAY<VECTOR<T,1> >& particle_distances)
{
    PHYSBAM_ASSERT(!rigid_deformable_collisions.tetrahedron_collision_bodies.m);
    return;
}
//#####################################################################
// Function Get_Point_Surface_Element_Pairs_Helper
//#####################################################################
template<class T> void
Get_Point_Surface_Element_Pairs_Helper(const RIGID_DEFORMABLE_COLLISIONS<VECTOR<T,2> >& rigid_deformable_collisions,const int particle_index,ARRAY<VECTOR<int,2> >& elements,
    ARRAY<VECTOR<T,2> >& weights,ARRAY<VECTOR<T,2> >& particle_distances)
{
    PHYSBAM_ASSERT(!rigid_deformable_collisions.tetrahedron_collision_bodies.m);
    return;
}
//#####################################################################
// Function Get_Point_Surface_Element_Pairs_Helper
//#####################################################################
template<class T> void
Get_Point_Surface_Element_Pairs_Helper(const RIGID_DEFORMABLE_COLLISIONS<VECTOR<T,3> >& rigid_deformable_collisions,const int particle_index,ARRAY<VECTOR<int,3> >& elements,
    ARRAY<VECTOR<T,3> >& weights,ARRAY<VECTOR<T,3> >& particle_distances)
{
    // TODO: compute once, and then update locally after push from particle (Skip_Collision_Check?)
    //tetrahedralized_volume.hierarchy->Update_Boxes(collision_thickness);
    typedef VECTOR<T,3> TV;
    const ARRAY<TETRAHEDRON_COLLISION_BODY<T>*>& tetrahedron_candidates=rigid_deformable_collisions.particle_tetrahedron_candidates.Get(particle_index);
    PARTICLES<TV>& particles=rigid_deformable_collisions.solid_body_collection.deformable_body_collection.particles;
    ARRAY<int> particles_to_ignore;particles_to_ignore.Append(particle_index);
    particles_to_ignore.Append_Elements(rigid_deformable_collisions.solid_body_collection.deformable_body_collection.soft_bindings.Parents(particle_index));
    for(int i=1;i<=tetrahedron_candidates.m;i++){TETRAHEDRON_COLLISION_BODY<T>& collision_body=*tetrahedron_candidates(i);
        TV w;int t=collision_body.Get_Tetrahedron_Near_Point(particles.X(particle_index),w,particles_to_ignore);if(!t) continue;
        TV surface_weights;int surface_triangle=collision_body.Get_Surface_Triangle(t,w,surface_weights,true);if(!surface_triangle) continue;
        VECTOR<int,3> element(collision_body.triangulated_surface.mesh.elements(surface_triangle));elements.Append(element);weights.Append(surface_weights);
        TV distance=-particles.X(particle_index);for(int j=1;j<=TV::dimension;j++) distance+=particles.X(element(j))*surface_weights(j);
        particle_distances.Append(distance);}
}
//#####################################################################
// Function Get_Objects_Intersecting_Particle
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_COLLISIONS<TV>::
Get_Objects_Intersecting_Particle(const int particle_index,ARRAY<ELEMENT>& triangles,ARRAY<TV>& weights,ARRAY<TV>& particle_distances,ARRAY<int>& bodies,
    ARRAY<TV>& body_distances)
{

    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    const TV& location=deformable_body_collection.particles.X(particle_index);

    // check all candidate rigid bodies
    if(ARRAY<RIGID_COLLISION_GEOMETRY<TV>*>* rigid_bodies_candidates=particle_rigid_body_candidates.Get_Pointer(particle_index)){
        for(int i=1;i<=rigid_bodies_candidates->m;i++){
            const RIGID_GEOMETRY<TV>& rigid_geometry=(*rigid_bodies_candidates)(i)->rigid_geometry;
            if(rigid_geometry.Implicit_Geometry_Lazy_Inside(location,solids_parameters.rigid_body_collision_parameters.collision_body_thickness)){
                bodies.Append(rigid_geometry.particle_index);
                T phi_value;TV normal=rigid_geometry.Implicit_Geometry_Normal(location,phi_value);
                body_distances.Append(-phi_value*normal);}}}

    Get_Point_Surface_Element_Pairs_Helper(*this,particle_index,triangles,weights,particle_distances);
}
//#####################################################################
// Function Apply_Impulse
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_COLLISIONS<TV>::
Apply_Impulse(const int particle,RIGID_BODY<TV>& rigid_body,const TV& impulse)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    // adjust the particle
    const TV particle_delta_V=particles.one_over_effective_mass(particle)*impulse;
    particles.V(particle)+=particle_delta_V;
    // apply impulse to parents if particle is soft bound
    const int soft_binding_index=soft_bindings.Soft_Binding(particle);
    if(soft_binding_index && soft_bindings.use_impulses_for_collisions(soft_binding_index)){
        const int parent=soft_bindings.bindings(soft_binding_index).y;
        BINDING<TV>* hard_binding=soft_bindings.binding_list.Binding(parent);
        if(hard_binding) hard_binding->Apply_Velocity_Change_To_Parents_Based_On_Embedding(particle_delta_V,0);
        else particles.V(parent)=particles.V(particle);} // TODO: add delta instead of setting?

    // adjust the rigid body
    rigid_body.Apply_Impulse_To_Body(particles.X(particle),-impulse);
}
//#####################################################################
// Function Apply_Displacement_To_Particle
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_COLLISIONS<TV>::
Apply_Displacement_To_Particle(const int particle_index,const TV& particle_delta_X)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    particles.X(particle_index)+=particle_delta_X;

    const int soft_binding_index=soft_bindings.Soft_Binding(particle_index);
    if(soft_binding_index && soft_bindings.use_impulses_for_collisions(soft_binding_index)){
        const int parent=soft_bindings.bindings(soft_binding_index).y;BINDING<TV>* hard_binding=soft_bindings.binding_list.Binding(parent);
        if(hard_binding) hard_binding->Apply_Displacement_To_Parents_Based_On_Embedding(particle_delta_X,0);
        else particles.X(parent)=particles.X(particle_index);} // TODO: add delta instead of setting?
}
//#####################################################################
// Function Get_Rigid_And_Tetrahedron_Collision_Bodies
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_COLLISIONS<TV>::
Get_Rigid_And_Tetrahedron_Collision_Bodies()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    tetrahedron_collision_bodies.Remove_All();rigid_collision_bodies.Remove_All();
    for(COLLISION_GEOMETRY_ID i(1);i<=deformable_body_collection.collisions.collision_body_list.bodies.m;i++) if(deformable_body_collection.collisions.collision_body_list.Is_Active(i)){
        COLLISION_GEOMETRY<TV>* collision_body=deformable_body_collection.collisions.collision_body_list.bodies(i);
        if(dynamic_cast<TETRAHEDRON_COLLISION_BODY<T>*>(collision_body)) tetrahedron_collision_bodies.Append(collision_body);
        else if(RIGID_COLLISION_GEOMETRY<TV>* body=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(collision_body))
            if(solid_body_collection.rigid_body_collection.Is_Active(body->rigid_geometry.particle_index))
                rigid_collision_bodies.Append(collision_body);}
}
//#####################################################################
// Function Apply_Rigid_Deformable_Collision_Impulse
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_COLLISIONS<TV>::
Apply_Rigid_Deformable_Collision_Impulse(RIGID_BODY<TV>& rigid_body,const int particle,const TV& location,const TV& normal,const TV& relative_velocity,const T coefficient_of_restitution,
    const T coefficient_of_friction,const bool clamp_friction_magnitude,TV& impulse,bool allow_pull,bool apply_impulse)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;

    // TODO: can this be made more efficient for the case of a static/kinematic rigid body?

    T relative_normal_velocity=TV::Dot_Product(relative_velocity,normal);
    if(relative_normal_velocity>0 && !allow_pull) relative_normal_velocity=0;

    if(!coefficient_of_friction || (allow_pull && relative_normal_velocity>=0)){ // frictionless case
        T nT_impulse1_n=particles.one_over_mass(particle),nT_impulse2_n=0;
        if(!rigid_body.Has_Infinite_Inertia()){
            T_SPIN r2xn=TV::Cross_Product(location-rigid_body.X(),normal);
            nT_impulse2_n=1/rigid_body.Mass()+Dot_Product(r2xn,rigid_body.World_Space_Inertia_Tensor_Inverse()*r2xn);}
        impulse=-(1+coefficient_of_restitution)*relative_normal_velocity/(nT_impulse1_n+nT_impulse2_n)*normal;}
    else{ // friction case
        T_SYMMETRIC_MATRIX impulse_factor=rigid_body.Impulse_Factor(location)+particles.one_over_mass(particle),impulse_factor_inverse=impulse_factor.Inverse();
        // see if friction stops sliding
        impulse=-impulse_factor_inverse*(coefficient_of_restitution*relative_normal_velocity*normal+relative_velocity); // sticking impulse
        T normal_component=TV::Dot_Product(impulse,normal);
        if((impulse-normal_component*normal).Magnitude()>coefficient_of_friction*normal_component){ // sticking impulse was not admissible
            // friction does not stop sliding
            TV relative_tangential_velocity=relative_velocity.Projected_Orthogonal_To_Unit_Direction(normal);
            TV tangential_direction=relative_tangential_velocity;T relative_tangential_velocity_magnitude=tangential_direction.Normalize();
            tangential_direction=tangential_direction.Projected_Orthogonal_To_Unit_Direction(normal); // used to combat error when the "robust" part of normalize is used
            TV impulse_factor_times_direction=impulse_factor*(normal-coefficient_of_friction*tangential_direction);
            assert(TV::Dot_Product(impulse_factor_times_direction,normal));
            TV delta=-(1+coefficient_of_restitution)*relative_normal_velocity/TV::Dot_Product(impulse_factor_times_direction,normal)*impulse_factor_times_direction;
            // clamp friction magnitude: should be (true) in the elastic collision case, (false) in the inelastic collision case
            if(clamp_friction_magnitude){ // should only clamp friction magnitude in the elastic case!
                TV new_relative_velocity=relative_velocity+delta;
                TV new_normal_velocity=new_relative_velocity.Projected_On_Unit_Direction(normal);
                TV new_tangential_velocity=new_relative_velocity-new_normal_velocity;T new_tangential_velocity_magnitude=new_tangential_velocity.Magnitude();
                if(new_tangential_velocity_magnitude>relative_tangential_velocity_magnitude)
                    delta=new_normal_velocity+(relative_tangential_velocity_magnitude/new_tangential_velocity_magnitude)*new_tangential_velocity-relative_velocity;}
            impulse=impulse_factor_inverse*delta;}}

    if(apply_impulse) Apply_Impulse(particle,rigid_body,impulse);
}
//#####################################################################
// Function Update_Rigid_Deformable_Collision_Pair
//#####################################################################
template<class TV> bool RIGID_DEFORMABLE_COLLISIONS<TV>::
Update_Rigid_Deformable_Collision_Pair(RIGID_BODY<TV>& rigid_body,const int particle_index,const T dt,const T time,ARRAY<TV>& X_save,ARRAY<TV>& V_save,ARRAY<TV>& rigid_X_save,
    ARRAY<ROTATION<TV> >& rigid_rotation_save,ARRAY<typename TV::SPIN>& rigid_angular_momentum_difference,ARRAY<TV>& rigid_velocity_difference)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;

    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    // do not collide with hard bound particles or if both rigid body and particle have infinite inertia
    if(solid_body_collection.deformable_body_collection.binding_list.Binding_Index_From_Particle_Index(particle_index) || (rigid_body.Has_Infinite_Inertia() && particles.mass(particle_index)==FLT_MAX)) return false;

    // sync soft bound particles before processing
    const int soft_binding_index=soft_bindings.Soft_Binding(particle_index);
    BINDING<TV>* hard_binding=0;int parent=0;
    const bool impulse_based_binding=soft_binding_index && soft_bindings.use_impulses_for_collisions(soft_binding_index);
    if(impulse_based_binding){
        parent=soft_bindings.bindings(soft_binding_index).y;hard_binding=soft_bindings.binding_list.Binding(parent);
        if(hard_binding){hard_binding->Clamp_To_Embedded_Position();hard_binding->Clamp_To_Embedded_Velocity();}
        if(soft_bindings.use_gauss_seidel_for_impulse_based_collisions){particles.X(particle_index)=particles.X(parent);particles.V(particle_index)=particles.V(parent);}}

    T depth;if(!rigid_body.Implicit_Geometry_Lazy_Inside_And_Value(particles.X(particle_index),depth,solids_parameters.rigid_body_collision_parameters.collision_body_thickness)) return false;
    particles_collided_with_rigid_body.Get_Or_Insert(rigid_body.particle_index).Set(particle_index);
    depth-=solids_parameters.rigid_body_collision_parameters.collision_body_thickness;
    TV collision_normal=rigid_body.Implicit_Geometry_Normal(particles.X(particle_index));
    TV body_V=rigid_body.Pointwise_Object_Velocity(particles.X(particle_index));
    TV collision_relative_velocity=particles.V(particle_index)-body_V;
    //if(fix_position_gap) collision_relative_velocity+=V_save(particle_index);
    if(TV::Dot_Product(collision_relative_velocity,collision_normal)>=0) return false; // do nothing in case of separating pair
    T normal_collision_relative_velocity=0,normal_collision_relative_position=0;
    if(fix_position_gap){
        normal_collision_relative_velocity=TV::Dot_Product(collision_relative_velocity,collision_normal);
        normal_collision_relative_position=TV::Dot_Product(particles.X(particle_index)-X_save(particle_index)+body_V*dt,collision_normal);
        T normal_velocity_scale;
        if(normal_collision_relative_position) normal_velocity_scale=normal_collision_relative_velocity<-normal_relative_velocity_threshold?depth/normal_collision_relative_position:1;
        else normal_velocity_scale=1;
        normal_velocity_scale=clamp(normal_velocity_scale,(T)0,(T)1);
        collision_relative_velocity+=(normal_velocity_scale-1)*normal_collision_relative_velocity*collision_normal;}
    TV impulse;
    Apply_Rigid_Deformable_Collision_Impulse(rigid_body,particle_index,particles.X(particle_index),collision_normal,collision_relative_velocity,0,rigid_body.coefficient_of_friction,true,impulse);

    T dt_recomputed;
    if(fix_position_gap){
        dt_recomputed=0;
        if(normal_collision_relative_velocity) dt_recomputed=normal_collision_relative_position/normal_collision_relative_velocity;}
    else dt_recomputed=dt;

    // update the particle's position
    Apply_Displacement_To_Particle(particle_index,dt_recomputed*particles.one_over_mass(particle_index)*impulse);

    // update the rigid body's position.  Euler_Step_Position_With_New_Velocity will call skip_collision_check.Set_Last_Moved()
    if(!rigid_body.Has_Infinite_Inertia()){
        rigid_body_collisions.collision_callbacks.Restore_Position(rigid_body.particle_index);
        rigid_body_collisions.collision_callbacks.Euler_Step_Position_With_New_Velocity(rigid_body.particle_index,dt_recomputed,time);
        rigid_body_collisions.skip_collision_check.Set_Last_Moved(rigid_body.particle_index);}

    return true;
}
//#####################################################################
// Function Add_Elastic_Collisions
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_COLLISIONS<TV>::
Add_Elastic_Collisions(const T dt,const T time,ARRAY<ROTATION<TV> >& rigid_rotation_save,ARRAY<typename TV::SPIN>& rigid_angular_momentum_difference,
    ARRAY<TV>& rigid_velocity_difference,ARRAY<TV>& rigid_X_save,ARRAY<TWIST<TV> >& rigid_velocity_save,ARRAY<typename TV::SPIN> &rigid_angular_momentum_save,ARRAY<TV>& X_save,ARRAY<TV>& V_save)
{
    rigid_body_collisions.rigid_body_particle_intersections.Remove_All();
    particles_collided_with_rigid_body.Remove_All();
    Get_Rigid_And_Tetrahedron_Collision_Bodies();

    // TODO: Update deformable_body_collection_collisions.collision_body_list.spatial_partition whenever a collision body moves
    if(solids_parameters.deformable_object_collision_parameters.use_spatial_partition_for_levelset_collision_objects)
        solid_body_collection.collision_body_list.Update_Spatial_Partition(solids_parameters.deformable_object_collision_parameters.spatial_partition_voxel_size_heuristic,
            solids_parameters.deformable_object_collision_parameters.spatial_partition_number_of_cells,solids_parameters.deformable_object_collision_parameters.spatial_partition_voxel_size_scale_factor);

    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;

    const bool kinematic_rigid_bodies_only=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m==0;

    V_save-=deformable_body_collection.particles.V;
    if(!kinematic_rigid_bodies_only)
        rigid_body_collisions.rigid_body_particle_intersections.Remove_All();
    particles_collided_with_rigid_body.Remove_All();
    Get_Rigid_And_Tetrahedron_Collision_Bodies();

    ARRAY<VECTOR<int,2> > pairs;
    if(!kinematic_rigid_bodies_only){pairs.Preallocate(10);rigid_body_collisions.skip_collision_check.Reset();rigid_body_collisions.pairs_processed_by_collisions.Remove_All();}

    // TODO: Update deformable_body_collection_collisions.collision_body_list.spatial_partition whenever a collision body moves
    if(solids_parameters.deformable_object_collision_parameters.use_spatial_partition_for_levelset_collision_objects)
        solid_body_collection.collision_body_list.Update_Spatial_Partition(solids_parameters.deformable_object_collision_parameters.spatial_partition_voxel_size_heuristic,
            solids_parameters.deformable_object_collision_parameters.spatial_partition_number_of_cells,solids_parameters.deformable_object_collision_parameters.spatial_partition_voxel_size_scale_factor);

    bool need_another_iteration=true;
    for(int i=1;i<=solids_parameters.rigid_body_collision_parameters.collision_iterations && need_another_iteration;i++){need_another_iteration=false;
        // rigid/rigid collisions
        if(!kinematic_rigid_bodies_only && solids_parameters.rigid_body_collision_parameters.perform_collisions){
            // Calls skip_collision_check.Skip_Pair and if necessary skip_collision_check.Set_Last_Checked
            rigid_body_collisions.Get_Bounding_Box_Collision_Pairs(dt,time,pairs,i==solids_parameters.rigid_body_collision_parameters.collision_iterations,
                i==1,solids_parameters.rigid_body_collision_parameters.collision_bounding_box_thickness);
            for(int j=1;j<=pairs.m;j++){
                rigid_body_collisions.added_bodies(1).Remove_All();rigid_body_collisions.added_bodies(2).Remove_All();
                int id_1=pairs(j)(1),id_2=pairs(j)(2);
                if(rigid_body_collisions.Update_Collision_Pair(id_1,id_2,dt,time,false)){
                    rigid_body_collisions.spatial_partition->Update_Body(solid_body_collection.rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(id_1));
                    rigid_body_collisions.spatial_partition->Update_Body(solid_body_collection.rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(id_2));
                    need_another_iteration=true;}
                rigid_body_collisions.Clean_Up_Fractured_Items_From_Lists(pairs,j,false);}}

        // rigid/deformable collisions
        // TODO: support use_protectors and disable_multiple_levelset_collisions
        if(solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions && rigid_collision_bodies.m){
            deformable_body_collection.collisions.Compute_Candidate_Nodes_For_Collision_Body_Collisions(rigid_collision_bodies);
            for(COLLISION_GEOMETRY_ID collision_body_id(1);collision_body_id<=rigid_collision_bodies.m;collision_body_id++){
                RIGID_BODY<TV>& rigid_body=dynamic_cast<RIGID_BODY<TV>&>(dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>&>(*rigid_collision_bodies(collision_body_id)).rigid_geometry);
                for(int j=1;j<=deformable_body_collection.collisions.collision_body_candidate_nodes(collision_body_id).m;j++){
                    int k=deformable_body_collection.collisions.collision_body_candidate_nodes(collision_body_id)(j);
                    if(Update_Rigid_Deformable_Collision_Pair(rigid_body,k,dt,time,X_save,V_save,rigid_X_save,rigid_rotation_save,rigid_angular_momentum_difference,
                            rigid_velocity_difference) && !rigid_body.Has_Infinite_Inertia()){
                        need_another_iteration=true;
                        if(solids_parameters.deformable_object_collision_parameters.use_spatial_partition_for_levelset_collision_objects)
                            solid_body_collection.collision_body_list.spatial_partition->Update_Body(collision_body_id);}}}}

        // deformable/deformable collisions
        // TODO: is this missing from contact?
        if(solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions){
            int interactions=deformable_body_collection.collisions.Adjust_Nodes_For_Collision_Body_Collisions(solid_body_collection.deformable_body_collection.binding_list,
                solid_body_collection.deformable_body_collection.soft_bindings,deformable_body_collection.particles.X,dt,&tetrahedron_collision_bodies);
            if(interactions!=0 && !kinematic_rigid_bodies_only) need_another_iteration=true;}

        if(!need_another_iteration && i<solids_parameters.rigid_body_collision_parameters.collision_iterations){
            i=solids_parameters.rigid_body_collision_parameters.collision_iterations-1;
            need_another_iteration=true;}} // force it to the last iteration (so it picks up contact pairs for rigid/rigid one time at least)
    V_save+=deformable_body_collection.particles.V;
}
//#####################################################################
// Function Process_Deformable_Contact_With_Kinematic_Rigid_Body
//#####################################################################
template<class TV> bool RIGID_DEFORMABLE_COLLISIONS<TV>::
Process_Deformable_Contact_With_Kinematic_Rigid_Body(RIGID_BODY<TV>& rigid_body,const T dt,const T time,ARRAY<TV>& rigid_X_save,ARRAY<ROTATION<TV> >& rigid_rotation_save,ARRAY<TV>& X_save)
{
    PHYSBAM_ASSERT(rigid_body.Has_Infinite_Inertia());
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;

    // TODO: deal with soft bound particles
    ARRAY<int> particles;Get_Particles_Contacting_Rigid_Body(rigid_body,particles,false);
    if(!particles.m) return false;

    PRECOMPUTE_CONTACT_PROJECTION precompute(rigid_body,false);
    ARRAY_VIEW<TV> V(deformable_body_collection.particles.V);
    exchange(rigid_X_save(rigid_body.particle_index),rigid_body.X());
    exchange(rigid_rotation_save(rigid_body.particle_index),rigid_body.Rotation());
    rigid_body.Update_Angular_Velocity();
    // NOTE: positions and normals are now at time n
    ARRAY<TV> relative_velocities;
    for(int i=1;i<=particles.m;i++){const int p=particles(i);
        const TV &X=X_save(p),V_rel=rigid_body.Pointwise_Object_Velocity(X)-V(p);
        // TODO: double check time
        if(TV::Dot_Product(V_rel,rigid_body.implicit_object->Extended_Normal(X))>=0){precompute.particles.Append(p);relative_velocities.Append(V_rel);}}

    // TODO: note we are recomputing normals inside the following call.  we could reuse them from the above computation.
    Initialize_Rigid_Deformable_Contact_Projection(precompute,X_save);
    // TODO: note we are recomputing V_rel inside the following call.  we could reuse them from the above computation.
    Apply_Rigid_Deformable_Contact_Projection(X_save,V,rigid_body.V(),rigid_body.Angular_Velocity(),precompute);

    // apply friction
    T mu=rigid_body.coefficient_of_friction;
    for(int i=1;i<=precompute.particles.m;i++){int p=precompute.particles(i);
        T delta_VN=TV::Dot_Product(relative_velocities(i),precompute.N(i));
        TV VT_direction=relative_velocities(i).Projected_Orthogonal_To_Unit_Direction(precompute.N(i));T VT_relative_magnitude=VT_direction.Normalize();
        V(p)+=min(mu*delta_VN,VT_relative_magnitude)*VT_direction;}

    // update positions
    for(int i=1;i<=particles.m;i++){const int p=particles(i);deformable_body_collection.particles.X(p)=X_save(p)+dt*V(p);}
    exchange(rigid_X_save(rigid_body.particle_index),rigid_body.X());
    exchange(rigid_rotation_save(rigid_body.particle_index),rigid_body.Rotation());
    rigid_body.Update_Angular_Velocity();

    return true;
}
//#####################################################################
// Function Update_Rigid_Deformable_Contact_Pair
//#####################################################################
template<class TV> bool RIGID_DEFORMABLE_COLLISIONS<TV>::
Update_Rigid_Deformable_Contact_Pair(RIGID_BODY<TV>& rigid_body,const int particle_index,const T dt,const T time,const T epsilon_scale,ARRAY<TV>& X_save,ARRAY<TV>& rigid_X_save,
    ARRAY<ROTATION<TV> >& rigid_rotation_save,const T collision_body_thickness,const bool process_contact_unconditionally)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    ARRAY<TV>& X=X_save;

    // do not collide with hard bound particles or if both rigid body and particle have infinite inertia
    if(solid_body_collection.deformable_body_collection.binding_list.Binding_Index_From_Particle_Index(particle_index) || (rigid_body.Has_Infinite_Inertia() && particles.mass(particle_index)==FLT_MAX)) return false;

    exchange(rigid_X_save(rigid_body.particle_index),rigid_body.X());
    exchange(rigid_rotation_save(rigid_body.particle_index),rigid_body.Rotation());
    rigid_body.Update_Angular_Velocity();
    // NOTE: positions and normals are now at time n

    // sync soft bound particles before processing
    const int soft_binding_index=soft_bindings.Soft_Binding(particle_index);
    BINDING<TV>* hard_binding=0;int parent=0;
    const bool impulse_based_binding=soft_binding_index && soft_bindings.use_impulses_for_collisions(soft_binding_index);
    if(impulse_based_binding){
        parent=soft_bindings.bindings(soft_binding_index).y;hard_binding=soft_bindings.binding_list.Binding(parent);
        // TODO: the first position clamp is probably not necessary if hard bindings were satisfied in X_save
        if(hard_binding){hard_binding->Clamp_To_Embedded_Position(X);hard_binding->Clamp_To_Embedded_Position();hard_binding->Clamp_To_Embedded_Velocity();}
        if(soft_bindings.use_gauss_seidel_for_impulse_based_collisions){
            X(particle_index)=X(parent);particles.X(particle_index)=particles.X(parent);particles.V(particle_index)=particles.V(parent);}}

    TV relative_velocity=particles.V(particle_index)-rigid_body.Pointwise_Object_Velocity(X(particle_index));
    TV normal=rigid_body.implicit_object->Extended_Normal(X(particle_index));
    if(!process_contact_unconditionally && TV::Dot_Product(relative_velocity,normal)>=0){
        exchange(rigid_X_save(rigid_body.particle_index),rigid_body.X());
        exchange(rigid_rotation_save(rigid_body.particle_index),rigid_body.Rotation());
        rigid_body.Update_Angular_Velocity();
        return false;} // do nothing in case of separating pair

    T depth;if(!process_contact_unconditionally && !rigid_body.Implicit_Geometry_Lazy_Inside_And_Value(particles.X(particle_index),depth,collision_body_thickness)){
        exchange(rigid_X_save(rigid_body.particle_index),rigid_body.X());
        exchange(rigid_rotation_save(rigid_body.particle_index),rigid_body.Rotation());
        rigid_body.Update_Angular_Velocity();
        return false;}
    depth-=collision_body_thickness;
    T normal_relative_velocity=0,normal_relative_position=0;
    if(fix_position_gap && !process_contact_unconditionally){
        normal_relative_velocity=TV::Dot_Product(relative_velocity,normal);
        normal_relative_position=TV::Dot_Product(particles.X(particle_index)-X_save(particle_index),normal);
        T normal_velocity_scale;
        if(normal_relative_position) normal_velocity_scale=normal_relative_velocity<-normal_relative_velocity_threshold?depth/normal_relative_position:1;
        else normal_velocity_scale=1;
        normal_velocity_scale=clamp(normal_velocity_scale,(T)0,(T)1);
        relative_velocity+=(normal_velocity_scale-1)*normal_relative_velocity*normal;}
    TV impulse;Apply_Rigid_Deformable_Collision_Impulse(rigid_body,particle_index,X(particle_index),normal,relative_velocity,-1+epsilon_scale,rigid_body.coefficient_of_friction,true,impulse,process_contact_unconditionally);

    T dt_recomputed;
    if(fix_position_gap && !process_contact_unconditionally){
        dt_recomputed=0;
        if(!normal_relative_velocity) dt_recomputed=normal_relative_position/normal_relative_velocity;}
    else dt_recomputed=dt;

    // update the particle's position
    Apply_Displacement_To_Particle(particle_index,dt_recomputed*particles.one_over_mass(particle_index)*impulse);

    // update the rigid body's position
    if(!rigid_body.Has_Infinite_Inertia()){
        rigid_body_collisions.collision_callbacks.Save_Position(rigid_body.particle_index);
        rigid_body_collisions.Euler_Step_Position(rigid_body.particle_index,dt_recomputed,time);}
    else{
        exchange(rigid_X_save(rigid_body.particle_index),rigid_body.X());
        exchange(rigid_rotation_save(rigid_body.particle_index),rigid_body.Rotation());
        rigid_body.Update_Angular_Velocity();}

    return true;
}
//#####################################################################
// Function Process_Precomputed_Contact_With_Kinematic_Rigid_Bodies
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_COLLISIONS<TV>::
Process_Precomputed_Contact_With_Rigid_Bodies()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    // TODO: This is very crude and should be fixed.

    bool has_infinite=false,has_finite=false;
    for(int j=1;j<=precompute_contact_projections.m;j++){if(precompute_contact_projections(j)->rigid_body.Has_Infinite_Inertia()) has_infinite=true;else has_finite=true;}

    if(has_finite)
        for(int iteration=1;iteration<=solids_parameters.rigid_body_collision_parameters.contact_iterations*rigid_body_collisions.contact_level_iterations;iteration++){
            for(int j=1;j<=precompute_contact_projections.m;j++){PRECOMPUTE_CONTACT_PROJECTION& precompute=*precompute_contact_projections(j);
                if(precompute.rigid_body.Has_Infinite_Inertia()) continue;
                TV impulse;
                for(int i=1;i<=precompute.particles.m;i++){const int p=precompute.particles(i);
                    TV relative_velocity=deformable_body_collection.particles.V(p)-precompute.rigid_body.Pointwise_Object_Velocity(deformable_body_collection.particles.X(p));
                    Apply_Rigid_Deformable_Collision_Impulse((RIGID_BODY<TV>&)precompute.rigid_body,p,deformable_body_collection.particles.X(p),precompute.N(i),relative_velocity,0,
                        precompute.rigid_body.coefficient_of_friction,false,impulse,true);}}}
    if(has_infinite) Process_Precomputed_Contact_With_Kinematic_Rigid_Bodies();
}
//#####################################################################
// Function Process_Precomputed_Contact_With_Kinematic_Rigid_Bodies
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_COLLISIONS<TV>::
Process_Precomputed_Contact_With_Kinematic_Rigid_Bodies()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    for(int j=1;j<=precompute_contact_projections.m;j++){PRECOMPUTE_CONTACT_PROJECTION& precompute=*precompute_contact_projections(j);
        if(!precompute.rigid_body.Has_Infinite_Inertia()) continue;
        ARRAY<TV> relative_velocities(precompute.particles.m);
        for(int i=1;i<=precompute.particles.m;i++){const int p=precompute.particles(i);
            relative_velocities(i)=precompute.rigid_body.Pointwise_Object_Velocity(deformable_body_collection.particles.X(p))-deformable_body_collection.particles.V(p);}
        const TWIST<TV>& twist=precompute.rigid_body.Twist();
        Apply_Rigid_Deformable_Contact_Projection(deformable_body_collection.particles.X,deformable_body_collection.particles.V,const_cast<TWIST<TV>&>(twist).linear,const_cast<TWIST<TV>&>(twist).angular,precompute);
        // apply friction
        for(int i=1;i<=precompute.particles.m;i++){int p=precompute.particles(i);
            T delta_VN=TV::Dot_Product(relative_velocities(i),precompute.N(i));
            if(delta_VN<=0) continue;
            TV VT_direction=relative_velocities(i).Projected_Orthogonal_To_Unit_Direction(precompute.N(i));
            T VT_relative_magnitude=VT_direction.Normalize();
            deformable_body_collection.particles.V(p)+=min(precompute.rigid_body.coefficient_of_friction*delta_VN,VT_relative_magnitude)*VT_direction;}}
}
// TODO: start here checking skip_collision_check
//#####################################################################
// Function Process_Contact
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_COLLISIONS<TV>::
Process_Contact(const T dt,const T time,ARTICULATED_RIGID_BODY<TV>* articulated_rigid_body,const bool use_saved_pairs,const bool use_existing_contact,
    ARRAY<TV>& rigid_X_save,ARRAY<ROTATION<TV> >& rigid_rotation_save,ARRAY<TV>& rigid_velocity_difference,ARRAY<typename TV::SPIN>& rigid_angular_momentum_difference,ARRAY<TV>& X_save,
    const T collision_body_thickness)
{
    RIGID_BODY_SKIP_COLLISION_CHECK& skip_collision_check=rigid_body_collisions.skip_collision_check;
    if(!use_saved_pairs) particles_contacting_rigid_body.Remove_All();
    skip_collision_check.Reset();

    bool do_dr=!use_existing_contact || !solids_parameters.rigid_body_collision_parameters.use_persistant_contact;

    // TODO: handle soft bound particles
    solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Positions();
    solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities(); // TODO: necessary?
    solid_body_collection.deformable_body_collection.soft_bindings.Clamp_Particles_To_Embedded_Positions(true);
    solid_body_collection.deformable_body_collection.soft_bindings.Clamp_Particles_To_Embedded_Velocities(true);

    // TODO: reconsider construction of contact graph

    bool need_another_iteration=true;int iteration=0;T epsilon_scale=1;
    ARRAY<ARRAY<VECTOR<int,2> > >& contact_pairs_for_level=use_saved_pairs?rigid_body_collisions.saved_contact_pairs_for_level:rigid_body_collisions.precomputed_contact_pairs_for_level;
    while(need_another_iteration && ++iteration<=solids_parameters.rigid_body_collision_parameters.contact_iterations){need_another_iteration=false;
        if(solids_parameters.rigid_body_collision_parameters.use_epsilon_scaling) epsilon_scale=(T)iteration/solids_parameters.rigid_body_collision_parameters.contact_iterations;
        for(int level=1;level<=rigid_body_collisions.contact_graph.Number_Of_Levels();level++){
            ARRAY<VECTOR<int,2> >& pairs=contact_pairs_for_level(level);
            const ARRAY<int>& rigid_bodies_in_level=rigid_body_collisions.contact_graph.directed_graph.Nodes_In_Level(level);
            bool need_another_level_iteration=true;int level_iteration=0;
            while(need_another_level_iteration && ++level_iteration<=rigid_body_collisions.contact_level_iterations){need_another_level_iteration=false;
                if(solids_parameters.rigid_body_collision_parameters.use_epsilon_scaling_for_level)
                    epsilon_scale=(T)iteration*level_iteration/(solids_parameters.rigid_body_collision_parameters.contact_iterations*rigid_body_collisions.contact_level_iterations);

                // rigid/rigid
                if(solids_parameters.rigid_body_collision_parameters.perform_contact) for(int i=1;i<=pairs.m;i++){int id_1=pairs(i)(1),id_2=pairs(i)(2);
                    if(skip_collision_check.Skip_Pair(id_1,id_2)) continue;
                    if(rigid_body_collisions.prune_contact_using_velocity){
                        TRIPLE<int,int,int>& pair_scale=rigid_body_collisions.pairs_scale.Get(pairs(i).Sorted());
                        if((rigid_body_collisions.contact_level_iterations-level_iteration)%pair_scale.y || (solids_parameters.rigid_body_collision_parameters.contact_iterations-iteration)%pair_scale.x) continue;
                        if(SOLVE_CONTACT::Update_Contact_Pair(rigid_body_collisions,rigid_body_collisions.collision_callbacks,rigid_body_collisions.analytic_contact_registry,id_1,id_2,
                                solids_parameters.rigid_body_evolution_parameters.correct_contact_energy,(int)(rigid_body_collisions.contact_pair_iterations/pair_scale.z),epsilon_scale,dt,time,false)){
                            if(!use_saved_pairs) rigid_body_collisions.saved_contact_pairs_for_level(level).Append_Unique(pairs(i)); // TODO: Potentially inefficient
                            need_another_level_iteration=true;need_another_iteration=true;}}
                    else{
                        if(SOLVE_CONTACT::Update_Contact_Pair(rigid_body_collisions,rigid_body_collisions.collision_callbacks,rigid_body_collisions.analytic_contact_registry,id_1,id_2,
                                solids_parameters.rigid_body_evolution_parameters.correct_contact_energy,rigid_body_collisions.contact_pair_iterations,epsilon_scale,dt,time,false)){
                            if(!use_saved_pairs) rigid_body_collisions.saved_contact_pairs_for_level(level).Append_Unique(pairs(i)); // TODO: Potentially inefficient
                            need_another_level_iteration=true;need_another_iteration=true;}}}

                // prestabilization
                if(articulated_rigid_body)
                    for(int i=1;i<=articulated_rigid_body->contact_level_iterations;i++) for(int j=1;j<=articulated_rigid_body->process_list(level).m;j++){
                        rigid_body_collisions.Apply_Prestabilization_To_Joint(dt,time,*articulated_rigid_body,articulated_rigid_body->process_list(level)(j),epsilon_scale);
                        need_another_level_iteration=need_another_iteration=true;}

                // rigid/deformable
                if(solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions && do_dr)
                    for(int i=1;i<=rigid_bodies_in_level.m;i++){RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(rigid_bodies_in_level(i));
                        if(rigid_body.Has_Infinite_Inertia()) continue; // processed below outside of level iterations
                        // TODO: this needs to use the same used_saved_pairs stuff as below to get the contact pairs right
                        ARRAY<int> particles;Get_Particles_Contacting_Rigid_Body(rigid_body,particles,true);if(!particles.m) continue;
                        for(int k=1;k<=particles.m;k++)
                            if(Update_Rigid_Deformable_Contact_Pair(rigid_body,particles(k),dt,time,epsilon_scale,X_save,rigid_X_save,rigid_rotation_save,collision_body_thickness))
                                need_another_level_iteration=need_another_iteration=true;}}}

                if(articulated_rigid_body && articulated_rigid_body->use_krylov_prestab){
                    articulated_rigid_body->Apply_Poststabilization_With_CG(dt,true,solids_parameters.implicit_solve_parameters.test_system,
                        solids_parameters.implicit_solve_parameters.print_matrix);
                    solid_body_collection.rigid_body_collection.Update_Angular_Momentum();
                    rigid_body_collisions.collision_callbacks.Restore_Positions();
                    for(int i=1;i<=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m;i++)
                        rigid_body_collisions.collision_callbacks.Euler_Step_Position(solid_body_collection.rigid_body_collection.simulated_rigid_body_particles(i),dt,time);}}

    // rigid/deformable (kinematic & static)
    if(do_dr) if(solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions)
        for(int i=1;i<=solid_body_collection.rigid_body_collection.static_and_kinematic_rigid_bodies.m;i++){
            RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(solid_body_collection.rigid_body_collection.static_and_kinematic_rigid_bodies(i));
            if(use_saved_pairs){
                const HASHTABLE<int>* contact_particles=particles_contacting_rigid_body.Get_Pointer(rigid_body.particle_index);
                if(!contact_particles) continue;
                for(HASHTABLE<int>::ITERATOR i(*contact_particles);i.Valid();i.Next()){int p=i.Key();
                    Update_Rigid_Deformable_Contact_Pair(rigid_body,p,dt,time,1,X_save,rigid_X_save,rigid_rotation_save,collision_body_thickness,true);}}
            else{
                ARRAY<int> particles;
                Get_Particles_Contacting_Rigid_Body(rigid_body,particles,true);
                for(int k=1;k<=particles.m;k++){
                    Update_Rigid_Deformable_Contact_Pair(rigid_body,particles(k),dt,time,1,X_save,rigid_X_save,rigid_rotation_save,collision_body_thickness);}}}
    if(rigid_body_collisions.prune_stacks_from_contact) rigid_body_collisions.Apply_Stacking_Contact();
}
//#####################################################################
// Function Get_Particles_Contacting_Rigid_Body
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_COLLISIONS<TV>::
Get_Particles_Contacting_Rigid_Body(const RIGID_BODY<TV>& rigid_body,ARRAY<int>& particles,const bool include_soft_bound)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    particles.Remove_All();
    const HASHTABLE<int>* candidate_particles=particles_collided_with_rigid_body.Get_Pointer(rigid_body.particle_index);
    if(!candidate_particles) return;

    // NOTE: use time n normals, but get intersections based on time n+1
    for(HASHTABLE<int>::ITERATOR i(*candidate_particles);i.Valid();i.Next()){int p=i.Key();
        if(!solid_body_collection.deformable_body_collection.binding_list.Binding_Index_From_Particle_Index(p) && (include_soft_bound || !soft_bindings.Particle_Is_Bound(p))
          && rigid_body.Implicit_Geometry_Lazy_Inside(deformable_body_collection.particles.X(p),solids_parameters.rigid_body_collision_parameters.collision_body_thickness)){
            particles.Append(p);
            particles_contacting_rigid_body.Get_Or_Insert(rigid_body.particle_index).Set(p);}}
}
//#####################################################################
// Function Set_Collision_Velocities
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_COLLISIONS<TV>::
Set_Collision_Velocities(ARRAY_VIEW<TV> V,ARRAY_VIEW<TWIST<TV> > twist,ARRAY<TV>& X_save,ARRAY<TV>& rigid_X_save,ARRAY_VIEW<const ROTATION<TV> >& rigid_rotation_save,
    ARRAY<TWIST<TV> >& rigid_velocity_save,ARRAY<typename TV::SPIN>& rigid_angular_momentum_save,ARRAY<TV>& V_save)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    // TODO: rigid/rigid
    ARRAY_VIEW<const TV> X(solids_parameters.use_trapezoidal_rule_for_velocities?ARRAY_VIEW<const TV>(X_save)
        :ARRAY_VIEW<const TV>(deformable_body_collection.particles.X)); // n+1 positions
    ARRAY_VIEW<const TV> rigid_X(solids_parameters.use_trapezoidal_rule_for_velocities?ARRAY_VIEW<const TV>(rigid_X_save)
        :ARRAY_VIEW<const TV>(solid_body_collection.rigid_body_collection.rigid_body_particle.X)); // n+1 frames
    ARRAY_VIEW<const ROTATION<TV> > rigid_rotation(solids_parameters.use_trapezoidal_rule_for_velocities?ARRAY_VIEW<const ROTATION<TV> >(rigid_rotation_save)
        :ARRAY_VIEW<const ROTATION<TV> >(solid_body_collection.rigid_body_collection.rigid_body_particle.rotation)); // n+1 frames
    ARRAY<PRECOMPUTE_CONTACT_PROJECTION*> precompute_boundary_list;
    for(int body=1;body<=precompute_contact_projections.m;body++){
        PRECOMPUTE_CONTACT_PROJECTION& precompute=*precompute_contact_projections(body),*precompute_boundary=0;const int rigid_p=precompute.rigid_body.particle_index;
        T_WORLD_SPACE_INERTIA_TENSOR inverse_inertia;
        if(solids_parameters.use_trapezoidal_rule_for_velocities) inverse_inertia=precompute.rigid_body.World_Space_Inertia_Tensor_Inverse();
        for(int i=1;i<=precompute.particles.m;i++){const int p=precompute.particles(i);
            TV V_rel=RIGID_BODY<TV>::Pointwise_Object_Velocity(twist(rigid_p),rigid_X(rigid_p),X(p))-V(p);
            TV V_rel_target;
            if(solids_parameters.use_trapezoidal_rule_for_velocities){
                // X_save and frame_save should be at n+1 for trapezoidal rule
                TWIST<TV> twist_n(rigid_velocity_save(rigid_p).linear,inverse_inertia*rigid_angular_momentum_save(rigid_p));
                // want N^T V^{n+1}_rel = N^T(2 V^{n+1/2}_rel - V^n_rel) >= 0
                V_rel_target=(T).5*(RIGID_BODY<TV>::Pointwise_Object_Velocity(twist_n,rigid_X(rigid_p),X(p))-V_save(p));}
            if(TV::Dot_Product(V_rel-V_rel_target,precompute.N(i))>=0){ // TODO: double check normal time should be n+1 for both trap and b.e.
                if(!precompute_boundary) precompute_boundary=new PRECOMPUTE_CONTACT_PROJECTION(precompute.rigid_body,true);
                precompute_boundary->particles.Append(p);
                precompute_boundary->V_rel_target.Append(V_rel_target);}
            else{
                precompute.particles.Remove_Index_Lazy(i);
                precompute.V_rel_target.Remove_Index_Lazy(i);
                precompute.N_over_NT_K_N.Remove_Index_Lazy(i);
                precompute.r.Remove_Index_Lazy(i);
                precompute.N.Remove_Index_Lazy(i);
                precompute.rN.Remove_Index_Lazy(i);
                i--;}}
        if(precompute_boundary){
            Initialize_Rigid_Deformable_Contact_Projection(*precompute_boundary,X);
            precompute_boundary_list.Append(precompute_boundary);}}

    for(int i=1;i<=solids_parameters.rigid_body_collision_parameters.contact_project_iterations;i++) for(int body=1;body<=precompute_boundary_list.m;body++)
        /*if(1 || i==1 || !precompute_boundary_list(body)->rigid_body.Has_Infinite_Inertia())*/{ // only process static/kinematic bodies once
            PRECOMPUTE_CONTACT_PROJECTION& precompute=*precompute_boundary_list(body);
            Apply_Rigid_Deformable_Contact_Projection(X,V,twist(precompute.rigid_body.particle_index).linear,twist(precompute.rigid_body.particle_index).angular,precompute);}
    precompute_boundary_list.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Project_Contact_Pairs
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_COLLISIONS<TV>::
Project_Contact_Pairs(ARRAY_VIEW<TV> V,ARRAY_VIEW<TWIST<TV> > twist)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    // TODO: rigid/rigid
    const int iterations=solids_parameters.rigid_body_collision_parameters.contact_project_iterations;
    for(int k=1;k<=iterations;k++) for(int i=1;i<=precompute_contact_projections.m;i++) if(k==1 || !precompute_contact_projections(i)->rigid_body.Has_Infinite_Inertia()){
        PRECOMPUTE_CONTACT_PROJECTION& precompute=*precompute_contact_projections(i);
        Apply_Rigid_Deformable_Contact_Projection(deformable_body_collection.particles.X,V,twist(precompute.rigid_body.particle_index).linear,twist(precompute.rigid_body.particle_index).angular,precompute);}
    for(int k=1;k<=iterations;k++) for(int i=precompute_contact_projections.m;i>=1;i--) if(k==iterations || !precompute_contact_projections(i)->rigid_body.Has_Infinite_Inertia()){
        PRECOMPUTE_CONTACT_PROJECTION& precompute=*precompute_contact_projections(i);
        Apply_Rigid_Deformable_Contact_Projection(deformable_body_collection.particles.X,V,twist(precompute.rigid_body.particle_index).linear,twist(precompute.rigid_body.particle_index).angular,precompute);}
}
//#####################################################################
// Function Initialize_All_Contact_Projections
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_COLLISIONS<TV>::
Initialize_All_Contact_Projections()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    contact_joints.Remove_All();

    if(solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions && deformable_body_collection.collisions.collisions_on){
        // rigid/deformable
        precompute_contact_projections.Delete_Pointers_And_Clean_Memory();
        for(int i=1;i<=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles(i);
            if(solid_body_collection.rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Contains(p)){
                RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(p);
                if(rigid_body.Has_Infinite_Inertia()) continue;
                if(HASHTABLE<int>* particles=particles_contacting_rigid_body.Get_Pointer(rigid_body.particle_index)) if(particles->Size()){
                    PRECOMPUTE_CONTACT_PROJECTION* precompute=new PRECOMPUTE_CONTACT_PROJECTION(rigid_body,true);
                    particles->Get_Keys(precompute->particles);
                    soft_bindings.Remove_Soft_Bound_Particles(precompute->particles); // TODO: fix
                    if(!precompute->particles.m){delete precompute;continue;}
                    Initialize_Rigid_Deformable_Contact_Projection(*precompute,deformable_body_collection.particles.X);
                    precompute_contact_projections.Append(precompute);}}}

        // set up static and kinematic bodies
        for(int i=1;i<=solid_body_collection.rigid_body_collection.static_and_kinematic_rigid_bodies.m;i++){int p=solid_body_collection.rigid_body_collection.static_and_kinematic_rigid_bodies(i);
            if(solid_body_collection.rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Contains(p)){
                RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(p);
                PRECOMPUTE_CONTACT_PROJECTION* precompute=0;
                if(HASHTABLE<int>* particles=particles_contacting_rigid_body.Get_Pointer(rigid_body.particle_index))
                    for(HASHTABLE<int>::ITERATOR iter(*particles);iter.Valid();iter.Next()){
                        if(!soft_bindings.Particle_Is_Bound(iter.Key())){ // skip soft bound particles. TODO: fix.
                            if(!precompute) precompute=new PRECOMPUTE_CONTACT_PROJECTION(rigid_body,true);
                            precompute->particles.Append(iter.Key());}}
                if(precompute){
                    Initialize_Rigid_Deformable_Contact_Projection(*precompute,deformable_body_collection.particles.X);
                    precompute_contact_projections.Append(precompute);}}}}

    // rigid/rigid
    rigid_body_collisions.Initialize_All_Contact_Projections(solids_parameters.rigid_body_collision_parameters.enforce_rigid_rigid_contact_in_cg);
}
//#####################################################################
// Function Initialize_Rigid_Deformable_Contact_Projection
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_COLLISIONS<TV>::
Initialize_Rigid_Deformable_Contact_Projection(PRECOMPUTE_CONTACT_PROJECTION& precompute,ARRAY_VIEW<const TV> X)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    // precompute.particles is the input list of particles to collide against
    if(!precompute.particles.m) return;

    precompute.N_over_NT_K_N.Resize(precompute.particles.m);
    precompute.r.Resize(precompute.particles.m);
    precompute.N.Resize(precompute.particles.m);
    precompute.rN.Resize(precompute.particles.m);
    precompute.V_rel_target.Resize(precompute.particles.m); // may be prefilled

    const bool has_infinite_inertia=precompute.rigid_body.Has_Infinite_Inertia();

    MATRIX<T,TV::dimension> L;MATRIX<T,T_SPIN::dimension,TV::dimension> rL;MATRIX<T,T_SPIN::dimension> rLr;
    for(int i=1;i<=precompute.particles.m;i++){const int p=precompute.particles(i);
        precompute.N(i)=precompute.rigid_body.implicit_object->Extended_Normal(X(p));
        const T one_over_NT_K_N=deformable_body_collection.particles.mass(p);precompute.r(i)=X(p)-precompute.rigid_body.X();precompute.N_over_NT_K_N(i)=precompute.N(i)*one_over_NT_K_N;
        if(has_infinite_inertia) continue;
        precompute.rN(i)=TV::Cross_Product(precompute.r(i),precompute.N(i));
        L+=MATRIX<T,TV::dimension>::Outer_Product(precompute.N(i),precompute.N_over_NT_K_N(i));
        rL+=MATRIX<T,T_SPIN::dimension,TV::dimension>::Outer_Product(precompute.rN(i),precompute.N_over_NT_K_N(i));
        rLr+=MATRIX<T,T_SPIN::dimension>::Outer_Product(precompute.rN(i),precompute.rN(i)*one_over_NT_K_N);}

    if(!has_infinite_inertia){
        typename MATRIX_POLICY<T_SPIN>::SYMMETRIC_MATRIX inertia=precompute.rigid_body.World_Space_Inertia_Tensor();
        precompute.A.Set_Submatrix(1,1,L+precompute.rigid_body.Mass());
        precompute.A.Set_Submatrix(TV::dimension+1,1,rL);
        precompute.A.Set_Submatrix(1,TV::dimension+1,rL.Transposed());
        precompute.A.Set_Submatrix(TV::dimension+1,TV::dimension+1,rLr+inertia);
        if(precompute.A_inverted) precompute.A.In_Place_Cholesky_Inverse();}
}
//#####################################################################
// Function Apply_Rigid_Deformable_Contact_Projection
//#####################################################################
template<class TV> bool RIGID_DEFORMABLE_COLLISIONS<TV>::
Apply_Rigid_Deformable_Contact_Projection(ARRAY_VIEW<const TV> X,ARRAY_VIEW<TV> V,TV& rigid_V,T_SPIN& rigid_angular_velocity,PRECOMPUTE_CONTACT_PROJECTION& precompute)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    if(!precompute.particles.m) return false;

    const bool has_infinite_inertia=precompute.rigid_body.Has_Infinite_Inertia();

    // construct the rhs
    TV Lv;T_SPIN rLv;ARRAY<T> N_V_over_NT_K_N(precompute.particles.m);
    for(int i=1;i<=precompute.particles.m;i++){const int p=precompute.particles(i);
        const TV V_rel=rigid_V+TV::Cross_Product(rigid_angular_velocity,X(p)-precompute.rigid_body.X())-V(p);
        N_V_over_NT_K_N(i)=TV::Dot_Product(precompute.N_over_NT_K_N(i),precompute.V_rel_target(i)-V_rel);
        if(has_infinite_inertia) continue;
        Lv+=precompute.N(i)*N_V_over_NT_K_N(i);
        rLv+=precompute.rN(i)*N_V_over_NT_K_N(i);}

    // compute the rigid body impulse and apply to the body
    TV m_inverse_j;T_SPIN m_inverse_j_tau;
    if(!has_infinite_inertia){
        VECTOR<T,TV::dimension+T_SPIN::dimension> b(Lv,rLv);
        VECTOR<T,TV::dimension+T_SPIN::dimension> m_inverse_j_full;
        if(precompute.A_inverted) m_inverse_j_full=precompute.A*b;
        else m_inverse_j_full=precompute.A.In_Place_Cholesky_Solve(b);
        m_inverse_j_full.Get_Subvector(1,m_inverse_j);m_inverse_j_full.Get_Subvector(TV::dimension+1,m_inverse_j_tau);
        rigid_V+=m_inverse_j;rigid_angular_velocity+=m_inverse_j_tau;}

    // compute particle impulses and apply to particles
    for(int i=1;i<=precompute.particles.m;i++){const int p=precompute.particles(i);
        T jc=-N_V_over_NT_K_N(i);
        if(!has_infinite_inertia){
            TV delta_pointwise_object_velocity=m_inverse_j+TV::Cross_Product(m_inverse_j_tau,precompute.r(i));
            jc+=TV::Dot_Product(precompute.N_over_NT_K_N(i),delta_pointwise_object_velocity);}
        V(p)+=deformable_body_collection.particles.one_over_mass(p)*jc*precompute.N(i);}

    return true;
}
//#####################################################################
// Function Push_Out_From_Rigid_Body
//#####################################################################
template<class TV> bool RIGID_DEFORMABLE_COLLISIONS<TV>::
Push_Out_From_Rigid_Body(RIGID_BODY<TV>& rigid_body,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const T move_fraction)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    RIGID_BODY<TV>& parent_rigid_body=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Get_Parent(rigid_body);
    T threshold=rigid_body.Length_Scale_Squared()*(T)1e-2;

    PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    
    ARRAY<int> particle_interactions,rigid_body_interactions;ARRAY<TV> particle_distances,rigid_body_collision_locations,rigid_body_distances;
    if(solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions)
        Get_Particles_Intersecting_Rigid_Body(rigid_body,particle_interactions,particle_distances);
    if(rigid_body_collection.simulated_rigid_body_particles.m)
        Get_Rigid_Bodies_Intersecting_Rigid_Body(rigid_body.particle_index,rigid_body_interactions,rigid_body_collision_locations,rigid_body_distances,particle_intersections);
    particle_distances*=move_fraction;rigid_body_distances*=move_fraction;

    HASHTABLE<int> affected_particles;
    affected_particles.Set_All(particle_interactions);
    for(int i=1;i<=particle_interactions.m;i++){int p=particle_interactions(i);
        const int soft_binding_index=solid_body_collection.deformable_body_collection.soft_bindings.Soft_Binding(p);
        BINDING<TV>* hard_binding=0;int parent=0;
        const bool impulse_based_binding=soft_binding_index && solid_body_collection.deformable_body_collection.soft_bindings.use_impulses_for_collisions(soft_binding_index);
        bool discard=false;
        if(impulse_based_binding){
            parent=solid_body_collection.deformable_body_collection.soft_bindings.bindings(soft_binding_index).y;
            hard_binding=solid_body_collection.deformable_body_collection.soft_bindings.binding_list.Binding(parent);
            if(hard_binding){
                ARRAY<int> parents=hard_binding->Parents();
                for(int i=1;i<=parents.m;i++) if(affected_particles.Contains(parents(i))) discard=true;
                affected_particles.Set_All(parents);}}
        if(discard){
            particle_interactions.Remove_Index_Lazy(i);
            particle_distances.Remove_Index_Lazy(i);
            i--;}}

    // Sync soft bindings
    for(int i=1;i<=particle_interactions.m;i++){int p=particle_interactions(i);
        const int soft_binding_index=solid_body_collection.deformable_body_collection.soft_bindings.Soft_Binding(p);
        BINDING<TV>* hard_binding=0;int parent=0;
        const bool impulse_based_binding=soft_binding_index && solid_body_collection.deformable_body_collection.soft_bindings.use_impulses_for_collisions(soft_binding_index);
        if(impulse_based_binding){
            parent=solid_body_collection.deformable_body_collection.soft_bindings.bindings(soft_binding_index).y;
            hard_binding=solid_body_collection.deformable_body_collection.soft_bindings.binding_list.Binding(parent);
            if(hard_binding) hard_binding->Clamp_To_Embedded_Position();
            if(solid_body_collection.deformable_body_collection.soft_bindings.use_gauss_seidel_for_impulse_based_collisions) particles.X(p)=particles.X(parent);}}

    if(particle_interactions.m==0 && rigid_body_interactions.m==0) return false;

    ARRAY<T_SYMMETRIC_MATRIX> K_inverse(rigid_body_interactions.m);
    for(int i=1;i<=rigid_body_interactions.m;i++){
        RIGID_BODY<TV>& parent_other_rigid_body=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Get_Parent(rigid_body_collection.Rigid_Body(rigid_body_interactions(i)));
        if(!parent_other_rigid_body.Has_Infinite_Inertia() || rigid_body_collisions.use_static_body_masses) K_inverse(i)=parent_other_rigid_body.Impulse_Factor(rigid_body_collision_locations(i)).Inverse();
        else K_inverse(i)=T_SYMMETRIC_MATRIX::Identity_Matrix();}

    TV velocity;T_SPIN angular_velocity;
    if(!parent_rigid_body.Has_Infinite_Inertia()){
        // TODO: static deformable particles
        TV ms[2];VECTOR<T,T_SPIN::dimension> mrs[2];MATRIX<T,T_SPIN::dimension,TV::dimension> mr[2];MATRIX<T,T_SPIN::dimension> mrr[2];
        T m=0; // non static particles
        for(int i=1;i<=particle_interactions.m;i++){int p=particle_interactions(i);
            T mass=particles.mass(p);TV radius=particles.X(p)-parent_rigid_body.X();
            m+=mass;ms[0]+=mass*particle_distances(i);mrs[0]+=mass*TV::Cross_Product(radius,particle_distances(i));
            MATRIX<T,T_SPIN::dimension,TV::dimension> ri=MATRIX<T,T_SPIN::dimension,TV::dimension>::Cross_Product_Matrix(radius);
            mr[0]+=mass*ri;mrr[0]+=mass*ri.Times_Transpose(ri);}

        T_SYMMETRIC_MATRIX K_inverse_sum[2]; // non static bodies

        int number_of_static_bodies=0;TV centroid;
        for(int i=1;i<=rigid_body_interactions.m;i++){
            RIGID_BODY<TV>& parent_other_rigid_body=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Get_Parent(rigid_body_collection.Rigid_Body(rigid_body_interactions(i)));
            int index=parent_other_rigid_body.Has_Infinite_Inertia()?1:0;
            if(index){centroid+=rigid_body_collision_locations(i)-parent_rigid_body.X();number_of_static_bodies++;}
            TV K_inverse_s=K_inverse(i)*rigid_body_distances(i),radius=rigid_body_collision_locations(i)-parent_rigid_body.X();
            K_inverse_sum[index]+=K_inverse(i);ms[index]+=K_inverse_s;mrs[index]+=TV::Cross_Product(radius,K_inverse_s);
            MATRIX<T,T_SPIN::dimension,TV::dimension> r_K_inverse=K_inverse(i).Cross_Product_Matrix_Times(radius);
            mr[index]+=r_K_inverse;mrr[index]+=r_K_inverse.Times_Cross_Product_Matrix_Transpose(radius);}

        if(number_of_static_bodies) centroid/=(T)number_of_static_bodies;

        // finish setting up the terms for the dynamic bodies, including diagonal terms and dynamic particle masses
        K_inverse_sum[0]+=parent_rigid_body.Mass()+m;
        mrr[0]+=parent_rigid_body.World_Space_Inertia_Tensor();

        // determine how static bodies affect the system
        typename MATRIX_POLICY<TV>::DIAGONAL_MATRIX eigenvalues;MATRIX<T,TV::dimension> eigenvectors;
        int equation_type=0; // indicate which type of system we will be solving below
        if(number_of_static_bodies==0) equation_type=0; // 0 dof specified by static bodies
        else if(number_of_static_bodies==1) equation_type=1;
        else{ // need to figure out if the static bodies are effectively coincident or colinear
            T_SYMMETRIC_MATRIX R_R_transpose;
            for(int i=1;i<=rigid_body_interactions.m;i++){
                RIGID_BODY<TV>& parent_other_rigid_body=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Get_Parent(rigid_body_collection.Rigid_Body(rigid_body_interactions(i)));
                if(parent_other_rigid_body.Has_Infinite_Inertia()) R_R_transpose+=T_SYMMETRIC_MATRIX::Outer_Product(rigid_body_collision_locations(i)-parent_rigid_body.X()-centroid);}
            R_R_transpose.Solve_Eigenproblem(eigenvalues,eigenvectors);
            int rank=TV::Componentwise_Greater_Equal(eigenvalues.To_Vector(),TV::All_Ones_Vector()*threshold).Number_True();
            if(rank==0) equation_type=1; // 3/6 dof (2/3 dof in 2d) specified by static bodies
            else if(rank==1 && TV::dimension==3) equation_type=2; // 5 dof specified by static bodies (not possible in 2d)
            else equation_type=3;} // full constrained by static bodies

        MATRIX<T,TV::dimension+T_SPIN::dimension> M[2];
        for(int i=0;i<=1;i++) if(equation_type!=(i?0:3)){
            M[i].Set_Submatrix(1,1,K_inverse_sum[i]);
                M[i].Set_Submatrix(TV::dimension+1,1,mr[i]);
            M[i].Set_Submatrix(1,TV::dimension+1,mr[i].Transposed());
            M[i].Set_Submatrix(TV::dimension+1,TV::dimension+1,mrr[i]);}

        // set up our system
        MATRIX<T,TV::dimension+T_SPIN::dimension> A;VECTOR<T,TV::dimension+T_SPIN::dimension> b;
        if(equation_type==0){A=M[0];b=VECTOR<T,TV::dimension+T_SPIN::dimension>(ms[0],mrs[0]);}
        else if(equation_type==3){A=M[1];b=VECTOR<T,TV::dimension+T_SPIN::dimension>(ms[1],mrs[1]);}
        else if(equation_type==1){ // one outer body static
            MATRIX<T,T_SPIN::dimension,TV::dimension+T_SPIN::dimension> Z;MATRIX<T,TV::dimension,TV::dimension+T_SPIN::dimension> Y;
            Z.Set_Submatrix(1,1,-MATRIX_POLICY<TV>::CROSS_PRODUCT_MATRIX::Cross_Product_Matrix(centroid));Z.Set_Submatrix(1,TV::dimension+1,MATRIX<T,T_SPIN::dimension>::Identity_Matrix());
            Y.Set_Submatrix(1,1,MATRIX<T,TV::dimension>::Identity_Matrix());Y.Set_Submatrix(1,T_SPIN::dimension+1,-MATRIX_POLICY<TV>::CROSS_PRODUCT_MATRIX::Cross_Product_Matrix(centroid));
            A.Set_Submatrix(1,1,Z*M[0]);
            A.Set_Submatrix(T_SPIN::dimension+1,1,Y*M[1]);
            b.Set_Subvector(1,Z*VECTOR<T,TV::dimension+T_SPIN::dimension>(ms[0],mrs[0]));
            b.Set_Subvector(T_SPIN::dimension+1,Y*VECTOR<T,TV::dimension+T_SPIN::dimension>(ms[1],mrs[1]));}
        else{
            assert(equation_type==2 && TV::dimension==3);
            int u_index=eigenvalues.To_Vector().Arg_Max();
            MATRIX<T,TV::dimension,TV::dimension+T_SPIN::dimension> Z_helper;
            Z_helper.Set_Submatrix(1,1,-MATRIX_POLICY<TV>::CROSS_PRODUCT_MATRIX::Cross_Product_Matrix(centroid));Z_helper.Set_Submatrix(1,TV::dimension+1,MATRIX<T,T_SPIN::dimension>::Identity_Matrix());
            Z_helper=eigenvectors.Transpose_Times(Z_helper);
            MATRIX<T,TV::dimension,TV::dimension+T_SPIN::dimension> Y;
            Y.Set_Submatrix(1,1,MATRIX<T,TV::dimension>::Identity_Matrix());Y.Set_Submatrix(1,T_SPIN::dimension+1,MATRIX_POLICY<TV>::CROSS_PRODUCT_MATRIX::Cross_Product_Matrix(centroid));
            VECTOR<T,TV::dimension+T_SPIN::dimension> b_helper[2]={VECTOR<T,TV::dimension+T_SPIN::dimension>(ms[0],mrs[0]),VECTOR<T,TV::dimension+T_SPIN::dimension>(ms[1],mrs[1])};
            A.Set_Submatrix(TV::dimension*(TV::dimension==3)+1,1,Y*M[1]); // Don't compile index out of bounds in 2D, even though this case cannot happen in 2D.
            b.Set_Subvector(TV::dimension+1,Y*b_helper[1]);
            MATRIX<T,1,TV::dimension+T_SPIN::dimension> A_row;
            for(int i=1;i<=TV::dimension;i++){
                Z_helper.Get_Submatrix(i,1,A_row);
                int matrix_index=i==u_index?0:1;
                A.Set_Submatrix(i,1,A_row*M[matrix_index]);
                b(i)=(A_row*b_helper[matrix_index]).x;}}

        // compute impulse as delta twist
        VECTOR<T,TV::dimension+T_SPIN::dimension> delta_twist;
        if(equation_type==0 || equation_type==3) delta_twist=A.In_Place_Cholesky_Solve(b);
        else delta_twist=A.Solve_Linear_System(b);
        delta_twist.Get_Subvector(1,velocity);delta_twist.Get_Subvector(TV::dimension+1,angular_velocity);}

    // apply push to particles
    for(int i=1;i<=particle_interactions.m;i++){int p=particle_interactions(i);
        Apply_Displacement_To_Particle(p,-particle_distances(i)+velocity+TV::Cross_Product(angular_velocity,particles.X(p)-parent_rigid_body.X()));}

    // apply push to other rigid bodies
    for(int i=1;i<=rigid_body_interactions.m;i++) if(!rigid_body_collection.Rigid_Body(rigid_body_interactions(i)).Has_Infinite_Inertia()){
        RIGID_BODY<TV>& other_rigid_body=rigid_body_collection.Rigid_Body(rigid_body_interactions(i));
        RIGID_BODY<TV>& parent_other_rigid_body=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Get_Parent(other_rigid_body);
        TV impulse=K_inverse(i)*(-rigid_body_distances(i)+velocity+TV::Cross_Product(angular_velocity,rigid_body_collision_locations(i)-parent_rigid_body.X()));
        T_SPIN other_angular_velocity=parent_other_rigid_body.World_Space_Inertia_Tensor_Inverse()*TV::Cross_Product(rigid_body_collision_locations(i)-parent_other_rigid_body.X(),impulse);
        parent_other_rigid_body.X()+=impulse/parent_other_rigid_body.Mass();
        if(solids_parameters.rigid_body_collision_parameters.use_push_out_rotation) parent_other_rigid_body.Rotation()=ROTATION<TV>::From_Rotation_Vector(other_angular_velocity)*parent_other_rigid_body.Rotation();parent_other_rigid_body.Rotation().Normalize();
        parent_other_rigid_body.Update_Angular_Velocity();
        parent_other_rigid_body.Update_Bounding_Box();
        rigid_body_collisions.skip_collision_check.Set_Last_Moved(other_rigid_body.particle_index);}

    // apply push to the rigid body
    if(!parent_rigid_body.Has_Infinite_Inertia()){
        parent_rigid_body.X()+=velocity;
        if(solids_parameters.rigid_body_collision_parameters.use_push_out_rotation) parent_rigid_body.Rotation()=ROTATION<TV>::From_Rotation_Vector(angular_velocity)*parent_rigid_body.Rotation();parent_rigid_body.Rotation().Normalize();
        parent_rigid_body.Update_Angular_Velocity();
        parent_rigid_body.Update_Bounding_Box();
        rigid_body_collisions.skip_collision_check.Set_Last_Moved(rigid_body.particle_index);}

    return true;
}
//#####################################################################
// Function Compute_Particle_Candidates_Fill_Hash
//#####################################################################
template<class TV,class T_BODY> void Compute_Particle_Candidates_Fill_Hash(HASHTABLE<int,ARRAY<T_BODY*> >& particle_collision_candidates,
    const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>& collision_bodies,const ARRAY<ARRAY<int>,COLLISION_GEOMETRY_ID>& collision_body_candidate_nodes)
{
    particle_collision_candidates.Remove_All();
    for(COLLISION_GEOMETRY_ID body(1);body<=collision_bodies.m;body++) if(collision_bodies(body))
        for(int j=1;j<=collision_body_candidate_nodes(collision_bodies(body)->collision_geometry_id).m;j++){int k=collision_body_candidate_nodes(body)(j);
            particle_collision_candidates.Get_Or_Insert(k).Append(dynamic_cast<T_BODY*>(collision_bodies(body)));}
}
//#####################################################################
// Function Compute_Particle_Candidates_Helper
//#####################################################################
template<class T> void Compute_Particle_Candidates_Helper(RIGID_DEFORMABLE_COLLISIONS<VECTOR<T,1> >& rigid_deformable_collisions,
    const ARRAY<COLLISION_GEOMETRY<VECTOR<T,1> >*,COLLISION_GEOMETRY_ID>& tetrahedron_collision_bodies,const ARRAY<COLLISION_GEOMETRY<VECTOR<T,1> >*,COLLISION_GEOMETRY_ID>& rigid_collision_bodies)
{}
template<class T> void Compute_Particle_Candidates_Helper(RIGID_DEFORMABLE_COLLISIONS<VECTOR<T,2> >& rigid_deformable_collisions,
    const ARRAY<COLLISION_GEOMETRY<VECTOR<T,2> >*,COLLISION_GEOMETRY_ID>& tetrahedron_collision_bodies,const ARRAY<COLLISION_GEOMETRY<VECTOR<T,2> >*,COLLISION_GEOMETRY_ID>& rigid_collision_bodies)
{}
template<class T> void Compute_Particle_Candidates_Helper(RIGID_DEFORMABLE_COLLISIONS<VECTOR<T,3> >& rigid_deformable_collisions,
    const ARRAY<COLLISION_GEOMETRY<VECTOR<T,3> >*,COLLISION_GEOMETRY_ID>& tetrahedron_collision_bodies,const ARRAY<COLLISION_GEOMETRY<VECTOR<T,3> >*,COLLISION_GEOMETRY_ID>& rigid_collision_bodies)
{
    typedef VECTOR<T,3> TV;
    DEFORMABLE_OBJECT_COLLISIONS<TV>& deformable_object_collisions=rigid_deformable_collisions.solid_body_collection.deformable_body_collection.collisions;
    // Push out from particles only if there are tetrahedron collision bodies; if only rigid bodies and particles, just process around rigid bodies
    deformable_object_collisions.Compute_Candidate_Nodes_For_Collision_Body_Collisions(deformable_object_collisions.collision_body_list.bodies);
    Compute_Particle_Candidates_Fill_Hash(rigid_deformable_collisions.particle_tetrahedron_candidates,tetrahedron_collision_bodies,deformable_object_collisions.collision_body_candidate_nodes);
    Compute_Particle_Candidates_Fill_Hash(rigid_deformable_collisions.particle_rigid_body_candidates,rigid_collision_bodies,deformable_object_collisions.collision_body_candidate_nodes);
}
//#####################################################################
// Function Push_Out_From_Particle
//#####################################################################
template<class TV> bool RIGID_DEFORMABLE_COLLISIONS<TV>::
Push_Out_From_Particle(const int particle)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    // Sync soft binding
    const int soft_binding_index=solid_body_collection.deformable_body_collection.soft_bindings.Soft_Binding(particle);
    BINDING<TV>* hard_binding=0;int parent=0;
    const bool impulse_based_binding=soft_binding_index && solid_body_collection.deformable_body_collection.soft_bindings.use_impulses_for_collisions(soft_binding_index);
    if(impulse_based_binding){
        parent=solid_body_collection.deformable_body_collection.soft_bindings.bindings(soft_binding_index).y;
        hard_binding=solid_body_collection.deformable_body_collection.soft_bindings.binding_list.Binding(parent);
        if(hard_binding) hard_binding->Clamp_To_Embedded_Position();
        if(solid_body_collection.deformable_body_collection.soft_bindings.use_gauss_seidel_for_impulse_based_collisions) particles.X(particle)=particles.X(parent);}

    // TODO: Note that we have not synced the bindings for these.
    ARRAY<ELEMENT> elements;ARRAY<TV> weights;ARRAY<int> bodies;ARRAY<TV> particle_distances,rigid_body_distances;
    Get_Objects_Intersecting_Particle(particle,elements,weights,particle_distances,bodies,rigid_body_distances);
    if(elements.m==0 && bodies.m==0) return false;

    ARRAY<LINEAR_BINDING<TV,ELEMENT::dimension>*> bindings(elements.m);TV ms;T m=0;
    for(int i=1;i<=elements.m;i++){
        bindings(i)=new LINEAR_BINDING<TV,ELEMENT::dimension>(particles,0,elements(i),weights(i));
        T effective_mass=1/bindings(i)->One_Over_Effective_Mass();
        ms+=effective_mass*particle_distances(i);m+=effective_mass;}

    TV& particle_location=particles.X(particle);
    ARRAY<T_SYMMETRIC_MATRIX> K_inverse(bodies.m);T_SYMMETRIC_MATRIX K_inverse_sum;
    bool have_static_bodies=false;for(int i=1;i<=bodies.m;i++) if(rigid_body_collection.Rigid_Body(bodies(i)).Has_Infinite_Inertia()){have_static_bodies=true;break;}
    for(int i=1;i<=bodies.m;i++) if(rigid_body_collection.Rigid_Body(bodies(i)).Has_Infinite_Inertia()==have_static_bodies){ // process only the static bodies if we have any
        RIGID_BODY<TV>& other_rigid_body=rigid_body_collection.Rigid_Body(bodies(i));
        RIGID_BODY<TV>& parent_other_rigid_body=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Get_Parent(other_rigid_body);
        if(!parent_other_rigid_body.Has_Infinite_Inertia() || rigid_body_collisions.use_static_body_masses) K_inverse(i)=parent_other_rigid_body.Impulse_Factor(particle_location).Inverse();
        else K_inverse(i)=T_SYMMETRIC_MATRIX::Identity_Matrix();
        K_inverse_sum+=K_inverse(i);ms+=K_inverse(i)*rigid_body_distances(i);}

    TV velocity;
    if(have_static_bodies) velocity=K_inverse_sum.Solve_Linear_System(ms);
    else velocity=(K_inverse_sum+m+particles.mass(particle)).Solve_Linear_System(ms);

    // apply push to elements
    for(int i=1;i<=elements.m;i++){
        TV impulse=(velocity-particle_distances(i))/bindings(i)->One_Over_Effective_Mass();
        for(int j=1;j<=ELEMENT::dimension;j++)
            Apply_Displacement_To_Particle(bindings(i)->parents[j],particles.one_over_mass(bindings(i)->parents[j])*bindings(i)->weights[j]*impulse);}

    // apply push to rigid bodies
    for(int i=1;i<=bodies.m;i++) if(!rigid_body_collection.Rigid_Body(bodies(i)).Has_Infinite_Inertia()){
        RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(bodies(i));
        RIGID_BODY<TV>& parent_rigid_body=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Get_Parent(rigid_body);
        TV impulse=K_inverse(i)*(-rigid_body_distances(i)+velocity);
        T_SPIN angular_velocity=parent_rigid_body.World_Space_Inertia_Tensor_Inverse()*TV::Cross_Product(particle_location-parent_rigid_body.X(),impulse);
        rigid_body.X()+=impulse/parent_rigid_body.Mass();
        rigid_body.Rotation()=ROTATION<TV>::From_Rotation_Vector(angular_velocity)*rigid_body.Rotation();rigid_body.Rotation().Normalize();
        rigid_body_collisions.skip_collision_check.Set_Last_Moved(rigid_body.particle_index);
        parent_rigid_body.Update_Angular_Velocity();}

    // apply push to the particle
    Apply_Displacement_To_Particle(particle,velocity);

    bindings.Delete_Pointers_And_Clean_Memory();

    return true;
}
//#####################################################################
// Function Process_Push_Out
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_COLLISIONS<TV>::
Process_Push_Out()
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    const bool kinematic_rigid_bodies_only=rigid_body_collection.simulated_rigid_body_particles.m==0;

    // TODO: handle soft bound particles
    solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Positions();
    solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities(); // TODO: necessary?
    solid_body_collection.deformable_body_collection.soft_bindings.Clamp_Particles_To_Embedded_Positions(true);
    solid_body_collection.deformable_body_collection.soft_bindings.Clamp_Particles_To_Embedded_Velocities(true);

    rigid_body_collisions.skip_collision_check.Reset();
    ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> > particle_intersections;
    particle_intersections.Preallocate(100);

    bool need_another_iteration=true;int iteration=0;
    while(need_another_iteration && ++iteration<=rigid_body_collisions.push_out_iterations){need_another_iteration=false;
        if(!kinematic_rigid_bodies_only){
            if(rigid_body_collisions.use_freezing_with_push_out) rigid_body_collisions.Clear_Temporarily_Static();
            for(int level=1;level<=rigid_body_collisions.contact_graph.Number_Of_Levels();level++){
                const ARRAY<int>& rigid_bodies_in_level=rigid_body_collisions.contact_graph.directed_graph.Nodes_In_Level(level);
                int level_iteration=0;bool need_more_level_iterations=true;T move_fraction=1;
                while(need_more_level_iterations && ++level_iteration<=rigid_body_collisions.push_out_level_iterations){need_more_level_iterations=false;
                    if(rigid_body_collisions.use_gradual_push_out) move_fraction=(T)level_iteration/rigid_body_collisions.push_out_level_iterations;
                    for(int i=1;i<=rigid_bodies_in_level.m;i++){
                        bool need_more_pair_iterations=true;int pair_iteration=0;
                        while(need_more_pair_iterations && ++pair_iteration<=rigid_body_collisions.push_out_pair_iterations){
                            need_more_pair_iterations=false;
                            if(Push_Out_From_Rigid_Body(rigid_body_collection.Rigid_Body(rigid_bodies_in_level(i)),particle_intersections,move_fraction)){
                                need_more_pair_iterations=need_more_level_iterations=need_another_iteration=true;
                                rigid_body_collisions.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Positions();}}}}
                if(rigid_body_collisions.use_freezing_with_push_out) rigid_body_collisions.Set_Level_Temporarily_Static(level);}
            if(rigid_body_collisions.use_freezing_with_push_out) rigid_body_collisions.Clear_Temporarily_Static();}

        // Push out from particles only if there are tetrahedron collision bodies; if only rigid bodies and particles, just process around rigid bodies
        if(tetrahedron_collision_bodies.m){
            Compute_Particle_Candidates_Helper(*this,tetrahedron_collision_bodies,rigid_collision_bodies);
            for(typename HASHTABLE<int,ARRAY<TETRAHEDRON_COLLISION_BODY<T>*> >::ITERATOR iterator(particle_tetrahedron_candidates);iterator.Valid();iterator.Next())
                if(Push_Out_From_Particle(iterator.Key()) && !kinematic_rigid_bodies_only){
                    need_another_iteration=true;rigid_body_collisions.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Positions();}}}

    // Process static/kinematic rigid bodies
    if(kinematic_rigid_bodies_only && solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions)
        for(int i=1;i<=solid_body_collection.rigid_body_collection.static_and_kinematic_rigid_bodies.m;i++){
            RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(solid_body_collection.rigid_body_collection.static_and_kinematic_rigid_bodies(i));
            if(Push_Out_From_Rigid_Body(rigid_body,particle_intersections,(T)1)) rigid_body_collisions.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Positions();}

    // TODO: handle soft bound particles
    //if(interactions) soft_bindings.Adjust_Parents_For_Changes_In_Surface_Children(particle_on_surface);

    solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Positions();
    solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities(); // TODO: necessary?
    solid_body_collection.deformable_body_collection.soft_bindings.Clamp_Particles_To_Embedded_Positions(true);
    solid_body_collection.deformable_body_collection.soft_bindings.Clamp_Particles_To_Embedded_Velocities(true);
}
//#####################################################################
template class RIGID_DEFORMABLE_COLLISIONS<VECTOR<float,1> >;
template class RIGID_DEFORMABLE_COLLISIONS<VECTOR<float,2> >;
template class RIGID_DEFORMABLE_COLLISIONS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_DEFORMABLE_COLLISIONS<VECTOR<double,1> >;
template class RIGID_DEFORMABLE_COLLISIONS<VECTOR<double,2> >;
template class RIGID_DEFORMABLE_COLLISIONS<VECTOR<double,3> >;
#endif

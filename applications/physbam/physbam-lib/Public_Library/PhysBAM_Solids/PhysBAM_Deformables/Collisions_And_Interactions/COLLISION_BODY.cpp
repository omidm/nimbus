//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Math_Tools/sign.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_PARTICLE_STATE.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COLLISION_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Function Adjust_Nodes_For_Collisions
//#####################################################################
template<class TV> int COLLISION_BODY_HELPER<TV>::
Adjust_Nodes_For_Collisions(COLLISION_GEOMETRY<TV>& body,ARRAY_VIEW<const TV> X_old,PARTICLES<TV>& collision_particles,SOFT_BINDINGS<TV>& soft_bindings,
    const ARRAY<int>& nodes_to_check,const ARRAY<bool>& particle_on_surface,const T collision_tolerance,ARRAY<COLLISION_PARTICLE_STATE<TV> >& collision_particle_state,
    ARRAY<COLLISION_GEOMETRY_ID>& particle_to_collision_geometry_id,const T max_relative_velocity,const T dt,const HASHTABLE<int,T> *friction_table,const HASHTABLE<int,T> *thickness_table)
{
    int interactions=0;T depth,one_over_dt=1/dt;ARRAY_VIEW<TV> X(collision_particles.X),V(collision_particles.V);
    for(int pp=1;pp<=nodes_to_check.m;pp++){
        int p=nodes_to_check(pp),soft_binding=soft_bindings.Soft_Binding(p);
        T thickness=thickness_table?thickness_table->Get_Default(p,0):0;
        COLLISION_PARTICLE_STATE<TV>& collision=collision_particle_state(p);
        if(soft_binding && soft_bindings.use_impulses_for_collisions(soft_binding)){
            int parent=soft_bindings.bindings(soft_binding).y;
            BINDING<TV>* hard_binding=soft_bindings.Hard_Binding(p);
            if(hard_binding) hard_binding->Clamp_To_Embedded_Position(); // sync hard binding position, allow (position only) drift in soft binding
            if(soft_bindings.use_gauss_seidel_for_impulse_based_collisions) X(p)=X(parent); // clamp soft binding position
            if(body.Implicit_Geometry_Lazy_Inside_And_Value(X(p),depth,thickness)){
                depth-=thickness;depth=max((T)0,-depth)+collision_tolerance;collision.enforce=true;interactions++;particle_to_collision_geometry_id(p)=body.collision_geometry_id;
                if(hard_binding) hard_binding->Clamp_To_Embedded_Velocity(); // sync hard binding velocity
                V(p)=V(parent); // sync soft binding velocity as well
                Adjust_Point_For_Collision(body,X_old(p),X(p),V(p),collision_particles.mass(p),depth,dt,one_over_dt,max_relative_velocity,collision,
                    friction_table?friction_table->Get_Default(p,-1):-1);
                // PASS IN ALL FALSE ARRAY FOR PARTICLE_ON_SURFACE (I.E., DON'T SKIP ANY NODES) TO BETTER PRESERVE MOMENTUM - BUT TETS ROLL STRANGELY
                if(hard_binding){
                    TV delta_X=X(p)-X(parent),delta_V=V(p)-V(parent);
                    hard_binding->Apply_Displacement_To_Parents_Based_On_Embedding(delta_X,&particle_on_surface); // skip parent particles on surface
                    hard_binding->Apply_Velocity_Change_To_Parents_Based_On_Embedding(delta_V,&particle_on_surface);}
                else if(!particle_on_surface(parent)){X(parent)=X(p);V(parent)=V(p);}}}
        else if(body.Implicit_Geometry_Lazy_Inside_And_Value(X(p),depth,thickness)){
            depth-=thickness;depth=max((T)0,-depth)+collision_tolerance;collision.enforce=true;interactions++;particle_to_collision_geometry_id(p)=body.collision_geometry_id;
            Adjust_Point_For_Collision(body,X_old(p),X(p),V(p),collision_particles.mass(p),depth,dt,one_over_dt,max_relative_velocity,collision,
                friction_table?friction_table->Get_Default(p,-1):-1);}}
    return interactions;
}
//#####################################################################
// Function Adjust_Point_For_Collision
//#####################################################################
template<class TV> void COLLISION_BODY_HELPER<TV>::
Adjust_Point_For_Collision(COLLISION_GEOMETRY<TV>& body,const TV& X_old,TV& X,TV& V,const T point_mass,const T penetration_depth,const T dt,const T one_over_dt,const T max_relative_velocity,
    COLLISION_PARTICLE_STATE<TV>& collision,T local_coefficient_of_friction)
{
    TV normal=body.Implicit_Geometry_Normal(X),V_body=body.Pointwise_Object_Velocity(X);T VN_body=TV::Dot_Product(V_body,normal);
    T clamped_penetration_depth=penetration_depth;
    if(penetration_depth*one_over_dt>max_relative_velocity) clamped_penetration_depth=dt*max_relative_velocity;
    TV V_half=one_over_dt*(X-X_old),V_half_relative=V_half-V_body;T VN_half_relative=TV::Dot_Product(V_half_relative,normal);
    T depth_old=clamped_penetration_depth+dt*VN_half_relative,theta;
    TV X_collision;
    if(depth_old<0){ // starts outside object, find fraction of timestep where collision occurs
        theta=LEVELSET_UTILITIES<T>::Theta(depth_old,clamped_penetration_depth);
        X_collision=(1-theta)*X_old+theta*X;}
    else{ // starts inside object, project outside
        theta=0;X_collision=X_old+depth_old*normal;}
    // TODO: maybe recompute normal at X_collision shifted to time n+1

    // compute for V velocities
    T VN=TV::Dot_Product(V,normal);TV VT=V-VN*normal;
    T VN_final=max(VN,VN_body);
    // compute for V_half velocities
    T VN_half=TV::Dot_Product(V_half,normal);TV VT_half=V_half-VN_half*normal;
    T VN_final_half=max(VN_half,VN_body);

    RIGID_GEOMETRY<TV>& rigid_geometry=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>&>(body).rigid_geometry;
    // Negative friction means look it up.
    if(local_coefficient_of_friction<0) local_coefficient_of_friction=rigid_geometry.coefficient_of_friction;
    if(local_coefficient_of_friction){
        // friction for V velocities
        {T normal_force=VN_final-VN;
        TV VT_body=V_body-VN_body*normal,VT_relative=VT-VT_body;T VT_relative_magnitude=VT_relative.Magnitude();
        T friction_magnitude=1;if(local_coefficient_of_friction*normal_force<VT_relative_magnitude) friction_magnitude=local_coefficient_of_friction*normal_force/VT_relative_magnitude;
        VT=VT_body+VT_relative*(1-friction_magnitude);}
        // friction for V_half velocities
        {T normal_force=VN_final_half-VN_half;
        TV VT_body=V_body-VN_body*normal,VT_relative=VT_half-VT_body;T VT_relative_magnitude=VT_relative.Magnitude();
        T friction_magnitude=1;if(local_coefficient_of_friction*normal_force<VT_relative_magnitude) friction_magnitude=local_coefficient_of_friction*normal_force/VT_relative_magnitude;
        VT_half=VT_body+VT_relative*(1-friction_magnitude);}}
    // compute final V and V_half velocities
    TV V_save=V;
    V=VN_final*normal+VT;
    V_half=VN_final_half*normal+VT_half;
    // impulse accumulation
    if(body.impulse_accumulator){
        TV impulse_on_particle=point_mass*(V-V_save);
        body.impulse_accumulator->Add_Impulse(X,-TWIST<TV>(impulse_on_particle,TV::Cross_Product(X-rigid_geometry.X(),impulse_on_particle)));}
    // update position
    X=X_collision+(1-theta)*dt*(VN_body*normal+VT_half);
    // set collision state
    collision.normal=body.Implicit_Geometry_Normal(X);collision.VT_body=body.Pointwise_Object_Velocity(X).Projected_Orthogonal_To_Unit_Direction(collision.normal);
    collision.VN=TV::Dot_Product(V_half,collision.normal);collision.friction=local_coefficient_of_friction;
}
//#####################################################################
// Function Adjust_Point_For_Collision
//#####################################################################
template<class TV> void COLLISION_BODY_HELPER<TV>::
Adjust_Point_For_Collision(COLLISION_GEOMETRY<TV>& body,TV& X,const T penetration_depth)
{
    X+=penetration_depth*body.Implicit_Geometry_Normal(X); // push the point out of the rigid body
}
//#####################################################################
template class COLLISION_BODY_HELPER<VECTOR<float,1> >;
template class COLLISION_BODY_HELPER<VECTOR<float,2> >;
template class COLLISION_BODY_HELPER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COLLISION_BODY_HELPER<VECTOR<double,1> >;
template class COLLISION_BODY_HELPER<VECTOR<double,2> >;
template class COLLISION_BODY_HELPER<VECTOR<double,3> >;
#endif

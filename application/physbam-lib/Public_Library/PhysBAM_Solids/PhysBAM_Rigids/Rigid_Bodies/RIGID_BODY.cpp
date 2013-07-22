//#####################################################################
// Copyright 2003-2009, Zhaosheng Bao, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Frank Losasso, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/EXTERNAL_ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Arrays/PROJECTED_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Vectors/Dot_Product.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Computations/LEVELSET_DYADIC_SIGNED_DISTANCE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/TRIANGULATED_SURFACE_SIGNED_DISTANCE_UNIFORM.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Dyadic/DYADIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Particles/RIGIDS_PARTICLES_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_BODY<TV>::
RIGID_BODY(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,bool create_collision_geometry)
    :BASE(rigid_body_collection_input.rigid_geometry_collection,create_collision_geometry,0),rigid_body_collection(rigid_body_collection_input),
    is_temporarily_static(false),fracture_threshold(FLT_MAX),thin_shell(false),CFL_initialized(false)
{
    EXTERNAL_ARRAY_COLLECTION* array_collection=dynamic_cast<EXTERNAL_ARRAY_COLLECTION*>(rigid_body_collection.rigid_body_particle.array_collection);
    if(!array_collection || array_collection->Owns_Element(ATTRIBUTE_ID_RIGID_MASS)) Set_Rigid_Mass(RIGID_BODY_MASS<TV>::Identity_Mass());
    Angular_Momentum()=T_SPIN();
    rigid_body_collection.rigid_body_particle.kinematic(particle_index)=false;
    Set_Coefficient_Of_Restitution();
    Set_Coefficient_Of_Rolling_Friction();
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_BODY<TV>::
RIGID_BODY(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,bool create_collision_geometry,int index)
    :BASE(rigid_body_collection_input.rigid_geometry_collection,create_collision_geometry,index),rigid_body_collection(rigid_body_collection_input),
    is_temporarily_static(false),fracture_threshold(FLT_MAX),thin_shell(false),CFL_initialized(false)
{
    EXTERNAL_ARRAY_COLLECTION* array_collection=dynamic_cast<EXTERNAL_ARRAY_COLLECTION*>(rigid_body_collection.rigid_body_particle.array_collection);
    if(!array_collection || array_collection->Owns_Element(ATTRIBUTE_ID_RIGID_MASS)) Set_Rigid_Mass(RIGID_BODY_MASS<TV>::Identity_Mass());
    Angular_Momentum()=T_SPIN();
    rigid_body_collection.rigid_body_particle.kinematic(particle_index)=false;
    Set_Coefficient_Of_Restitution();
    Set_Coefficient_Of_Rolling_Friction();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_BODY<TV>::
~RIGID_BODY()
{
}
//#####################################################################
// Function Print_Pairs
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Print_Pairs(const RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const ARRAY<VECTOR<int,2> >& pairs)
{
    {std::stringstream ss;ss<<"{";LOG::filecout(ss.str());}
    for(int i=1;i<=pairs.m;i++){
        {std::stringstream ss;ss<<"(\""<<rigid_body_collection.Rigid_Body(pairs(i)(1)).name<<"\", \""<<rigid_body_collection.Rigid_Body(pairs(i)(2)).name<<"\")";LOG::filecout(ss.str());}
        if(i<pairs.m) {std::stringstream ss;ss<<", ";LOG::filecout(ss.str());}}
    {std::stringstream ss;ss<<"}";LOG::filecout(ss.str());}
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Initialize_CFL()
{
    // Assumes bounding box is up to date
    const RANGE<TV>& box=BASE::Object_Space_Bounding_Box();
    bounding_box_radius=TV::Componentwise_Max(abs(box.Minimum_Corner()),abs(box.Maximum_Corner())).Magnitude();
    CFL_initialized=true;
}
//#####################################################################
// Function CFL
//#####################################################################
// set max_distance_per_time_step and/or max_rotation_per_time_step to zero to ignore those checks
template<class TV> typename TV::SCALAR RIGID_BODY<TV>::
CFL(const T max_distance_per_time_step,const T max_rotation_per_time_step,const bool verbose)
{
    if(is_static) return FLT_MAX;
    if(!CFL_initialized) Initialize_CFL();

    T angular_speed=Twist().angular.Magnitude(),linear_speed=Twist().linear.Magnitude()+bounding_box_radius*angular_speed;
    T linear_dt=max_distance_per_time_step!=0?Robust_Divide(max_distance_per_time_step,linear_speed):FLT_MAX;
    //angular speed check avoids a divide by zero warning on Windows
    T angular_dt=(max_rotation_per_time_step!=0 && angular_speed!=0)?Robust_Divide(max_rotation_per_time_step,angular_speed):FLT_MAX;
    T dt=min(linear_dt,angular_dt);
    return dt;
}
//#####################################################################
// Function Interpolate_Between_States
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Interpolate_Between_States(const RIGID_BODY_STATE<TV>& state1,const RIGID_BODY_STATE<TV>& state2,const T time,RIGID_BODY_STATE<TV>& interpolated_state)
{
    T alpha=(time-state1.time)/(state2.time-state1.time);alpha=clamp(alpha,(T)0,(T)1);
    BASE::Interpolate_Between_States(state1,state2,time,interpolated_state);
    interpolated_state.angular_momentum=(1-alpha)*state1.angular_momentum+alpha*state2.angular_momentum;
    Update_Angular_Velocity(interpolated_state);
}
//#####################################################################
// Function Compute_Velocity_Between_States
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Compute_Velocity_Between_States(const RIGID_BODY_STATE<TV>& state1,const RIGID_BODY_STATE<TV>& state2,RIGID_BODY_STATE<TV>& result_state)
{
    BASE::Compute_Velocity_Between_States(state1,state2,result_state);
    Update_Angular_Momentum(result_state); // Assumes result_state has a valid orientation
}
//#####################################################################
// Function Apply_Impulse_To_Body
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Apply_Impulse_To_Body(const TV& location,const TV& impulse,const T_SPIN& angular_impulse,const bool half_impulse_for_accumulator)
{
    if(Has_Infinite_Inertia()) return;
    T_SPIN total_angular_impulse=TV::Cross_Product(location-X(),impulse)+angular_impulse;
    V()+=impulse/Mass();
    Angular_Momentum()+=total_angular_impulse;
    Update_Angular_Velocity();
    if(this->impulse_accumulator) this->impulse_accumulator->Add_Impulse(location,(half_impulse_for_accumulator?(T).5:(T)1)*TWIST<TV>(impulse,total_angular_impulse));
}
//#####################################################################
// Function Apply_Impulse
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Apply_Impulse(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const TV& location,const TV& impulse,const T_SPIN& angular_impulse,const bool half_impulse_for_accumulator)
{
    body1.Apply_Impulse_To_Body(location,impulse,angular_impulse,half_impulse_for_accumulator);
    body2.Apply_Impulse_To_Body(location,-impulse,-angular_impulse,half_impulse_for_accumulator);
}
//#####################################################################
// Function Apply_Clamped_Impulse
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Compute_Clamped_Impulse(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const TV& location,TWIST<TV>& impulse,const ROTATION<TV>& saved_rotation_1,const ROTATION<TV>& saved_rotation_2)
{
    RIGID_BODY<TV>* bodies[2]={&body1,&body2};
    const ROTATION<TV>* saved_rotation[2]={&saved_rotation_1,&saved_rotation_2};
    T_SPIN jr[2],I_inverse_jr[2];
    T impulse_ratio_num[2]={0},impulse_ratio_denom[2]={0},impulse_linear_magnitude_squared=impulse.linear.Magnitude_Squared();
    TV velocity_adjustment;
    for(int i=0;i<2;i++) if(bodies[i]->rigid_body_collection.rigid_body_particle.kinematic(bodies[i]->particle_index)) velocity_adjustment=bodies[i]->Pointwise_Object_Velocity(location);
    for(int i=0;i<2;i++) if(!bodies[i]->Has_Infinite_Inertia()){
        jr[i]=TV::Cross_Product(location-bodies[i]->X(),impulse.linear)+impulse.angular;
        I_inverse_jr[i]=bodies[i]->Rigid_Mass().World_Space_Inertia_Tensor_Inverse_Times(*saved_rotation[i],jr[i]);
        impulse_ratio_num[i]=Dot_Product(impulse.linear,bodies[i]->Twist().linear-velocity_adjustment)+Dot_Product(I_inverse_jr[i],bodies[i]->Angular_Momentum());
        impulse_ratio_denom[i]=impulse_linear_magnitude_squared/bodies[i]->Mass()+Dot_Product(jr[i],I_inverse_jr[i]);}
    T total_impulse_ratio_num=-2*(impulse_ratio_num[0]-impulse_ratio_num[1]);
    if(total_impulse_ratio_num>0){
        T total_impulse_ratio_denom=impulse_ratio_denom[0]+impulse_ratio_denom[1];
        T scale=total_impulse_ratio_num<total_impulse_ratio_denom?total_impulse_ratio_num/total_impulse_ratio_denom:1;
        impulse*=scale;}
    else impulse=TWIST<TV>();
}
//#####################################################################
// Function Compute_Collision_Impulse
//#####################################################################
// clamp friction magnitude: should be (true) in the elastic collision case, (false) in the inelastic collision case
template<class TV> TWIST<TV> RIGID_BODY<TV>::
Compute_Collision_Impulse(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const ROTATION<TV>& saved_rotation_1,const ROTATION<TV>& saved_rotation_2,const TV& location,const TV& normal,
    const TV& relative_velocity,const T coefficient_of_restitution,const T coefficient_of_friction,const bool clamp_friction_magnitude,const bool rolling_friction,const bool clamp_energy)
{
    TWIST<TV> impulse;
    if(body1.Has_Infinite_Inertia() && body2.Has_Infinite_Inertia()) return TWIST<TV>();
    if(!coefficient_of_friction){ // frictionless case
        T nT_impulse1_n=0;
        if(!body1.Has_Infinite_Inertia()){
            T_SPIN r1xn=TV::Cross_Product(location-body1.X(),normal);
            nT_impulse1_n=1/body1.Mass()+Dot_Product(r1xn,body1.World_Space_Inertia_Tensor_Inverse()*r1xn);}
        T nT_impulse2_n=0;
        if(!body2.Has_Infinite_Inertia()){
            T_SPIN r2xn=TV::Cross_Product(location-body2.X(),normal);
            nT_impulse2_n=1/body2.Mass()+Dot_Product(r2xn,body2.World_Space_Inertia_Tensor_Inverse()*r2xn);}
        impulse.linear=-(1+coefficient_of_restitution)*TV::Dot_Product(relative_velocity,normal)*normal/(nT_impulse1_n+nT_impulse2_n);}
    else{ // friction case
        T relative_normal_velocity=min((T)0,TV::Dot_Product(relative_velocity,normal));
        T_SYMMETRIC_MATRIX impulse_factor=Impulse_Factor(body1,body2,location),impulse_factor_inverse=impulse_factor.Inverse();
        // see if friction stops sliding
        TV sticking_impulse=impulse_factor_inverse*(-coefficient_of_restitution*relative_normal_velocity*normal-relative_velocity);
        T normal_component=TV::Dot_Product(sticking_impulse,normal);
        if((sticking_impulse-normal_component*normal).Magnitude()<=coefficient_of_friction*normal_component){
            impulse.linear=sticking_impulse;
            if(rolling_friction) impulse+=Apply_Rolling_Friction(body1,body2,location,normal,normal_component);}
        // friction does not stop sliding
        else{
            TV relative_tangential_velocity=relative_velocity.Projected_Orthogonal_To_Unit_Direction(normal);
            TV tangential_direction=relative_tangential_velocity;T relative_tangential_velocity_magnitude=tangential_direction.Normalize();
            TV impulse_factor_times_direction=impulse_factor*(normal-coefficient_of_friction*tangential_direction);
            PHYSBAM_ASSERT(TV::Dot_Product(impulse_factor_times_direction,normal));
            TV delta=-(1+coefficient_of_restitution)*relative_normal_velocity/TV::Dot_Product(impulse_factor_times_direction,normal)*impulse_factor_times_direction;
            if(clamp_friction_magnitude && coefficient_of_restitution){ // should only clamp friction magnitude in the elastic case!
                TV new_relative_velocity=relative_velocity+delta;
                TV new_normal_velocity=new_relative_velocity.Projected_On_Unit_Direction(normal);
                TV new_tangential_velocity=new_relative_velocity-new_normal_velocity;T new_tangential_velocity_magnitude=new_tangential_velocity.Magnitude();
                if(new_tangential_velocity_magnitude > relative_tangential_velocity_magnitude)
                    delta=new_normal_velocity+(relative_tangential_velocity_magnitude/new_tangential_velocity_magnitude)*new_tangential_velocity-relative_velocity;}
            impulse.linear=impulse_factor_inverse*delta;}}
    if(clamp_energy) Compute_Clamped_Impulse(body1,body2,location,impulse,saved_rotation_1,saved_rotation_2);
    return impulse;
}
//#####################################################################
// Function Apply_Collision_Impulse
//#####################################################################
// clamp friction magnitude: should be (true) in the elastic collision case, (false) in the inelastic collision case
template<class TV> void RIGID_BODY<TV>::
Apply_Collision_Impulse(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const ROTATION<TV>& saved_rotation_1,const ROTATION<TV>& saved_rotation_2,const TV& location,const TV& normal,
    const TV& relative_velocity,const T coefficient_of_restitution,const T coefficient_of_friction,const bool clamp_friction_magnitude,const bool rolling_friction,const bool clamp_energy,
    const bool half_impulse_for_accumulator)
{
    TWIST<TV> impulse=Compute_Collision_Impulse(body1,body2,saved_rotation_1,saved_rotation_2,location,normal,relative_velocity,coefficient_of_restitution,
        coefficient_of_friction,clamp_friction_magnitude,rolling_friction,clamp_energy);
    Apply_Impulse(body1,body2,location,impulse.linear,impulse.angular,half_impulse_for_accumulator);
}
//#####################################################################
// Function Apply_Sticking_And_Angular_Sticking_Impulse
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Apply_Sticking_And_Angular_Sticking_Impulse(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const TV& location,const TWIST<TV>& delta_relative_twist,const MATRIX_MXN<T>& angular_constraint_matrix,
    const MATRIX_MXN<T>& prismatic_constraint_matrix)
{
    //body1.Update_Angular_Velocity();body2.Update_Angular_Velocity();
    TWIST<TV> impulse;

    // faster version for fully constrained prismatic and angular impulse
    if(angular_constraint_matrix.Columns()==T_SPIN::dimension && prismatic_constraint_matrix.Columns()==TV::dimension)
        impulse=Find_Impulse_And_Angular_Impulse(body1,body2,location,delta_relative_twist);
    else if(angular_constraint_matrix.Columns()==0 && prismatic_constraint_matrix.Columns()==TV::dimension)
        impulse.linear=Impulse_Factor(body1,body2,location).Inverse()*delta_relative_twist.linear;
    else impulse=Find_Impulse_And_Angular_Impulse(body1,body2,location,delta_relative_twist,angular_constraint_matrix,prismatic_constraint_matrix);
    Apply_Impulse(body1,body2,location,impulse.linear,impulse.angular);
}
//#####################################################################
// Function Apply_Rolling_Friction
//#####################################################################
// location is a point in world space about which body is rolling
template<class T,class TV> static TWIST<TV> Apply_Rolling_Friction_Helper(RIGID_BODY<VECTOR<T,1> >& body1,RIGID_BODY<TV>& body2,const TV& location,const TV& normal,const T normal_impulse)
{
    return TWIST<TV>(); // TODO: implement
}
template<class T,class TV> static TWIST<TV> Apply_Rolling_Friction_Helper(RIGID_BODY<VECTOR<T,2> >& body1,RIGID_BODY<TV>& body2,const TV& location,const TV& normal,const T normal_impulse)
{
    return TWIST<TV>(); // TODO: implement
}
template<class T,class TV> static TWIST<TV> Apply_Rolling_Friction_Helper(RIGID_BODY<VECTOR<T,3> >& body1,RIGID_BODY<TV>& body2,const TV& location,const TV& normal,const T normal_impulse)
{
    PHYSBAM_ASSERT(!body1.Has_Infinite_Inertia() || !body2.Has_Infinite_Inertia());PHYSBAM_ASSERT(normal_impulse>=0);
    T coefficient_of_rolling_friction=RIGID_BODY<TV>::Coefficient_Of_Rolling_Friction(body1,body2);if(!coefficient_of_rolling_friction) return TWIST<TV>();
    body1.Update_Angular_Velocity();body2.Update_Angular_Velocity();
    TV relative_angular_velocity=RIGID_BODY<TV>::Relative_Angular_Velocity(body1,body2);
    T normal_component=TV::Dot_Product(relative_angular_velocity,normal),normal_magnitude=abs(normal_component);
    TV tangential_component=relative_angular_velocity-normal_component*normal;T tangential_magnitude=tangential_component.Magnitude();
    TV tangential_direction;if(tangential_magnitude!=0) tangential_direction=tangential_component/tangential_magnitude;
    normal_magnitude-=coefficient_of_rolling_friction*normal_impulse;tangential_magnitude-=coefficient_of_rolling_friction*normal_impulse;
    TV new_relative_angular_velocity=sign(normal_component)*max((T)0,normal_magnitude)*normal+max((T)0,tangential_magnitude)*tangential_direction;
    return RIGID_BODY<TV>::Find_Impulse_And_Angular_Impulse(body1,body2,location,TWIST<TV>(TV(),new_relative_angular_velocity-relative_angular_velocity));
}
template<class TV> TWIST<TV> RIGID_BODY<TV>::
Apply_Rolling_Friction(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const TV& location,const TV& normal,const T normal_impulse)
{
    return Apply_Rolling_Friction_Helper(body1,body2,location,normal,normal_impulse);
}
//#####################################################################
// Function Find_Impulse_And_Angular_Impulse
//#####################################################################
template<class T,class TV> TWIST<TV>
Find_Impulse_And_Angular_Impulse_Helper(const RIGID_BODY<VECTOR<T,1> >& body1,const RIGID_BODY<VECTOR<T,1> >& body2,const TV& location,const TWIST<TV>& delta_relative_twist_at_location,
    const MATRIX_MXN<T>& angular_constraint_matrix,const MATRIX_MXN<T>& prismatic_constraint_matrix)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
template<class T,class TV> typename ENABLE_IF<(TV::m>1),TWIST<TV> >::TYPE
Find_Impulse_And_Angular_Impulse_Helper(const RIGID_BODY<TV>& body1,const RIGID_BODY<TV>& body2,const TV& location,const TWIST<TV>& delta_relative_twist_at_location,
    const MATRIX_MXN<T>& angular_constraint_matrix,const MATRIX_MXN<T>& prismatic_constraint_matrix)
{
    // compute blocks of constrained matrix
    typedef typename RIGID_BODY_POLICY<TV>::WORLD_SPACE_INERTIA_TENSOR T_WORLD_SPACE_INERTIA_TENSOR;
    T_WORLD_SPACE_INERTIA_TENSOR I_inverse_1,I_inverse_2;T m1_inv_plus_m2_inv=0;
    if(!body1.Has_Infinite_Inertia()){I_inverse_1=body1.World_Space_Inertia_Tensor_Inverse();m1_inv_plus_m2_inv=1/body1.Mass();}
    if(!body2.Has_Infinite_Inertia()){I_inverse_2=body2.World_Space_Inertia_Tensor_Inverse();m1_inv_plus_m2_inv+=1/body2.Mass();}

    // fill in NXN constrained matrix C
    TV r1=location-body1.X(),r2=location-body2.X();
    int angular_constrained_axes=angular_constraint_matrix.Columns(),prismatic_constrained_axes=prismatic_constraint_matrix.Columns();
    MATRIX_MXN<T> r_cross_P_1=prismatic_constraint_matrix.Cross_Product_Matrix_Times(r1),r_cross_P_2=prismatic_constraint_matrix.Cross_Product_Matrix_Times(r2);
    MATRIX_MXN<T> P_T_r_cross_T_I_inverse_1=r_cross_P_1.Transpose_Times(I_inverse_1),P_T_r_cross_T_I_inverse_2=r_cross_P_2.Transpose_Times(I_inverse_2);
    MATRIX_MXN<T> C(prismatic_constrained_axes+angular_constrained_axes);
    C.Set_Submatrix(1,1,P_T_r_cross_T_I_inverse_1*r_cross_P_1+P_T_r_cross_T_I_inverse_2*r_cross_P_2+m1_inv_plus_m2_inv*prismatic_constraint_matrix.Transpose_Times(prismatic_constraint_matrix));
    if(angular_constrained_axes){
        C.Set_Submatrix(prismatic_constrained_axes+1,prismatic_constrained_axes+1,angular_constraint_matrix.Transpose_Times((I_inverse_1+I_inverse_2)*angular_constraint_matrix));
        MATRIX_MXN<T> c12=(P_T_r_cross_T_I_inverse_1+P_T_r_cross_T_I_inverse_2)*angular_constraint_matrix;
        C.Set_Submatrix(1,prismatic_constrained_axes+1,c12);C.Set_Submatrix(prismatic_constrained_axes+1,1,c12.Transposed());}

    VECTOR_ND<T> b(prismatic_constrained_axes+angular_constrained_axes);
    b.Set_Subvector(1,prismatic_constraint_matrix.Transpose_Times(delta_relative_twist_at_location.linear));
    if(angular_constrained_axes) b.Set_Subvector(prismatic_constrained_axes+1,angular_constraint_matrix.Transpose_Times(delta_relative_twist_at_location.angular));

    TWIST<TV> impulse;VECTOR_ND<T> x=C.Cholesky_Solve(b),linear(prismatic_constrained_axes),angular(angular_constrained_axes);
    x.Get_Subvector(1,linear);impulse.linear=prismatic_constraint_matrix*linear;
    if(angular_constrained_axes){x.Get_Subvector(prismatic_constrained_axes+1,angular);impulse.angular=angular_constraint_matrix*angular;}
    return impulse;
}
template<class TV> TWIST<TV> RIGID_BODY<TV>::
Find_Impulse_And_Angular_Impulse(const RIGID_BODY<TV>& body1,const RIGID_BODY<TV>& body2,const TV& location,const TWIST<TV>& delta_relative_twist_at_location,
    const MATRIX_MXN<T>& angular_constraint_matrix,const MATRIX_MXN<T>& prismatic_constraint_matrix)
{
    PHYSBAM_ASSERT(!body1.Has_Infinite_Inertia() || !body2.Has_Infinite_Inertia());
    return Find_Impulse_And_Angular_Impulse_Helper(body1,body2,location,delta_relative_twist_at_location,angular_constraint_matrix,prismatic_constraint_matrix);
}
//#####################################################################
// Function Find_Impulse_And_Angular_Impulse
//#####################################################################
template<class T,class TV> static TWIST<TV> Find_Impulse_And_Angular_Impulse_Helper(const RIGID_BODY<VECTOR<T,1> >& body1,const RIGID_BODY<TV>& body2,const TV& location,
    const TWIST<TV>& delta_relative_twist_at_location)
{
    return delta_relative_twist_at_location;
}
template<class T,class TV> static TWIST<TV> Find_Impulse_And_Angular_Impulse_Helper(const RIGID_BODY<VECTOR<T,2> >& body1,const RIGID_BODY<TV>& body2,const TV& location,
    const TWIST<TV>& delta_relative_twist_at_location)
{
    MATRIX<T,1> I_inverse_1=body1.World_Space_Inertia_Tensor_Inverse(),I_inverse_2=body2.World_Space_Inertia_Tensor_Inverse();
    SYMMETRIC_MATRIX<T,2> c11=I_inverse_1.Conjugate_With_Cross_Product_Matrix(location-body1.X())+
        I_inverse_2.Conjugate_With_Cross_Product_Matrix(location-body2.X())+1/body1.Mass()+1/body2.Mass();
    MATRIX<T,1,2> c12=I_inverse_1.Times_Cross_Product_Matrix(location-body1.X())+I_inverse_2.Times_Cross_Product_Matrix(location-body2.X());
    MATRIX<T,1> c22=I_inverse_1+I_inverse_2;
    SYMMETRIC_MATRIX<T,3> A(c11.x11,c11.x21,c12(1,1),c11.x22,c12(1,2),c22.x11);

    VECTOR<T,3> b(delta_relative_twist_at_location.linear.x,delta_relative_twist_at_location.linear.y,delta_relative_twist_at_location.angular.x);
    VECTOR<T,3> all_impulses=A.Inverse()*b;
    return TWIST<TV>(VECTOR<T,2>(all_impulses.x,all_impulses.y),VECTOR<T,1>(all_impulses.z));
}
template<class T,class TV> static TWIST<TV> Find_Impulse_And_Angular_Impulse_Helper(const RIGID_BODY<VECTOR<T,3> >& body1,const RIGID_BODY<TV>& body2,const TV& location,
    const TWIST<TV>& delta_relative_twist_at_location)
{
    SYMMETRIC_MATRIX<T,3> I_inverse_1=body1.World_Space_Inertia_Tensor_Inverse(),I_inverse_2=body2.World_Space_Inertia_Tensor_Inverse();
    MATRIX<T,3> r_cross_1=MATRIX<T,3>::Cross_Product_Matrix(location-body1.X()),r_cross_I_inverse_1=r_cross_1*I_inverse_1,
                r_cross_2=MATRIX<T,3>::Cross_Product_Matrix(location-body2.X()),r_cross_I_inverse_2=r_cross_2*I_inverse_2;
    SYMMETRIC_MATRIX<T,3> c11=SYMMETRIC_MATRIX<T,3>::Times_Transpose_With_Symmetric_Result(r_cross_I_inverse_1,r_cross_1)+
                              SYMMETRIC_MATRIX<T,3>::Times_Transpose_With_Symmetric_Result(r_cross_I_inverse_2,r_cross_2)+1/body1.Mass()+1/body2.Mass();
    MATRIX<T,3> c12=-(r_cross_I_inverse_1+r_cross_I_inverse_2);
    SYMMETRIC_MATRIX<T,3> c22=I_inverse_1+I_inverse_2,c22_inverse=c22.Inverse();
    MATRIX<T,3> c12_c22_inverse=c12*c22_inverse;
    SYMMETRIC_MATRIX<T,3> A=c11-SYMMETRIC_MATRIX<T,3>::Times_Transpose_With_Symmetric_Result(c12_c22_inverse,c12);
    TV b=delta_relative_twist_at_location.linear-c12_c22_inverse*delta_relative_twist_at_location.angular;
    TV linear_impulse=A.Inverse()*b;
    TV angular_impulse=c22_inverse*(delta_relative_twist_at_location.angular-c12.Transpose_Times(linear_impulse));
    return TWIST<TV>(linear_impulse,angular_impulse);
}
template<class TV> TWIST<TV> RIGID_BODY<TV>::
Find_Impulse_And_Angular_Impulse(const RIGID_BODY<TV>& body1,const RIGID_BODY<TV>& body2,const TV& location,const TWIST<TV>& delta_relative_twist_at_location)
{
    PHYSBAM_ASSERT(!body1.Has_Infinite_Inertia() || !body2.Has_Infinite_Inertia());

    if(body1.Has_Infinite_Inertia() || body2.Has_Infinite_Inertia()){
        const RIGID_BODY<TV>& body=(!body1.Has_Infinite_Inertia())?body1:body2;TV r(location-body.X());
        TV impulse_linear=body.Mass()*(delta_relative_twist_at_location.linear+TV::Cross_Product(r,delta_relative_twist_at_location.angular));
        T_SPIN impulse_angular=body.World_Space_Inertia_Tensor()*delta_relative_twist_at_location.angular-TV::Cross_Product(r,impulse_linear);
        return TWIST<TV>(impulse_linear,impulse_angular);}

    // c11*impulse.linear+c12*impulse.angular=delta_relative_twist_at_location.linear
    // c21*impulse.linear+c22*impulse.angular=delta_relative_twist_at_location.angular
    // Note: c21=c12^T
    return Find_Impulse_And_Angular_Impulse_Helper(body1,body2,location,delta_relative_twist_at_location);
}
//#####################################################################
// Function Apply_Push
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Apply_Push(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const TV& location,const TV& normal,const T distance,const bool half_impulse_for_accumulator)
{
    if(body1.Has_Infinite_Inertia() && body2.Has_Infinite_Inertia()) return;
    TV impulse=Impulse_Factor(body1,body2,location).Inverse()*(distance*normal);
    RIGID_BODY<TV>* body[]={&body1,&body2};
    for(int i=0;i<=1;i++) if(!body[i]->Has_Infinite_Inertia()){
        T sign=(T)1-2*i;
        TV velocity=impulse/(body[i]->Mass()*sign);
        T_SPIN angular_velocity=body[i]->World_Space_Inertia_Tensor_Inverse()*TV::Cross_Product(location-body[i]->X(),sign*impulse);
        body[i]->X()+=velocity;
        body[i]->Rotation()=ROTATION<TV>::From_Rotation_Vector(angular_velocity)*body[i]->Rotation();body[i]->Rotation().Normalize();
        body[i]->Update_Angular_Velocity();
        if(body[i]->impulse_accumulator) body[i]->impulse_accumulator->Add_Impulse(location,(half_impulse_for_accumulator?(T).5:(T)1)*TWIST<TV>(velocity,angular_velocity));} //NOTE: These are not impulses but thats ok as we just want to keep track of the changes
}
//#####################################################################
// Function Volume
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_BODY<TV>::
Volume() const
{
    PHYSBAM_ASSERT(simplicial_object);return simplicial_object->Volumetric_Volume();
}
//#####################################################################
// Function Volumetric_Density
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_BODY<TV>::
Volumetric_Density() const
{
    PHYSBAM_ASSERT(simplicial_object);return Mass()/Volume();
}
//#####################################################################
// Function Diagonalize_Inertia_Tensor
//#####################################################################
template<class TV,class T_INERTIA> void
Diagonalize_Inertia_Tensor_Helper(const T_INERTIA& inertia_input,T_INERTIA& inertia,ROTATION<TV>& rotation)
{
    inertia=inertia_input;
}
//#####################################################################
// Function Diagonalize_Inertia_Tensor
//#####################################################################
template<class T,class T_WORLD_INERTIA,class T_INERTIA> void
Diagonalize_Inertia_Tensor_Helper(const T_WORLD_INERTIA& inertia_input,T_INERTIA& inertia,ROTATION<VECTOR<T,3> >& rotation)
{
    MATRIX<T,3> rotation_matrix;
    inertia_input.Solve_Eigenproblem(inertia,rotation_matrix);
    rotation=ROTATION<VECTOR<T,3> >(rotation_matrix);
}
//#####################################################################
// Function Diagonalize_Inertia_Tensor
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Diagonalize_Inertia_Tensor(const T_WORLD_SPACE_INERTIA_TENSOR& inertia_tensor_at_center_of_mass)
{
    Diagonalize_Inertia_Tensor_Helper(inertia_tensor_at_center_of_mass,Inertia_Tensor(),Rotation());
    Update_Angular_Velocity();
}
//#####################################################################
// Function Initialize_From_Tetrahedralized_Volume_And_Triangulated_Surface
//#####################################################################
template<class TV> template<class T2> void RIGID_BODY<TV>::
Initialize_From_Tetrahedralized_Volume_And_Triangulated_Surface(TETRAHEDRALIZED_VOLUME<T2>& tetrahedralized_volume,TRIANGULATED_SURFACE<T>& triangulated_surface,const T cell_size,
    const int subdivision_loops,const bool (*create_levelset_test)(TETRAHEDRALIZED_VOLUME<T>&),const bool use_implicit_surface_maker,const int levels_of_octree,const T shrink_levelset_amount)
{
    Add_Structure(tetrahedralized_volume);
    TRIANGULATED_SURFACE<T>* triangulated_surface_condensed=TRIANGULATED_SURFACE<T>::Create();
    triangulated_surface_condensed->mesh.Initialize_Mesh(triangulated_surface.mesh);
    triangulated_surface_condensed->particles.array_collection->Initialize(*triangulated_surface.particles.array_collection); 
    triangulated_surface_condensed->Discard_Valence_Zero_Particles_And_Renumber();
    for(int i=1;i<=subdivision_loops;i++) triangulated_surface_condensed->Linearly_Subdivide();
    Add_Structure(*triangulated_surface_condensed);
    if(!create_levelset_test || (*create_levelset_test)(tetrahedralized_volume)){
        LEVELSET_IMPLICIT_OBJECT<TV>* levelset=LEVELSET_IMPLICIT_OBJECT<TV>::Create();
        levelset->levelset.grid=GRID<TV>::Create_Grid_Given_Cell_Size(*triangulated_surface_condensed->bounding_box,cell_size,false,2); 
        if(use_implicit_surface_maker){
            LEVELSET_MAKER_UNIFORM<T> levelset_maker;
            levelset_maker.Compute_Signed_Distance_Function(); 
            levelset_maker.Compute_Level_Set(*triangulated_surface_condensed,levelset->levelset.grid,levelset->levelset.phi);
            levelset->levelset.phi+=shrink_levelset_amount;}
        else SIGNED_DISTANCE::Calculate(*triangulated_surface_condensed,levelset->levelset.grid,levelset->levelset.phi);
        IMPLICIT_OBJECT<TV>* implicit_surface=0;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
        if(levels_of_octree){
            int m=levelset->levelset.grid.counts.x>>(levels_of_octree-1),n=levelset->levelset.grid.counts.y>>(levels_of_octree-1),mn=levelset->levelset.grid.counts.z>>(levels_of_octree-1);
            if(m%2==0)m++;if(n%2==0)n++;if(mn%2==0)mn++;
            if(m<=1) m=3;if(n<=1) n=3;if(mn<=1) mn=3;
            if(m>3&&n>3&&mn>3){ //create the levelset only if the body is sufficiently large
                DYADIC_IMPLICIT_OBJECT<TV>* octree_levelset=DYADIC_IMPLICIT_OBJECT<TV>::Create();
                VECTOR<T,3> margin=levelset->levelset.grid.domain.Edge_Lengths()-levelset->levelset.grid.dX.Min()*VECTOR<T,3>(levelset->levelset.grid.Domain_Indices().Edge_Lengths());
                RANGE<TV> box=levelset->levelset.grid.domain+RANGE<TV>(-margin,margin);
                octree_levelset->levelset.grid.Initialize(GRID<TV>(m,n,mn,box),levels_of_octree,4,true,false);
                SIGNED_DISTANCE::Calculate(levelset->levelset,octree_levelset->levelset.grid,octree_levelset->levelset.phi);
                implicit_surface=octree_levelset;}
            delete levelset;}
        else
#endif
            implicit_surface=levelset;
        if(implicit_surface) Add_Structure(*implicit_surface);}
}
//#####################################################################
// Function Effective_Inertia_Inverse
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Effective_Inertia_Inverse(MATRIX<T,dimension>& extended_mass_inverse,const TV& location) const
{
    if(Has_Infinite_Inertia()){extended_mass_inverse=MATRIX<T,TV::m+T_SPIN::m>();return;}

    TV r=location-X();
    MATRIX<T,T_SPIN::m> Ii=World_Space_Inertia_Tensor_Inverse();
    MATRIX<T,T_SPIN::m,TV::m> Ii_rs=Ii.Times_Cross_Product_Matrix(r);
    MATRIX<T,TV::m> mi_rst_Ii_rs=Ii_rs.Cross_Product_Matrix_Transpose_Times(r)+1/Mass();
    extended_mass_inverse.Set_Submatrix(1,1,mi_rst_Ii_rs);
    extended_mass_inverse.Set_Submatrix(1,1+TV::m,Ii_rs.Transposed());
    extended_mass_inverse.Set_Submatrix(1+TV::m,1,Ii_rs);
    extended_mass_inverse.Set_Submatrix(1+TV::m,1+TV::m,Ii);
}
//#####################################################################
// Function Effective_Inertia_At_Point
//#####################################################################
template<class TV> void RIGID_BODY<TV>::
Effective_Inertia(MATRIX<T,dimension>& extended_mass_inverse,const TV& location) const
{
    PHYSBAM_ASSERT(!Has_Infinite_Inertia());

    TV r=location-X();
    T m=Mass();
    MATRIX<T,TV::m> A=MATRIX<T,TV::m>()+m;
    MATRIX<T,T_SPIN::m,TV::m> B=-m*MATRIX<T,T_SPIN::m,TV::m>::Cross_Product_Matrix(r);
    MATRIX<T,T_SPIN::m> C=World_Space_Inertia_Tensor_Inverse()-B.Times_Cross_Product_Matrix_Transpose(r);
    extended_mass_inverse.Set_Submatrix(1,1,A);
    extended_mass_inverse.Set_Submatrix(1,1+TV::m,B.Transposed());
    extended_mass_inverse.Set_Submatrix(1+TV::m,1,B);
    extended_mass_inverse.Set_Submatrix(1+TV::m,1+TV::m,C);
}
//#####################################################################
// Function Effective_Inertia_Inverse_Times
//#####################################################################
template<class TV> TWIST<TV> RIGID_BODY<TV>::
Effective_Inertia_Inverse_Times(const TWIST<TV>& wrench,const TV& location) const
{
    if(Has_Infinite_Inertia()) return TWIST<TV>();
    return Scatter(Inertia_Inverse_Times(Gather(wrench,location)),location);
}
//#####################################################################
// Function Effective_Inertia_Times
//#####################################################################
template<class TV> TWIST<TV> RIGID_BODY<TV>::
Effective_Inertia_Times(const TWIST<TV>& twist,const TV& location) const
{
    PHYSBAM_ASSERT(!Has_Infinite_Inertia());

    TV r=location-X();
    T_SPIN torque=twist.angular-TV::Cross_Product(r,twist.linear);
    T_SPIN omega=World_Space_Inertia_Tensor_Times(torque);
    return TWIST<TV>(twist.linear*Mass()-TV::Cross_Product(omega,r),omega);
}
template<class TV> void RIGID_BODY<TV>::
Gather_Matrix(MATRIX<T,dimension>& gather,const TV& location) const
{
    gather.Set_Identity_Matrix();
    gather.Set_Submatrix(1+TV::m,1,MATRIX<T,T_SPIN::m,TV::m>::Cross_Product_Matrix(location-X()));
}
template<class TV> void RIGID_BODY<TV>::
Scatter_Matrix(MATRIX<T,dimension>& scatter,const TV& location) const
{
    scatter.Set_Identity_Matrix();
    scatter.Set_Submatrix(1,1+TV::m,MATRIX<T,T_SPIN::m,TV::m>::Cross_Product_Matrix(location-X()).Transposed());
}
//#####################################################################
template class RIGID_BODY<VECTOR<float,1> >;
template class RIGID_BODY<VECTOR<float,2> >;
template class RIGID_BODY<VECTOR<float,3> >;
template void RIGID_BODY<VECTOR<float,3> >::Initialize_From_Tetrahedralized_Volume_And_Triangulated_Surface<float>(TETRAHEDRALIZED_VOLUME<float>&,
    TRIANGULATED_SURFACE<float>&,float,int,bool const (*)(TETRAHEDRALIZED_VOLUME<float>&),bool,int,float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_BODY<VECTOR<double,1> >;
template class RIGID_BODY<VECTOR<double,2> >;
template class RIGID_BODY<VECTOR<double,3> >;
template void RIGID_BODY<VECTOR<double,3> >::Initialize_From_Tetrahedralized_Volume_And_Triangulated_Surface<double>(TETRAHEDRALIZED_VOLUME<double>&,
    TRIANGULATED_SURFACE<double>&,double,int,bool const (*)(TETRAHEDRALIZED_VOLUME<double>&),bool,int,double);
#endif
}

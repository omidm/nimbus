//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TETRAHEDRON_COLLISION_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> static TRIANGULATED_SURFACE<T>& Triangulated_Surface_Helper(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,TRIANGULATED_SURFACE<T>* triangulated_surface)
{
    if(triangulated_surface) return *triangulated_surface;
    if(!tetrahedralized_volume.triangulated_surface) tetrahedralized_volume.Initialize_Triangulated_Surface();
    return *tetrahedralized_volume.triangulated_surface;
}
template<class T> TETRAHEDRON_COLLISION_BODY<T>::
TETRAHEDRON_COLLISION_BODY(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume_input,TRIANGULATED_SURFACE<T>& undeformed_triangulated_surface_input,
    IMPLICIT_OBJECT<TV>& implicit_surface_input,TRIANGULATED_SURFACE<T>* triangulated_surface_input)
    :particles(dynamic_cast<PARTICLES<TV>&>(tetrahedralized_volume_input.particles)),undeformed_particles(dynamic_cast<PARTICLES<TV>&>(undeformed_triangulated_surface_input.particles)),tetrahedralized_volume(tetrahedralized_volume_input),
    undeformed_triangulated_surface(undeformed_triangulated_surface_input),triangulated_surface(Triangulated_Surface_Helper(tetrahedralized_volume_input,triangulated_surface_input)),
    implicit_surface(implicit_surface_input)
{
    Set_Max_Min_Barycentric_Weight_Tolerance();Set_Min_Tet_Volume_Tolerance();
    Set_Relaxation_Factor();Set_Friction_Coefficient();Set_Normal_Interpolation_Scale_Factor();Set_Self_Collision_Normal_Angle_Tolerance();
    tetrahedralized_volume.Initialize_Hierarchy();triangulated_surface.Initialize_Hierarchy();
    undeformed_triangulated_surface.Update_Triangle_List();undeformed_triangulated_surface.Initialize_Hierarchy();
    triangulated_surface.avoid_normal_interpolation_across_sharp_edges=false;
    if(&particles==&undeformed_particles) PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Surface_Normal
//#####################################################################
template<class T> inline VECTOR<T,3> TETRAHEDRON_COLLISION_BODY<T>::
Surface_Normal(const int triangle,const TV& weights) const
{
    int i,j,k;triangulated_surface.mesh.elements(triangle).Get(i,j,k);
    TV n1=(*triangulated_surface.vertex_normals)(i),n2=(*triangulated_surface.vertex_normals)(j),n3=(*triangulated_surface.vertex_normals)(k);
    return (weights.x*n1+weights.y*n2+weights.z*n3).Normalized();
}
//#####################################################################
// Function Implicit_Geometry_Lazy_Inside
//#####################################################################
template<class T> bool TETRAHEDRON_COLLISION_BODY<T>::
Implicit_Geometry_Lazy_Inside(const TV& location,T contour_value) const
{
    assert(contour_value<=0);
    TV tet_weights,tri_weights;
    int tet_nearest_point=Get_Tetrahedron_Near_Point(location,tet_weights);
    if(!tet_nearest_point) return false;
    int surface_triangle=Get_Surface_Triangle(tet_nearest_point,tet_weights,tri_weights,true);
    if(!surface_triangle) return false;else if(!contour_value) return true;
    int i,j,k;triangulated_surface.mesh.elements(surface_triangle).Get(i,j,k);
    TV surface_point=particles.X(i)*tri_weights.x+particles.X(j)*tri_weights.y+particles.X(k)*tri_weights.z;
    return (location-surface_point).Magnitude()+contour_value>=0;
}
//#####################################################################
// Function Implicit_Geometry_Lazy_Inside_And_Value
//#####################################################################
template<class T> bool TETRAHEDRON_COLLISION_BODY<T>::
Implicit_Geometry_Lazy_Inside_And_Value(const TV& location,T& phi,T contour_value) const
{
    assert(contour_value==0);
    TV tet_weights;int tet_nearest_point=Get_Tetrahedron_Near_Point(location,tet_weights);
    if(tet_nearest_point){
        TV tri_weights;int triangle_nearest_point=Get_Surface_Triangle(tet_nearest_point,tet_weights,tri_weights,true);
        if(triangle_nearest_point){
            int i,j,k;triangulated_surface.mesh.elements(triangle_nearest_point).Get(i,j,k);
            phi=-((tri_weights.x*particles.X(i)+tri_weights.y*particles.X(j)+tri_weights.z*particles.X(k)-location).Magnitude());
            return true;}}
    return false;
}
//#####################################################################
// Function Implicit_Geometry_Lazy_Outside_Extended_Levelset_And_Value
//#####################################################################
template<class T> bool TETRAHEDRON_COLLISION_BODY<T>::
Implicit_Geometry_Lazy_Outside_Extended_Levelset_And_Value(const TV& location,T& phi_value,T contour_value) const
{
    assert(contour_value==0);
    TRIANGLE_MESH& triangle_mesh=triangulated_surface.mesh;
    TV tet_weights;int tet_nearest_point=Get_Tetrahedron_Near_Point(location,tet_weights);
    if(tet_nearest_point){
        TV tri_weights;int triangle_nearest_point=Get_Surface_Triangle(tet_nearest_point,tet_weights,tri_weights,false,true);
        if(triangle_nearest_point){
            int i,j,k;triangle_mesh.elements(triangle_nearest_point).Get(i,j,k);
            phi_value=(tri_weights.x*particles.X(i)+tri_weights.y*particles.X(j)+tri_weights.z*particles.X(k)-location).Magnitude();
            return true;}}
    else{ // outside for certain
        ARRAY<int> intersection_list;triangulated_surface.hierarchy->Intersection_List(location,intersection_list,collision_thickness);
        T closest_distance=FLT_MAX;
        for(int t=1;t<=intersection_list.m;t++){
            int i,j,k,tri=intersection_list(t);triangle_mesh.elements(tri).Get(i,j,k);
            T distance=TRIANGLE_3D<T>(particles.X(i),particles.X(j),particles.X(k)).Distance_To_Triangle(location);if(distance<closest_distance) closest_distance=distance;}
        phi_value=closest_distance;return true;}
    return false;
}
//#####################################################################
// Function Implicit_Geometry_Normal
//#####################################################################
template<class T> VECTOR<T,3> TETRAHEDRON_COLLISION_BODY<T>::
Implicit_Geometry_Normal(const TV& location,const int aggregate) const
{
    TRIANGLE_MESH& triangle_mesh=triangulated_surface.mesh;
    TV tet_weights;int tet_nearest_point=Get_Tetrahedron_Near_Point(location,tet_weights);
    if(tet_nearest_point){
        bool inside;TV tri_weights;int triangle_nearest_point=Get_Surface_Triangle(tet_nearest_point,tet_weights,tri_weights,false,false,&inside);
        assert(triangle_nearest_point); // with false/false it should always find one
        int i,j,k;triangle_mesh.elements(triangle_nearest_point).Get(i,j,k);
        TV normal=tri_weights.x*particles.X(i)+tri_weights.y*particles.X(j)+tri_weights.z*particles.X(k)-location;
        T unsigned_phi=normal.Normalize();
        if(!inside) normal*=-1; // always an outward normal
        if(unsigned_phi>normal_interpolation_scale_factor) return normal;
        else return Depth_Interpolated_Normal(unsigned_phi,normal,Surface_Normal(triangle_nearest_point,tri_weights));}
    else{ // outside for certain
        ARRAY<int> intersection_list;triangulated_surface.hierarchy->Intersection_List(location,intersection_list,collision_thickness);
        T closest_distance_squared=FLT_MAX;TV closest_normal;int closest_triangle=0;TV closest_tri_weights;
        for(int t=1;t<=intersection_list.m;t++){
            int i,j,k,tri=intersection_list(t);triangle_mesh.elements(tri).Get(i,j,k);
            TV tri_weights,closest_point=TRIANGLE_3D<T>(particles.X(i),particles.X(j),particles.X(k)).Closest_Point(location,tri_weights),normal=location-closest_point;
            T distance_squared=normal.Magnitude_Squared();
            if(distance_squared<closest_distance_squared){closest_distance_squared=distance_squared;closest_normal=normal;closest_triangle=tri;closest_tri_weights=tri_weights;}}
        if(closest_distance_squared==FLT_MAX) return TV();
        T unsigned_phi=closest_normal.Normalize();
        if(unsigned_phi>normal_interpolation_scale_factor) return closest_normal;
        else return Depth_Interpolated_Normal(unsigned_phi,closest_normal,Surface_Normal(closest_triangle,closest_tri_weights));}
}
//#####################################################################
// Function Implicit_Geometry_Normal
//#####################################################################
template<class T> VECTOR<T,3> TETRAHEDRON_COLLISION_BODY<T>::
Implicit_Geometry_Normal(const TV& location,T& phi_value,const int aggregate,const int location_particle_index) const
{
    TRIANGLE_MESH& triangle_mesh=triangulated_surface.mesh;
    ARRAY<int> particles_to_ignore(1);particles_to_ignore(1)=location_particle_index;
    TV tet_weights;int tet_nearest_point=Get_Tetrahedron_Near_Point(location,tet_weights,particles_to_ignore);
    if(tet_nearest_point){
        bool inside;TV tri_weights;int triangle_nearest_point=Get_Surface_Triangle(tet_nearest_point,tet_weights,tri_weights,false,false,&inside);
        assert(triangle_nearest_point); // with false/false it should always find one
        int i,j,k;triangle_mesh.elements(triangle_nearest_point).Get(i,j,k);
        TV normal=tri_weights.x*particles.X(i)+tri_weights.y*particles.X(j)+tri_weights.z*particles.X(k)-location;
        T unsigned_phi=normal.Normalize();
        if(inside) phi_value=-unsigned_phi;else{phi_value=unsigned_phi;normal*=-1;} // always an outward normal
        if(unsigned_phi>normal_interpolation_scale_factor) return normal;
        else return Depth_Interpolated_Normal(unsigned_phi,normal,Surface_Normal(triangle_nearest_point,tri_weights));}
    else{ // outside for certain
        ARRAY<int> intersection_list;triangulated_surface.hierarchy->Intersection_List(location,intersection_list,collision_thickness);
        T closest_distance_squared=FLT_MAX;TV closest_normal;int closest_triangle=0;TV closest_tri_weights;
        for(int t=1;t<=intersection_list.m;t++){
            int i,j,k,tri=intersection_list(t);triangle_mesh.elements(tri).Get(i,j,k);
            TRIANGLE_3D<T> triangle(particles.X(i),particles.X(j),particles.X(k));
            TV tri_weights,closest_point=triangle.Closest_Point(location,tri_weights),normal=location-closest_point;
            if(location_particle_index && (location_particle_index==i || location_particle_index==j || location_particle_index==k || //TODO: fix for embedded
                TV::Dot_Product((*triangulated_surface.vertex_normals)(location_particle_index),triangle.normal)>=self_collision_normal_angle_tolerance)) continue;
            T distance_squared=normal.Magnitude_Squared();
            if(distance_squared<closest_distance_squared){closest_distance_squared=distance_squared;closest_normal=normal;closest_triangle=tri;closest_tri_weights=tri_weights;}}
        if(closest_distance_squared==FLT_MAX){phi_value=FLT_MAX;return TV();}
        T unsigned_phi=closest_normal.Normalize();phi_value=unsigned_phi;
        if(unsigned_phi>normal_interpolation_scale_factor) return closest_normal;
        else return Depth_Interpolated_Normal(unsigned_phi,closest_normal,Surface_Normal(closest_triangle,closest_tri_weights));}
}
//#####################################################################
// Function Implicit_Geometry_Extended_Normal
//#####################################################################
template<class T> VECTOR<T,3> TETRAHEDRON_COLLISION_BODY<T>::
Implicit_Geometry_Extended_Normal(const TV& location,T& phi_value,const int aggregate,const int location_particle_index) const
{
    return Implicit_Geometry_Normal(location,phi_value,aggregate,location_particle_index);
}
//#####################################################################
// Function Get_Tetrahedron_Near_Point
//#####################################################################
template<class T> int TETRAHEDRON_COLLISION_BODY<T>::
Get_Tetrahedron_Near_Point(const TV& point,TV& weights,const ARRAY<int>& particles_to_ignore) const
{
    TETRAHEDRON_MESH& tetrahedron_mesh=tetrahedralized_volume.mesh;
    ARRAY<int> intersection_list;tetrahedralized_volume.hierarchy->Intersection_List(point,intersection_list);
    T max_min_weight=-(T)FLT_MAX;int closest=0;
    for(int p=1;p<=intersection_list.m;p++){
        int t=intersection_list(p);VECTOR<int,4> nodes=tetrahedron_mesh.elements(t);int i,j,k,l;nodes.Get(i,j,k,l);
        for(int q=1;q<=particles_to_ignore.m;q++) if(nodes.Contains(particles_to_ignore(q))) goto CONTINUE;
        {TV w=TETRAHEDRON<T>::First_Three_Barycentric_Coordinates(point,particles.X(i),particles.X(j),particles.X(k),particles.X(l));T min_weight=min(w.x,w.y,w.z,1-w.x-w.y-w.z);
        if(min_weight>max_min_weight){
            weights=w;if(min_weight>0 && TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Signed_Volume()>min_tet_volume_tolerance) return t;
            max_min_weight=min_weight;closest=t;}}
        CONTINUE:;}
    if(closest && max_min_weight>max_min_barycentric_weight_tolerance){
        int i,j,k,l;tetrahedron_mesh.elements(closest).Get(i,j,k,l);
        if(TETRAHEDRON<T>(particles.X(i),particles.X(j),particles.X(k),particles.X(l)).Signed_Volume()>min_tet_volume_tolerance) return closest;}
    return 0;
}
//#####################################################################
// Function Get_Surface_Triangle
//#####################################################################
template<class T> int TETRAHEDRON_COLLISION_BODY<T>::
Get_Surface_Triangle(const int tetrahedron_index,const TV& tetrahedron_weights,TV& surface_weights,const bool omit_outside_points,const bool omit_inside_points,bool* inside) const
{
    int i,j,k,l;tetrahedralized_volume.mesh.elements(tetrahedron_index).Get(i,j,k,l);
    TV location=tetrahedron_weights.x*undeformed_particles.X(i)+tetrahedron_weights.y*undeformed_particles.X(j)
        +tetrahedron_weights.z*undeformed_particles.X(k)+(1-tetrahedron_weights.x-tetrahedron_weights.y-tetrahedron_weights.z)*undeformed_particles.X(l);
    TV surface_location=implicit_surface.Closest_Point_On_Boundary(location);ARRAY<int> nearby_surface_triangles;
    for(T thickness=collision_thickness;!nearby_surface_triangles.m;thickness*=(T)2)
        undeformed_triangulated_surface.hierarchy->Intersection_List(surface_location,nearby_surface_triangles,thickness);
    int closest_triangle=0;T closest_distance_squared=(T)FLT_MAX;TV closest_projected_point;
    for(int t=1;t<=nearby_surface_triangles.m;t++){
        TRIANGLE_3D<T>& surface_triangle=(*undeformed_triangulated_surface.triangle_list)(nearby_surface_triangles(t));
        TV weights,projected_point=surface_triangle.Closest_Point(surface_location,weights);
        T distance_squared=(surface_location-projected_point).Magnitude_Squared();
        if(distance_squared<closest_distance_squared){
            closest_distance_squared=distance_squared;closest_triangle=nearby_surface_triangles(t);surface_weights=weights;closest_projected_point=projected_point;}}
    if(inside || omit_outside_points || omit_inside_points){
        bool inside_temp=implicit_surface.Extended_Phi(location)<=implicit_surface.Extended_Phi(closest_projected_point);
        if(inside) *inside=inside_temp;
        if((omit_outside_points && !inside_temp) || (omit_inside_points && inside_temp)) return 0;}
    return closest_triangle;
}
//#####################################################################
// Function Adjust_Point_Face_Collision_Position_And_Velocity
//#####################################################################
template<class T> void TETRAHEDRON_COLLISION_BODY<T>::
Adjust_Point_Face_Collision_Position_And_Velocity(const int triangle_index,TV& X,TV& V,SOFT_BINDINGS<TV>& soft_bindings,const T one_over_mass,const T dt,const TV& weights,TV& position_change)
{
    // compute desired point adjustment
    int i,j,k;triangulated_surface.mesh.elements(triangle_index).Get(i,j,k);
    int ei=soft_bindings.Driftless_Particle(i),ej=soft_bindings.Driftless_Particle(j),ek=soft_bindings.Driftless_Particle(k);
    TV surface_point=weights.x*particles.X(ei)+weights.y*particles.X(ej)+weights.z*particles.X(ek);
    TV outward_direction=surface_point-X;position_change=(T).5*outward_direction;

    // apply point/face collision
    TV relative_velocity=V-(weights.x*particles.V(ei)+weights.y*particles.V(ej)+weights.z*particles.V(ek));
    TV N=TRIANGLE_3D<T>::Normal(particles.X(ei),particles.X(ej),particles.X(ek));
    T relative_speed=TV::Dot_Product(relative_velocity,N);
    if(relative_speed<0){
        T scalar_impulse=-relative_speed;TV impulse=-scalar_impulse*N;
        if(friction_coefficient){
            TV relative_tangential_velocity=relative_velocity-relative_speed*N;T relative_tangential_velocity_magnitude=relative_tangential_velocity.Magnitude();
            T friction_based_velocity_change=friction_coefficient*scalar_impulse;
            scalar_impulse=1;if(friction_based_velocity_change<relative_tangential_velocity_magnitude) scalar_impulse=friction_based_velocity_change/relative_tangential_velocity_magnitude;
            impulse+=scalar_impulse*relative_tangential_velocity;}
        T one_over_mass_1=particles.one_over_effective_mass(ei),one_over_mass_2=particles.one_over_effective_mass(ej),one_over_mass_3=particles.one_over_effective_mass(ek);
        impulse/=(one_over_mass+sqr(weights.x)*one_over_mass_1+sqr(weights.y)*one_over_mass_2+sqr(weights.z)*one_over_mass_3);
        V-=one_over_mass*impulse;particles.V(i)+=weights.x*one_over_mass_1*impulse;particles.V(j)+=weights.y*one_over_mass_2*impulse;particles.V(k)+=weights.z*one_over_mass_3*impulse;}
}
//#####################################################################
// Function Adjust_Nodes_For_Collisions
//#####################################################################
template<class T> int TETRAHEDRON_COLLISION_BODY<T>::
Adjust_Nodes_For_Collisions(ARRAY_VIEW<const TV> X_old,PARTICLES<TV>& collision_particles,SOFT_BINDINGS<TV>& soft_bindings,const ARRAY<int>& nodes_to_check,
    const ARRAY<bool>& particle_on_surface,const T collision_tolerance,ARRAY<COLLISION_PARTICLE_STATE<TV> >& collision_particle_state,
    ARRAY<COLLISION_GEOMETRY_ID>& particle_to_collision_body_id,const T max_relative_velocity,const T dt,const HASHTABLE<int,T> *friction_table,const HASHTABLE<int,T> *thickness_table)
{
    assert(!friction_table && !thickness_table);
    int interactions=0;ARRAY_VIEW<TV> X(collision_particles.X),V(collision_particles.V);
    tetrahedralized_volume.hierarchy->Update_Boxes(collision_thickness);
    ARRAY<VECTOR<int,2> > interaction_pair;ARRAY<TV> weights,position_change;
    for(int pp=1;pp<=nodes_to_check.m;pp++){int p=nodes_to_check(pp);
        ARRAY<int> particles_to_ignore;
        particles_to_ignore.Append(p);particles_to_ignore.Append_Elements(soft_bindings.Parents(p));
        TV w;int t=Get_Tetrahedron_Near_Point(X(p),w,particles_to_ignore);
        if(t){interaction_pair.Append(VECTOR<int,2>(p,t));weights.Append(w);}}
    for(int k=1;k<=interaction_pair.m;k++){
        int p,t;interaction_pair(k).Get(p,t);
        TV surface_weights;int surface_triangle=Get_Surface_Triangle(t,weights(k),surface_weights,true);
        if(surface_triangle){
            interactions++;interaction_pair(k)(2)=surface_triangle;TV change;
            // TODO: make this function take max_relative_velocity
            Adjust_Point_Face_Collision_Position_And_Velocity(surface_triangle,X(p),V(p),soft_bindings,collision_particles.one_over_effective_mass(p),dt,surface_weights,change);
            position_change.Append(change);}
        else interaction_pair(k)(2)=0;}
    int count=0;
    for(int pair=1;pair<=interaction_pair.m;pair++){
        int p,t;interaction_pair(pair).Get(p,t);
        if(t) X(p)+=relaxation_factor*position_change(++count);}
    if(interactions) soft_bindings.Adjust_Parents_For_Changes_In_Surface_Children(particle_on_surface);
    return interactions;
}
//#####################################################################
// Function Update_Bounding_Box
//#####################################################################
template<class T> void TETRAHEDRON_COLLISION_BODY<T>::
Update_Bounding_Box()
{
    tetrahedralized_volume.Update_Bounding_Box();
}
//#####################################################################
// Function Axis_Aligned_Bounding_Box
//#####################################################################
template<class T> const RANGE<VECTOR<T,3> >& TETRAHEDRON_COLLISION_BODY<T>::
Axis_Aligned_Bounding_Box() const
{
    return *tetrahedralized_volume.bounding_box;
}
//#####################################################################
// Function void Read_State
//#####################################################################
template<class T> void TETRAHEDRON_COLLISION_BODY<T>::
Read_State(TYPED_ISTREAM& input,const int state_index)
{PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function void Write_State
//#####################################################################
template<class T> void TETRAHEDRON_COLLISION_BODY<T>::
Write_State(TYPED_OSTREAM& output,const int state_index) const
{PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
template class TETRAHEDRON_COLLISION_BODY<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TETRAHEDRON_COLLISION_BODY<double>;
#endif

//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Avi Robinson-Mosher, Andrew Selle, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_OBJECT_FLUID_COLLISIONS
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Computations/RIGID_GEOMETRY_RASTERIZATION_DYADIC.h>
#include <PhysBAM_Geometry/Grids_RLE_Computations/RIGID_GEOMETRY_RASTERIZATION_RLE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/RIGID_GEOMETRY_RASTERIZATION_UNIFORM.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/SEGMENTED_CURVE_INSIDE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TETRAHEDRALIZED_VOLUME_INSIDE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGULATED_AREA_INSIDE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_POINT_SIMPLICES_1D_INTERSECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_SEGMENTED_CURVE_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_TRIANGULATED_SURFACE_INTERSECTION.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
namespace PhysBAM{
//##################################################################### 
// Function Initialize_For_Thin_Shells_Fluid_Coupling
//##################################################################### 
template<class T> static void Initialize_For_Thin_Shells_Fluid_Coupling_Helper(POINT_SIMPLICES_1D<T>& point_simplices)
{}
template<class T> static void Initialize_For_Thin_Shells_Fluid_Coupling_Helper(SEGMENTED_CURVE_2D<T>& segmented_curve)
{
    segmented_curve.Update_Bounding_Box();
    segmented_curve.Update_Segment_List();
    segmented_curve.Initialize_Hierarchy();
}
template<class T> static void Initialize_For_Thin_Shells_Fluid_Coupling_Helper(TRIANGULATED_SURFACE<T>& triangulated_surface)
{
    triangulated_surface.Update_Bounding_Box();
    triangulated_surface.Update_Triangle_List();
    triangulated_surface.Initialize_Hierarchy();
}
template<class TV> void DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
Initialize_For_Thin_Shells_Fluid_Coupling()
{
    Initialize_For_Thin_Shells_Fluid_Coupling_Helper(object);
}
//##################################################################### 
// Function Pointwise_Object_Pseudo_Velocity
//##################################################################### 
template<class TV> TV DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
Pointwise_Object_Pseudo_Velocity(const int simplex_id,const TV& location,const int state1,const int state2) const
{
    VECTOR<int,TV::dimension> nodes=object.mesh.elements(simplex_id);
    typename IF<TV::m==1,ONE,TV>::TYPE weights=T_SIMPLEX::Clamped_Barycentric_Coordinates(location,object.particles.X.Subset(nodes));
    TV dX=T_SIMPLEX::Point_From_Barycentric_Coordinates(weights,(saved_states(state2).x->X-saved_states(state1).x->X).Subset(nodes));
    return dX/(saved_states(state2).y-saved_states(state1).y);
}
//##################################################################### 
// Function Simplex_Crossover
//##################################################################### 
template<class T> POINT_SIMPLEX_COLLISION_TYPE Simplex_Crossover_Helper(const VECTOR<T,1>& start_X,const VECTOR<T,1>& end_X,const T dt,T& hit_time,ONE& weights,
    const T& collision_thickness,POINT_SIMPLEX_1D<T> initial_simplex,POINT_SIMPLEX_1D<T> final_simplex) {PHYSBAM_NOT_IMPLEMENTED();}
template<class T> POINT_SIMPLEX_COLLISION_TYPE Simplex_Crossover_Helper(const VECTOR<T,2>& start_X,const VECTOR<T,2>& end_X,const T dt,T& hit_time,T& weights,
    const T& collision_thickness,typename BASIC_GEOMETRY_POLICY<VECTOR<T,2> >::SEGMENT initial_simplex,typename BASIC_GEOMETRY_POLICY<VECTOR<T,2> >::SEGMENT final_simplex)
{
    typedef typename BASIC_GEOMETRY_POLICY<VECTOR<T,2> >::SEGMENT T_SIMPLEX;
    VECTOR<T,2> normal;T relative_speed; // unused
    return T_SIMPLEX::Robust_Point_Segment_Collision(initial_simplex,final_simplex,start_X,end_X,dt,collision_thickness,hit_time,normal,weights,relative_speed);
}
template<class T> POINT_SIMPLEX_COLLISION_TYPE Simplex_Crossover_Helper(const VECTOR<T,3>& start_X,const VECTOR<T,3>& end_X,const T dt,T& hit_time,VECTOR<T,3>& weights,
    const T& collision_thickness,typename BASIC_GEOMETRY_POLICY<VECTOR<T,3> >::TRIANGLE initial_simplex,typename BASIC_GEOMETRY_POLICY<VECTOR<T,3> >::TRIANGLE final_simplex)
{
    typedef typename BASIC_GEOMETRY_POLICY<VECTOR<T,3> >::TRIANGLE T_SIMPLEX;
    VECTOR<T,3> normal;T relative_speed; // unused
    return T_SIMPLEX::Robust_Point_Triangle_Collision(initial_simplex,final_simplex,start_X,end_X,dt,collision_thickness,hit_time,normal,weights,relative_speed);
}
template<class TV> POINT_SIMPLEX_COLLISION_TYPE DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,T_WEIGHTS& weights,const int simplex) const 
{
    T_SIMPLEX initial_simplex=World_Space_Simplex(simplex),final_simplex=World_Space_Simplex(simplex,*saved_states(1).x);
    return Simplex_Crossover_Helper(start_X,end_X,dt,hit_time,weights,collision_thickness,initial_simplex,final_simplex);
}
//##################################################################### 
// Function World_Space_Simplex
//##################################################################### 
template<class TV> typename DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::T_SIMPLEX DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
World_Space_Simplex(const int simplex_id,const bool use_saved_state) const
{
    if(use_saved_state){assert(saved_states(1).x);return World_Space_Simplex(simplex_id,*saved_states(1).x);}
    else return T_SIMPLEX(object.particles.X.Subset(object.mesh.elements(simplex_id)));
}
//##################################################################### 
// Function World_Space_Simplex
//##################################################################### 
template<class TV> typename DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::T_SIMPLEX DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
World_Space_Simplex(const int simplex_id,const GEOMETRY_PARTICLES<TV>& state) const
{
    return T_SIMPLEX(state.X.Subset(object.mesh.elements(simplex_id)));
}
//##################################################################### 
// Function Earliest_Simplex_Crossover
//##################################################################### 
template<class TV> bool DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
Earliest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,T_WEIGHTS& weights,int& simplex_id) const
{
    if(!saved_states(1).x){LOG::cerr<<"No saved_state"<<std::endl;PHYSBAM_FATAL_ERROR();}
    T min_time=FLT_MAX;T current_hit_time;T_WEIGHTS current_weights;
    if(object.hierarchy){
        ARRAY<int> triangles_to_check;
        if(start_X==end_X) object.hierarchy->Intersection_List(start_X,triangles_to_check,(T)5*collision_thickness);
        else object.hierarchy->Intersection_List(RANGE<TV>::Bounding_Box(start_X,end_X),triangles_to_check,(T)5*collision_thickness);
        for(int i=1;i<=triangles_to_check.m;i++){int t=triangles_to_check(i);
            POINT_SIMPLEX_COLLISION_TYPE collision_type=Simplex_Crossover(start_X,end_X,dt,current_hit_time,current_weights,t);
            if(collision_type!=POINT_SIMPLEX_NO_COLLISION && current_hit_time<min_time){
                min_time=hit_time=current_hit_time;weights=current_weights;simplex_id=t;}}}
    else{
        LOG::cerr<<"Earliest_Simplex_Crossover: no triangle hierarchy"<<std::endl;
        for(int t=1;t<=object.mesh.elements.m;t++){   
            POINT_SIMPLEX_COLLISION_TYPE collision_type=Simplex_Crossover(start_X,end_X,dt,current_hit_time,current_weights,t);
            if(collision_type!=POINT_SIMPLEX_NO_COLLISION && current_hit_time<min_time){
                min_time=hit_time=current_hit_time;weights=current_weights;simplex_id=t;}}}
    return min_time<FLT_MAX;
}
//##################################################################### 
// Function Latest_Simplex_Crossover
//##################################################################### 
template<class TV> bool DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
Latest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,T_WEIGHTS& weights,int& simplex_id,POINT_SIMPLEX_COLLISION_TYPE& returned_collision_type) const
{   
    if(!saved_states(1).x){LOG::cerr<<"No saved_state"<<std::endl;PHYSBAM_FATAL_ERROR();}
    returned_collision_type=POINT_SIMPLEX_NO_COLLISION;
    T max_time=-FLT_MAX;T current_hit_time;T_WEIGHTS current_weights;
    if(object.hierarchy){
        ARRAY<int> triangles_to_check;
        if(start_X==end_X) object.hierarchy->Intersection_List(start_X,triangles_to_check,(T)5*collision_thickness);
        else object.hierarchy->Intersection_List(RANGE<TV>::Bounding_Box(start_X,end_X),triangles_to_check,(T)5*collision_thickness);
        for(int i=1;i<=triangles_to_check.m;i++){int t=triangles_to_check(i);
            POINT_SIMPLEX_COLLISION_TYPE collision_type=Simplex_Crossover(start_X,end_X,dt,current_hit_time,current_weights,t);
            if(collision_type!=POINT_SIMPLEX_NO_COLLISION && current_hit_time>max_time){
                max_time=hit_time=current_hit_time;weights=current_weights;simplex_id=t;returned_collision_type=collision_type;}}}
    else{
        LOG::cerr<<"Latest_Simplex_Crossover: no triangle hierarchy"<<std::endl;
        for(int t=1;t<=object.mesh.elements.m;t++){   
            POINT_SIMPLEX_COLLISION_TYPE collision_type=Simplex_Crossover(start_X,end_X,dt,current_hit_time,current_weights,t);
            if(collision_type!=POINT_SIMPLEX_NO_COLLISION && current_hit_time>max_time){
                max_time=hit_time=current_hit_time;weights=current_weights;simplex_id=t;returned_collision_type=collision_type;}}}
    return max_time>-FLT_MAX;
}
//#####################################################################
// Function Any_Simplex_Crossover
//#####################################################################
template<class TV> bool DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
Any_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt) const
{
    if(!saved_states(1).x){LOG::cerr<<"No saved_state"<<std::endl;PHYSBAM_FATAL_ERROR();}
    T hit_time;T_WEIGHTS weights;
    if(object.hierarchy){
        ARRAY<int> triangles_to_check;
        if(start_X==end_X) object.hierarchy->Intersection_List(start_X,triangles_to_check,(T)5*collision_thickness);
        else object.hierarchy->Intersection_List(RANGE<TV>::Bounding_Box(start_X,end_X),triangles_to_check,(T)5*collision_thickness);
        for(int i=1;i<=triangles_to_check.m;i++){int t=triangles_to_check(i);
            POINT_SIMPLEX_COLLISION_TYPE collision_type=Simplex_Crossover(start_X,end_X,dt,hit_time,weights,t);
            if(collision_type!=POINT_SIMPLEX_NO_COLLISION) return true;}}
    else{
        LOG::cerr<<"Any_Simplex_Crossover: no triangle hierarchy"<<std::endl;
        for(int t=1;t<=object.mesh.elements.m;t++){   
            POINT_SIMPLEX_COLLISION_TYPE collision_type=Simplex_Crossover(start_X,end_X,dt,hit_time,weights,t);
            if(collision_type!=POINT_SIMPLEX_NO_COLLISION) return true;}}
    return false;
}
//#####################################################################
// Function Get_Simplex_Bounding_Boxes
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
Get_Simplex_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,const bool with_body_motion,const T extra_thickness,const T body_thickness_factor) const
{
    if(with_body_motion && !saved_states(1).x){LOG::cerr<<"No saved_state"<<std::endl;PHYSBAM_FATAL_ERROR();}
    for(int t=1;t<=object.mesh.elements.m;t++){
        VECTOR<int,TV::dimension> nodes=object.mesh.elements(t);
        RANGE<TV> box=RANGE<TV>::Bounding_Box(object.particles.X.Subset(nodes));
        if(with_body_motion) box.Enlarge_Nonempty_Box_To_Include_Points(saved_states(1).x->X.Subset(nodes));
        box.Change_Size(extra_thickness+body_thickness_factor*collision_thickness);
        bounding_boxes.Append(box);}
}
//##################################################################### 
// Function Update_Intersection_Acceleration_Structures
//##################################################################### 
template<class T> static void Update_Intersection_Acceleration_Structures_Helper(DEFORMABLE_OBJECT_FLUID_COLLISIONS<VECTOR<T,1> >& collisions,const bool,const int,const int)
{}
template<class T> static void Update_Intersection_Acceleration_Structures_Helper(DEFORMABLE_OBJECT_FLUID_COLLISIONS<VECTOR<T,2> >& collisions,const bool use_swept_hierarchy,
    const int state1,const int state2)
{
    collisions.object.Update_Segment_List();
    if(!collisions.object.hierarchy) collisions.object.Initialize_Hierarchy(false);
    if(use_swept_hierarchy) collisions.object.hierarchy->Update_Boxes(collisions.saved_states(state1).x->X,collisions.saved_states(state2).x->X);
    else collisions.object.hierarchy->Update_Boxes();
    collisions.object.Update_Bounding_Box();
    if(collisions.volume_object){
        if(!collisions.volume_object->hierarchy) collisions.volume_object->Initialize_Hierarchy(false);
        if(use_swept_hierarchy) collisions.volume_object->hierarchy->Update_Boxes(collisions.saved_states(state1).x->X,collisions.saved_states(state2).x->X);
        else collisions.volume_object->hierarchy->Update_Boxes();
        collisions.volume_object->Update_Bounding_Box();}
}
template<class T> static void Update_Intersection_Acceleration_Structures_Helper(DEFORMABLE_OBJECT_FLUID_COLLISIONS<VECTOR<T,3> >& collisions,const bool use_swept_hierarchy,
    const int state1,const int state2)
{
    collisions.object.Update_Triangle_List();
    if(!collisions.object.hierarchy) collisions.object.Initialize_Hierarchy(false);
    if(use_swept_hierarchy) collisions.object.hierarchy->Update_Boxes(collisions.saved_states(state1).x->X,collisions.saved_states(state2).x->X);
    else collisions.object.hierarchy->Update_Boxes();
    collisions.object.Update_Bounding_Box();
}
template<class TV> void DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
Update_Intersection_Acceleration_Structures(const bool use_swept_hierarchy,const int state1,const int state2)
{
    Update_Intersection_Acceleration_Structures_Helper(*this,use_swept_hierarchy,state1,state2);
}
template<class TV> void DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
//#####################################################################
// Function Update_Bounding_Box
//#####################################################################
Update_Bounding_Box()
{
    if(volume_object) volume_object->Update_Bounding_Box();
    else object.Update_Bounding_Box();
}
//#####################################################################
// Function Axis_Aligned_Bounding_Box
//#####################################################################
template<class TV> const RANGE<TV>& DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
Axis_Aligned_Bounding_Box() const
{ 
    if(volume_object) return *volume_object->bounding_box;
    else return *object.bounding_box;
}
//##################################################################### 
// Function Pointwise_Object_Velocity
//##################################################################### 
template<class TV> TV DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
Pointwise_Object_Velocity(const int simplex_id,const TV& location) const
{
    if(volume_object){
        typedef VECTOR<T,TV::m+1> T_VOLUME_WEIGHTS;
        VECTOR<int,TV::dimension+1> nodes=(*volume_object).mesh.elements(simplex_id);
        T_VOLUME_WEIGHTS weights=T_VOLUME_SIMPLEX::Barycentric_Coordinates(location,particles.X.Subset(nodes));
        return T_VOLUME_SIMPLEX::Point_From_Barycentric_Coordinates(weights,particles.V.Subset(nodes));}
    else{
        VECTOR<int,TV::dimension> nodes=object.mesh.elements(simplex_id);
        typename IF<TV::m==1,ONE,TV>::TYPE weights=T_SIMPLEX::Clamped_Barycentric_Coordinates(location,particles.X.Subset(nodes));
        return T_SIMPLEX::Point_From_Barycentric_Coordinates(weights,particles.V.Subset(nodes));}
}
//##################################################################### 
// Function Simplex_Intersection
//##################################################################### 
template<class TV> bool DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
Simplex_Intersection(RAY<TV>& ray) const
{
    return INTERSECTION::Intersects(ray,object,collision_thickness);
}
//##################################################################### 
// Function Simplex_Closest_Non_Intersecting_Point
//##################################################################### 
template<class TV> bool DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
Simplex_Closest_Non_Intersecting_Point(RAY<TV>& ray) const
{
    return INTERSECTION::Closest_Non_Intersecting_Point(ray,object,collision_thickness);
}
//##################################################################### 
// Function Inside_Any_Simplex
//##################################################################### 
template<class TV> bool DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
Inside_Any_Simplex(const TV& location,int& simplex_id) const
{
    if(volume_object){
        return TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Inside_Any_Simplex(*volume_object,location,simplex_id,(T).5*collision_thickness);}
    else{
        return object.Inside_Any_Simplex(location,simplex_id,(T).5*collision_thickness);}
}
//##################################################################### 
// Function Inside
//##################################################################### 
template<class TV> bool DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
Inside(const TV& location,const T thickness_over_two) const
{
    return object.Inside(location,thickness_over_two);
}
//##################################################################### 
// Function Implicit_Geometry_Lazy_Inside
//##################################################################### 
template<class TV> bool DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
Implicit_Geometry_Lazy_Inside(const TV& location,T contour_value) const
{
    assert(contour_value==0); // TODO(jontg): As implemented, this function does not with non-default arguments
    return object.Inside(location,0);
}
//##################################################################### 
// Function Implicit_Geometry_Extended_Value
//##################################################################### 
template<class TV> typename TV::SCALAR DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
Implicit_Geometry_Extended_Value(const TV& location) const
{
    return object.Calculate_Signed_Distance(location);
}
//##################################################################### 
// Function Simplex_Closest_Point_On_Boundary
//##################################################################### 
template<class TV> TV DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
Simplex_Closest_Point_On_Boundary(const TV& location,const T max_distance,const T thickness_over_2,int* simplex_id,T* distance) const
{
    // TODO Is there a better way to force the use of collision_thickness?
    if(thickness_over_2>=(T)0)
        return object.Closest_Point_On_Boundary(location,max_distance,thickness_over_2,simplex_id,distance);
    else
        return object.Closest_Point_On_Boundary(location,max_distance,collision_thickness,simplex_id,distance);
}
//##################################################################### 
// Function Simplex_World_Space_Point_From_Barycentric_Coordinates
//##################################################################### 
template<class TV> TV DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
Simplex_World_Space_Point_From_Barycentric_Coordinates(const int simplex_id,const T_WEIGHTS& weights) const
{
    return T_SIMPLEX::Point_From_Barycentric_Coordinates(weights,particles.X.Subset(object.mesh.elements(simplex_id)));
}
//#####################################################################
// Function Number_Of_Simplices
//#####################################################################
template<class TV> int DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>::
Number_Of_Simplices() const
{
    return object.mesh.elements.m;
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template class DEFORMABLE_OBJECT_FLUID_COLLISIONS<VECTOR<T,d> >;
INSTANTIATION_HELPER(float,1)
INSTANTIATION_HELPER(float,2)
INSTANTIATION_HELPER(float,3)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1)
INSTANTIATION_HELPER(double,2)
INSTANTIATION_HELPER(double,3)
#endif
}

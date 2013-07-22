//#####################################################################
// Copyright 2003-2008, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Data_Structures/STACK.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOUNDED_HORIZONTAL_PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Collision_Detection/COLLISION_GEOMETRY_SPATIAL_PARTITION.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/MULTIBODY_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_CONTACT_GRAPH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_SKIP_COLLISION_CHECK.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_STRUCTURE_INTERACTION_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGIDS_COLLISION_CALLBACKS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/SOLVE_CONTACT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Fracture/FRACTURE_PATTERN.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/NORMAL_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_CLUSTER_CONSTITUENT_ID.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <cfloat>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_BODY_COLLISIONS<TV>::
RIGID_BODY_COLLISIONS(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,RIGID_BODY_COLLISION_PARAMETERS<TV>& parameters_input,
    RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks_input,RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>& rigids_example_forces_and_velocities_input)
    :parameters(parameters_input),collision_callbacks(collision_callbacks_input),verbose(false),prune_stacks_from_contact(false),prune_contact_using_velocity(false),collision_manager(0),
    skip_collision_check(*new RIGID_BODY_SKIP_COLLISION_CHECK),rigid_body_collection(rigid_body_collection_input),rigids_example_forces_and_velocities(rigids_example_forces_and_velocities_input),
    spatial_partition(new COLLISION_GEOMETRY_SPATIAL_PARTITION<COLLISION_GEOMETRY<TV>,ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>(rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies)),
    intersections(*new RIGID_BODY_INTERSECTIONS<TV>(rigid_body_collection)),triangle_collisions(0),triangle_collision_geometry(0),contact_graph(*new RIGID_BODY_CONTACT_GRAPH<TV>(rigid_body_collection.rigid_body_particle)),
    store_collision_intersections_for_projection(false),use_static_body_masses(false),use_parent_normal(false),
    rigid_body_cluster_bindings(rigid_body_collection.rigid_body_cluster_bindings),fracture_pattern(0),mpi_rigids(0)
{
    Set_Collision_Pair_Iterations();
    Set_Contact_Level_Iterations();Set_Contact_Pair_Iterations();
    Set_Shock_Propagation_Iterations();Set_Shock_Propagation_Level_Iterations();Set_Shock_Propagation_Pair_Iterations();
    Set_Push_Out_Iterations();Set_Push_Out_Level_Iterations();Set_Push_Out_Pair_Iterations();
    Use_Freezing_With_Push_Out();Use_Gradual_Push_Out();Set_Desired_Separation_Distance();
    Use_Rolling_Friction(false);object_indices.Preallocate(50);

    if(parameters.use_analytic_collisions) Register_Analytic_Collisions();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_BODY_COLLISIONS<TV>::
~RIGID_BODY_COLLISIONS()
{
    delete &contact_graph;
    delete &intersections;
    delete spatial_partition;
    delete &skip_collision_check;
    delete fracture_pattern;
    delete triangle_collisions;
    delete triangle_collision_geometry;
}
//#####################################################################
// Function Initialize_Data_Structures
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Initialize_Data_Structures(const bool reset)
{
    contact_graph.Initialize();skip_collision_check.Initialize(rigid_body_collection.rigid_body_particle.array_collection->Size(),reset);
    if(parameters.use_ccd){delete triangle_collisions;delete triangle_collision_geometry;
        triangle_collision_geometry=new RIGID_TRIANGLE_COLLISIONS_GEOMETRY<TV>();
        for(int n=1;n<=rigid_body_collection.rigid_body_particle.array_collection->Size();n++)
            triangle_collision_geometry->structures.Append(rigid_body_collection.Rigid_Body(n).simplicial_object);
        triangle_collision_geometry->Build_Collision_Geometry();
        triangle_collisions=new RIGID_TRIANGLE_COLLISIONS<TV>(*triangle_collision_geometry);
        triangle_collisions->Update_Local_Hierachies();}
}
//#####################################################################
// Function Adjust_Bounding_Boxes
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Adjust_Bounding_Boxes(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const T false_thickness,const T extra_padding_distance)
{
    for(COLLISION_GEOMETRY_ID i(1);i<=rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies.Size();i++){
        RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_geometry=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies(i));
        if(rigid_collision_geometry){RIGID_BODY<TV>* rigid_body=dynamic_cast<RIGID_BODY<TV>*>(&rigid_collision_geometry->rigid_geometry);
            if(rigid_body) if(rigid_body_collection.Is_Active(rigid_body->particle_index)){
                if(rigid_body->simplicial_object){ // only modify triangulated surface bounding box
                    assert(rigid_body->simplicial_object->bounding_box);
                    RANGE<TV>& box=*rigid_body->simplicial_object->bounding_box;
                    bool already_adjusted=false; // multiple rigid bodies can share a triangulated surface - don't adjust bounding box multiple times
                    for(COLLISION_GEOMETRY_ID j(1);j<=i-1;j++) if(rigid_body_collection.rigid_geometry_collection.collision_body_list->Is_Active(j)){
                        int rigid_body_id_j=rigid_body_collection.rigid_geometry_collection.collision_body_list->collision_geometry_id_to_geometry_id.Get(j);
                        if(rigid_body_id_j && rigid_body_collection.Is_Active(rigid_body_id_j) && rigid_body_collection.Rigid_Body(rigid_body_id_j).simplicial_object &&
                            rigid_body->simplicial_object==rigid_body_collection.Rigid_Body(rigid_body_id_j).simplicial_object){
                            already_adjusted=true;break;}}
                    if(!already_adjusted){
                        TV delta,edge_lengths=box.Edge_Lengths();for(int d=1;d<=TV::dimension;d++) if(edge_lengths[d]==0) delta[d]=false_thickness;
                        box.Change_Size(delta+extra_padding_distance);}
                    rigid_body->Update_Bounding_Box();}}}}
}
//#####################################################################
// Function Get_Deepest_Intersection_Point
//#####################################################################
template<class TV> bool RIGID_BODY_COLLISIONS<TV>::
Get_Deepest_Intersection_Point(const int id_1,const int id_2,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,
    T& smallest_value,int& smallest_index,TV& collision_location,TV& collision_normal,TV& collision_relative_velocity,const bool ignore_separating,const T desired_separation_distance,
    bool& ignored_separating) const
{
    ignored_separating=false;
    particle_intersections.Remove_All();
    intersections.Append_All_Intersections(id_1,id_2,particle_intersections,desired_separation_distance);
    if(!particle_intersections.m){skip_collision_check.Set_Last_Checked(id_1,id_2);return false;}
    smallest_value=FLT_MAX;smallest_index=0;
    for(int i=1;i<=particle_intersections.m;i++){
        const RIGID_BODY_PARTICLE_INTERSECTION<TV>& intersection=particle_intersections(i);
        T phi=(*rigid_body_collection.Rigid_Body(intersection.levelset_body).implicit_object)(rigid_body_collection.Rigid_Body(intersection.particle_body).World_Space_Point(intersection.particle_location));
        if(phi<smallest_value){
            RIGID_BODY<TV> &body1=rigid_body_collection.Rigid_Body(intersection.particle_body),&body2=rigid_body_collection.Rigid_Body(intersection.levelset_body);
            RIGID_BODY<TV> &body_parent1=rigid_body_cluster_bindings.Get_Parent(body1),&body_parent2=rigid_body_cluster_bindings.Get_Parent(body2);
            if(ignore_separating){
                TV location=body1.World_Space_Point(intersection.particle_location);
                TV normal=use_parent_normal?body_parent2.Implicit_Geometry_Normal(location):body2.Implicit_Geometry_Normal(location);
                TV relative_velocity=RIGID_GEOMETRY<TV>::Relative_Velocity_At_Geometry1_Particle(body_parent1,body_parent2,location,intersection.particle_index);
                if(TV::Dot_Product(relative_velocity,normal)>=0){ignored_separating=true;continue;}
                collision_location=location;collision_normal=normal;collision_relative_velocity=relative_velocity;}
            smallest_value=phi;smallest_index=i;}}

    if(smallest_index && !ignore_separating){ // these quantities are already computed if we ignore_separating
        const RIGID_BODY_PARTICLE_INTERSECTION<TV>& intersection=particle_intersections(smallest_index);
        RIGID_BODY<TV> &body1=rigid_body_collection.Rigid_Body(intersection.particle_body),&body2=rigid_body_collection.Rigid_Body(intersection.levelset_body);
        collision_location=body1.World_Space_Point(intersection.particle_location);collision_normal=body2.Implicit_Geometry_Normal(collision_location);collision_relative_velocity=TV();}

    return smallest_index!=0;
}
//#####################################################################
// Function Get_First_Intersection_Point
//#####################################################################
template<class T,class T_ARRAY> void Create_Edges(const T_ARRAY& X1,const T_ARRAY& X2,const VECTOR<int,2>& nodes,POINT_2D<T>& point1,POINT_2D<T>& point2)
{
    point1=POINT_2D<T>(X1(nodes[1]));
    point2=POINT_2D<T>(X2(nodes[2]));
}
template<class T,class T_ARRAY> void Create_Edges(const T_ARRAY& X1,const T_ARRAY& X2,const VECTOR<int,4>& nodes,SEGMENT_3D<T>& segment1,SEGMENT_3D<T>& segment2)
{
    segment1=SEGMENT_3D<T>(X1(nodes[1]),X1(nodes[2]));
    segment2=SEGMENT_3D<T>(X2(nodes[3]),X2(nodes[4]));
}
template<> bool RIGID_BODY_COLLISIONS<VECTOR<float,1> >::
Get_First_Intersection_Point(const int id_1,const int id_2,T& smallest_value,int& smallest_index,VECTOR<float,1>& collision_location,VECTOR<float,1>& collision_normal,VECTOR<float,1>& collision_relative_velocity,
    const bool ignore_separating,const T desired_separation_distance,bool& ignored_separating,const T dt,const bool swap)
{PHYSBAM_FATAL_ERROR();}
template<> bool RIGID_BODY_COLLISIONS<VECTOR<double,1> >::
Get_First_Intersection_Point(const int id_1,const int id_2,T& smallest_value,int& smallest_index,VECTOR<double,1>& collision_location,VECTOR<double,1>& collision_normal,VECTOR<double,1>& collision_relative_velocity,
    const bool ignore_separating,const T desired_separation_distance,bool& ignored_separating,const T dt,const bool swap)
{PHYSBAM_FATAL_ERROR();}
template<class TV> bool RIGID_BODY_COLLISIONS<TV>::
Get_First_Intersection_Point(const int id_1,const int id_2,T& smallest_value,int& smallest_index,TV& collision_location,TV& collision_normal,TV& collision_relative_velocity,
    const bool ignore_separating,const T desired_separation_distance,bool& ignored_separating,const T dt,const bool swap)
{
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX T_FACE;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX_FACE T_EDGE;
    
    PHYSBAM_ASSERT(triangle_collisions);
    RIGID_BODY_PARTICLES<TV>& particle=rigid_body_collection.rigid_body_particle;
    triangle_collision_geometry->interacting_structure_pairs.Remove_All();triangle_collision_geometry->interacting_structure_pairs.Append(VECTOR<int,2>(id_1,id_2));
    VECTOR<FRAME<TV>,2> frame1,frame2;frame1(2)=rigid_body_collection.Rigid_Body(id_1).Frame();frame2(2)=rigid_body_collection.Rigid_Body(id_2).Frame();
    triangle_collision_geometry->structure_geometries(id_1)->Save_Current_State(particle.X(id_1),particle.rotation(id_1),particle.V(id_1));
    triangle_collision_geometry->structure_geometries(id_2)->Save_Current_State(particle.X(id_2),particle.rotation(id_2),particle.V(id_2));
    if(swap) collision_callbacks.Swap_States(id_1,id_2);
    else{collision_callbacks.Restore_Position(id_1);collision_callbacks.Restore_Position(id_2);}
    frame1(1)=rigid_body_collection.Rigid_Body(id_1).Frame();frame2(1)=rigid_body_collection.Rigid_Body(id_2).Frame();
    triangle_collision_geometry->structure_geometries(id_1)->Save_Self_Collision_Free_State(particle.X(id_1),particle.rotation(id_1),particle.V(id_1));
    triangle_collision_geometry->structure_geometries(id_2)->Save_Self_Collision_Free_State(particle.X(id_2),particle.rotation(id_2),particle.V(id_2));
    triangle_collision_geometry->structure_geometries(id_1)->Save_Pseudo_State(dt);
    triangle_collision_geometry->structure_geometries(id_2)->Save_Pseudo_State(dt);
    triangle_collisions->Compute_Pairs(parameters.collision_thickness,&frame1,&frame2);
    smallest_value=2*dt;
    for(int pair_index=1;pair_index<=triangle_collisions->point_face_pairs_internal.m;pair_index++){int structure_1=id_1,structure_2=id_2;
        if(pair_index>triangle_collisions->swap_index){structure_1=id_2;structure_2=id_1;}
        RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_geom1=*triangle_collision_geometry->structure_geometries(structure_1);
        RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_geom2=*triangle_collision_geometry->structure_geometries(structure_2);
        const VECTOR<int,TV::dimension+1>& nodes=triangle_collisions->point_face_pairs_internal(pair_index);
        T relative_speed;T collision_time;TV normal,weights;T_FACE face(structure_geom2.X_self_collision_free.Subset(nodes.Remove_Index(1)));
        if(face.Point_Face_Collision(structure_geom1.X_self_collision_free(nodes[1]),structure_geom1.V(nodes[1]),structure_geom2.V.Subset(nodes.Remove_Index(1)),
            dt,parameters.collision_thickness,collision_time,normal,weights,relative_speed,false)){
            if(collision_time<smallest_value){
                smallest_value=collision_time;
                collision_relative_velocity=structure_geom1.V(nodes[1])-(weights(1)*structure_geom2.V(nodes[2])+weights(2)*structure_geom2.V(nodes[3])+weights(3)*structure_geom2.V(nodes[4]));
                if(structure_1!=id_1) collision_relative_velocity*=-1;
                if(TV::Dot_Product(collision_relative_velocity,normal)>0) collision_normal=-normal; else collision_normal=normal;
                collision_location=structure_geom1.X_self_collision_free(nodes[1])+smallest_value*structure_geom1.V(nodes[1]);}}}
    for(int pair_index=1;pair_index<=triangle_collisions->edge_edge_pairs_internal.m;pair_index++){
        RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_geom1=*triangle_collision_geometry->structure_geometries(id_1);
        RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_geom2=*triangle_collision_geometry->structure_geometries(id_2);
        const VECTOR<int,2*TV::dimension-2>& nodes=triangle_collisions->edge_edge_pairs_internal(pair_index);
        T relative_speed;T collision_time;TV normal;VECTOR<T,2> weights;T_EDGE edge1,edge2;Create_Edges(structure_geom1.X_self_collision_free,structure_geom2.X_self_collision_free,nodes,edge1,edge2);
        ARRAY<TV> V_edges(2*TV::dimension-2);V_edges(1)=structure_geom1.V(nodes[1]);V_edges(2)=structure_geom1.V(nodes[2]);V_edges(3)=structure_geom2.V(nodes[3]);V_edges(4)=structure_geom2.V(nodes[4]);
        if(edge1.Edge_Edge_Collision(edge2,V_edges,dt,parameters.collision_thickness,collision_time,normal,weights,relative_speed,false,triangle_collision_geometry->small_number,false)){
            if(collision_time<smallest_value){
                smallest_value=collision_time;
                collision_relative_velocity=(((T)1.0-weights(1))*structure_geom1.V(nodes[1])+weights(1)*structure_geom1.V(nodes[2]))-(((T)1.0-weights(2))*structure_geom2.V(nodes[3])+weights(2)*structure_geom2.V(nodes[4]));
                if(TV::Dot_Product(collision_relative_velocity,normal)>0) collision_normal=-normal; else collision_normal=normal;
                collision_location=(1-weights(1))*structure_geom1.X_self_collision_free(nodes[1])+weights(1)*structure_geom1.X_self_collision_free(nodes[2])+
                    smallest_value*((1-weights(1))*structure_geom1.V(nodes[1])+weights(1)*structure_geom1.V(nodes[2]));}}}
    assert(smallest_value>=0);
    return smallest_value<dt;
}
//#####################################################################
// Function Update_Collision_Pair
//#####################################################################
template<class TV> bool RIGID_BODY_COLLISIONS<TV>::
Update_Collision_Pair(const int id_1,const int id_2,const T dt,const T time,const bool mpi_one_ghost)
{
    RIGID_BODY<TV>& body1=rigid_body_collection.Rigid_Body(id_1),&body2=rigid_body_collection.Rigid_Body(id_2);
    VECTOR<std::string,2> key(typeid(*body1.implicit_object->object_space_implicit_object).name(),typeid(*body2.implicit_object->object_space_implicit_object).name());
    if(key.x==typeid(MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>).name()||key.y==typeid(MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>).name()){
        if(!body1.simplicial_object && !body2.simplicial_object) return Update_Analytic_Multibody_Collision(body1,body2,dt,time,mpi_one_ghost);}
    if(bool (**collision_function)(RIGID_BODY_COLLISIONS&,const int,const int,IMPLICIT_OBJECT<TV>*,IMPLICIT_OBJECT<TV>*,const T,const T,const bool)=analytic_collision_registry.Get_Pointer(key.Sorted()))
        return (*collision_function)(*this,id_1,id_2,body1.implicit_object->object_space_implicit_object,body2.implicit_object->object_space_implicit_object,dt,time,mpi_one_ghost);
    if(parameters.use_ccd) return Update_Surface_Collision_Pair(id_1,id_2,dt,time,mpi_one_ghost);
    return Update_Levelset_Collision_Pair(id_1,id_2,dt,time,mpi_one_ghost);
}
//#####################################################################
// Function Update_Collision_Pair_Helper
//#####################################################################
template<class TV> bool RIGID_BODY_COLLISIONS<TV>::
Update_Collision_Pair_Helper(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const T dt,const T time,const TV& collision_location,const TV& collision_normal,const TV& collision_relative_velocity,
    const bool mpi_one_ghost)
{
    fractured_bodies=VECTOR<int,2>();
    RIGID_BODY<TV> &body_parent1=rigid_body_cluster_bindings.Get_Parent(body1),&body_parent2=rigid_body_cluster_bindings.Get_Parent(body2);
    // TODO: ok to store object space location of intersection point?
    if(store_collision_intersections_for_projection)
        rigid_body_particle_intersections.Set(Tuple(body1.particle_index,body2.particle_index,collision_location));
    assert(Either_Body_Collides_With_The_Other(body1.particle_index,body2.particle_index));
    // TODO: put an assert where we compute joints that collisions are two-way (later make joint respect one-way collisions)
    if(!Body_Collides_With_The_Other(body1.particle_index,body2.particle_index)) body_parent1.is_temporarily_static=true;
    if(!Body_Collides_With_The_Other(body2.particle_index,body1.particle_index)) body_parent2.is_temporarily_static=true;
    TWIST<TV> impulse=collision_callbacks.Compute_Collision_Impulse(body_parent1,body_parent2,collision_location,collision_normal,collision_relative_velocity,
        RIGID_BODY<TV>::Coefficient_Of_Restitution(body_parent1,body_parent2),RIGID_BODY<TV>::Coefficient_Of_Friction(body_parent1,body_parent2),true,rolling_friction,true);

    if(parameters.use_fracture_pattern){
        if(!fracture_pattern) fracture_pattern=new FRACTURE_PATTERN<T>;
        if(!body1.Has_Infinite_Inertia() && impulse.linear.Magnitude()>=body1.fracture_threshold){
            TV object_space_collision_location=body1.Object_Space_Point(collision_location);
            collision_callbacks.Begin_Fracture(body1.particle_index);
            fracture_pattern->Intersect_With_Rigid_Body(body1,body1.World_Space_Point(object_space_collision_location),added_bodies(1),parameters.allow_refracturing,
                parameters.use_fracture_particle_optimization);
            fractured_bodies(1)=body1.particle_index;
            rigid_body_collection.rigid_body_particle.Remove_Body(body1.particle_index);
            rigid_body_collection.rigid_geometry_collection.Destroy_Unreferenced_Geometry();
            collision_callbacks.End_Fracture(body1.particle_index,added_bodies(1));}
        else body1.Apply_Impulse_To_Body(collision_location,impulse.linear,impulse.angular,mpi_one_ghost);
        if(!body2.Has_Infinite_Inertia() && impulse.linear.Magnitude()>=body2.fracture_threshold){
            TV object_space_collision_location=body2.Object_Space_Point(collision_location);
            collision_callbacks.Begin_Fracture(body2.particle_index);
            fracture_pattern->Intersect_With_Rigid_Body(body2,body2.World_Space_Point(object_space_collision_location),added_bodies(2),parameters.allow_refracturing,
                parameters.use_fracture_particle_optimization);
            fractured_bodies(2)=body2.particle_index;
            rigid_body_collection.rigid_body_particle.Remove_Body(body2.particle_index);
            rigid_body_collection.rigid_geometry_collection.Destroy_Unreferenced_Geometry();
            collision_callbacks.End_Fracture(body2.particle_index,added_bodies(2));}
        else body2.Apply_Impulse_To_Body(collision_location,-impulse.linear,-impulse.angular,mpi_one_ghost);
        if(fractured_bodies(1) || fractured_bodies(2)){
            rigid_body_collection.Update_Simulated_Particles();
            Initialize_Data_Structures(false);
            for(int i=1;i<=added_bodies(1).m;i++){
                collision_callbacks.Euler_Step_Position_With_New_Velocity(added_bodies(1)(i),dt,time);
                skip_collision_check.Set_Last_Moved(added_bodies(1)(i));}
            for(int i=1;i<=added_bodies(2).m;i++){
                collision_callbacks.Euler_Step_Position_With_New_Velocity(added_bodies(2)(i),dt,time);
                skip_collision_check.Set_Last_Moved(added_bodies(2)(i));}}}
    else RIGID_BODY<TV>::Apply_Impulse(body_parent1,body_parent2,collision_location,impulse.linear,impulse.angular,mpi_one_ghost);
    if(!fractured_bodies(1)){
        body_parent1.is_temporarily_static=false;
        int parent_id_1=body_parent1.particle_index;
        if(!body_parent1.Has_Infinite_Inertia()){
            collision_callbacks.Euler_Step_Position_With_New_Velocity(parent_id_1,dt,time);
            skip_collision_check.Set_Last_Moved(parent_id_1);}}
    if(!fractured_bodies(2)){
        body_parent2.is_temporarily_static=false;
        int parent_id_2=body_parent2.particle_index;
        if(!body_parent2.Has_Infinite_Inertia()){
            collision_callbacks.Euler_Step_Position_With_New_Velocity(parent_id_2,dt,time);
            skip_collision_check.Set_Last_Moved(parent_id_2);}}
    return fractured_bodies(1)||fractured_bodies(2);
}
//#####################################################################
// Function Update_Analytic_Multibody_Collision
//#####################################################################
template<class TV> bool RIGID_BODY_COLLISIONS<TV>::
Update_Analytic_Multibody_Collision(const int id_1,const int id_2,MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>& multibody,IMPLICIT_OBJECT<TV>& levelset,const typename TV::SCALAR dt,const typename TV::SCALAR time,const bool mpi_one_ghost)
{
    bool return_val=false;
    IMPLICIT_OBJECT<TV>* levelset_base=&levelset;
    if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* levelset_transformed=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(levelset_base)) 
        levelset_base=levelset_transformed->object_space_implicit_object;
    for(int i=1;i<=multibody.levelsets->m;i++){
        VECTOR<std::string,2> key(typeid(*dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >&>(*(*multibody.levelsets)(i)).object_space_implicit_object).name(),typeid(*levelset_base).name());
        if(bool (**collision_function)(RIGID_BODY_COLLISIONS&,const int,const int,IMPLICIT_OBJECT<TV>*,IMPLICIT_OBJECT<TV>*,const T,const T,const bool)=analytic_collision_registry.Get_Pointer(key.Sorted()))
            return_val|=(*collision_function)(*this,id_1,id_2,(*multibody.levelsets)(i),&levelset,dt,time,mpi_one_ghost);
        else PHYSBAM_FATAL_ERROR();}
    return return_val;
}
//#####################################################################
// Function Update_Analytic_Multibody_Collision
//#####################################################################
template<class TV> bool RIGID_BODY_COLLISIONS<TV>::
Update_Analytic_Multibody_Collision(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const typename TV::SCALAR dt,const typename TV::SCALAR time,const bool mpi_one_ghost)
{
    bool return_val=false;
    MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>* multibody1=dynamic_cast<MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>*>(body1.implicit_object->object_space_implicit_object);
    MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>* multibody2=dynamic_cast<MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>*>(body2.implicit_object->object_space_implicit_object);
    if(multibody1 && multibody2) for(int i=1;i<=multibody2->levelsets->m;i++) return_val|=Update_Analytic_Multibody_Collision(body1.particle_index,body2.particle_index,*multibody1,*(*multibody2->levelsets)(i),dt,time,mpi_one_ghost);
    else if(multibody1) return_val=Update_Analytic_Multibody_Collision(body1.particle_index,body2.particle_index,*multibody1,*body2.implicit_object->object_space_implicit_object,dt,time,mpi_one_ghost);
    else return_val=Update_Analytic_Multibody_Collision(body1.particle_index,body2.particle_index,*multibody2,*body1.implicit_object->object_space_implicit_object,dt,time,mpi_one_ghost);
    return return_val;
}
//#####################################################################
// Function Update_Box_Box_Collision
//#####################################################################
template<class TV> bool
Update_Box_Box_Collision(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,const int i1,const int i2,IMPLICIT_OBJECT<TV>* object1,
    IMPLICIT_OBJECT<TV>* object2,const typename TV::SCALAR dt,const typename TV::SCALAR time,const bool mpi_one_ghost)
{
    typedef typename TV::SCALAR T;
    RIGID_BODY<TV>& body1=rigid_body_collisions.rigid_body_collection.Rigid_Body(i1);
    RIGID_BODY<TV>& body2=rigid_body_collisions.rigid_body_collection.Rigid_Body(i2);
    FRAME<TV> transform1,transform2;
    if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* object_transformed=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(object1)){
        transform1=*object_transformed->transform;object1=object_transformed->object_space_implicit_object;}
    if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* object_transformed=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(object2)){
        transform2=*object_transformed->transform;object2=object_transformed->object_space_implicit_object;}
    BOX<TV>& box1=dynamic_cast<ANALYTIC_IMPLICIT_OBJECT<BOX<TV> >&>(*object1).analytic;
    BOX<TV>& box2=dynamic_cast<ANALYTIC_IMPLICIT_OBJECT<BOX<TV> >&>(*object2).analytic;
    ORIENTED_BOX<TV> box1_transformed(box1,body1.Frame());
    ORIENTED_BOX<TV> box2_transformed(box2,body2.Frame());

    TV collision_normal=body1.X()-body2.X();collision_normal.Normalize();
    if(!box1_transformed.Intersection(box2_transformed)){rigid_body_collisions.skip_collision_check.Set_Last_Checked(i1,i2);return false;}
    rigid_body_collisions.pairs_processed_by_collisions.Set(VECTOR<int,2>(i1,i2).Sorted());
    if(TV::Dot_Product(collision_normal,body1.Twist().linear-body2.Twist().linear)>=0) return false;
    TV collision_location=(body1.X()-body2.X())*.5+body1.X();
    rigid_body_collisions.Update_Collision_Pair_Helper(body1,body2,dt,time,collision_location,collision_normal,
        body1.Pointwise_Object_Velocity(collision_location)-body2.Pointwise_Object_Velocity(collision_location),mpi_one_ghost);
    return true;
}
//#####################################################################
// Function Update_Sphere_Sphere_Collision
//#####################################################################
template<class TV> bool
Update_Sphere_Sphere_Collision(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,const int i1,const int i2,IMPLICIT_OBJECT<TV>* object1,
    IMPLICIT_OBJECT<TV>* object2,const typename TV::SCALAR dt,const typename TV::SCALAR time,const bool mpi_one_ghost)
{
    typedef typename TV::SCALAR T;
    RIGID_BODY<TV>& body1=rigid_body_collisions.rigid_body_collection.Rigid_Body(i1);
    RIGID_BODY<TV>& body2=rigid_body_collisions.rigid_body_collection.Rigid_Body(i2);
    FRAME<TV> transform1,transform2;
    if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* object_transformed=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(object1)){
        transform1=*object_transformed->transform;object1=object_transformed->object_space_implicit_object;}
    if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* object_transformed=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(object2)){
        transform2=*object_transformed->transform;object2=object_transformed->object_space_implicit_object;}
    SPHERE<TV>& sphere1=dynamic_cast<ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >&>(*object1).analytic;
    SPHERE<TV>& sphere2=dynamic_cast<ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >&>(*object2).analytic;

    TV sphere1_center=(body1.Frame()*transform1).t,sphere2_center=(body2.Frame()*transform2).t;
    TV collision_normal=body1.X()-body2.X();collision_normal.Normalize();
    T d=(sphere1_center-sphere2_center).Magnitude(),r1=sphere1.radius,r2=sphere2.radius;
    if(d>r1+r2){rigid_body_collisions.skip_collision_check.Set_Last_Checked(i1,i2);return false;}
    rigid_body_collisions.pairs_processed_by_collisions.Set(VECTOR<int,2>(i1,i2).Sorted());
    if(TV::Dot_Product(collision_normal,body1.Twist().linear-body2.Twist().linear)>=0) return false;
    TV collision_location=sphere1_center+(T).5*(d+min(d,r2)-min(d,r1))*collision_normal;
    rigid_body_collisions.Update_Collision_Pair_Helper(body1,body2,dt,time,collision_location,collision_normal,
        body1.Pointwise_Object_Velocity(collision_location)-body2.Pointwise_Object_Velocity(collision_location),mpi_one_ghost);
    return true;
}
//#####################################################################
// Function Update_Box_Plane_Collision
//#####################################################################
template<class TV> bool
Update_Box_Plane_Collision(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,const int i1,const int i2,IMPLICIT_OBJECT<TV>* object1,
    IMPLICIT_OBJECT<TV>* object2,const typename TV::SCALAR dt,const typename TV::SCALAR time,const bool mpi_one_ghost)
{
    typedef typename TV::SCALAR T;
    FRAME<TV> transform;
    RIGID_BODY<TV>* body1=&rigid_body_collisions.rigid_body_collection.Rigid_Body(i1);
    RIGID_BODY<TV>* body2=&rigid_body_collisions.rigid_body_collection.Rigid_Body(i2);
    if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* object_transformed=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(object2)){
        transform=*object_transformed->transform;object2=object_transformed->object_space_implicit_object;}
    ANALYTIC_IMPLICIT_OBJECT<BOX<TV> >* implicit_box=dynamic_cast<ANALYTIC_IMPLICIT_OBJECT<BOX<TV> >*>(object2);
    if(!implicit_box){
        if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* object_transformed=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(object1)){
            transform=*object_transformed->transform;object1=object_transformed->object_space_implicit_object;}
        exchange(object1,object2);
        exchange(body1,body2);
        implicit_box=dynamic_cast<ANALYTIC_IMPLICIT_OBJECT<BOX<TV> >*>(object2);}
    BOX<TV>& box=implicit_box->analytic;

    ARRAY<TV> points;bool intersect=false;
    for(int i=1;i<=8;i++){
        TV point=box.min_corner;
        if(i>4) point(1)=box.max_corner(1);
        if(i%4==0||i%4==3) point(2)=box.max_corner(2);
        if(i%2==0) point(3)=box.max_corner(3);
        TV transformed_point=body1->Frame().Inverse()*body2->Frame()*point;
        if(transformed_point(2)<0){
            intersect=true;
            transformed_point(2)*=.5;
            points.Append(transformed_point);}}

    TV collision_normal=-body1->Rotation().Rotated_Axis(2);
    if(!intersect){rigid_body_collisions.skip_collision_check.Set_Last_Checked(i1,i2);return false;}
    rigid_body_collisions.pairs_processed_by_collisions.Set(VECTOR<int,2>(i1,i2).Sorted());
    if(TV::Dot_Product(body1->Twist().linear-body2->Twist().linear,collision_normal)>=0) return false;

    TV collision_location;for(int i=1;i<=points.m;i++) collision_location+=points(i);collision_location/=(T)points.m;collision_location=body1->Frame()*collision_location;
    rigid_body_collisions.Update_Collision_Pair_Helper(*body1,*body2,dt,time,collision_location,collision_normal,
        body1->Pointwise_Object_Velocity(collision_location)-body2->Pointwise_Object_Velocity(collision_location),mpi_one_ghost);
    return true;
}
//#####################################################################
// Function Update_Sphere_Plane_Collision
//#####################################################################
template<class TV> bool
Update_Sphere_Plane_Collision(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,const int i1,const int i2,IMPLICIT_OBJECT<TV>* object1,
    IMPLICIT_OBJECT<TV>* object2,const typename TV::SCALAR dt,const typename TV::SCALAR time,const bool mpi_one_ghost)
{
    typedef typename TV::SCALAR T;
    FRAME<TV> transform;
    RIGID_BODY<TV>* body1=&rigid_body_collisions.rigid_body_collection.Rigid_Body(i1);
    RIGID_BODY<TV>* body2=&rigid_body_collisions.rigid_body_collection.Rigid_Body(i2);
    if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* object_transformed=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(object2)){
        transform=*object_transformed->transform;object2=object_transformed->object_space_implicit_object;}
    ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >* implicit_sphere=dynamic_cast<ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >*>(object2);
    if(!implicit_sphere){
        if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* object_transformed=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(object1)){
            transform=*object_transformed->transform;object1=object_transformed->object_space_implicit_object;}
        exchange(object1,object2);
        exchange(body1,body2);
        implicit_sphere=dynamic_cast<ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >*>(object2);}
    SPHERE<TV>& sphere=implicit_sphere->analytic;

    TV sphere_center=(body2->Frame()*transform).t;
    TV collision_normal=-body1->Rotation().Rotated_Axis(2);
    T separation=TV::Dot_Product(body1->X()-sphere_center,collision_normal);
    if(separation>=sphere.radius){rigid_body_collisions.skip_collision_check.Set_Last_Checked(i1,i2);return false;}
    rigid_body_collisions.pairs_processed_by_collisions.Set(VECTOR<int,2>(i1,i2).Sorted());
    if(TV::Dot_Product(body1->Twist().linear-body2->Twist().linear,collision_normal)>=0) return false;

    T collision_depth=(T).5*(-sphere.radius+min(-separation,sphere.radius));
    TV collision_location=sphere_center-collision_depth*collision_normal;
    rigid_body_collisions.Update_Collision_Pair_Helper(*body1,*body2,dt,time,collision_location,collision_normal,
        body1->Pointwise_Object_Velocity(collision_location)-body2->Pointwise_Object_Velocity(collision_location),mpi_one_ghost);
    return true;
}
//#####################################################################
// Function Update_Surface_Collision_Pair
//#####################################################################
template<class TV> bool RIGID_BODY_COLLISIONS<TV>::
Update_Surface_Collision_Pair(const int id_1,const int id_2,const T dt,const T time,const bool mpi_one_ghost)
{
    RIGID_BODY<TV>& body1=rigid_body_collection.Rigid_Body(id_1);RIGID_BODY<TV>& body2=rigid_body_collection.Rigid_Body(id_2);int i;
    for(i=1;i<=collision_pair_iterations;i++){
        T smallest_value;int smallest_index;TV collision_location,collision_normal,collision_relative_velocity;bool ignored_separating=false;
        bool found_intersection=Get_First_Intersection_Point(id_1,id_2,smallest_value,smallest_index,collision_location,collision_normal,
            collision_relative_velocity,true,0,ignored_separating,dt,true);
        collision_callbacks.Swap_States(id_1,id_2);
        if(found_intersection || ignored_separating){
            pairs_processed_by_collisions.Set(VECTOR<int,2>(id_1,id_2).Sorted());}
        if(Update_Collision_Pair_Helper(body1,body2,dt,time,collision_location,collision_normal,collision_relative_velocity,mpi_one_ghost)) return false;}
    return i>1; // if i==1, we didn't process any collisions
}
//#####################################################################
// Function Update_Levelset_Collision_Pair
//#####################################################################
template<class TV> bool RIGID_BODY_COLLISIONS<TV>::
Update_Levelset_Collision_Pair(const int id_1,const int id_2,const T dt,const T time,const bool mpi_one_ghost)
{
    // get all points in levelset, and iterate until each is (at some point) no longer interpenetrating or moving out of the levelset
    ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> > particle_intersections;int i;
    particle_intersections.Preallocate(100);
    for(i=1;i<=collision_pair_iterations;i++){
        T smallest_value;int smallest_index;TV collision_location,collision_normal,collision_relative_velocity;bool ignored_separating;
        bool found_intersection=Get_Deepest_Intersection_Point(id_1,id_2,particle_intersections,smallest_value,smallest_index,collision_location,collision_normal,
            collision_relative_velocity,true,0,ignored_separating);
        if(found_intersection || ignored_separating)
            pairs_processed_by_collisions.Set(VECTOR<int,2>(id_1,id_2).Sorted());
        if(!found_intersection){
            RIGID_BODY<TV>& body_1=rigid_body_collection.Rigid_Body(id_1);
            RIGID_BODY_STATE<TV> saved_state;body_1.Save_State(saved_state);
            collision_callbacks.Restore_Position(id_1);Euler_Step_Position(id_1,dt,time);
            bool ignored_separating_local=false;
            found_intersection=Get_Deepest_Intersection_Point(id_1,id_2,particle_intersections,smallest_value,smallest_index,collision_location,collision_normal,
                collision_relative_velocity,true,0,ignored_separating_local);
            ignored_separating|=ignored_separating_local;
            body_1.Restore_State(saved_state);
            body_1.Update_Bounding_Box();
            if(found_intersection || ignored_separating) pairs_processed_by_collisions.Set(VECTOR<int,2>(id_1,id_2).Sorted());
            if(!found_intersection){skip_collision_check.Set_Last_Checked(id_1,id_2);break;}
            else{
                RIGID_BODY_PARTICLE_INTERSECTION<TV>& intersection=particle_intersections(smallest_index);
                collision_location=rigid_body_collection.Rigid_Body(intersection.particle_body).World_Space_Point(intersection.particle_location);
                if(use_parent_normal){
                    RIGID_BODY<TV>& levelset_body_parent=rigid_body_cluster_bindings.Get_Parent(rigid_body_collection.Rigid_Body(intersection.levelset_body));
                    collision_normal=levelset_body_parent.Implicit_Geometry_Normal(collision_location);}
                    else collision_normal=rigid_body_cluster_bindings.Get_Parent(rigid_body_collection.Rigid_Body(intersection.levelset_body)).Implicit_Geometry_Normal(collision_location);}}
        const RIGID_BODY_PARTICLE_INTERSECTION<TV>& intersection=particle_intersections(smallest_index);
        RIGID_BODY<TV> &body1=rigid_body_collection.Rigid_Body(intersection.particle_body),&body2=rigid_body_collection.Rigid_Body(intersection.levelset_body);
        if(Update_Collision_Pair_Helper(body1,body2,dt,time,collision_location,collision_normal,collision_relative_velocity,mpi_one_ghost)) return false;}
    return i>1; // if i==1, we didn't process any collisions
}
//#####################################################################
// Function Get_Contact_Pairs
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Get_Contact_Pairs(const T dt,const T time,ARRAY<VECTOR<int,2> >& pairs)
{
    LOG::SCOPE("get contact pairs",parameters.threadid);

    typedef typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER T_CLUSTER;
    typedef typename HASHTABLE<int,T_CLUSTER*>::ITERATOR T_REVERSE_BINDING_ITERATOR;
    // advance all bodies using old velocities (and save new states)
    for(int i=1;i<=rigid_body_collection.dynamic_rigid_body_particles.m;i++){
        int id=rigid_body_collection.dynamic_rigid_body_particles(i);
        if(rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Contains(id))
            Euler_Step_Position(id,dt,time);}
    
    // reinitialize in candidate position with old velocities 
    if(parameters.use_projected_gauss_seidel)
        spatial_partition->Set_Collision_Body_Thickness(parameters.contact_proximity);
    else
        spatial_partition->Set_Collision_Body_Thickness(0);
    
    ARRAY<int> colliding_rigid_body_particles;
    if(rigid_body_cluster_bindings.collide_constituent_bodies) 
        for(T_REVERSE_BINDING_ITERATOR i(rigid_body_cluster_bindings.reverse_bindings);i.Valid();i.Next()){
            const T_CLUSTER& bindings=*i.Data();
            for(RIGID_CLUSTER_CONSTITUENT_ID j(1);j<=bindings.children.Size();j++)
                if(rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Contains(bindings.children(j)))
                    colliding_rigid_body_particles.Append(bindings.children(j));}
    for(int n=1;n<=rigid_body_collection.dynamic_rigid_body_particles.m;n++)
        if(rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Contains(rigid_body_collection.dynamic_rigid_body_particles(n)) &&
            (!mpi_rigids || mpi_rigids->Is_Real_Body(rigid_body_collection.dynamic_rigid_body_particles(n))))
            colliding_rigid_body_particles.Append(rigid_body_collection.dynamic_rigid_body_particles(n));
    for(int n=1;n<=colliding_rigid_body_particles.m;n++){int p=colliding_rigid_body_particles(n);
        if(!rigid_body_collection.rigid_geometry_collection.collision_body_list->Has_Collision_Body(p)) continue;
        COLLISION_GEOMETRY_ID id=rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(p);
        RIGID_BODY_STATE<TV> saved_state;rigid_body_collection.Rigid_Body(p).Save_State(saved_state);
        if(parameters.use_ccd) rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies(id)->Save_State(1,time);
        collision_callbacks.Reevolve_Body_With_Saved_State(p,dt,time); // Get X_n+1 
        skip_collision_check.Set_Last_Moved(p);
        spatial_partition->Get_Potential_Collisions_Using_Current_Position(id,object_indices,get_potential_collisions_already_added,false);
        for(int k=1;k<=object_indices.m;k++){
            int j=rigid_body_collection.rigid_geometry_collection.collision_body_list->collision_geometry_id_to_geometry_id.Get(object_indices(k));
            if(!j) continue;
            if(Either_Body_Collides_With_The_Other(p,j) && ((parameters.use_ccd && 
                rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies(id)->Axis_Aligned_Bounding_Box().Intersection(rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies(object_indices(k))->Axis_Aligned_Bounding_Box()))|| 
                    intersections.Bounding_Boxes_Intersect(p,j,parameters.collision_bounding_box_thickness+parameters.use_projected_gauss_seidel?parameters.contact_proximity:0))){ // this should *not* be ids, but should remain indices
                int particle_body,levelset_body;
                if((parameters.use_analytic_collisions || parameters.use_projected_gauss_seidel || intersections.Find_Any_Intersection(p,j,particle_body,levelset_body)))
                    pairs.Append_Unique(VECTOR<int,2>(j,p));}}
        rigid_body_collection.Rigid_Body(p).Restore_State(saved_state);
        if(parameters.use_ccd) rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies(id)->Remove_States();
        rigid_body_collection.Rigid_Body(p).Update_Bounding_Box();} // restore body to its saved position
}
//#####################################################################
// Function Compute_Contact_Graph
//#####################################################################
// Computed by going through each body, advancing it using its new velocity while advancing the rest using the old velocities,
// and seeing which bodies it intersects. Assumes velocities haven't been updated to new velocities yet!
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Compute_Contact_Graph(const T dt,const T time,ARTICULATED_RIGID_BODY<TV>* articulated_rigid_body)
{
    LOG::SCOPE scope("rigid body compute contact graph",parameters.threadid);
    contact_graph.Reset(); // Call in case number of rigid bodies has changed

    ARRAY<VECTOR<int,2> > edge_pairs;
    if(parameters.perform_collisions && !parameters.use_ccd)
        edge_pairs=contact_pairs_from_collisions;
    else{
        Get_Contact_Pairs(dt,time,edge_pairs);
        if(mpi_rigids) mpi_rigids->Exchange_Bounding_Box_Collision_Pairs(rigid_body_collection,edge_pairs,true);}

    if(prune_stacks_from_contact){
        ARRAY<int> body_stack(rigid_body_collection.rigid_body_particle.array_collection->Size());
        HASHTABLE<PAIR<int,int> > stack_static_bodies;
        for(int i=1;i<=contact_stack.m;i++){
            INDIRECT_ARRAY<ARRAY<int>,ARRAY<int>&> contact_subset=body_stack.Subset(contact_stack(i));
            ARRAYS_COMPUTATIONS::Fill(contact_subset,i);
            for(int j=1;j<=contact_stack(i).m;j++) if(rigid_body_collection.Rigid_Body(contact_stack(i)(j)).Has_Infinite_Inertia()) stack_static_bodies.Set(PAIR<int,int>(i,contact_stack(i)(j)));}
        for(COLLISION_GEOMETRY_ID i(1);i<=rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies.Size();i++){
            int rigid_body_id=rigid_body_collection.rigid_geometry_collection.collision_body_list->collision_geometry_id_to_geometry_id.Get(i);
            if(rigid_body_id) if(rigid_body_collection.Is_Active(rigid_body_id) && rigid_body_collection.Rigid_Body(rigid_body_id).Has_Infinite_Inertia()) body_stack(rigid_body_id)=-2;}

        for(int i=1;i<=edge_pairs.m;i++){
            VECTOR<int,2> state(body_stack.Subset(edge_pairs(i)));
            if(state(1)==state(2)) continue;
            if(state(1)==-2 && stack_static_bodies.Contains(PAIR<int,int>(state(2),edge_pairs(i)(1)))) continue;
            if(state(2)==-2 && stack_static_bodies.Contains(PAIR<int,int>(state(1),edge_pairs(i)(2)))) continue;
            for(int j=1;j<=2;j++){int& stack_value=body_stack(edge_pairs(i)(j));if(stack_value!=-2) stack_value=-1;}}
        for(int i=contact_stack.m;i>=1;i--) if(body_stack.Subset(contact_stack(i)).Find(-1)){
            for(int j=1;j<=contact_stack(i).m;j++) if(body_stack(contact_stack(i)(j))!=-2) body_stack(contact_stack(i)(j))=-1;
            contact_stack.Remove_Index_Lazy(i);}
        for(int i=edge_pairs.m;i>=1;i--){
            VECTOR<int,2> state(body_stack.Subset(edge_pairs(i)));
            if(!state.Find(-1) && !state.Find(0)) edge_pairs.Remove_Index_Lazy(i);}}
    if(prune_contact_using_velocity)
        for(int i=1;i<=edge_pairs.m;i++)
            if(!pairs_scale.Contains(edge_pairs(i).Sorted()))
                pairs_scale.Set(edge_pairs(i).Sorted(),TRIPLE<int,int,int>(parameters.contact_iterations,contact_level_iterations,1));
    for(int i=1;i<=edge_pairs.m;i++) if(rigid_body_collection.Is_Active(edge_pairs(i)(1)) && rigid_body_collection.Is_Active(edge_pairs(i)(2))) contact_graph.Add_Edge(edge_pairs(i)(1),edge_pairs(i)(2));

    collision_callbacks.Restore_Positions(); // Get X_n

    contact_graph.directed_graph.Generate_Levels();

    precomputed_contact_pairs_for_level.Resize(contact_graph.Number_Of_Levels());saved_contact_pairs_for_level.Resize(contact_graph.Number_Of_Levels());
    for(int level=1;level<=contact_graph.Number_Of_Levels();level++){
        ARRAY<VECTOR<int,2> >& pairs=precomputed_contact_pairs_for_level(level);pairs.Resize(0);saved_contact_pairs_for_level(level).Resize(0);
        for(int i=1;i<=contact_graph.directed_graph.Nodes_In_Level(level).m;i++){
            int body_id=contact_graph.directed_graph.Nodes_In_Level(level)(i);
            if(!rigid_body_collection.Is_Active(body_id)){
                PHYSBAM_ASSERT(contact_graph.directed_graph.Nodes_In_Level(level).m==1);
                contact_graph.directed_graph.Nodes_In_Level(level).Remove_All();
                continue;}
            for(int j=1;j<=contact_graph.directed_graph.Parents(body_id).m;j++){
                int other_body_id=contact_graph.directed_graph.Parents(body_id)(j);
                int other_body_level=contact_graph.directed_graph.Level_Of_Node(other_body_id);assert(other_body_level<=level);
                // only add the pair once if they are both parents of each other
                if(other_body_level!=level || !contact_graph.directed_graph.Parents(other_body_id).Contains(body_id) || body_id<other_body_id)
                    if(!parameters.perform_collisions || parameters.use_ccd || pairs_processed_by_collisions.Contains(VECTOR<int,2>(body_id,other_body_id).Sorted()))
                        pairs.Append(VECTOR<int,2>(body_id,other_body_id));}}}

    if(articulated_rigid_body) articulated_rigid_body->Generate_Process_List_Using_Contact_Graph(contact_graph);
}
//#####################################################################
// Function Process_Push_Out_Legacy
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Process_Push_Out_Legacy()
{
    LOG::SCOPE scope("rigid body push out legacy",parameters.threadid);
    skip_collision_check.Reset();
    ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> > particle_intersections;particle_intersections.Preallocate(100);
    bool need_another_iteration=true;int iteration=0;
    VECTOR<COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR<TV>*,2> stored_accumulator;
    while(need_another_iteration && ++iteration<=push_out_iterations){
        if(mpi_rigids){
            mpi_rigids->Clear_Impulse_Accumulators(rigid_body_collection);
            mpi_rigid_X_save.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
            mpi_rigid_rotation_save.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
            for(int p=1;p<=rigid_body_collection.rigid_body_particle.array_collection->Size();p++) {
                mpi_rigid_X_save(p)=rigid_body_collection.rigid_body_particle.X(p);
                mpi_rigid_rotation_save(p)=rigid_body_collection.rigid_body_particle.rotation(p);}}

        need_another_iteration=false;
        if(use_freezing_with_push_out)
            for(COLLISION_GEOMETRY_ID i(1);i<=rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies.Size();i++){
                RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_geometry=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies(i));
                if(rigid_collision_geometry){RIGID_BODY<TV>* rigid_body=dynamic_cast<RIGID_BODY<TV>*>(&rigid_collision_geometry->rigid_geometry);
                    if(rigid_body) if(rigid_body_collection.Is_Active(rigid_body->particle_index)) rigid_body_collection.Rigid_Body(rigid_body->particle_index).is_temporarily_static=false;}}
        for(int level=1;level<=contact_graph.Number_Of_Levels();level++){
            ARRAY<VECTOR<int,2> > pairs=precomputed_contact_pairs_for_level(level);
            int iterations=0;bool need_more_iterations=true;T move_fraction=1;
            while(need_more_iterations && ++iterations<=push_out_level_iterations){need_more_iterations=false;
                if(use_gradual_push_out) move_fraction=(T)iterations/push_out_level_iterations;
                for(int i=1;i<=pairs.m;i++){int id_1=pairs(i)(1),id_2=pairs(i)(2);
                    if(parameters.use_projected_gauss_seidel && !pairs_processed_by_contact.Contains(pairs(i).Sorted()))
                        continue;
                    int parent_id_1=rigid_body_cluster_bindings.Get_Parent(rigid_body_collection.Rigid_Body(id_1)).particle_index;
                    int parent_id_2=rigid_body_cluster_bindings.Get_Parent(rigid_body_collection.Rigid_Body(id_2)).particle_index;
                    if(skip_collision_check.Skip_Pair(id_1,id_2)) continue;
                    bool need_more_iterations2=true;int iterations2=0;
                    while(need_more_iterations2 && ++iterations2<=push_out_pair_iterations){need_more_iterations2=false;
                        particle_intersections.Remove_All();
                        intersections.Append_All_Intersections(id_1,id_2,particle_intersections,desired_separation_distance);
                        if(!particle_intersections.m){skip_collision_check.Set_Last_Checked(id_1,id_2);continue;}
                        T smallest_value=FLT_MAX;int smallest_index=0;
                        TV collision_location,collision_normal;T collision_push_distance=0;
                        RIGID_BODY<TV> *collision_body1=0,*collision_body2=0;
                        for(int i=1;i<=particle_intersections.m;i++){
                            const RIGID_BODY_PARTICLE_INTERSECTION<TV>& intersection=particle_intersections(i);
                            T phi=(*rigid_body_collection.Rigid_Body(intersection.levelset_body).implicit_object)(rigid_body_collection.Rigid_Body(intersection.particle_body).World_Space_Point(intersection.particle_location));
                            if(phi<smallest_value){
                                RIGID_BODY<TV> &body1=rigid_body_collection.Rigid_Body(intersection.particle_body),&body2=rigid_body_collection.Rigid_Body(intersection.levelset_body);
                                RIGID_BODY<TV> &body_parent1=rigid_body_cluster_bindings.Get_Parent(body1),&body_parent2=rigid_body_cluster_bindings.Get_Parent(body2);
                                TV location=body1.World_Space_Point(intersection.particle_location);
                                TV normal=use_parent_normal?body_parent2.Implicit_Geometry_Normal(location):body2.Implicit_Geometry_Normal(location);
                                T push_distance=move_fraction*max((T)0,desired_separation_distance-phi);
                                collision_location=location;collision_normal=normal;collision_push_distance=push_distance;
                                collision_body1=&body_parent1;collision_body2=&body_parent2;
                                smallest_value=phi;smallest_index=i;}}
                        RIGID_BODY<TV>::Apply_Push(*collision_body1,*collision_body2,collision_location,collision_normal,collision_push_distance,
                            (mpi_rigids && (mpi_rigids->Is_Dynamic_Ghost_Body(rigid_body_collection.Rigid_Body(id_1)) || 
                                mpi_rigids->Is_Dynamic_Ghost_Body(rigid_body_collection.Rigid_Body(id_2)))));
                        rigid_body_collection.Rigid_Body(particle_intersections(smallest_index).particle_body).Update_Bounding_Box();
                        rigid_body_collection.Rigid_Body(particle_intersections(smallest_index).levelset_body).Update_Bounding_Box();
                        skip_collision_check.Set_Last_Moved(id_1);skip_collision_check.Set_Last_Moved(id_2);
                        if(parent_id_1!=id_1)
                            rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Positions(rigid_body_cluster_bindings.Get_Parent(rigid_body_collection.Rigid_Body(id_1)).particle_index);
                        if(parent_id_2!=id_2)
                            rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Positions(rigid_body_cluster_bindings.Get_Parent(rigid_body_collection.Rigid_Body(id_2)).particle_index);
                        need_more_iterations2=true;need_another_iteration=true;need_more_iterations=true;}}}
            if(use_freezing_with_push_out) for(int i=1;i<=contact_graph.directed_graph.Nodes_In_Level(level).m;i++)
                rigid_body_collection.Rigid_Body(contact_graph.directed_graph.Nodes_In_Level(level)(i)).is_temporarily_static=true;}

        if(mpi_rigids){
            int need_another_iteration_int=(int)need_another_iteration;
            need_another_iteration=mpi_rigids->Reduce_Max(need_another_iteration_int)>0?true:false;
            mpi_rigids->Exchange_All_Pushes(rigid_body_collection,mpi_rigid_X_save,mpi_rigid_rotation_save,*this);}}
    if(use_freezing_with_push_out) for(COLLISION_GEOMETRY_ID i(1);i<=rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies.Size();i++){
        RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_geometry=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies(i));
        if(rigid_collision_geometry){RIGID_BODY<TV>* rigid_body=dynamic_cast<RIGID_BODY<TV>*>(&rigid_collision_geometry->rigid_geometry);
            if(rigid_body) if(rigid_body_collection.Is_Active(rigid_body->particle_index)) rigid_body_collection.Rigid_Body(rigid_body->particle_index).is_temporarily_static=false;}}
}
//#####################################################################
// Function Process_Push_Out_Projected_Gauss_Seidel
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Process_Push_Out_Projected_Gauss_Seidel()
{
    PHYSBAM_ASSERT(parameters.use_projected_gauss_seidel);
    SOLVE_CONTACT::Push_Out(collision_callbacks,rigid_body_collection,parameters,pairs_processed_by_contact);
}
//#####################################################################
// Function Process_Push_Out
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Process_Push_Out(const bool perform_collision_body_collisions,const T residual_push_out_depth)
{
    LOG::SCOPE scope("rigid body push out new",parameters.threadid);
    const bool kinematic_rigid_bodies_only=rigid_body_collection.simulated_rigid_body_particles.m==0;

    skip_collision_check.Reset();
    ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> > particle_intersections;
    particle_intersections.Preallocate(100);

    bool need_another_iteration=true;int iteration=0;
    while(need_another_iteration && ++iteration<=push_out_iterations){need_another_iteration=false;
        if(!kinematic_rigid_bodies_only){
            if(use_freezing_with_push_out) Clear_Temporarily_Static();
            for(int level=1;level<=contact_graph.Number_Of_Levels();level++){
                const ARRAY<int>& rigid_bodies_in_level=contact_graph.directed_graph.Nodes_In_Level(level);
                int level_iteration=0;bool need_more_level_iterations=true;T move_fraction=1;
                while(need_more_level_iterations && ++level_iteration<=push_out_level_iterations){need_more_level_iterations=false;
                    if(use_gradual_push_out) move_fraction=(T)level_iteration/push_out_level_iterations;
                    for(int i=1;i<=rigid_bodies_in_level.m;i++){
                        bool need_more_pair_iterations=true;int pair_iteration=0;
                        while(need_more_pair_iterations && ++pair_iteration<=push_out_pair_iterations){
                            need_more_pair_iterations=false;
                            if(Push_Out_From_Rigid_Body(rigid_body_collection.Rigid_Body(rigid_bodies_in_level(i)),particle_intersections,move_fraction,residual_push_out_depth)){
                                need_more_pair_iterations=need_more_level_iterations=need_another_iteration=true;
                                rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Positions();}}}}
                if(use_freezing_with_push_out) Set_Level_Temporarily_Static(level);}
            if(use_freezing_with_push_out) Clear_Temporarily_Static();}}

    // Process static/kinematic rigid bodies
    if(kinematic_rigid_bodies_only && perform_collision_body_collisions)
        for(int i=1;i<=rigid_body_collection.static_and_kinematic_rigid_bodies.m;i++){
            RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(rigid_body_collection.static_and_kinematic_rigid_bodies(i));
            if(Push_Out_From_Rigid_Body(rigid_body,particle_intersections,(T)1,residual_push_out_depth)) rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Positions();}
}
//#####################################################################
// Function Get_Rigid_Bodies_Intersecting_Rigid_Body
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Get_Rigid_Bodies_Intersecting_Rigid_Body(const int particle_index,ARRAY<int>& rigid_bodies,ARRAY<TV>& collision_locations,ARRAY<TV>& body_distances,
    ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const T residual_push_out_depth) const
{
    const RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(particle_index);
    const int rigid_body_id=rigid_body.particle_index;
    const int rigid_body_level=contact_graph.directed_graph.Level_Of_Node(rigid_body_id);

    const ARRAY<VECTOR<int,2> >& level_contact_pairs=precomputed_contact_pairs_for_level(rigid_body_level);
    for(int i=1;i<=level_contact_pairs.m;i++) if(level_contact_pairs(i)(1)==rigid_body_id || level_contact_pairs(i)(2)==rigid_body_id){
        const int other_body_id=level_contact_pairs(i)(1)==rigid_body_id?level_contact_pairs(i)(2):level_contact_pairs(i)(1);
        const RIGID_BODY<TV>& other_rigid_body=rigid_body_collection.Rigid_Body(other_body_id);
        if(skip_collision_check.Skip_Pair(rigid_body_id,other_body_id)
            || (collision_manager && !collision_manager->Either_Body_Collides_With_The_Other(other_rigid_body.particle_index,rigid_body.particle_index))) continue;
        T smallest_value;int smallest_index;TV collision_location,collision_normal,collision_relative_velocity;bool ignored_separating;
        if(Get_Deepest_Intersection_Point(rigid_body_id,other_body_id,particle_intersections,
                smallest_value,smallest_index,collision_location,collision_normal,collision_relative_velocity,false,
                desired_separation_distance-residual_push_out_depth,ignored_separating)){
            // reverse the direction of separation distance to be relative to rigid_body_id
            if(particle_intersections(smallest_index).levelset_body!=rigid_body_id) smallest_value=-smallest_value;
            rigid_bodies.Append(other_rigid_body.particle_index);collision_locations.Append(collision_location);
            smallest_value-=residual_push_out_depth;
            body_distances.Append(smallest_value*collision_normal);}
        else skip_collision_check.Set_Last_Checked(rigid_body_id,other_body_id);} // set last checked even if ignored_separating==true
}
//}################################################################
// Function Push_Out_From_Rigid_Body
//#####################################################################
template<class TV> bool RIGID_BODY_COLLISIONS<TV>::
Push_Out_From_Rigid_Body(RIGID_BODY<TV>& rigid_body,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const T move_fraction,const T residual_push_out_depth)
{
    RIGID_BODY<TV>& parent_rigid_body=rigid_body_collection.rigid_body_cluster_bindings.Get_Parent(rigid_body);
    T threshold=rigid_body.Length_Scale_Squared()*(T)1e-2;

    ARRAY<int> particle_interactions,rigid_body_interactions;ARRAY<TV> particle_distances,rigid_body_collision_locations,rigid_body_distances;
    if(rigid_body_collection.simulated_rigid_body_particles.m)
        Get_Rigid_Bodies_Intersecting_Rigid_Body(rigid_body.particle_index,rigid_body_interactions,rigid_body_collision_locations,rigid_body_distances,particle_intersections,residual_push_out_depth);
    particle_distances*=move_fraction;rigid_body_distances*=move_fraction;

    if(rigid_body_interactions.m==0) return false;

    ARRAY<SYMMETRIC_MATRIX<T,TV::dimension> > K_inverse(rigid_body_interactions.m);
    for(int i=1;i<=rigid_body_interactions.m;i++){
        RIGID_BODY<TV>& parent_other_rigid_body=rigid_body_collection.rigid_body_cluster_bindings.Get_Parent(rigid_body_collection.Rigid_Body(rigid_body_interactions(i)));
        if(!parent_other_rigid_body.Has_Infinite_Inertia() || use_static_body_masses) K_inverse(i)=parent_other_rigid_body.Impulse_Factor(rigid_body_collision_locations(i)).Inverse();
        else K_inverse(i)=SYMMETRIC_MATRIX<T,TV::dimension>::Identity_Matrix();}

    TV velocity;T_SPIN angular_velocity;
    if(!parent_rigid_body.Has_Infinite_Inertia()){
        // TODO: static deformable particles
        TV ms[2];VECTOR<T,T_SPIN::dimension> mrs[2];MATRIX<T,T_SPIN::dimension,TV::dimension> mr[2];MATRIX<T,T_SPIN::dimension> mrr[2];
        T m=0; // non static particles
        
        SYMMETRIC_MATRIX<T,TV::dimension> K_inverse_sum[2]; // non static bodies

        int number_of_static_bodies=0;TV centroid;
        for(int i=1;i<=rigid_body_interactions.m;i++){
            if(parameters.use_projected_gauss_seidel && !pairs_processed_by_contact.Contains(VECTOR<int,2>(rigid_body.particle_index,rigid_body_interactions(i)).Sorted()))
                continue;
            RIGID_BODY<TV>& parent_other_rigid_body=rigid_body_collection.rigid_body_cluster_bindings.Get_Parent(rigid_body_collection.Rigid_Body(rigid_body_interactions(i)));
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
            SYMMETRIC_MATRIX<T,TV::dimension> R_R_transpose;
            for(int i=1;i<=rigid_body_interactions.m;i++){
                RIGID_BODY<TV>& parent_other_rigid_body=rigid_body_collection.rigid_body_cluster_bindings.Get_Parent(rigid_body_collection.Rigid_Body(rigid_body_interactions(i)));
                if(parent_other_rigid_body.Has_Infinite_Inertia()) R_R_transpose+=SYMMETRIC_MATRIX<T,TV::dimension>::Outer_Product(rigid_body_collision_locations(i)-parent_rigid_body.X()-centroid);}
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

    // apply push to other rigid bodies
    for(int i=1;i<=rigid_body_interactions.m;i++) if(!rigid_body_collection.Rigid_Body(rigid_body_interactions(i)).Has_Infinite_Inertia()){
        RIGID_BODY<TV>& other_rigid_body=rigid_body_collection.Rigid_Body(rigid_body_interactions(i));
        RIGID_BODY<TV>& parent_other_rigid_body=rigid_body_collection.rigid_body_cluster_bindings.Get_Parent(other_rigid_body);
        TV impulse=K_inverse(i)*(-rigid_body_distances(i)+velocity+TV::Cross_Product(angular_velocity,rigid_body_collision_locations(i)-parent_rigid_body.X()));
        T_SPIN other_angular_velocity=parent_other_rigid_body.World_Space_Inertia_Tensor_Inverse()*TV::Cross_Product(rigid_body_collision_locations(i)-parent_other_rigid_body.X(),impulse);
        parent_other_rigid_body.X()+=impulse/parent_other_rigid_body.Mass();
        if(parameters.use_push_out_rotation)
            parent_other_rigid_body.Rotation()=ROTATION<TV>::From_Rotation_Vector(other_angular_velocity)*parent_other_rigid_body.Rotation();parent_other_rigid_body.Rotation().Normalize();
        parent_other_rigid_body.Update_Angular_Velocity();
        parent_other_rigid_body.Update_Bounding_Box();
        skip_collision_check.Set_Last_Moved(other_rigid_body.particle_index);}

    // apply push to the rigid body
    if(!parent_rigid_body.Has_Infinite_Inertia()){
        parent_rigid_body.X()+=velocity;
        if(parameters.use_push_out_rotation)
            parent_rigid_body.Rotation()=ROTATION<TV>::From_Rotation_Vector(angular_velocity)*parent_rigid_body.Rotation();parent_rigid_body.Rotation().Normalize();
        parent_rigid_body.Update_Angular_Velocity();
        parent_rigid_body.Update_Bounding_Box();
        skip_collision_check.Set_Last_Moved(rigid_body.particle_index);}

    return true;
}
    
template<> bool RIGID_BODY_COLLISIONS<VECTOR<float,1> >::
Push_Out_From_Rigid_Body(RIGID_BODY<VECTOR<float,1> >& rigid_body,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<VECTOR<float,1> > >& particle_intersections,const float move_fraction,
    const float residual_push_out_depth)
{
    PHYSBAM_NOT_IMPLEMENTED();
};
template<> bool RIGID_BODY_COLLISIONS<VECTOR<double,1> >::
Push_Out_From_Rigid_Body(RIGID_BODY<VECTOR<double,1> >& rigid_body,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<VECTOR<double,1> > >& particle_intersections,const double move_fraction,
    const double residual_push_out_depth)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Get_Rigid_Body_Depth
//#####################################################################
template<class TV> int  RIGID_BODY_COLLISIONS<TV>::
Get_Rigid_Body_Depth(int i,ARRAY<int,int>& depths)
{
    if(depths(i)!=-1) return depths(i);
    else{
        depths(i)=0;
        for(int child=1;child<=contact_graph.directed_graph.Children(i).m;child++)
            depths(i)+=Get_Rigid_Body_Depth(contact_graph.directed_graph.Children(i)(child),depths);
        depths(i)++;
        return depths(i);}
}
//#####################################################################
// Function Compute_Contact_Frequency
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Compute_Contact_Frequency()
{
    T linear_threshold_squared=(T)1e-3,angular_threshold_squared=(T)1e-3;
    int stack_depth=2;

    ARRAY<ARRAY<int>,int> adj(contact_graph.directed_graph.Number_Nodes());
    for(int i(1);i<=adj.m;i++) for(int j=1;j<=contact_graph.directed_graph.Parents(i).m;j++){
        adj(i).Append(contact_graph.directed_graph.Parents(i)(j));adj(contact_graph.directed_graph.Parents(i)(j)).Append(i);}

    ARRAY<int,int> depths(adj.m);ARRAYS_COMPUTATIONS::Fill(depths,-1);
    for(int i(1);i<=adj.m;i++) if(depths(i)==-1) Get_Rigid_Body_Depth(i,depths);

    pairs_scale.Remove_All();
    for(int i(1);i<=adj.m;i++){
        for(int j=1;j<=adj(i).m;j++){
            VECTOR<int,2> sorted_pair=VECTOR<int,2>(i,adj(i)(j)).Sorted();
            if(!pairs_scale.Contains(sorted_pair)){
                TWIST<TV> relative_twist=rigid_body_collection.Rigid_Body(i).Twist()-rigid_body_collection.Rigid_Body(adj(i)(j)).Twist();
                if(relative_twist.linear.Magnitude_Squared() <= linear_threshold_squared && relative_twist.angular.Magnitude_Squared() <= angular_threshold_squared){
                    T scale=(T).5*((relative_twist.linear.Magnitude_Squared()/linear_threshold_squared)+
                        (relative_twist.angular.Magnitude_Squared()/angular_threshold_squared));
                    if(depths(i) > stack_depth || depths(adj(i)(j)) > stack_depth) scale=0;
                    pairs_scale.Set(sorted_pair,TRIPLE<int,int,int>((int)(scale*(parameters.contact_iterations-1)+1),(int)(scale*(contact_level_iterations-1)+1),
                            (int)((1-scale)*(contact_pair_iterations-1)+1)));}
                else pairs_scale.Set(sorted_pair,TRIPLE<int,int,int>(parameters.contact_iterations,contact_level_iterations,1));}}}
}
//#####################################################################
// Function Construct_Stacks
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Construct_Stacks()
{
    T linear_threshold_squared=(T)1e-4,angular_threshold_squared=(T)1e-4;

    ARRAY<ARRAY<int>,int> adj(contact_graph.directed_graph.Number_Nodes());
    for(int i(1);i<=adj.m;i++) for(int j=1;j<=contact_graph.directed_graph.Parents(i).m;j++){
        adj(i).Append(contact_graph.directed_graph.Parents(i)(j));adj(contact_graph.directed_graph.Parents(i)(j)).Append(i);}

    ARRAY<int,int> visited(rigid_body_collection.rigid_body_particle.array_collection->Size());

    STACK<int> stack;
    for(int i(1);i<=visited.m;i++) if(!visited(i)){
        ARRAY<int> list;
        visited(i)=i;
        stack.Push(i);
        bool has_infinite_inertia_body=false;
        while(!stack.Empty()){
            int j=stack.Pop();
            list.Append(j);
            if(rigid_body_collection.Rigid_Body(j).Has_Infinite_Inertia()){has_infinite_inertia_body=true;continue;}
            for(int k=1;k<=adj(j).m;k++) if(visited(adj(j)(k))!=i){
                visited(adj(j)(k))=i;
                stack.Push(adj(j)(k));}}
        if(!has_infinite_inertia_body || list.m<=1) continue;
        for(int j=2;j<=list.m;j++) if(rigid_body_collection.Rigid_Body(list(j)).Has_Infinite_Inertia()){exchange(list(j),list(1));break;}
        const TWIST<TV>& base_twist=rigid_body_collection.Rigid_Body(list(1)).Twist();
        bool same_velocity=true;
        for(int j=2;j<=list.m;j++){
            TWIST<TV> relative_twist=rigid_body_collection.Rigid_Body(list(j)).Twist()-base_twist;
            if(relative_twist.linear.Magnitude_Squared() > linear_threshold_squared || relative_twist.angular.Magnitude_Squared() > angular_threshold_squared){
                same_velocity=false;break;}}
        if(!same_velocity) continue;
        contact_stack.Append(ARRAY<int>());
        for(int j=1;j<=list.m;j++) contact_stack.Last().Append(list(j));}
}
//#####################################################################
// Function Check_For_Any_Interpenetration
//#####################################################################
template<class TV> bool RIGID_BODY_COLLISIONS<TV>::
Check_For_Any_Interpenetration()
{
    spatial_partition->Set_Collision_Body_Thickness(0);
    for(COLLISION_GEOMETRY_ID i(1);i<=rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies.Size();i++){
        RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_geometry=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies(i));
        if(rigid_collision_geometry){RIGID_BODY<TV>* rigid_body=dynamic_cast<RIGID_BODY<TV>*>(&rigid_collision_geometry->rigid_geometry);
            if(rigid_body) if(rigid_body_collection.Is_Active(rigid_body->particle_index) && rigid_body_collection.Rigid_Body(rigid_body->particle_index).Is_Simulated()){
            int particle_body(0),levelset_body(0);
            spatial_partition->Get_Potential_Collisions(i,object_indices,get_potential_collisions_already_added,false);
            for(int t=1;t<=object_indices.m;t++) if(object_indices(t)){
                int id=rigid_body_collection.rigid_geometry_collection.collision_body_list->collision_geometry_id_to_geometry_id.Get(object_indices(t));
                if(id && Either_Body_Collides_With_The_Other(rigid_body->particle_index,id) &&
                    intersections.Intersection_Check(rigid_body->particle_index,id,particle_body,levelset_body)){
                    {std::stringstream ss;ss<<"!!!! Interpenetration detected: id:"<<particle_body<<" \""<<rigid_body_collection.Rigid_Body(particle_body).name<<"\" point inside id: "<<levelset_body<<" \""
                             <<rigid_body_collection.Rigid_Body(levelset_body).name<<"\" levelset"<<std::endl;LOG::filecout(ss.str());}
                    return true;}}}}}
    return false;
}
//#####################################################################
// Function Print_Interpenetration_Statistics
//#####################################################################
template<class TV> ARRAY<VECTOR<int,2> > RIGID_BODY_COLLISIONS<TV>::
Find_All_Bounding_Box_Pairs(const T thickness)
{
    spatial_partition->Reinitialize();
    ARRAY<VECTOR<int,2> > pairs;
    for(COLLISION_GEOMETRY_ID i(1);i<=rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies.Size();i++){
        RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_geometry=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies(i));
        if(rigid_collision_geometry){RIGID_BODY<TV>* rigid_body=dynamic_cast<RIGID_BODY<TV>*>(&rigid_collision_geometry->rigid_geometry);if(rigid_body){
                spatial_partition->Get_Potential_Collisions(i,object_indices,get_potential_collisions_already_added);
            for(int t=1;t<=object_indices.m;t++){int j=rigid_body_collection.rigid_geometry_collection.collision_body_list->collision_geometry_id_to_geometry_id.Get(object_indices(t));
                if(j && Either_Body_Collides_With_The_Other(rigid_body->particle_index,j) && intersections.Bounding_Boxes_Intersect(rigid_body->particle_index,j,thickness))
                    pairs.Append(VECTOR<int,2>(j,rigid_body->particle_index));}}}}
    return pairs;
}
//#####################################################################
// Function Print_Interpenetration_Statistics
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Print_Interpenetration_Statistics()
{
    spatial_partition->Set_Collision_Body_Thickness(0);

    int num_interpenetration_points=0;int worst_id_1(0),worst_id_2(0);
    T smallest_phi=FLT_MAX,avg_phi=0;
    ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> > particle_intersections;
    particle_intersections.Preallocate(1000);

    for(COLLISION_GEOMETRY_ID i(1);i<=rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies.Size();i++){
        RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_geometry=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies(i));
        if(rigid_collision_geometry){RIGID_BODY<TV>* rigid_body=dynamic_cast<RIGID_BODY<TV>*>(&rigid_collision_geometry->rigid_geometry);
            if(rigid_body) if(rigid_body_collection.Is_Active(rigid_body->particle_index) && rigid_body_collection.Rigid_Body(rigid_body->particle_index).Is_Simulated()){
                    spatial_partition->Get_Potential_Collisions(i,object_indices,get_potential_collisions_already_added);
            for(int t=1;t<=object_indices.m;t++){int j=rigid_body_collection.rigid_geometry_collection.collision_body_list->collision_geometry_id_to_geometry_id.Get(object_indices(t));
                if(j && Either_Body_Collides_With_The_Other(rigid_body->particle_index,j) && intersections.Bounding_Boxes_Intersect(rigid_body->particle_index,j))
                    intersections.Append_All_Intersections(rigid_body->particle_index,j,particle_intersections);}}}}
    num_interpenetration_points=particle_intersections.m;
    for(int i=1;i<=particle_intersections.m;i++){
        const RIGID_BODY_PARTICLE_INTERSECTION<TV>& intersection=particle_intersections(i);
        TV location=rigid_body_collection.Rigid_Body(intersection.particle_body).World_Space_Point(intersection.particle_location);
        T phi=(*rigid_body_collection.Rigid_Body(intersection.levelset_body).implicit_object)(location);
        if(phi<smallest_phi){worst_id_1=intersection.particle_body;worst_id_2=intersection.levelset_body;smallest_phi=phi;}
        avg_phi+=phi;}
    if(num_interpenetration_points>0){
        avg_phi/=num_interpenetration_points;
        {std::stringstream ss;ss<<"Interpenetration statistics: num interpenetrating points = "<<num_interpenetration_points<<", smallest phi = "
           <<smallest_phi<<"(\""<<rigid_body_collection.Rigid_Body(worst_id_1).name<<"\" and \""<<rigid_body_collection.Rigid_Body(worst_id_2).name<<"\")"
           <<", avg phi = "<<avg_phi<<std::endl;LOG::filecout(ss.str());}}
    else {std::stringstream ss;ss<<"No interpenetration"<<std::endl;LOG::filecout(ss.str());}
}
//#####################################################################
// Function Get_Bounding_Box_Collision_Pairs_Of_Body
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Get_Bounding_Box_Collision_Pairs_Of_Body(ARRAY<VECTOR<int,2> >& pairs,int id,const T thickness)
{
    spatial_partition->Get_Potential_Collisions(rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(id),object_indices,get_potential_collisions_already_added,false);
    for(int t=1;t<=object_indices.m;t++){int j=rigid_body_collection.rigid_geometry_collection.collision_body_list->collision_geometry_id_to_geometry_id.Get(object_indices(t));
        if(!j || !Either_Body_Collides_With_The_Other(id,j)) continue;
        if(!rigid_body_collection.Rigid_Body(j).Has_Infinite_Inertia()){
            if(j<id) continue;} // to avoid processing the same pair of dynamic bodies twice
        if(!skip_collision_check.Skip_Pair(id,j)){
            if(intersections.Bounding_Boxes_Intersect(id,j,thickness)) pairs.Append(VECTOR<int,2>(id,j));
            else skip_collision_check.Set_Last_Checked(id,j);}}
}
//#####################################################################
// Function Get_Bounding_Box_Collision_Pairs
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Get_Bounding_Box_Collision_Pairs(const T dt,const T time,ARRAY<VECTOR<int,2> >& pairs,const bool add_contact_pairs,const bool reinitialize_spatial_partition,const T thickness)
{
    typedef typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER T_CLUSTER;
    typedef typename HASHTABLE<int,T_CLUSTER*>::ITERATOR T_REVERSE_BINDING_ITERATOR;
    pairs.Remove_All();

    if(reinitialize_spatial_partition)
        spatial_partition->Set_Collision_Body_Thickness(0);

    ARRAY<int> colliding_rigid_body_particles;
    if(rigid_body_cluster_bindings.collide_constituent_bodies) 
        for(T_REVERSE_BINDING_ITERATOR i(rigid_body_cluster_bindings.reverse_bindings);i.Valid();i.Next()){
            const T_CLUSTER& bindings=*i.Data();
            for(RIGID_CLUSTER_CONSTITUENT_ID j(1);j<=bindings.children.Size();j++)
                if(rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Contains(bindings.children(j)))
                    colliding_rigid_body_particles.Append(bindings.children(j));}
    for(int n=1;n<=rigid_body_collection.dynamic_rigid_body_particles.m;n++){
        int p=rigid_body_collection.dynamic_rigid_body_particles(n);
        if(rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Contains(p) && (!mpi_rigids || mpi_rigids->Is_Real_Body(p)))
            colliding_rigid_body_particles.Append(p);}
    
    for(int n=1;n<=colliding_rigid_body_particles.m;n++){int p=colliding_rigid_body_particles(n);
        Get_Bounding_Box_Collision_Pairs_Of_Body(pairs,p,thickness);}
    if(mpi_rigids) mpi_rigids->Exchange_Bounding_Box_Collision_Pairs(rigid_body_collection,pairs,false);

    if(add_contact_pairs){
        collision_callbacks.Restore_Positions();
        contact_pairs_from_collisions.Remove_All();
        Get_Contact_Pairs(dt,time,contact_pairs_from_collisions);
        if(mpi_rigids) mpi_rigids->Exchange_Bounding_Box_Collision_Pairs(rigid_body_collection,contact_pairs_from_collisions,true);
        pairs.Append_Elements(contact_pairs_from_collisions);
        for(int n=1;n<=rigid_body_collection.dynamic_rigid_body_particles.m;n++)
            collision_callbacks.Euler_Step_Position_With_New_Velocity(rigid_body_collection.dynamic_rigid_body_particles(n),dt,time);}
}
//#####################################################################
// Function Register_Analytic_Collisions
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Register_Analytic_Collisions()
{
    const char* sphere=typeid(ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >).name();
    const char* box=typeid(ANALYTIC_IMPLICIT_OBJECT<BOX<TV> >).name();
    const char* plane=typeid(ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<TV> >).name();
    analytic_collision_registry.Set(VECTOR<std::string,2>(sphere,sphere),Update_Sphere_Sphere_Collision<TV>);
    analytic_collision_registry.Set(VECTOR<std::string,2>(sphere,plane).Sorted(),Update_Sphere_Plane_Collision<TV>);
    analytic_collision_registry.Set(VECTOR<std::string,2>(box,box),Update_Box_Box_Collision<TV>);
    analytic_collision_registry.Set(VECTOR<std::string,2>(box,plane).Sorted(),Update_Box_Plane_Collision<TV>);
    SOLVE_CONTACT::Register_Analytic_Contacts<TV>(analytic_contact_registry);
}
//#####################################################################
// Function Apply_Stacking_Contact
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Apply_Stacking_Contact()
{
    for(int i=1;i<=contact_stack.m;i++){
        TV& static_V=rigid_body_collection.rigid_body_particle.V(contact_stack(i)(1));
        T_SPIN& static_angular_velocity=rigid_body_collection.rigid_body_particle.angular_velocity(contact_stack(i)(1));
        for(int j=2;j<=contact_stack(i).m;j++){
            RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(contact_stack(i)(j));
            collision_callbacks.Restore_Position(rigid_body.particle_index);
            rigid_body.V()=static_V;rigid_body.Angular_Velocity()=static_angular_velocity;rigid_body.Update_Angular_Momentum();}}
}
//#####################################################################
// Function Add_Elastic_Collisions
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Add_Elastic_Collisions(const T dt,const T time)
{
    skip_collision_check.Reset();pairs_processed_by_collisions.Remove_All();

    LOG::SCOPE scope("rigid body collisions",parameters.threadid);
    skip_collision_check.Reset();pairs_processed_by_collisions.Remove_All();
    rigid_body_particle_intersections.Remove_All();
    VECTOR<COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR<TV>*,2> stored_accumulator;
    ARRAY<VECTOR<int,2> > pairs;pairs.Preallocate(10);bool need_another_iteration=true;
    for(int i=1;i<=parameters.collision_iterations && need_another_iteration;i++){
        if(mpi_rigids){
            mpi_rigids->Clear_Impulse_Accumulators(rigid_body_collection);
            mpi_rigid_velocity_save.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
            mpi_rigid_angular_momentum_save.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
            for(int p=1;p<=rigid_body_collection.rigid_body_particle.array_collection->Size();p++) {
                mpi_rigid_velocity_save(p).linear=rigid_body_collection.rigid_body_particle.V(p);
                mpi_rigid_velocity_save(p).angular=rigid_body_collection.rigid_body_particle.angular_velocity(p);
                mpi_rigid_angular_momentum_save(p)=rigid_body_collection.rigid_body_particle.angular_momentum(p);}}

        need_another_iteration=false;Get_Bounding_Box_Collision_Pairs(dt,time,pairs,i==parameters.collision_iterations,i==1);
        for(int j=1;j<=pairs.m;j++){
            int id_1=pairs(j)(1),id_2=pairs(j)(2);
            bool mpi_one_ghost=mpi_rigids && (mpi_rigids->Is_Dynamic_Ghost_Body(rigid_body_collection.Rigid_Body(id_1)) || mpi_rigids->Is_Dynamic_Ghost_Body(rigid_body_collection.Rigid_Body(id_2)));
            if(Update_Collision_Pair(id_1,id_2,dt,time,mpi_one_ghost)){
                spatial_partition->Update_Body(rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(id_1));
                spatial_partition->Update_Body(rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(id_2));
                need_another_iteration=true;}}

        if(mpi_rigids){
            int need_another_iteration_int=(int)need_another_iteration;
            need_another_iteration=mpi_rigids->Reduce_Max(need_another_iteration_int)>0?true:false;}

        // force it to the last iteration (so it picks up contact pairs one time at least)
        if(!need_another_iteration && i<parameters.collision_iterations){
            i=parameters.collision_iterations-1;
            need_another_iteration=true;}

        if(mpi_rigids)
            mpi_rigids->Exchange_All_Impulses(rigid_body_collection,mpi_rigid_velocity_save,mpi_rigid_angular_momentum_save,*this,true,dt,time);}
}
//#####################################################################
// Function Process_Contact_Using_Graph
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Process_Contact_Using_Graph(const T dt,const T time,ARTICULATED_RIGID_BODY<TV>* articulated_rigid_body,const bool correct_contact_energy,const bool use_saved_pairs)
{
    LOG::SCOPE scope("rigid body contact",parameters.threadid);

    SOLVE_CONTACT::Solve(*this,collision_callbacks,rigid_body_collection,parameters,correct_contact_energy,use_saved_pairs,dt,time,mpi_rigids,mpi_rigid_velocity_save,mpi_rigid_angular_momentum_save);
}
//#####################################################################
// Function Shock_Propagation_Using_Graph
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Shock_Propagation_Using_Graph(const T dt,const T time,ARTICULATED_RIGID_BODY<TV>* articulated_rigid_body,const bool use_saved_pairs)
{
    LOG::SCOPE scope("shock propagation",parameters.threadid);
    skip_collision_check.Reset();bool need_another_iteration=true;int iteration=0;T epsilon_scale=1;
    ARRAY<ARRAY<VECTOR<int,2> > >& contact_pairs_for_level=use_saved_pairs?saved_contact_pairs_for_level:precomputed_contact_pairs_for_level;
    VECTOR<COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR<TV>*,2> stored_accumulator;
    while(need_another_iteration && ++iteration<=shock_propagation_iterations){
        if(mpi_rigids){
            mpi_rigids->Clear_Impulse_Accumulators(rigid_body_collection);
            mpi_rigid_velocity_save.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
            mpi_rigid_angular_momentum_save.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
            for(int p=1;p<=rigid_body_collection.rigid_body_particle.array_collection->Size();p++) {
                mpi_rigid_velocity_save(p).linear=rigid_body_collection.rigid_body_particle.V(p);
                mpi_rigid_velocity_save(p).angular=rigid_body_collection.rigid_body_particle.angular_velocity(p);
                mpi_rigid_angular_momentum_save(p)=rigid_body_collection.rigid_body_particle.angular_momentum(p);}}
        
        need_another_iteration=false;
        Clear_Temporarily_Static();
        for(int level=1;level<=contact_graph.Number_Of_Levels();level++){
            ARRAY<VECTOR<int,2> >& pairs=contact_pairs_for_level(level);
            bool need_another_level_iteration=true;int level_iteration=0;
            while(need_another_level_iteration && ++level_iteration<=shock_propagation_level_iterations){need_another_level_iteration=false;
                if(parameters.use_epsilon_scaling_for_level) epsilon_scale=(T)level_iteration/shock_propagation_level_iterations;
                for(int i=1;i<=pairs.m;i++){int id_1=pairs(i)(1),id_2=pairs(i)(2);
                    if(skip_collision_check.Skip_Pair(id_1,id_2)) continue;
                    if(!parameters.use_projected_gauss_seidel || pairs_processed_by_contact.Contains(pairs(i).Sorted()))
                        if(SOLVE_CONTACT::Update_Contact_Pair(*this,collision_callbacks,analytic_contact_registry,id_1,id_2,false,shock_propagation_pair_iterations,epsilon_scale,dt,time,
                                (mpi_rigids && (mpi_rigids->Is_Dynamic_Ghost_Body(rigid_body_collection.Rigid_Body(id_1)) ||
                                    mpi_rigids->Is_Dynamic_Ghost_Body(rigid_body_collection.Rigid_Body(id_2)))))){
                            if(!use_saved_pairs) saved_contact_pairs_for_level(level).Append_Unique(pairs(i));need_another_level_iteration=true;need_another_iteration=true;}}
                if(articulated_rigid_body && articulated_rigid_body->use_shock_propagation)
                    for(int i=1;i<=articulated_rigid_body->shock_propagation_level_iterations;i++) for(int j=1;j<=articulated_rigid_body->process_list(level).m;j++){
                        JOINT_ID joint_id=articulated_rigid_body->process_list(level)(j);
                        RIGID_BODY<TV> &parent=*articulated_rigid_body->Parent(joint_id),&child=*articulated_rigid_body->Child(joint_id);
                        bool old_temporarily_static_parent=parent.is_temporarily_static,old_temporarily_static_child=child.is_temporarily_static;
                        parent.is_temporarily_static=false;child.is_temporarily_static=false;
                        Apply_Prestabilization_To_Joint(dt,time,*articulated_rigid_body,joint_id,epsilon_scale);
                        parent.is_temporarily_static=old_temporarily_static_parent;child.is_temporarily_static=old_temporarily_static_child;
                        need_another_level_iteration=need_another_iteration=true;}}
            // note that the bodies involved in joints get frozen too, since both bodies in a joint must be at least this level or lower to get processed
            Set_Level_Temporarily_Static(level);}

        if(mpi_rigids){
            int need_another_iteration_int=(int)need_another_iteration;
            need_another_iteration=mpi_rigids->Reduce_Max(need_another_iteration_int)>0?true:false;
            mpi_rigids->Exchange_All_Impulses(rigid_body_collection,mpi_rigid_velocity_save,mpi_rigid_angular_momentum_save,*this,false,dt,time);}}

    Clear_Temporarily_Static();

    // enforce articulation constraints once more
    if(articulated_rigid_body && !articulated_rigid_body->use_shock_propagation && articulated_rigid_body->do_final_pass){
        articulated_rigid_body->Store_Velocities_And_Momenta();
        for(int i=1;i<=articulated_rigid_body->shock_propagation_level_iterations;i++)
            for(int level=1;level<=contact_graph.Number_Of_Levels();level++) for(int j=1;j<=articulated_rigid_body->process_list(level).m;j++)
               Apply_Prestabilization_To_Joint(dt,time,*articulated_rigid_body,articulated_rigid_body->process_list(level)(j),epsilon_scale);
        articulated_rigid_body->Restore_Velocities_And_Momenta();}
    if(prune_stacks_from_contact) Apply_Stacking_Contact();
}
//#####################################################################
// Function Apply_Prestabilization_To_Joint
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Apply_Prestabilization_To_Joint(const T dt,const T time,ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body,const JOINT_ID joint_id,const T epsilon_scale)
{
    int parent_id=articulated_rigid_body.Parent_Id(joint_id),child_id=articulated_rigid_body.Child_Id(joint_id);
    RIGID_BODY<TV> &parent=rigid_body_collection.Rigid_Body(parent_id),&child=rigid_body_collection.Rigid_Body(child_id);
    // revert to the saved positions & save the proposed positions in rigid_frame_save - restore rigid_frame_save below
    collision_callbacks.Exchange_Frame(parent_id);
    collision_callbacks.Exchange_Frame(child_id);
    parent.Update_Angular_Velocity();child.Update_Angular_Velocity(); // needed for relative velocity
    articulated_rigid_body.Apply_Prestabilization_To_Joint(joint_id,dt,articulated_rigid_body.use_epsilon_scale?epsilon_scale:(T)1);
    collision_callbacks.Save_Position(parent_id);collision_callbacks.Save_Position(child_id); // fix saved values & re-evolve bodies
    Euler_Step_Position(parent_id,dt,time);Euler_Step_Position(child_id,dt,time);
}
//#####################################################################
// Function Either_Body_Collides_With_The_Other
//#####################################################################
template<class TV> bool RIGID_BODY_COLLISIONS<TV>::
Either_Body_Collides_With_The_Other(const int rigid_body_id_1,const int rigid_body_id_2) const
{
    return rigid_body_cluster_bindings.Valid_Cluster_Collision(rigid_body_id_1,rigid_body_id_2)
       && (!collision_manager || collision_manager->Either_Body_Collides_With_The_Other(rigid_body_id_1,rigid_body_id_2));
}
//#####################################################################
// Function Body_Collides_With_The_Other
//#####################################################################
template<class TV> bool RIGID_BODY_COLLISIONS<TV>::
Body_Collides_With_The_Other(const int rigid_body_id_1,const int rigid_body_id_2) const
{
    return rigid_body_cluster_bindings.Valid_Cluster_Collision(rigid_body_id_1,rigid_body_id_2)
       && (!collision_manager || collision_manager->Body_Collides_With_The_Other(rigid_body_id_1,rigid_body_id_2));
}
//#####################################################################
// Function Euler_Step_Position
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Euler_Step_Position(const int id,const T dt,const T time)
{
    collision_callbacks.Euler_Step_Position(id,dt,time);
    rigid_body_collection.Rigid_Body(id).Update_Bounding_Box();skip_collision_check.Set_Last_Moved(id);
}
//#####################################################################
// Function Set_Level_Temporarily_Static
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Set_Level_Temporarily_Static(const int level)
{
    const ARRAY<int>& nodes=contact_graph.directed_graph.Nodes_In_Level(level);
    for(int i=1;i<=nodes.m;i++) if(rigid_body_collection.Is_Active(nodes(i))) rigid_body_collection.Rigid_Body(nodes(i)).is_temporarily_static=true;
}
//#####################################################################
// Function Clear_Temporarily_Static
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Clear_Temporarily_Static()
{
    for(COLLISION_GEOMETRY_ID i(1);i<=rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies.Size();i++){
        RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_geometry=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies(i));
        if(rigid_collision_geometry){RIGID_BODY<TV>* rigid_body=dynamic_cast<RIGID_BODY<TV>*>(&rigid_collision_geometry->rigid_geometry);
            if(rigid_body) if(rigid_body_collection.Is_Active(rigid_body->particle_index)) rigid_body_collection.Rigid_Body(rigid_body->particle_index).is_temporarily_static=false;}}
}
//#####################################################################
// Function Clean_Up_Fractured_Items_From_Lists
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Clean_Up_Fractured_Items_From_Lists(ARRAY<VECTOR<int,2> >& pairs,const int current_pair,const bool called_from_contact)
{
    if(added_bodies(1).m+added_bodies(2).m>0){
        // Go through the remaining pairs, and remove any pairs that contain either of the removed fractured bodies
        for(int k=pairs.m;k>=current_pair+1;k--){
            if(pairs(k)(1)==fractured_bodies(1) || pairs(k)(2)==fractured_bodies(1) || 
                pairs(k)(1)==fractured_bodies(2) || pairs(k)(2)==fractured_bodies(2))
                pairs.Remove_Index_Lazy(k);}
        // Go through the pairs processed by collisions, and add the cartesian product
        ARRAY<VECTOR<int,2> > keys;pairs_processed_by_collisions.Get_Keys(keys);
        for(int k=1;k<=keys.m;k++){
            ARRAY<int> first_ids,second_ids;
            if(keys(k)(1)==fractured_bodies(1)) first_ids=added_bodies(1);
            else if(keys(k)(1)==fractured_bodies(2)) first_ids=added_bodies(2);
            if(keys(k)(2)==fractured_bodies(1)) second_ids=added_bodies(1);
            else if(keys(k)(2)==fractured_bodies(2)) second_ids=added_bodies(2);
            if(!first_ids.m && !second_ids.m) continue;
            pairs_processed_by_collisions.Delete(keys(k));
            if(!first_ids.m && second_ids.m) 
                for(int body=1;body<=second_ids.m;body++) pairs_processed_by_collisions.Set(VECTOR<int,2>(keys(k)(1),second_ids(body)).Sorted());
            else if(first_ids.m && !second_ids.m)
                for(int body=1;body<=first_ids.m;body++) pairs_processed_by_collisions.Set(VECTOR<int,2>(keys(k)(2),first_ids(body)).Sorted());
            else
                for(int body1=1;body1<=first_ids.m;body1++)
                    for(int body2=1;body2<=second_ids.m;body2++)
                        pairs_processed_by_collisions.Set(VECTOR<int,2>(first_ids(body1),second_ids(body2)).Sorted());}
        if(!called_from_contact){ // TODO: add new contact pairs to the list as well
            // Go through added bodies, and collide them against all other bodies and add them to the pairs list
            spatial_partition->Set_Collision_Body_Thickness(0);
            for(int k=1;k<=added_bodies.m;k++) for(int body=1;body<=added_bodies(k).m;body++)
                Get_Bounding_Box_Collision_Pairs_Of_Body(pairs,added_bodies(k)(body),
                    parameters.collision_bounding_box_thickness);}}
}
//#####################################################################
// Function Initialize_All_Contact_Projections
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Initialize_All_Contact_Projections(const bool enforce_rigid_rigid_contact_in_cg)
{
    // rigid/rigid
    if(rigid_body_collection.dynamic_rigid_body_particles.m && enforce_rigid_rigid_contact_in_cg){
        for(typename HASHTABLE<TRIPLE<int,int,TV> >::ITERATOR iterator(rigid_body_particle_intersections);iterator.Valid();iterator.Next()){
            const TRIPLE<int,int,TV>& intersection=iterator.Key();
            const RIGID_BODY<TV> &particle_body=rigid_body_collection.Rigid_Body(intersection.x),
                &levelset_body=rigid_body_collection.Rigid_Body(intersection.y);
            if(rigid_body_collection.Is_Active(particle_body.particle_index) && rigid_body_collection.Is_Active(levelset_body.particle_index)) Create_Contact_Joint(particle_body,levelset_body,intersection.z);}}
}
//#####################################################################
// Function Remove_Contact_Joints
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Remove_Contact_Joints()
{
    for(int i=1;i<=contact_joints.m;i++) rigid_body_collection.articulated_rigid_body.joint_mesh.Remove_Articulation(contact_joints(i));
    contact_joints.Remove_All();
}
//#####################################################################
// Function Create_Contact_Joint
//#####################################################################
template<class TV> void RIGID_BODY_COLLISIONS<TV>::
Create_Contact_Joint(const RIGID_BODY<TV>& parent,const RIGID_BODY<TV>& child,const TV& location)
{
    NORMAL_JOINT<TV>* joint=new NORMAL_JOINT<TV>;
    rigid_body_collection.articulated_rigid_body.joint_mesh.Add_Articulation(parent.particle_index,child.particle_index,joint);
    T phi;const TV normal=child.Implicit_Geometry_Extended_Normal(location,phi);
    FRAME<TV> J(location,ROTATION<TV>::From_Rotated_Vector(TV::Axis_Vector(1),normal));
    joint->Set_Parent_To_Joint_Frame(J.Inverse_Times(parent.Frame())); // TODO: make setting up a joint from a frame part of the joint infrastructure
    joint->Set_Child_To_Joint_Frame(J.Inverse_Times(child.Frame()));

    contact_joints.Append(joint->id_number);
}
//#####################################################################
template class RIGID_BODY_COLLISIONS<VECTOR<float,1> >;
template class RIGID_BODY_COLLISIONS<VECTOR<float,2> >;
template class RIGID_BODY_COLLISIONS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_BODY_COLLISIONS<VECTOR<double,1> >;
template class RIGID_BODY_COLLISIONS<VECTOR<double,2> >;
template class RIGID_BODY_COLLISIONS<VECTOR<double,3> >;
#endif
}

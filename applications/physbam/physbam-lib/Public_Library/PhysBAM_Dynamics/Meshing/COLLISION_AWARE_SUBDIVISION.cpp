//#####################################################################
// Copyright 2003-2007, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGLE_SUBDIVISION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COLLISION_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/STRUCTURE_INTERACTION_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Dynamics/Meshing/COLLISION_AWARE_SUBDIVISION.h>
using namespace PhysBAM;
//#####################################################################
// Subdivision Small Number
//#####################################################################
template<class T>
inline T Subdivision_Small_Number()
{
    return (T)1e-5;
}
template<>
inline double Subdivision_Small_Number()
{
    return 1e-12;
}
//#####################################################################
// Function Push_Surface_Outside_Of_Collision_Bodies
//#####################################################################
// assumes the cloth is free of self-intersections, but may penetrate collision bodies in the scene
template<class T> bool COLLISION_AWARE_SUBDIVISION<T>::
Push_Surface_Outside_Of_Collision_Bodies(const T push_distance,const int max_number_of_attempts)
{
    int attempts=0,interactions=0;T depth;
    TRIANGLE_REPULSIONS<TV> triangle_repulsions(geometry);
    TRIANGLE_COLLISIONS<TV> triangle_collisions(geometry,triangle_repulsions.repulsion_thickness);
    geometry.Build_Collision_Geometry();triangle_collisions.Set_Collision_Thickness(collision_tolerance);
    geometry.Set_Small_Number((T)1e-12);
    while(!attempts || (interactions && attempts<max_number_of_attempts)){
        attempts++;interactions=0;
        for(int s=1;s<=geometry.structure_geometries.m;s++){
            STRUCTURE_INTERACTION_GEOMETRY<TV>& structure=*geometry.structure_geometries(s);
            GEOMETRY_PARTICLES<TV>& particles=structure.triangulated_surface->particles;
            ARRAY_VIEW<TV> X(particles.X),V(particles.V);
            V=X;
            for(int p=1;p<=particles.array_collection->Size();p++) for(COLLISION_GEOMETRY_ID r(1);r<=collision_body_list.bodies.m;r++)
                if(collision_body_list.Is_Active(r) && collision_body_list.bodies(r)->Implicit_Geometry_Lazy_Inside_And_Value(X(p),depth,push_distance)){
                    depth=push_distance-depth+collision_tolerance;interactions++;
                    COLLISION_BODY_HELPER<TV>::Adjust_Point_For_Collision(*collision_body_list.bodies(r),X(p),depth);}}
        if(interactions && geometry.Check_For_Intersection(false,collision_tolerance)){
            for(int s=1;s<=geometry.structure_geometries.m;s++){
                STRUCTURE_INTERACTION_GEOMETRY<TV>& structure=*geometry.structure_geometries(s);
                GEOMETRY_PARTICLES<TV>& particles=structure.triangulated_surface->particles;
                ARRAY<TV> new_V(particles.X-particles.V);
                particles.X=particles.V;
                particles.V=new_V;}
            geometry.Save_Self_Collision_Free_State();
            triangle_collisions.Stop_Nodes_Before_Self_Collision(1);}}
    return interactions==0;
}
//#####################################################################
// Function Subdivide
//#####################################################################
template<class T> void COLLISION_AWARE_SUBDIVISION<T>::
Subdivide(const int number_of_subdivisions,const bool push_outside_collision_bodies,const T push_distance,const int push_attempts,const bool verbose)
{
    if(push_outside_collision_bodies){
        bool sucessful=Push_Surface_Outside_Of_Collision_Bodies(push_distance,push_attempts);
        if(verbose && !sucessful) {std::stringstream ss;ss<<"unable to push surface ouside of collision bodies in "<<push_attempts<<" attempts"<<std::endl;LOG::filecout(ss.str());}}
    for(int k=1;k<=number_of_subdivisions;k++){
        if(verbose) {std::stringstream ss;ss<<"Subdivision "<<k<<"...";LOG::filecout(ss.str());}
        Subdivide();
        if(verbose) {std::stringstream ss;ss<<"done."<<std::endl;LOG::filecout(ss.str());}}
}
//#####################################################################
// Function Subdivide
//#####################################################################
template<class T> void COLLISION_AWARE_SUBDIVISION<T>::
Subdivide()
{
    GEOMETRY_PARTICLES<TV>& full_particles=geometry.deformable_body_collection.particles;
    ARRAY<TV> goal_X;ARRAY_VIEW<TV> X(full_particles.X),V(full_particles.V);
    for(int s=1;s<=geometry.structures.m;s++){
        TRIANGULATED_SURFACE<T>& triangulated_surface=dynamic_cast<TRIANGULATED_SURFACE<T>&>(*geometry.structures(s)); // MUST BE triangulated surface
        triangulated_surface.mesh.Set_Number_Nodes(full_particles.array_collection->Size());
        TRIANGLE_SUBDIVISION triangle_subdivision(triangulated_surface.mesh);TRIANGLE_MESH refined_mesh;triangle_subdivision.Refine_Mesh(refined_mesh);
        int number_new_particles=refined_mesh.number_nodes-full_particles.array_collection->Size();full_particles.array_collection->Add_Elements(number_new_particles);
        triangle_subdivision.Apply_Linear_Subdivision(X,X);goal_X.Resize(full_particles.array_collection->Size());triangle_subdivision.Apply_Loop_Subdivision(X,goal_X);
        triangulated_surface.mesh.Initialize_Mesh(refined_mesh);triangulated_surface.Clean_Memory();} // auxiliary structures are wrong and need to be deleted

    TRIANGLE_REPULSIONS<TV> triangle_repulsions(geometry); // TODO: maybe we could move repulsion thickness into the geometry class, then this line goes away
    TRIANGLE_COLLISIONS<TV> triangle_collisions(geometry,triangle_repulsions.repulsion_thickness);
    geometry.Build_Collision_Geometry();triangle_collisions.Set_Collision_Thickness(collision_tolerance);
    geometry.Set_Small_Number(Subdivision_Small_Number<T>());
    if(geometry.Check_For_Intersection(false,collision_tolerance)){
        V=goal_X-X;
        geometry.deformable_body_collection.particles.Euler_Step_Position(1);
        geometry.Save_Self_Collision_Free_State();
        triangle_collisions.Stop_Nodes_Before_Self_Collision(1);}
    else X=goal_X;
}
//#####################################################################
template class COLLISION_AWARE_SUBDIVISION<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COLLISION_AWARE_SUBDIVISION<double>;
#endif

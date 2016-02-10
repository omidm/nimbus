//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace SOLVE_CONTACT
//##################################################################### 
#ifndef __DISCRETE_CONTACT__
#define __DISCRETE_CONTACT__

#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/CONTACT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/PROJECTED_GAUSS_SEIDEL.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_STRUCTURE_INTERACTION_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_TRIANGLE_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
namespace PhysBAM{

namespace DISCRETE_CONTACT
{
//#####################################################################
// Function Compute_Contacts
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
template<class T>
void Compute_Contacts(RIGID_TRIANGLE_COLLISIONS<VECTOR<T,1> >* triangle_collisions,RIGID_BODY<VECTOR<T,1> >& body1,RIGID_BODY<VECTOR<T,1> >& body2,ARRAY<CONTACT<VECTOR<T,1> > >& contacts,const T contact_thickness,const T collision_thickness,const T dt)
{}
template<class T>
void Compute_Contacts(RIGID_TRIANGLE_COLLISIONS<VECTOR<T,2> >* triangle_collisions,RIGID_BODY<VECTOR<T,2> >& body1,RIGID_BODY<VECTOR<T,2> >& body2,ARRAY<CONTACT<VECTOR<T,2> > >& contacts,const T contact_thickness,const T collision_thickness,const T dt)
{}
template<class T>
void Compute_Contacts(RIGID_TRIANGLE_COLLISIONS<VECTOR<T,3> >* triangle_collisions,RIGID_BODY<VECTOR<T,3> >& body1,RIGID_BODY<VECTOR<T,3> >& body2,ARRAY<CONTACT<VECTOR<T,3> > >& contacts,const T contact_thickness,const T collision_thickness,const T dt)
{
    typedef VECTOR<T,3> TV;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX T_FACE;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX_FACE T_EDGE;

    triangle_collisions->geometry.interacting_structure_pairs.Remove_All();triangle_collisions->geometry.interacting_structure_pairs.Append(VECTOR<int,2>(body1.particle_index,body2.particle_index));
    FRAME<TV> frame1,frame2;frame1=body1.Frame();frame2=body2.Frame();
    triangle_collisions->geometry.structure_geometries(body1.particle_index)->Save_Current_State(frame1.t,frame1.r,body1.Twist().linear);
    triangle_collisions->geometry.structure_geometries(body2.particle_index)->Save_Current_State(frame2.t,frame2.r,body2.Twist().linear);
    triangle_collisions->geometry.structure_geometries(body1.particle_index)->Save_Self_Collision_Free_State(frame1.t,frame1.r,body1.Twist().linear);
    triangle_collisions->geometry.structure_geometries(body2.particle_index)->Save_Self_Collision_Free_State(frame2.t,frame2.r,body2.Twist().linear);
    triangle_collisions->Compute_Pairs(contact_thickness,frame1,frame2);
    VECTOR<HASHTABLE<int,TRIPLE<T,TV,TV> >,2> vertices;
    for(int pair_index=1;pair_index<=triangle_collisions->point_face_pairs_internal.m;pair_index++){RIGID_BODY<TV>* point_body=&body1,*face_body=&body2;int hash_id=1;
        if(pair_index>triangle_collisions->swap_index){point_body=&body2;face_body=&body1;hash_id=2;}
        RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& point_str=*triangle_collisions->geometry.structure_geometries(point_body->particle_index);
        RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& face_str=*triangle_collisions->geometry.structure_geometries(face_body->particle_index);
        const VECTOR<int,TV::dimension+1>& nodes=triangle_collisions->point_face_pairs_internal(pair_index);
        T_FACE face(face_str.X.Subset(nodes.Remove_Index(1)));
        TV weights;TV closest_point=face.Closest_Point(point_str.X(nodes[1]),weights);        
        TV normal=closest_point-point_str.X(nodes[1]);
        T distance=normal.Magnitude();
        if(distance>0){normal/=distance;
            TRIPLE<T,TV,TV> cur_data(distance,closest_point,normal);
            TRIPLE<T,TV,TV>& min_data=vertices(hash_id).Get_Or_Insert(nodes[1],TRIPLE<T,TV,TV>(cur_data));
            if(distance<min_data.x){min_data.x=distance;min_data.y=closest_point;min_data.z=normal;}}}
    for(int i=1;i<=2;i++){RIGID_BODY<TV>* point_body=&body1,*face_body=&body2;
        if(i==2){point_body=&body2;face_body=&body1;}
        for(typename HASHTABLE<int,TRIPLE<T,TV,TV> >::ITERATOR iterator(vertices(i));iterator.Valid();iterator.Next())
            if(iterator.Data().x>0) contacts.Append(CONTACT<TV>(*point_body,*face_body,iterator.Data().y,iterator.Data().z,min((T)0,iterator.Data().x-collision_thickness),dt));}
    for(int pair_index=1;pair_index<=triangle_collisions->edge_edge_pairs_internal.m;pair_index++){
        RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_geom1=*triangle_collisions->geometry.structure_geometries(body1.particle_index);
        RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_geom2=*triangle_collisions->geometry.structure_geometries(body2.particle_index);
        const VECTOR<int,2*TV::dimension-2>& nodes=triangle_collisions->edge_edge_pairs_internal(pair_index);
        T_EDGE edge1,edge2;Create_Edges(structure_geom1.X,structure_geom2.X,nodes,edge1,edge2);
        VECTOR<T,2> weights;bool boundary;TV normal=edge2.Shortest_Vector_Between_Segments(edge1,weights,&boundary);
        T distance=normal.Magnitude();
        if(!boundary && distance<contact_thickness && distance>0){normal/=distance;
            TV closest_point=weights(1)*edge2.x2+((T)1.0-weights(1))*edge2.x1;
            contacts.Append(CONTACT<TV>(body1,body2,closest_point,normal,min((T)0,distance-collision_thickness),dt));}}
}
//#####################################################################
// Function Compute_Contacts
//#####################################################################
template<class TV>
void Compute_Contacts(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,RIGID_TRIANGLE_COLLISIONS<TV>* triangle_collisions,ARRAY<VECTOR<int,2> >& pairs,ARRAY<CONTACT<TV> >& contacts,const typename TV::SCALAR contact_thickness,const typename TV::SCALAR collision_thickness,const typename TV::SCALAR dt)
{
    typedef typename TV::SCALAR T;
    LOG::SCOPE scope("DISCRETE_CONTACT::Compute_Contacts");

    for(int i=1;i<=pairs.m;i++){
        int id1=pairs(i)(1),id2=pairs(i)(2);
        RIGID_BODY<TV>& body1=rigid_body_collection.Rigid_Body(id1);
        RIGID_BODY<TV>& body2=rigid_body_collection.Rigid_Body(id2);
        if(!body1.Has_Infinite_Inertia() || !body2.Has_Infinite_Inertia())
            if(body1.Bounding_Boxes_Intersect(body2,contact_thickness))
                Compute_Contacts(triangle_collisions,body1,body2,contacts,contact_thickness,collision_thickness,dt);}
}
//#####################################################################
// Function Compute_Contacts
//#####################################################################
template<class TV>
void Compute_Contacts(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,RIGID_TRIANGLE_COLLISIONS<TV>* triangle_collisions,ARRAY<CONTACT<TV> >& contacts,const typename TV::SCALAR contact_thickness,const typename TV::SCALAR collision_thickness,const typename TV::SCALAR dt)
{
    typedef typename TV::SCALAR T;
    LOG::SCOPE scope("DISCRETE_CONTACT::Compute_Contacts");

    COLLISION_GEOMETRY_SPATIAL_PARTITION<COLLISION_GEOMETRY<TV>,ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID> spatial_partition(rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies,contact_thickness);
    OPERATION_HASH<COLLISION_GEOMETRY_ID> get_potential_collisions_already_added;
    spatial_partition.Compute_Voxel_Size(SPATIAL_PARTITION_MAX_BOX_SIZE,Value(rigid_body_collection.rigid_geometry_collection.collision_body_list->bodies.Size()),(T)1);
    spatial_partition.Reinitialize();
    for(int id1=1;id1<=rigid_body_collection.rigid_body_particle.array_collection->Size();id1++) if(rigid_body_collection.Is_Active(id1)){
        ARRAY<COLLISION_GEOMETRY_ID> candidates;
        spatial_partition.Get_Potential_Collisions(rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(id1),candidates,get_potential_collisions_already_added);
        for(int i=1;i<=candidates.Size();i++){
            int id2=rigid_body_collection.rigid_geometry_collection.collision_body_list->collision_geometry_id_to_geometry_id.Get(candidates(i));
            RIGID_BODY<TV>& body1=rigid_body_collection.Rigid_Body(id1);
            RIGID_BODY<TV>& body2=rigid_body_collection.Rigid_Body(id2);
            if(!body1.Has_Infinite_Inertia() || !body2.Has_Infinite_Inertia())
                if(body1.Bounding_Boxes_Intersect(body2,contact_thickness))
                    Compute_Contacts(triangle_collisions,body1,body2,contacts,contact_thickness,collision_thickness,dt);}}
}

}
}
#endif

//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_PARTITION.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/PARTICLE_LEVELSET_VISITOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_PARTICLE_INTERSECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/PARTICLES_IN_IMPLICIT_OBJECT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>

using namespace PhysBAM;
//#####################################################################
// Function Destructor
//#####################################################################
template<class TV> RIGID_BODY_INTERSECTIONS<TV>::
~RIGID_BODY_INTERSECTIONS()
{
    typedef typename HASHTABLE<const GEOMETRY_PARTICLES<TV>*,const PARTICLE_HIERARCHY<TV>*>::ITERATOR T_HIERARCHY_ITERATOR;
    for(T_HIERARCHY_ITERATOR i(particle_hierarchies);i.Valid();i.Next()) delete i.Data();
}
//#####################################################################
// Function Intersection_Check
//#####################################################################
template<class TV> bool RIGID_BODY_INTERSECTIONS<TV>::
Intersection_Check(const int id_1,const int id_2,int& particle_body,int& levelset_body,const T thickness) const
{
    if(!Bounding_Boxes_Intersect(id_1,id_2,thickness)) return false;
    else return Find_Any_Intersection(id_1,id_2,particle_body,levelset_body);
}
//#####################################################################
// Function Bounding_Boxes_Intersect
//#####################################################################
template<class TV> bool RIGID_BODY_INTERSECTIONS<TV>::
Bounding_Boxes_Intersect(const int id_1,const int id_2,const T thickness) const
{
    RIGID_BODY<TV> &body_1=rigid_body_collection.Rigid_Body(id_1),&body_2=rigid_body_collection.Rigid_Body(id_2);
    if(!body_1.Is_Simulated() && !body_2.Is_Simulated()) return false; // don't check when neither object is dynamic
    return body_1.Bounding_Boxes_Intersect(body_2,thickness);
}
//#####################################################################
// Function Find_Any_Intersection
//#####################################################################
template<class TV> bool RIGID_BODY_INTERSECTIONS<TV>::
Find_Any_Intersection(const int id_1,const int id_2,int& particle_body,int& levelset_body) const
{
    RIGID_BODY<TV> &body1=rigid_body_collection.Rigid_Body(id_1),&body2=rigid_body_collection.Rigid_Body(id_2);
    Initialize_Transformation_From_Body1_To_Body2_Coordinates(body1,body2);
    ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> > intersection_list;
    if(!use_triangle_hierarchy){
        if(body1.simplicial_object && body2.implicit_object){
            Particles_In_Levelset(id_1,id_2,intersection_list,0,true);
            if(intersection_list.m){particle_body=id_1;levelset_body=id_2;return true;}}
        if(body1.implicit_object && body2.simplicial_object){
            Flip_Transformation();Particles_In_Levelset(id_2,id_1,intersection_list,0,true);
            if(intersection_list.m){particle_body=id_2;levelset_body=id_1;return true;}}}
    else{ // using triangle hierarchy
        if(!use_edge_intersection){
            ARRAY<int> triangle_list;triangle_list.Preallocate(50);
            if(body1.simplicial_object && body2.implicit_object){
                PARTICLES_IN_IMPLICIT_OBJECT::Get_Interfering_Simplices<TV>(body1,body2,triangle_list,rotation,translation,use_triangle_hierarchy_center_phi_test);
                PARTICLES_IN_IMPLICIT_OBJECT::Intersections_Using_Hierarchy<TV>(body1,body2,triangle_list,intersection_list,0,true,rotation,translation);
                if(intersection_list.m){particle_body=id_1;levelset_body=id_2;return true;}}
            if(body1.implicit_object && body2.simplicial_object){
                Flip_Transformation();PARTICLES_IN_IMPLICIT_OBJECT::Get_Interfering_Simplices<TV>(body2,body1,triangle_list,rotation,translation,use_triangle_hierarchy_center_phi_test);
                PARTICLES_IN_IMPLICIT_OBJECT::Intersections_Using_Hierarchy<TV>(body2,body1,triangle_list,intersection_list,0,true,rotation,translation);
                if(intersection_list.m){particle_body=id_2;levelset_body=id_1;return true;}}}
        else{ // use edge intersection
            ARRAY<int> triangle_list1,triangle_list2;triangle_list1.Preallocate(50);triangle_list2.Preallocate(50);
            if(body1.simplicial_object) PARTICLES_IN_IMPLICIT_OBJECT::Get_Interfering_Simplices<TV>(body1,body2,triangle_list1,rotation,translation,use_triangle_hierarchy_center_phi_test);
            Flip_Transformation();
            if(body2.simplicial_object){PARTICLES_IN_IMPLICIT_OBJECT::Get_Interfering_Simplices<TV>(body2,body1,triangle_list2,rotation,translation,use_triangle_hierarchy_center_phi_test);
                if(body1.implicit_object){PARTICLES_IN_IMPLICIT_OBJECT::Intersections_Using_Hierarchy_And_Edges<TV>(body2,body1,triangle_list2,triangle_list1,intersection_list,0,true,rotation,translation);
                    if(intersection_list.m){particle_body=id_2;levelset_body=id_1;return true;}}}
            if(body1.simplicial_object && body2.implicit_object){
                Flip_Transformation();
                PARTICLES_IN_IMPLICIT_OBJECT::Intersections_Using_Hierarchy_And_Edges<TV>(body1,body2,triangle_list1,triangle_list2,intersection_list,0,true,rotation,translation);
                if(intersection_list.m){particle_body=id_1;levelset_body=id_2;return true;}}}}
    return false;
}
//#####################################################################
// Function Append_All_Intersections
//#####################################################################
template<class TV> void RIGID_BODY_INTERSECTIONS<TV>::
Append_All_Intersections(const int id_1,const int id_2,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const T contour_value) const
{
    RIGID_BODY<TV> &body1=rigid_body_collection.Rigid_Body(id_1),&body2=rigid_body_collection.Rigid_Body(id_2);
    PARTICLES_IN_IMPLICIT_OBJECT::Append_All_Intersections(body1,body2,particle_intersections,contour_value,use_triangle_hierarchy,use_edge_intersection,use_triangle_hierarchy_center_phi_test);
}
//#####################################################################
// Function Particle_Partition
//#####################################################################
template<class TV> const PARTICLE_PARTITION<TV>& RIGID_BODY_INTERSECTIONS<TV>::
Particle_Partition(const RIGID_BODY<TV>& rigid_body) const
{
    // TODO: regenerate if particle_partition_size changes
    PHYSBAM_ASSERT(rigid_body.simplicial_object);
    if(!rigid_body.simplicial_object->particle_partition)
        rigid_body.simplicial_object->Initialize_Particle_Partition(particle_partition_size);
    return *rigid_body.simplicial_object->particle_partition;
}
//#####################################################################
// Function Particle_Hierarchy
//#####################################################################
template<class TV> const PARTICLE_HIERARCHY<TV>& RIGID_BODY_INTERSECTIONS<TV>::
Particle_Hierarchy(const RIGID_BODY<TV>& rigid_body) const
{
    PHYSBAM_ASSERT(rigid_body.simplicial_object);
    const GEOMETRY_PARTICLES<TV>& particles=rigid_body.simplicial_object->particles;
    try{
        return *particle_hierarchies.Get(&particles);}
    catch(KEY_ERROR&){
        PARTICLE_HIERARCHY<TV>* particle_hierarchy(new PARTICLE_HIERARCHY<TV>(particles.X));
        particle_hierarchy->Update_Box_Radii();
        particle_hierarchies.Insert(&particles,particle_hierarchy);
        return *particle_hierarchy;}
}
//#####################################################################
// Function Simplex_Hierarchy
//#####################################################################
template<class TV> const typename RIGID_BODY_INTERSECTIONS<TV>::T_SIMPLEX_HIERARCHY& RIGID_BODY_INTERSECTIONS<TV>::
Simplex_Hierarchy(const RIGID_BODY<TV>& rigid_body) const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Particles_In_Levelset
//#####################################################################
template<class TV> void RIGID_BODY_INTERSECTIONS<TV>::
Particles_In_Levelset(const int particle_body_id,const int levelset_body_id,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,
    const T contour_value,const bool exit_early) const
{
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;
    RIGID_BODY<TV> &body1=rigid_body_collection.Rigid_Body(particle_body_id),&body2=rigid_body_collection.Rigid_Body(levelset_body_id);
    if(use_particle_hierarchy){
        PARTICLES_IN_IMPLICIT_OBJECT::Particles_In_Implicit_Object_Hierarchy<TV>(body1,body2,particle_intersections,contour_value,particle_hierarchies);}
    else if(!use_particle_partition){
        PARTICLES_IN_IMPLICIT_OBJECT::Particles_In_Implicit_Object<TV>(body1,body2,particle_intersections,contour_value,exit_early);}
    else{ // use particle partition
        PARTICLES_IN_IMPLICIT_OBJECT::Particles_In_Implicit_Object_Partition<TV>(body1,body2,particle_intersections,contour_value,use_particle_partition_center_phi_test,particle_partition_size,exit_early);
    }
}
//#####################################################################
// Function Oriented_Box2_In_Body1_Coordinates
//#####################################################################
template<class TV> typename BASIC_GEOMETRY_POLICY<TV>::ORIENTED_BOX RIGID_BODY_INTERSECTIONS<TV>::
Oriented_Box2_In_Body1_Coordinates(const RIGID_BODY<TV>& body1,const RIGID_BODY<TV>& body2) const
{
    return T_ORIENTED_BOX(body2.Object_Space_Bounding_Box(),body1.Frame().Inverse_Times(body2.Frame()));
}
//#####################################################################
// Function Initialize_Transformation_From_Body1_To_Body2_Coordinates
//#####################################################################
template<class TV> void RIGID_BODY_INTERSECTIONS<TV>::
Initialize_Transformation_From_Body1_To_Body2_Coordinates(const RIGID_BODY<TV>& body1,const RIGID_BODY<TV>& body2) const
{
    FRAME<TV> frame=body2.Frame().Inverse_Times(body1.Frame());
    rotation=frame.r.Rotation_Matrix();translation=frame.t;
    rotation_reverse=rotation.Transposed();translation_reverse=-(rotation_reverse*translation); // reverse is from body2 to body1 coordinates
}
//#####################################################################
template class RIGID_BODY_INTERSECTIONS<VECTOR<float,1> >;
template class RIGID_BODY_INTERSECTIONS<VECTOR<float,2> >;
template class RIGID_BODY_INTERSECTIONS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_BODY_INTERSECTIONS<VECTOR<double,1> >;
template class RIGID_BODY_INTERSECTIONS<VECTOR<double,2> >;
template class RIGID_BODY_INTERSECTIONS<VECTOR<double,3> >;
#endif

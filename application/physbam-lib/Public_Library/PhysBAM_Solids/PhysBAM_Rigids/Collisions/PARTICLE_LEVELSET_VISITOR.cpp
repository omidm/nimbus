//#####################################################################
// Copyright 2002-2008, Ronald Fedkiw, Eran Guendelman, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY_DEFINITION.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/PARTICLE_LEVELSET_VISITOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_PARTICLE_INTERSECTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLE_LEVELSET_VISITOR<TV>::
PARTICLE_LEVELSET_VISITOR(const PARTICLE_HIERARCHY<TV>& particle_hierarchy,const IMPLICIT_OBJECT<TV>& object_space_implicit_object,const MATRIX<T,d>& rotation,
    const VECTOR<T,d>& translation,const T contour_value,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const int particle_body_id,
    const int levelset_body_id)
    :box_hierarchy(particle_hierarchy.box_hierarchy),box_radius(particle_hierarchy.box_radius),X(particle_hierarchy.X),
    particles_in_group(particle_hierarchy.particles_in_group),object_space_implicit_object(object_space_implicit_object),rotation(rotation),translation(translation),
    contour_value(contour_value),particle_intersections(particle_intersections),particle_body_id(particle_body_id),levelset_body_id(levelset_body_id)
{
    PHYSBAM_ASSERT(particle_hierarchy.particles_per_group);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PARTICLE_LEVELSET_VISITOR<TV>::
~PARTICLE_LEVELSET_VISITOR()
{}
//#####################################################################
// Function Store
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_VISITOR<TV>::
Store(const int box) const
{
    TV center=rotation*box_hierarchy(box).Center()+translation; // in a perfect world, CSE would merge this with center from Cull
    ARRAY_VIEW<const int> group(particles_in_group(box));
    if(object_space_implicit_object.box.Inside(center,box_radius(box))){ // if box is entirely inside implicit object box, we can use Lazy_Inside
        for(int i=1;i<=group.Size();i++){int p=group(i);
            if(object_space_implicit_object.Lazy_Inside(rotation*X(p)+translation,contour_value))
                particle_intersections.Append(RIGID_BODY_PARTICLE_INTERSECTION<TV>(X(p),p,particle_body_id,levelset_body_id));}}
    else{
        for(int i=1;i<=group.Size();i++){int p=group(i);
            if(object_space_implicit_object.Lazy_Inside_Extended_Levelset(rotation*X(p)+translation,contour_value))
                particle_intersections.Append(RIGID_BODY_PARTICLE_INTERSECTION<TV>(X(p),p,particle_body_id,levelset_body_id));}}
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template struct PARTICLE_LEVELSET_VISITOR<VECTOR<T,d> >; \
    template void BOX_HIERARCHY<VECTOR<T,d> >::Intersection_List<PARTICLE_LEVELSET_VISITOR<VECTOR<T,d> > >(PARTICLE_LEVELSET_VISITOR<VECTOR<T,d> >&) const;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif

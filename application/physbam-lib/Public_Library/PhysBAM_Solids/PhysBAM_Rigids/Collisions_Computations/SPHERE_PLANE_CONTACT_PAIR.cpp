//#####################################################################
// Copyright 2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace CONTACT_PAIRS
//##################################################################### 
#ifndef __SPHERE_PLANE_CONTACT_PAIR__
#define __SPHERE_PLANE_CONTACT_PAIR__

#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_PARTICLE_INTERSECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_SKIP_COLLISION_CHECK.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGIDS_COLLISION_CALLBACKS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions_Computations/SOLVE_CONTACT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>

namespace PhysBAM{

template<class TV> class RIGID_BODY_COLLECTION;
template<class TV> class RIGIDS_COLLISION_CALLBACKS;

namespace CONTACT_PAIRS
{
template<class TV>
bool Update_Sphere_Plane_Contact_Pair(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,
    const int id_1,const int id_2,IMPLICIT_OBJECT<TV>* object1,IMPLICIT_OBJECT<TV>* object2,const bool correct_contact_energy,const int max_iterations,
    const typename TV::SCALAR epsilon_scale,const typename TV::SCALAR dt,const typename TV::SCALAR time,const bool mpi_one_ghost)
{
    typedef typename TV::SCALAR T;

    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=rigid_body_collisions.rigid_body_collection;

    int i1=id_1;
    int i2=id_2;
    FRAME<TV> transform;
    RIGID_BODY<TV>* body1=&rigid_body_collection.Rigid_Body(i1);
    RIGID_BODY<TV>* body2=&rigid_body_collection.Rigid_Body(i2);
    if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* object_transformed=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(object2)){
        transform=*object_transformed->transform;object2=object_transformed->object_space_implicit_object;}
    ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >* implicit_sphere=dynamic_cast<ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >*>(object2);
    if(!implicit_sphere){
        if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* object_transformed=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(object1)){
            transform=*object_transformed->transform;object1=object_transformed->object_space_implicit_object;}
        exchange(body1,body2);
        exchange(i1,i2);
        exchange(object1,object2);
        implicit_sphere=dynamic_cast<ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >*>(object2);}
    SPHERE<TV>& sphere=implicit_sphere->analytic;

    TV sphere_center=(body2->Frame()*transform).t;
    TV collision_normal=-body1->Rotation().Rotated_Axis(2);
    T separation=TV::Dot_Product(body1->X()-sphere_center,collision_normal);
    if(separation>=sphere.radius){rigid_body_collisions.skip_collision_check.Set_Last_Checked(i1,i2);return false;}
    if(TV::Dot_Product(body1->Twist().linear-body2->Twist().linear,collision_normal)>=0) return false;
    T collision_depth=(T).5*(-sphere.radius+min(-separation,sphere.radius));
    TV collision_location=sphere_center-collision_depth*collision_normal;
    TV collision_relative_velocity=body1->Pointwise_Object_Velocity(collision_location)-body2->Pointwise_Object_Velocity(collision_location);

    collision_callbacks.Swap_States(i1,i2);
    SOLVE_CONTACT::Update_Contact_Pair_Helper<TV>(rigid_body_collisions,collision_callbacks,i1,i2,dt,time,epsilon_scale,collision_location,collision_normal,collision_relative_velocity,
        correct_contact_energy,rigid_body_collisions.rolling_friction,mpi_one_ghost);
    return true;
}

#define INSTANTIATION_HELPER(T,d) \
    template bool Update_Sphere_Plane_Contact_Pair(RIGID_BODY_COLLISIONS<VECTOR<T,d> >& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks,\
        const int id_1,const int id_2,IMPLICIT_OBJECT<VECTOR<T,d> >* object1,IMPLICIT_OBJECT<VECTOR<T,d> >* object2,const bool correct_contact_energy,const int max_iterations, \
        const T epsilon_scale,const T dt,const T time,const bool mpi_one_ghost);

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif
}
}
#endif

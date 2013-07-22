//#####################################################################
// Copyright 2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace CONTACT_PAIRS
//##################################################################### 
#ifndef __SPHERE_SPHERE_CONTACT_PAIR__
#define __SPHERE_SPHERE_CONTACT_PAIR__

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
bool Update_Sphere_Sphere_Contact_Pair(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,
    const int id_1,const int id_2,IMPLICIT_OBJECT<TV>* object1,IMPLICIT_OBJECT<TV>* object2,const bool correct_contact_energy,const int max_iterations,
    const typename TV::SCALAR epsilon_scale,const typename TV::SCALAR dt,const typename TV::SCALAR time,const bool mpi_one_ghost)
{
    typedef typename TV::SCALAR T;

    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=rigid_body_collisions.rigid_body_collection;

    FRAME<TV> transform1,transform2;
    RIGID_BODY<TV>& body1=rigid_body_collection.Rigid_Body(id_1);
    RIGID_BODY<TV>& body2=rigid_body_collection.Rigid_Body(id_2);
    if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* object_transformed=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(object1)){
        transform1=*object_transformed->transform;object1=object_transformed->object_space_implicit_object;}
    if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* object_transformed=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(object2)){
        transform2=*object_transformed->transform;object2=object_transformed->object_space_implicit_object;}
    SPHERE<TV>& sphere1=dynamic_cast<ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >&>(*object1).analytic;
    SPHERE<TV>& sphere2=dynamic_cast<ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >&>(*object2).analytic;

    TV sphere1_center=(body1.Frame()*transform1).t,sphere2_center=(body2.Frame()*transform2).t;
    TV collision_normal=body1.X()-body2.X();collision_normal.Normalize();
    T d=(sphere1_center-sphere2_center).Magnitude(),r1=sphere1.radius,r2=sphere2.radius;
    if(d>r1+r2){rigid_body_collisions.skip_collision_check.Set_Last_Checked(id_1,id_2);return false;}
    if(TV::Dot_Product(collision_normal,body1.Twist().linear-body2.Twist().linear)>=0) return false;
    TV collision_location=sphere1_center+(T).5*(d+min(d,r2)-min(d,r1))*collision_normal;
    TV collision_relative_velocity=body1.Pointwise_Object_Velocity(collision_location)-body2.Pointwise_Object_Velocity(collision_location);

    collision_callbacks.Swap_States(id_1,id_2);
    SOLVE_CONTACT::Update_Contact_Pair_Helper<TV>(rigid_body_collisions,collision_callbacks,id_1,id_2,dt,time,epsilon_scale,collision_location,collision_normal,collision_relative_velocity,
        correct_contact_energy,rigid_body_collisions.rolling_friction,mpi_one_ghost);
    return true;
}

#define INSTANTIATION_HELPER(T,d) \
    template bool Update_Sphere_Sphere_Contact_Pair(RIGID_BODY_COLLISIONS<VECTOR<T,d> >& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks, \
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

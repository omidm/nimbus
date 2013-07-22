//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_STANDARD_TESTS
//#####################################################################
#ifndef __RIGIDS_STANDARD_TESTS__
#define __RIGIDS_STANDARD_TESTS__

#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_ID.h>
namespace PhysBAM{

template<class TV> class EXAMPLE;
template<class TV> class RIGID_BODY;
template<class TV> class RIGID_BODY_COLLECTION;

template<class TV>
class RIGIDS_STANDARD_TESTS
{
    typedef typename TV::SCALAR T;
public:
    std::string& data_directory;
    const STREAM_TYPE& stream_type;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;

    RIGIDS_STANDARD_TESTS(std::string& data_directory_input,const STREAM_TYPE& stream_type_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input);

//#####################################################################
    RIGID_BODY<TV>& Add_Rigid_Body(const std::string& name,const T scaling_factor,const T friction,const bool read_implicit=true,const bool always_read_object=false, const STREAM_TYPE* stream_type_read=0);
    RIGID_BODY<TV>& Add_Ground(const T friction=(T).3,const T height=0,const T coefficient_of_restitution=(T).5,const T scale=(T)1);
    RIGID_BODY<TV>& Add_Analytic_Box(const VECTOR<T,1>& scaling_factor);
    RIGID_BODY<TV>& Add_Analytic_Box(const VECTOR<T,2>& scaling_factor,int segments_per_side=1);
    RIGID_BODY<TV>& Add_Analytic_Box(const VECTOR<T,3>& scaling_factor);
    RIGID_BODY<TV>& Add_Analytic_Torus(const T inner_radius,const T outer_radius,int inner_resolution=8,int outer_resolution=16);
    RIGID_BODY<TV>& Add_Analytic_Cylinder(const T height,const T radius,int resolution_radius=16,int resolution_height=4);
    RIGID_BODY<TV>& Add_Analytic_Sphere(const T radius,const T density,int levels=4);
    RIGID_BODY<TV>& Add_Analytic_Ellipse(const VECTOR<T,2> radii,const T density,int levels=4);
    RIGID_BODY<TV>& Add_Analytic_Ellipsoid(const VECTOR<T,3> radii,const T density,int levels=4);
    void Make_Lathe_Chain(const FRAME<TV>& frame,const T scale=1,const T friction=(T).6,const T cor=(T).3);
    void Set_Joint_Frames(JOINT_ID id,const TV& location);
    JOINT_ID Connect_With_Point_Joint(RIGID_BODY<TV>& parent,RIGID_BODY<TV>& child,const TV& location);
//#####################################################################
};
}
#endif

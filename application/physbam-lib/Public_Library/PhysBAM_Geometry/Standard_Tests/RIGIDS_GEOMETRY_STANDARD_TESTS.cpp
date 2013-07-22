//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_GEOMETRY_STANDARD_TESTS
//#####################################################################
#include <PhysBAM_Tools/Ordinary_Differential_Equations/EXAMPLE.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOUNDED_HORIZONTAL_PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/RING.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TORUS.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Tessellation/CYLINDER_TESSELLATION.h>
#include <PhysBAM_Geometry/Tessellation/IMPLICIT_OBJECT_TESSELLATION.h>
#include <PhysBAM_Geometry/Tessellation/RANGE_TESSELLATION.h>
#include <PhysBAM_Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <PhysBAM_Geometry/Tessellation/TORUS_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Standard_Tests/RIGIDS_GEOMETRY_STANDARD_TESTS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGIDS_GEOMETRY_STANDARD_TESTS<TV>::
RIGIDS_GEOMETRY_STANDARD_TESTS(EXAMPLE<TV>& example_input,RIGID_GEOMETRY_COLLECTION<TV>& rigid_body_collection_input)
    :example(example_input),rigid_geometry_collection(rigid_body_collection_input)
{
}
//#####################################################################
// Function Add_Analytic_Sphere
//#####################################################################
namespace{
template<class T,class TV>
static IMPLICIT_OBJECT<TV>* Add_Analytic_Sphere(const TV&,const T scaling_factor)
{
    return new ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >(SPHERE<TV>(TV(),scaling_factor));
}
//#####################################################################
// Function Add_Analytic_Box
//#####################################################################
template<class T,class TV>
static IMPLICIT_OBJECT<TV>* Add_Analytic_Box(const TV&,const T scaling_factor)
{
    return new ANALYTIC_IMPLICIT_OBJECT<BOX<TV> >(BOX<TV>(TV()-scaling_factor/(typename TV::SCALAR)2.,TV()+scaling_factor/(typename TV::SCALAR)2.));
}
//#####################################################################
// Function Add_Analytic_Cylinder
//#####################################################################
template<class T>
static IMPLICIT_OBJECT<VECTOR<T,2> >* Add_Analytic_Cylinder(const VECTOR<T,2>&,const T scaling_factor)
{
    return new ANALYTIC_IMPLICIT_OBJECT<BOX<VECTOR<T,2> > >(BOX<VECTOR<T,2> >(VECTOR<T,2>(-scaling_factor/2,-scaling_factor/2),VECTOR<T,2>(scaling_factor/2,scaling_factor/2)));
}
template<class T>
static IMPLICIT_OBJECT<VECTOR<T,3> >* Add_Analytic_Cylinder(const VECTOR<T,3>&,const T scaling_factor)
{
    return new ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> >(CYLINDER<T>(VECTOR<T,3>(0,-scaling_factor,0),VECTOR<T,3>(0,scaling_factor,0),(T).1*scaling_factor));
}
//#####################################################################
// Function Add_Analytic_Ring
//#####################################################################
template<class T>
static IMPLICIT_OBJECT<VECTOR<T,2> >* Add_Analytic_Ring(const VECTOR<T,2>&,const T scaling_factor)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
template<class T>
static IMPLICIT_OBJECT<VECTOR<T,3> >* Add_Analytic_Ring(const VECTOR<T,3>&,const T scaling_factor)
{
    T half_height=(T).5*scaling_factor;
    return new ANALYTIC_IMPLICIT_OBJECT<RING<T> >(RING<T>(VECTOR<T,3>(0,-half_height,0),VECTOR<T,3>(0,half_height,0),(T)3*scaling_factor,(T)2*scaling_factor));
}
//#####################################################################
// Function Add_Analytic_Ground
//#####################################################################
template<class TV,class T>
static IMPLICIT_OBJECT<TV>* Add_Analytic_Ground(const TV&,const T scaling_factor)
{
    return new ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<TV> >(BOUNDED_HORIZONTAL_PLANE<TV>(100*scaling_factor));
}
}
//#####################################################################
// Function Add_Rigid_Body
//#####################################################################
template<class TV> RIGID_GEOMETRY<TV>& RIGIDS_GEOMETRY_STANDARD_TESTS<TV>::
Add_Rigid_Body(const std::string& name,const T scaling_factor,const T friction,const bool read_implicit,const bool always_read_object)
{
    std::string rigid_directory=example.data_directory+"/Rigid_Bodies"+(TV::m==3?"":"_2D"),basename=rigid_directory+"/"+name;
    bool read_simplicial=true;
    if(!always_read_object){
        IMPLICIT_OBJECT<TV>* implicit=0;
        if(name=="sphere") implicit=::Add_Analytic_Sphere(TV(),scaling_factor);
        if(name=="skinnycyllink") implicit=::Add_Analytic_Cylinder(TV(),scaling_factor);
        if(name=="Rings_Test/ring_revolve") implicit=Add_Analytic_Ring(TV(),scaling_factor);
        if(name=="ground") implicit=Add_Analytic_Ground(TV(),scaling_factor);
        if(name=="box") implicit=::Add_Analytic_Box(TV(),scaling_factor);
        if(implicit){
            if(TV::m==3){
                STRUCTURE<TV>* structure=choice<TV::m-1>(implicit,TESSELLATION::Generate_Triangles(*implicit));
                if(!rigid_geometry_collection.Register_Analytic_Replacement_Structure(basename+".tri",scaling_factor,structure))
                    delete structure;}
            else
                read_simplicial=false;
            if(!rigid_geometry_collection.Register_Analytic_Replacement_Structure(basename+(TV::m==3?".phi":".phi2d"),scaling_factor,implicit))
                delete implicit;}}
    int id=rigid_geometry_collection.Add_Rigid_Geometry(example.stream_type,basename,scaling_factor,read_simplicial,read_implicit);
    RIGID_GEOMETRY<TV>& rigid_body=rigid_geometry_collection.Rigid_Geometry(id);
    rigid_body.name=name;
    return rigid_body;
}
//#####################################################################
// Function Add_Ground
//#####################################################################
template<class TV> RIGID_GEOMETRY<TV>& RIGIDS_GEOMETRY_STANDARD_TESTS<TV>::
Add_Ground(const T friction,const T height,const T scale)
{
    PHYSBAM_ASSERT(scale>0);
    RIGID_GEOMETRY<TV>& ground=Add_Rigid_Body("ground",scale,friction);
    ground.X().y=height;
    ground.is_static=true;
    ground.rigid_geometry_collection.collision_body_list->Get_Collision_Geometry(ground.particle_index)->add_to_spatial_partition=false;
    ground.Set_Name("ground");
    return ground;
}
//#####################################################################
#define INSTANTIATION_HELPER_1D(T) \
    template RIGIDS_GEOMETRY_STANDARD_TESTS<VECTOR<T,1> >::RIGIDS_GEOMETRY_STANDARD_TESTS(EXAMPLE<VECTOR<T,1> >&,RIGID_GEOMETRY_COLLECTION<VECTOR<T,1> >&);

#define INSTANTIATION_HELPER_ALL(T,d) \
    template RIGIDS_GEOMETRY_STANDARD_TESTS<VECTOR<T,d> >::RIGIDS_GEOMETRY_STANDARD_TESTS(EXAMPLE<VECTOR<T,d> >&,RIGID_GEOMETRY_COLLECTION<VECTOR<T,d> >&); \
    template RIGID_GEOMETRY<VECTOR<T,d> >& RIGIDS_GEOMETRY_STANDARD_TESTS<VECTOR<T,d> >::Add_Rigid_Body(const std::string&,const T,const T,const bool,const bool); \
    template RIGID_GEOMETRY<VECTOR<T,d> >& RIGIDS_GEOMETRY_STANDARD_TESTS<VECTOR<T,d> >::Add_Ground(const T,const T,const T);

#define INSTANTIATION_HELPER(T) \
    INSTANTIATION_HELPER_1D(T) \
    INSTANTIATION_HELPER_ALL(T,2) \
    INSTANTIATION_HELPER_ALL(T,3)

INSTANTIATION_HELPER(float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double);
#endif


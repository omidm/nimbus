//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_STANDARD_TESTS
//#####################################################################
#include <PhysBAM_Geometry/Basic_Geometry/BOUNDED_HORIZONTAL_PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/RING.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/ELLIPSE.h>
#include <PhysBAM_Geometry/Basic_Geometry/ELLIPSOID.h>
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
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Standard_Tests/RIGIDS_STANDARD_TESTS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGIDS_STANDARD_TESTS<TV>::
RIGIDS_STANDARD_TESTS(std::string& data_directory_input,const STREAM_TYPE& stream_type_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
    :data_directory(data_directory_input),stream_type(stream_type_input),rigid_body_collection(rigid_body_collection_input)
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
    return new ANALYTIC_IMPLICIT_OBJECT<BOX<TV> >(BOX<TV>(TV()-scaling_factor/(T)2.,TV()+scaling_factor/(T)2.));
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
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Rigid_Body(const std::string& name,const T scaling_factor,const T friction,const bool read_implicit,const bool always_read_object,const STREAM_TYPE* stream_type_read)
{
    std::string rigid_directory=data_directory+"/Rigid_Bodies"+(TV::m==3?"":"_2D"),basename=rigid_directory+"/"+name;
    bool read_simplicial=true;
    if(!always_read_object){
        IMPLICIT_OBJECT<TV>* implicit=0;
        if(name=="sphere") implicit=::Add_Analytic_Sphere(TV(),scaling_factor);
        if(name=="skinnycyllink") implicit=::Add_Analytic_Cylinder(TV(),scaling_factor);
        if(name=="Rings_Test/ring_revolve") implicit=Add_Analytic_Ring(TV(),scaling_factor);
        if(name=="ground") implicit=Add_Analytic_Ground(TV(),scaling_factor);
        if(name=="box") implicit=::Add_Analytic_Box(TV(),scaling_factor);
        if(name=="circle") implicit=::Add_Analytic_Sphere(TV(),scaling_factor);
        if(implicit){
            if(TV::m==3){
                STRUCTURE<TV>* structure=choice<TV::m-1>(implicit,TESSELLATION::Generate_Triangles(*implicit));
                if(!rigid_body_collection.rigid_geometry_collection.Register_Analytic_Replacement_Structure(basename+".tri",scaling_factor,structure))
                    delete structure;}
            else
                read_simplicial=false;
            if(!rigid_body_collection.rigid_geometry_collection.Register_Analytic_Replacement_Structure(basename+(TV::m==3?".phi":".phi2d"),scaling_factor,implicit))
                delete implicit;}}
    int id;
    if(stream_type_read)
        id=rigid_body_collection.Add_Rigid_Body(*stream_type_read,basename,scaling_factor,read_simplicial,read_implicit);
    else
        id=rigid_body_collection.Add_Rigid_Body(stream_type,basename,scaling_factor,read_simplicial,read_implicit);
    RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(id);
    rigid_body.Set_Coefficient_Of_Friction(friction);
    rigid_body.name=name;
    return rigid_body;
}
//#####################################################################
// Function Add_Ground
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Ground(const T friction,const T height,const T coefficient_of_restitution,const T scale)
{
    PHYSBAM_ASSERT(scale>0);
    RIGID_BODY<TV>& ground=Add_Rigid_Body("ground",scale,friction);
    ground.X().y=height;
    ground.is_static=true;
    ground.rigid_body_collection.rigid_geometry_collection.collision_body_list->Get_Collision_Geometry(ground.particle_index)->add_to_spatial_partition=false;
    ground.Set_Coefficient_Of_Restitution(coefficient_of_restitution);
    ground.Set_Name("ground");
    return ground;
}
//#####################################################################
// Function Add_Analytic_Box
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Analytic_Box(const VECTOR<T,1>& scaling_factor)
{
    RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection,true);
    BOX<TV> box((T)-.5*scaling_factor,(T).5*scaling_factor);
    rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<BOX<TV> >(box));
    GEOMETRY_PARTICLES<TV>& particles=*new GEOMETRY_PARTICLES<TV>;
    POINT_SIMPLICES_1D<T>& simplicial_object=*POINT_SIMPLICES_1D<T>::Create(particles);
    particles.array_collection->Add_Elements(2);
    POINT_SIMPLEX_MESH& segment_mesh=simplicial_object.mesh;segment_mesh.number_nodes=2;segment_mesh.elements.Preallocate(1);
    particles.X(1)=VECTOR<T,1>(box.min_corner.x);particles.X(2)=VECTOR<T,1>(box.max_corner.x);
    simplicial_object.mesh.elements.Append(VECTOR<int,1>(1));simplicial_object.mesh.directions.Append(false);
    simplicial_object.mesh.elements.Append(VECTOR<int,1>(2));simplicial_object.mesh.directions.Append(true);
    rigid_body.Add_Structure(simplicial_object);
    simplicial_object.Update_Point_Simplex_List();
    assert(simplicial_object.point_simplex_list);
    rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
    return rigid_body;   
}
//#####################################################################
// Function Add_Analytic_Box
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Analytic_Box(const VECTOR<T,2>& scaling_factor,int segments_per_side)
{
    RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection,true);
    BOX<TV> box((T)-.5*scaling_factor,(T).5*scaling_factor);
    rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<BOX<TV> >(box));
    GEOMETRY_PARTICLES<TV>& particles=*new GEOMETRY_PARTICLES<TV>;
    SEGMENTED_CURVE_2D<T>& simplicial_object=*SEGMENTED_CURVE_2D<T>::Create(particles);
    particles.array_collection->Add_Elements(4*segments_per_side);
    SEGMENT_MESH& segment_mesh=simplicial_object.mesh;segment_mesh.number_nodes=4*segments_per_side;segment_mesh.elements.Preallocate(4*segments_per_side);
    int last_node=0;VECTOR<T,2> position=VECTOR<T,2>(box.min_corner.x,box.min_corner.y);
    for(int side=1;side<=4;++side) for(int vertex=1;vertex<=segments_per_side;++vertex){
        int current_node=(side-1)*segments_per_side+vertex;particles.X(current_node)=position;
        position(side%2+1)+=(T)(side<3?1:-1)*(side%2?((box.max_corner.y-box.min_corner.y)/(T)segments_per_side):((box.max_corner.x-box.min_corner.x)/(T)segments_per_side));
        if(last_node) simplicial_object.mesh.elements.Append(VECTOR<int,2>(last_node,current_node));
        last_node=current_node;}
    simplicial_object.mesh.elements.Append(VECTOR<int,2>(last_node,1));
    rigid_body.Add_Structure(simplicial_object);
    simplicial_object.Update_Segment_List();
    assert(simplicial_object.segment_list);
    rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
    return rigid_body;
}
//#####################################################################
// Function Add_Analytic_Box
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Analytic_Box(const VECTOR<T,3>& scaling_factor)
{
    RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection,true);
    BOX<TV> box((T)-.5*scaling_factor,(T).5*scaling_factor);
    rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<BOX<TV> >(box));
    rigid_body.Add_Structure(*TESSELLATION::Generate_Triangles(box));
    rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
    return rigid_body;
}
//#####################################################################
// Function Add_Analytic_Torus
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Analytic_Torus(const T inner_radius,const T outer_radius,int inner_resolution,int outer_resolution)
{
    RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection,true);
    TORUS<T> torus(TV(),TV(0,0,1),inner_radius,outer_radius);
    rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<TORUS<T> >(torus));
    rigid_body.Add_Structure(*TESSELLATION::Generate_Triangles(torus,inner_resolution,outer_resolution));
    rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
    return rigid_body;
}
//#####################################################################
// Function Add_Analytic_Cylinder
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Analytic_Cylinder(const T height,const T radius,int resolution_radius,int resolution_height)
{
    RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection,true);
    CYLINDER<T> cylinder(TV(0,0,-height/2),TV(0,0,height/2),radius);
    rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> >(cylinder));
    rigid_body.Add_Structure(*TESSELLATION::Generate_Triangles(cylinder,resolution_height,resolution_radius));
    rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
    return rigid_body;
}
//#####################################################################
// Function Add_Analytic_Sphere
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Analytic_Sphere(const T radius,const T density,int levels)
{
    RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection,true);
    SPHERE<TV> sphere(TV(),radius);
    rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >(sphere));
    rigid_body.Add_Structure(*TESSELLATION::Tessellate_Boundary(sphere,levels));
    rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
    T r2=radius*radius,m2d=r2*(T)pi*density;
    if(TV::m==2){
        rigid_body.Mass()=m2d;
        rigid_body.Inertia_Tensor()(1,1)=rigid_body.Mass()*r2*(T).25;}
    else if(TV::m==3){
        rigid_body.Mass()=m2d*radius*((T)4/3);
        rigid_body.Inertia_Tensor()*=0;
        rigid_body.Inertia_Tensor()+=rigid_body.Mass()*r2*(T).4;}
    return rigid_body;
}
//#####################################################################
// Function Add_Analytic_Ellipse
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Analytic_Ellipse(const VECTOR<T,2> radii,const T density,int levels)
{
    RIGID_BODY<VECTOR<T,2> >& rigid_body=*new RIGID_BODY<VECTOR<T,2> >(rigid_body_collection,true);
    ELLIPSE<T> ellipse(VECTOR<T,2>(),DIAGONAL_MATRIX<T,2>(radii(1),radii(2)));
    //SPHERE<TV> sphere(TV(),radii(1));
    //rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >(sphere));
    rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<ELLIPSE<T> >(ellipse));
    rigid_body.Add_Structure(*TESSELLATION::Tessellate_Boundary(ellipse,levels));
    rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
    //The codes computing mass below is incorrect and needs to be fixed
    T r2=radii(1)*radii(2),m2d=r2*(T)pi*density;
    rigid_body.Mass()=m2d;
    rigid_body.Inertia_Tensor()(1,1)=rigid_body.Mass()*r2*(T).25;
    return rigid_body;
}
//#####################################################################
// Function Add_Analytic_Ellipse
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Analytic_Ellipsoid(const VECTOR<T,3> radii,const T density,int levels)
{
    RIGID_BODY<VECTOR<T,3> >& rigid_body=*new RIGID_BODY<VECTOR<T,3> >(rigid_body_collection,true);
    ELLIPSOID<T> ellipsoid(VECTOR<T,3>(),DIAGONAL_MATRIX<T,3>(radii(1),radii(2),radii(3)));
    //SPHERE<TV> sphere(TV(),radii(1));
    //rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >(sphere));
    rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<ELLIPSOID<T> >(ellipsoid));
    rigid_body.Add_Structure(*TESSELLATION::Generate_Triangles(ellipsoid,levels));
    rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
    //The codes computing mass below is incorrect and needs to be fixed
    T r2=radii(1)*radii(2),m2d=r2*(T)pi*density;
    rigid_body.Mass()=m2d;
    rigid_body.Inertia_Tensor()(1,1)=rigid_body.Mass()*r2*(T).25;
    return rigid_body;
}
//#####################################################################
// Function Make_Lathe_Chain
//#####################################################################
template<class TV> void RIGIDS_STANDARD_TESTS<TV>::
Make_Lathe_Chain(const FRAME<TV>& frame,const T scale,const T friction,const T cor)
{
    PHYSBAM_ASSERT(scale>0);
    ARTICULATED_RIGID_BODY<TV>& arb=rigid_body_collection.articulated_rigid_body;
    RIGID_BODY<TV>* links[7];

    for(int i=1;i<=6;i++){
        RIGID_BODY<TV>& rigid_body=Add_Rigid_Body("ARB/lathe_object",scale,friction);
        rigid_body.Set_Coefficient_Of_Restitution(cor);
        links[i-1]=&rigid_body;
        switch(i){
            case 1:rigid_body.Set_Frame(frame*FRAME<TV>(TV(0,4*sin((T)pi/3),0)*scale,ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>(-(T)pi,TV(1,0,0))));break;
            case 2:rigid_body.Set_Frame(frame*FRAME<TV>(TV(2+2*cos((T)pi/3),4*sin((T)pi/3)*(T).5,0)*scale,ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>(-2*(T)pi/3,TV(1,0,0))));break;
            case 3:rigid_body.Set_Frame(frame*FRAME<TV>(TV(2+2*cos((T)pi/3),-4*sin((T)pi/3)*(T).5,0)*scale,ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>(-(T)pi/3,TV(1,0,0))));break;
            case 4:rigid_body.Set_Frame(frame*FRAME<TV>(TV(0,-4*sin((T)pi/3),0)*scale,ROTATION<TV>((T)pi/2,TV(0,1,0))));break;
            case 5:rigid_body.Set_Frame(frame*FRAME<TV>(TV(-2-2*cos((T)pi/3),-4*sin((T)pi/3)*(T).5,0)*scale,ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>((T)pi/3,TV(1,0,0))));break;
            case 6:rigid_body.Set_Frame(frame*FRAME<TV>(TV(-2-2*cos((T)pi/3),4*sin((T)pi/3)*(T).5,0)*scale,ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>(2*(T)pi/3,TV(1,0,0))));break;}}

    links[6]=links[0];
    for(int i=1;i<=6;i++){
        JOINT<TV>* joint=new POINT_JOINT<TV>();
        arb.joint_mesh.Add_Articulation(links[i-1]->particle_index,links[i]->particle_index,joint);
        joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,2*scale)));
        joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(0,0,-2*scale)));}
}
template<class TV> void RIGIDS_STANDARD_TESTS<TV>::
Set_Joint_Frames(JOINT_ID id,const TV& location)
{
    RIGID_BODY<TV>& parent=*rigid_body_collection.articulated_rigid_body.Parent(id);
    rigid_body_collection.articulated_rigid_body.joint_mesh(id)->Set_Joint_To_Parent_Frame(FRAME<TV>(-parent.Object_Space_Point(location)));
    rigid_body_collection.articulated_rigid_body.Set_Consistent_Child_Frame(id);
}
template<class TV> JOINT_ID RIGIDS_STANDARD_TESTS<TV>::
Connect_With_Point_Joint(RIGID_BODY<TV>& parent,RIGID_BODY<TV>& child,const TV& location)
{
    POINT_JOINT<TV>* joint=new POINT_JOINT<TV>;
    rigid_body_collection.articulated_rigid_body.joint_mesh.Add_Articulation(parent.particle_index,child.particle_index,joint);
    Set_Joint_Frames(joint->id_number,location);
    return joint->id_number;
}
//#####################################################################
#define INSTANTIATION_HELPER_1D(T) \
    template RIGID_BODY<VECTOR<T,1> >& RIGIDS_STANDARD_TESTS<VECTOR<T,1> >::Add_Analytic_Box(VECTOR<T,1> const&); \
    template RIGIDS_STANDARD_TESTS<VECTOR<T,1> >::RIGIDS_STANDARD_TESTS(std::string&,const STREAM_TYPE&,RIGID_BODY_COLLECTION<VECTOR<T,1> >&);

#define INSTANTIATION_HELPER_ALL(T,d) \
    template RIGIDS_STANDARD_TESTS<VECTOR<T,d> >::RIGIDS_STANDARD_TESTS(std::string&,const STREAM_TYPE&,RIGID_BODY_COLLECTION<VECTOR<T,d> >&); \
    template RIGID_BODY<VECTOR<T,d> >& RIGIDS_STANDARD_TESTS<VECTOR<T,d> >::Add_Rigid_Body(const std::string&,const T,const T,const bool,const bool,const STREAM_TYPE*); \
    template RIGID_BODY<VECTOR<T,d> >& RIGIDS_STANDARD_TESTS<VECTOR<T,d> >::Add_Ground(const T,const T,const T,const T);

#define INSTANTIATION_HELPER(T) \
    INSTANTIATION_HELPER_1D(T) \
    INSTANTIATION_HELPER_ALL(T,2) \
    INSTANTIATION_HELPER_ALL(T,3) \
    template RIGID_BODY<VECTOR<T,2> >& RIGIDS_STANDARD_TESTS<VECTOR<T,2> >::Add_Analytic_Box(VECTOR<T,2> const&,int); \
    template RIGID_BODY<VECTOR<T,3> >& RIGIDS_STANDARD_TESTS<VECTOR<T,3> >::Add_Analytic_Box(VECTOR<T,3> const&); \
    template RIGID_BODY<VECTOR<T,3> >& RIGIDS_STANDARD_TESTS<VECTOR<T,3> >::Add_Analytic_Torus(const T,const T,int,int); \
    template RIGID_BODY<VECTOR<T,3> >& RIGIDS_STANDARD_TESTS<VECTOR<T,3> >::Add_Analytic_Cylinder(const T,const T,int,int); \
    template void RIGIDS_STANDARD_TESTS<VECTOR<T,3> >::Make_Lathe_Chain(const FRAME<VECTOR<T,3> >&,const T,const T,const T);

INSTANTIATION_HELPER(float);
template JOINT_ID RIGIDS_STANDARD_TESTS<VECTOR<float,3> >::Connect_With_Point_Joint(RIGID_BODY<VECTOR<float,3> >&,RIGID_BODY<VECTOR<float,3> >&,VECTOR<float,3> const&);
template RIGID_BODY<VECTOR<float,2> >& RIGIDS_STANDARD_TESTS<VECTOR<float,2> >::Add_Analytic_Sphere(float,float,int);
template RIGID_BODY<VECTOR<float,3> >& RIGIDS_STANDARD_TESTS<VECTOR<float,3> >::Add_Analytic_Sphere(float,float,int);
template RIGID_BODY<VECTOR<float,2> >& RIGIDS_STANDARD_TESTS<VECTOR<float,2> >::Add_Analytic_Ellipse(const VECTOR<float,2>,const T,int);
template RIGID_BODY<VECTOR<float,3> >& RIGIDS_STANDARD_TESTS<VECTOR<float,3> >::Add_Analytic_Ellipsoid(const VECTOR<float,3>,const T,int);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double);
template JOINT_ID RIGIDS_STANDARD_TESTS<VECTOR<double,3> >::Connect_With_Point_Joint(RIGID_BODY<VECTOR<double,3> >&,RIGID_BODY<VECTOR<double,3> >&,VECTOR<double,3> const&);
template RIGID_BODY<VECTOR<double,2> >& RIGIDS_STANDARD_TESTS<VECTOR<double,2> >::Add_Analytic_Sphere(double,double,int);
template RIGID_BODY<VECTOR<double,3> >& RIGIDS_STANDARD_TESTS<VECTOR<double,3> >::Add_Analytic_Sphere(double,double,int);
template RIGID_BODY<VECTOR<double,2> >& RIGIDS_STANDARD_TESTS<VECTOR<double,2> >::Add_Analytic_Ellipse(const VECTOR<double,2>,const T,int);
template RIGID_BODY<VECTOR<double,3> >& RIGIDS_STANDARD_TESTS<VECTOR<double,3> >::Add_Analytic_Ellipsoid(const VECTOR<double,3>,const T,int);
#endif


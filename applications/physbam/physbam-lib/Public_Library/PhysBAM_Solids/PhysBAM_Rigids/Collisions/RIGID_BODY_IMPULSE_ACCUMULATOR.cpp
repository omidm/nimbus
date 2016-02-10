//#####################################################################
// Copyright 2004-2005, Zhaosheng Bao, Eran Guendelman, Sergey Koltakov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_IMPULSE_ACCUMULATOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
using namespace PhysBAM;
//#####################################################################
// Function Reset
//#####################################################################
template<class TV,int d> void RIGID_BODY_IMPULSE_ACCUMULATOR<TV,d>::
Reset()
{
    accumulated_impulse=TWIST<TV>();
    if(accumulated_node_impulses){
        accumulated_node_impulses->Resize(simplicial_object->particles.array_collection->Size(),false,false);
        ARRAYS_COMPUTATIONS::Fill(*accumulated_node_impulses,TV());}
}
//#####################################################################
// Function Initialize_Simplicial_Object
//#####################################################################
template<class TV,int d> void RIGID_BODY_IMPULSE_ACCUMULATOR<TV,d>::
Initialize_Simplicial_Object(T_SIMPLICIAL_OBJECT* simplicial_object_input,ARRAY<TV>* accumulated_node_impulses_input)
{
    simplicial_object=simplicial_object_input;
    accumulated_node_impulses=accumulated_node_impulses_input;
    {std::stringstream ss;ss<<" "<<simplicial_object->particles.array_collection->Size();LOG::filecout(ss.str());}
}
//#####################################################################
// Function Add_Impulse
//#####################################################################
template<class T_OBJECT,class T,class TV> static int
Inside_Helper(T_OBJECT& object,const TV& object_space_location,const T surface_roughness)
{
    //return object.Inside(object_space_location,surface_roughness);
    return 0;
}
template<class T,class TV> static int
Inside_Helper(TRIANGULATED_SURFACE<T>& triangulated_surface,const TV& object_space_location,const T surface_roughness)
{
    if(!triangulated_surface.triangle_list) triangulated_surface.Update_Triangle_List();
    if(!triangulated_surface.hierarchy) triangulated_surface.Initialize_Hierarchy();
    int simplex=0;triangulated_surface.Inside_Any_Triangle(object_space_location,simplex,surface_roughness);
    assert(simplex);return simplex;
}   
template<class T,class TV> static int
Inside_Helper(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,const TV& object_space_location,const T surface_roughness)
{
    if(!tetrahedralized_volume.tetrahedron_list) tetrahedralized_volume.Update_Tetrahedron_List();
    if(!tetrahedralized_volume.hierarchy) tetrahedralized_volume.Initialize_Hierarchy(true);
    return tetrahedralized_volume.Inside(object_space_location,surface_roughness);
}
template<class TV,int d> void RIGID_BODY_IMPULSE_ACCUMULATOR<TV,d>::
Add_Impulse(const TV& location,const TWIST<TV>& impulse)
{
    accumulated_impulse+=impulse;

    if(simplicial_object && accumulated_node_impulses){
        typedef typename BASIC_SIMPLEX_POLICY<TV,d>::SIMPLEX T_SIMPLEX;
        TV object_space_location=rigid_body.Object_Space_Point(location);
        TV object_space_impulse=rigid_body.Object_Space_Vector(impulse.linear);
        ARRAY_VIEW<TV> X(simplicial_object->particles.X);
        int simplex=Inside_Helper(*simplicial_object,object_space_location,surface_roughness);
        if(simplex){
            VECTOR<int,d+1>& nodes=simplicial_object->mesh.elements(simplex);
            //VECTOR<T,d+1> weights=T_SIMPLEX::Barycentric_Coordinates(object_space_location,X.Subset(nodes));
            VECTOR<T,d+1> weights;
            for(int i=1;i<=weights.m;i++) (*accumulated_node_impulses)(nodes[i])+=weights[i]*object_space_impulse;}}
}
//#####################################################################
template class RIGID_BODY_IMPULSE_ACCUMULATOR<VECTOR<float,1>,0>;
template class RIGID_BODY_IMPULSE_ACCUMULATOR<VECTOR<float,1>,1>;
template class RIGID_BODY_IMPULSE_ACCUMULATOR<VECTOR<float,2>,1>;
template class RIGID_BODY_IMPULSE_ACCUMULATOR<VECTOR<float,2>,2>;
template class RIGID_BODY_IMPULSE_ACCUMULATOR<VECTOR<float,3>,2>;
template class RIGID_BODY_IMPULSE_ACCUMULATOR<VECTOR<float,3>,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_BODY_IMPULSE_ACCUMULATOR<VECTOR<double,1>,0>;
template class RIGID_BODY_IMPULSE_ACCUMULATOR<VECTOR<double,1>,1>;
template class RIGID_BODY_IMPULSE_ACCUMULATOR<VECTOR<double,2>,1>;
template class RIGID_BODY_IMPULSE_ACCUMULATOR<VECTOR<double,2>,2>;
template class RIGID_BODY_IMPULSE_ACCUMULATOR<VECTOR<double,3>,2>;
template class RIGID_BODY_IMPULSE_ACCUMULATOR<VECTOR<double,3>,3>;
#endif

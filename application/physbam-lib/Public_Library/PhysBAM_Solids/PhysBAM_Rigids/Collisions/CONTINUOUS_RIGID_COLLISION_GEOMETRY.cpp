//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_FRAME.h>
#endif
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_POINT_SIMPLICES_1D_INTERSECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_SEGMENTED_CURVE_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_TRIANGULATED_SURFACE_INTERSECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/CONTINUOUS_RIGID_COLLISION_GEOMETRY.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> CONTINUOUS_RIGID_COLLISION_GEOMETRY<TV>::
CONTINUOUS_RIGID_COLLISION_GEOMETRY(RIGID_GEOMETRY<TV>& rigid_geometry_input)
    :BASE(rigid_geometry_input),enlarged_box_up_to_date(false)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> CONTINUOUS_RIGID_COLLISION_GEOMETRY<TV>::
~CONTINUOUS_RIGID_COLLISION_GEOMETRY()
{
}
//#####################################################################
// Function Update_Bounding_Box
//#####################################################################
template<class TV> void CONTINUOUS_RIGID_COLLISION_GEOMETRY<TV>::
Update_Bounding_Box()
{
    BASE::Update_Bounding_Box();
    if(saved_states.m){
        enlarged_box_up_to_date=true;
        enlarged_box=RANGE<TV>::Combine(rigid_geometry.Axis_Aligned_Bounding_Box(),T_ORIENTED_BOX(rigid_geometry.Object_Space_Bounding_Box(),saved_states(1).x).Axis_Aligned_Bounding_Box());}
    else enlarged_box_up_to_date=false;
}
//#####################################################################
// Function Axis_Aligned_Bounding_Box
//#####################################################################
template<class TV> const RANGE<TV>& CONTINUOUS_RIGID_COLLISION_GEOMETRY<TV>::
Axis_Aligned_Bounding_Box() const
{
    if(saved_states.m){assert(enlarged_box_up_to_date);return enlarged_box;}
    else return BASE::Axis_Aligned_Bounding_Box();
}
//#####################################################################
template class CONTINUOUS_RIGID_COLLISION_GEOMETRY<VECTOR<float,1> >;
template class CONTINUOUS_RIGID_COLLISION_GEOMETRY<VECTOR<float,2> >;
template class CONTINUOUS_RIGID_COLLISION_GEOMETRY<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class CONTINUOUS_RIGID_COLLISION_GEOMETRY<VECTOR<double,1> >;
template class CONTINUOUS_RIGID_COLLISION_GEOMETRY<VECTOR<double,2> >;
template class CONTINUOUS_RIGID_COLLISION_GEOMETRY<VECTOR<double,3> >;
#endif
}

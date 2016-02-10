//#####################################################################
// Copyright 2012, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/SEGMENT_3D_TRIANGLE_3D_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(const TRIANGLE_3D<T>& triangle1,const TRIANGLE_3D<T>& triangle2,const T thickness_over_two)
{
    if(Intersects(SEGMENT_3D<T>(triangle1.x1,triangle1.x2),triangle2)) return true;
    if(Intersects(SEGMENT_3D<T>(triangle1.x2,triangle1.x3),triangle2)) return true;
    if(Intersects(SEGMENT_3D<T>(triangle1.x3,triangle1.x1),triangle2)) return true;
    if(Intersects(SEGMENT_3D<T>(triangle2.x1,triangle2.x2),triangle1)) return true;
    if(Intersects(SEGMENT_3D<T>(triangle2.x2,triangle2.x3),triangle1)) return true;
    if(Intersects(SEGMENT_3D<T>(triangle2.x3,triangle2.x1),triangle1)) return true;
    return false;
}
//#####################################################################
template bool Intersects(const TRIANGLE_3D<float>&,const TRIANGLE_3D<float>&,const float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(const TRIANGLE_3D<double>&,const TRIANGLE_3D<double>&,const double);
#endif
};
};

//#####################################################################
// Copyright 2012, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#ifndef __TRIANGLE_3D_TRIANGLE_3D_INTERSECTION__
#define __TRIANGLE_3D_TRIANGLE_3D_INTERSECTION__
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
template<class T> bool Intersects(const TRIANGLE_3D<T>& segment,const TRIANGLE_3D<T>& triangle,const T thickness_over_two=0);
//#####################################################################
};
};
#endif

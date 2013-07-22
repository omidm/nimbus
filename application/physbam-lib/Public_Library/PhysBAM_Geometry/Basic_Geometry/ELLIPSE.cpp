//#####################################################################
// Copyright 2011, Mridul Aanajneya.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Basic_Geometry/ELLIPSE.h>
namespace PhysBAM{
//#####################################################################
// Function Approximate_Surface
//#####################################################################
template<class T> T ELLIPSE<T>::
Approximate_Signed_Distance(const TV& location) const       
{
    TV offset=orientation.Inverse_Rotate(location-center),scaled_offset=radii.Solve_Linear_System(offset);
    T scaled_magnitude=scaled_offset.Normalize();
    TV closest_offset=radii*scaled_offset;
    T distance=(closest_offset-offset).Magnitude();
    return scaled_magnitude<1?-distance:distance;
}
//#####################################################################
template class ELLIPSE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ELLIPSE<double>;
#endif
}

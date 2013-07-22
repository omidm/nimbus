//#####################################################################
// Copyright 2002-2005, Eran Guendelman, Geoffrey Irving, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REMOVED_PARTICLES_BLENDER_3D
//##################################################################### 
#ifndef __REMOVED_PARTICLES_BLENDER_3D__
#define __REMOVED_PARTICLES_BLENDER_3D__

#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Tools/Polynomials/CUBIC.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class TV> class RANGE;
template<class T> class ORIENTED_BOX;
template<class T> class ELLIPSOID;

template<class T>
class REMOVED_PARTICLES_BLENDER_3D
{
    typedef VECTOR<T,3> TV;
    T blending_parameter;
    CUBIC<T> kernel;
    T small_number;
public:
    T R;

    REMOVED_PARTICLES_BLENDER_3D(T blending_parameter_input);

    T C(const T r) const
    {if(r>R) return 0;
    return kernel(r/R);}

    T C_Prime(const T r) const
    {if(r>R) return 0;
    return kernel.Prime(r/R);}
    
    static T Get_Distance(const T one_over_radius_x_squared,const T one_over_radius_yz_squared,const TV& position,const TV& major_axis,const TV& location)
    {TV dX=location-position;T dot=TV::Dot_Product(dX,major_axis);
    return sqrt(one_over_radius_yz_squared*dX.Magnitude_Squared()+(one_over_radius_x_squared-one_over_radius_yz_squared)*sqr(dot));}
    
    static T Get_Distance(const TV& position,const SYMMETRIC_MATRIX<T,3>& metric_tensor,const TV& location)
    {TV dX=location-position;
    return sqrt(TV::Dot_Product(dX,metric_tensor*dX));}

//#####################################################################
    ORIENTED_BOX<TV> Get_Oriented_Bounding_Box(T radius_x,T radius_yz,const TV& position,const TV& major_axis) const;
    RANGE<TV> Get_Bounding_Box(T radius_x,T radius_yz,const TV& position,const TV& major_axis) const;
    RANGE<TV> Get_Bounding_Box(const ELLIPSOID<T>& ellipsoid) const;
//#####################################################################
};
}
#endif

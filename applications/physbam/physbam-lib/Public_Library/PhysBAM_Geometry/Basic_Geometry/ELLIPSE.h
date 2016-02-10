//#####################################################################
// Copyright 2011, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ELLIPSE
//#####################################################################
#ifndef __ELLIPSE__
#define __ELLIPSE__

#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
namespace PhysBAM{

template<class T>
class ELLIPSE
{
    typedef VECTOR<T,2> TV;
public:
    typedef TV VECTOR_T;
    TV center;
    DIAGONAL_MATRIX<T,2> radii;
    ROTATION<TV> orientation;    
    
    ELLIPSE()
        :radii(1,1)
    {}

    ELLIPSE(TV center_input,DIAGONAL_MATRIX<T,2> radii_input,ROTATION<TV> orientation_input=ROTATION<TV>())
        :center(center_input),radii(radii_input),orientation(orientation_input.Normalized())
    {}

    ~ELLIPSE()
    {}

    TV Normal(const TV& location) const //this is not correct, can't be used now
    {return (location-center).Normalized();}

    T Signed_Distance(const TV& location) const
    {return Approximate_Signed_Distance(location);}

    RANGE<TV> Bounding_Box() const
    {return RANGE<TV>(center(1)-radii(1),center(1)+radii(1),center(2)-radii(2),center(2)+radii(2));}

    VECTOR<T,1> Principal_Curvatures(const TV& X) const //this is not correct, can't be used now
    {return VECTOR<T,1>::All_Ones_Vector()/radii.x11;}
    
    static std::string Name()
    {return "ELLIPSE<T>";}

    T Approximate_Signed_Distance(const TV& location) const; // not the true signed distance, but has correct inside/outside
};
}
#endif

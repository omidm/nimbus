//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONTINUOUS_RIGID_COLLISION_GEOMETRY
//#####################################################################
#ifndef __CONTINUOUS_RIGID_COLLISION_GEOMETRY__
#define __CONTINUOUS_RIGID_COLLISION_GEOMETRY__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_1D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
namespace PhysBAM{

template<class TV> class RIGID_GEOMETRY;

template<class TV_input>
class CONTINUOUS_RIGID_COLLISION_GEOMETRY:public RIGID_COLLISION_GEOMETRY<TV_input>
{
private:
    typedef TV_input TV;
    typedef typename TV::SCALAR T;
    typedef RIGID_COLLISION_GEOMETRY<TV> BASE;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::ORIENTED_BOX T_ORIENTED_BOX;    

public:
    using BASE::saved_states;using BASE::rigid_geometry;
    bool enlarged_box_up_to_date;
    mutable RANGE<TV> enlarged_box;

    CONTINUOUS_RIGID_COLLISION_GEOMETRY(RIGID_GEOMETRY<TV>& rigid_geometry_input);
    virtual ~CONTINUOUS_RIGID_COLLISION_GEOMETRY();

//#####################################################################
    void Update_Bounding_Box();
    const RANGE<TV>& Axis_Aligned_Bounding_Box() const;
//#####################################################################
};
}
#endif

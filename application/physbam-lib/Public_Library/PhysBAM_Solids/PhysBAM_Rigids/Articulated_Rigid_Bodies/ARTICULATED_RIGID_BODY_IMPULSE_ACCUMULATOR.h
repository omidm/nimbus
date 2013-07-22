//#####################################################################
// Copyright 2008, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR
//#####################################################################
#ifndef __ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR__
#define __ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR.h>
namespace PhysBAM{

template<class TV> class JOINT;
template<class TV> class ARTICULATED_RIGID_BODY;
template<class TV>
class ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR:public COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR<TV>
{
    typedef typename TV::SCALAR T;
public:
    JOINT<TV>& joint;
    ARTICULATED_RIGID_BODY<TV>& arb;
    T surface_roughness;
    T weight;
    TWIST<TV> accumulated_impulse; // in world space

    ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR(JOINT<TV>& joint_input,ARTICULATED_RIGID_BODY<TV>& arb_input,const T surface_roughness_input=(T)1e-4)
        :joint(joint_input),arb(arb_input),surface_roughness(surface_roughness_input),weight(1)
    {
        Reset();
    }

    T Energy()
    {
        if(TV::m==1) return 0; //in 1D there is no rotational energy
        return weight*(accumulated_impulse.angular.Magnitude());
    }

//#####################################################################
    void Reset() PHYSBAM_OVERRIDE;
    void Add_Impulse(const TV& location,const TWIST<TV>& impulse) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

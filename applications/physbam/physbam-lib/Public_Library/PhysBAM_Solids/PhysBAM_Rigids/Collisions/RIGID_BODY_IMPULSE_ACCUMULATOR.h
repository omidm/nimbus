//#####################################################################
// Copyright 2004-2005, Zhaosheng Bao, Eran Guendelman, Sergey Koltakov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_IMPULSE_ACCUMULATOR
//#####################################################################
#ifndef __RIGID_BODY_IMPULSE_ACCUMULATOR__
#define __RIGID_BODY_IMPULSE_ACCUMULATOR__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY;
template<class TV> class TWIST;

template<class TV,int d>
class RIGID_BODY_IMPULSE_ACCUMULATOR:public COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d>::OBJECT T_SIMPLICIAL_OBJECT;
public:
    RIGID_BODY<TV>& rigid_body;
    T surface_roughness;
    TWIST<TV> accumulated_impulse; // in world space
    ARRAY<TV>* accumulated_node_impulses; // in object space
    T_SIMPLICIAL_OBJECT* simplicial_object; // in object space

    RIGID_BODY_IMPULSE_ACCUMULATOR(RIGID_BODY<TV>& rigid_body_input,const T surface_roughness_input=(T)1e-4)
        :rigid_body(rigid_body_input),surface_roughness(surface_roughness_input),accumulated_node_impulses(0),simplicial_object(0)
    {
        Reset();
    }

//#####################################################################
    void Reset() PHYSBAM_OVERRIDE;
    void Initialize_Simplicial_Object(T_SIMPLICIAL_OBJECT* simplicial_object_input,ARRAY<TV>* accumulated_node_impulses_input);
    void Add_Impulse(const TV& location,const TWIST<TV>& impulse) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

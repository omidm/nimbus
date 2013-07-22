//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLES_PARAMETERS
//#####################################################################
#ifndef __DEFORMABLES_PARAMETERS__
#define __DEFORMABLES_PARAMETERS__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Collisions/COLLISIONS_GEOMETRY_FORWARD.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_EVOLUTION_PARAMETERS.h>
namespace PhysBAM{

template<class TV> class TRIANGLE_COLLISION_PARAMETERS;
template<class TV> class IMPLICIT_SOLVE_PARAMETERS;
template<class TV> class DEFORMABLE_OBJECT_COLLISION_PARAMETERS;
//template<class TV> class DEFORMABLE_OBJECT_EVOLUTION_PARAMETERS;

template<class TV>
class DEFORMABLES_PARAMETERS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    bool verbose,verbose_dt;
    T cfl;
    T min_dt;

    bool enforce_repulsions_in_cg;
    bool use_post_cg_constraints;

    RIGID_GEOMETRY_EVOLUTION_PARAMETERS rigid_geometry_evolution_parameters;
    IMPLICIT_SOLVE_PARAMETERS<TV>& implicit_solve_parameters;
    bool use_trapezoidal_rule_for_velocities; // otherwise use backward euler from n->n+1

    TRIANGLE_COLLISION_PARAMETERS<TV>& triangle_collision_parameters;
    DEFORMABLE_OBJECT_COLLISION_PARAMETERS<TV>& deformable_object_collision_parameters;
    //DEFORMABLE_OBJECT_EVOLUTION_PARAMETERS<TV>& deformable_object_evolution_parameters;

    DEFORMABLES_PARAMETERS();
    virtual ~DEFORMABLES_PARAMETERS();
//#####################################################################
};
}
#endif

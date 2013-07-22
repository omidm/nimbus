//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_PARAMETERS
//#####################################################################
#ifndef __RIGIDS_PARAMETERS__
#define __RIGIDS_PARAMETERS__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Collisions/COLLISIONS_GEOMETRY_FORWARD.h>
namespace PhysBAM{

template<class TV> class TRIANGLE_COLLISION_PARAMETERS;
template<class TV> class IMPLICIT_SOLVE_PARAMETERS;
template<class TV> class RIGID_BODY_COLLISION_PARAMETERS;
template<class TV> class RIGID_BODY_EVOLUTION_PARAMETERS;

template<class TV>
class RIGIDS_PARAMETERS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    bool verbose,verbose_dt;
    T cfl;
    T min_dt;

    bool use_spatial_partition_for_levelset_collision_objects;
    SPATIAL_PARTITION_VOXEL_SIZE_HEURISTIC spatial_partition_voxel_size_heuristic;
    int spatial_partition_number_of_cells;
    T spatial_partition_voxel_size_scale_factor;
    int rigid_cluster_fracture_frequency;

    bool use_post_cg_constraints;

    IMPLICIT_SOLVE_PARAMETERS<TV>& implicit_solve_parameters;
    bool use_trapezoidal_rule_for_velocities; // otherwise use backward euler from n->n+1
    bool enforce_poststabilization_in_cg;

    int threadid;

    RIGID_BODY_COLLISION_PARAMETERS<TV>& rigid_body_collision_parameters;
    RIGID_BODY_EVOLUTION_PARAMETERS<TV>& rigid_body_evolution_parameters;

    RIGIDS_PARAMETERS();
    virtual ~RIGIDS_PARAMETERS();
//#####################################################################
};
}
#endif

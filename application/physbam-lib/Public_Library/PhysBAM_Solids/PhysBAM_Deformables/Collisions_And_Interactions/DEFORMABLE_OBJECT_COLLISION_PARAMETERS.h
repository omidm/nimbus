//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_OBJECT_COLLISION_PARAMETERS
//#####################################################################
#ifndef __DEFORMABLE_OBJECT_COLLISION_PARAMETERS__
#define __DEFORMABLE_OBJECT_COLLISION_PARAMETERS__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Collisions/COLLISIONS_GEOMETRY_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_OBJECT_FORWARD.h>
namespace PhysBAM{

template<class TV>
class DEFORMABLE_OBJECT_COLLISION_PARAMETERS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    T collision_tolerance;
    bool collide_with_interior;
    bool perform_collision_body_collisions;
    bool use_spatial_partition_for_levelset_collision_objects;
    bool disable_multiple_levelset_collisions;
    bool use_existing_contact;
    T maximum_levelset_collision_projection_velocity;
    SPATIAL_PARTITION_VOXEL_SIZE_HEURISTIC spatial_partition_voxel_size_heuristic;
    int spatial_partition_number_of_cells;
    T spatial_partition_voxel_size_scale_factor;

    DEFORMABLE_OBJECT_COLLISION_PARAMETERS();
    virtual ~DEFORMABLE_OBJECT_COLLISION_PARAMETERS();
//#####################################################################
};
}
#endif

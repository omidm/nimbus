//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_OBJECT_COLLISION_PARAMETERS
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Collision_Detection/COLLISION_GEOMETRY_SPATIAL_PARTITION.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/STRUCTURE_INTERACTION_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLE_OBJECT_COLLISION_PARAMETERS<TV>::
DEFORMABLE_OBJECT_COLLISION_PARAMETERS()
    :collision_tolerance((T)1e-6),collide_with_interior(false),perform_collision_body_collisions(true),use_spatial_partition_for_levelset_collision_objects(false),
    disable_multiple_levelset_collisions(true),use_existing_contact(false),maximum_levelset_collision_projection_velocity(FLT_MAX),
    spatial_partition_voxel_size_heuristic(SPATIAL_PARTITION_MAX_BOX_SIZE),spatial_partition_number_of_cells(100),spatial_partition_voxel_size_scale_factor((T)8)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DEFORMABLE_OBJECT_COLLISION_PARAMETERS<TV>::
~DEFORMABLE_OBJECT_COLLISION_PARAMETERS()
{
}
//#####################################################################
template class DEFORMABLE_OBJECT_COLLISION_PARAMETERS<VECTOR<float,1> >;
template class DEFORMABLE_OBJECT_COLLISION_PARAMETERS<VECTOR<float,2> >;
template class DEFORMABLE_OBJECT_COLLISION_PARAMETERS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DEFORMABLE_OBJECT_COLLISION_PARAMETERS<VECTOR<double,1> >;
template class DEFORMABLE_OBJECT_COLLISION_PARAMETERS<VECTOR<double,2> >;
template class DEFORMABLE_OBJECT_COLLISION_PARAMETERS<VECTOR<double,3> >;
#endif
}

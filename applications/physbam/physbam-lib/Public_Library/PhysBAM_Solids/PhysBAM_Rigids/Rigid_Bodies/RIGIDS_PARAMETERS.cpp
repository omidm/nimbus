//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_PARAMETERS
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR_1D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGIDS_PARAMETERS.h>
#include <climits>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGIDS_PARAMETERS<TV>::
RIGIDS_PARAMETERS()
    :verbose(true),verbose_dt(false),cfl((T).9),min_dt((T)0),use_spatial_partition_for_levelset_collision_objects(false),
    spatial_partition_voxel_size_heuristic(SPATIAL_PARTITION_MAX_BOX_SIZE),spatial_partition_number_of_cells(100),
    spatial_partition_voxel_size_scale_factor((T)8),rigid_cluster_fracture_frequency(INT_MAX),use_post_cg_constraints(true),
    implicit_solve_parameters(*new IMPLICIT_SOLVE_PARAMETERS<TV>),
    use_trapezoidal_rule_for_velocities(true),enforce_poststabilization_in_cg(true),
    threadid(1),
    rigid_body_collision_parameters(*new RIGID_BODY_COLLISION_PARAMETERS<TV>),
    rigid_body_evolution_parameters(*new RIGID_BODY_EVOLUTION_PARAMETERS<TV>)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGIDS_PARAMETERS<TV>::
~RIGIDS_PARAMETERS()
{
    delete &implicit_solve_parameters;
    delete &rigid_body_collision_parameters;
    delete &rigid_body_evolution_parameters;
}
//#####################################################################
template class RIGIDS_PARAMETERS<VECTOR<float,1> >;
template class RIGIDS_PARAMETERS<VECTOR<float,2> >;
template class RIGIDS_PARAMETERS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGIDS_PARAMETERS<VECTOR<double,1> >;
template class RIGIDS_PARAMETERS<VECTOR<double,2> >;
template class RIGIDS_PARAMETERS<VECTOR<double,3> >;
#endif
}

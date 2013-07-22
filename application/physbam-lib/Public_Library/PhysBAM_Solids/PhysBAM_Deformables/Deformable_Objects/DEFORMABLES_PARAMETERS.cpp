//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLES_PARAMETERS
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR_1D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLES_PARAMETERS.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLES_PARAMETERS<TV>::
DEFORMABLES_PARAMETERS()
    :verbose(true),verbose_dt(false),cfl((T).9),min_dt((T)0),enforce_repulsions_in_cg(true),use_post_cg_constraints(true),
    implicit_solve_parameters(*new IMPLICIT_SOLVE_PARAMETERS<TV>()),use_trapezoidal_rule_for_velocities(false),
    triangle_collision_parameters(*new TRIANGLE_COLLISION_PARAMETERS<TV>()),
    deformable_object_collision_parameters(*new DEFORMABLE_OBJECT_COLLISION_PARAMETERS<TV>())
    //deformable_object_evolution_parameters(*new DEFORMABLE_OBJECT_EVOLUTION_PARAMETERS<TV>())
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DEFORMABLES_PARAMETERS<TV>::
~DEFORMABLES_PARAMETERS()
{
    delete &implicit_solve_parameters;
    delete &triangle_collision_parameters;
    delete &deformable_object_collision_parameters;
    //delete &deformable_object_evolution_parameters;
}
//#####################################################################
template class DEFORMABLES_PARAMETERS<VECTOR<float,1> >;
template class DEFORMABLES_PARAMETERS<VECTOR<float,2> >;
template class DEFORMABLES_PARAMETERS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DEFORMABLES_PARAMETERS<VECTOR<double,1> >;
template class DEFORMABLES_PARAMETERS<VECTOR<double,2> >;
template class DEFORMABLES_PARAMETERS<VECTOR<double,3> >;
#endif
}

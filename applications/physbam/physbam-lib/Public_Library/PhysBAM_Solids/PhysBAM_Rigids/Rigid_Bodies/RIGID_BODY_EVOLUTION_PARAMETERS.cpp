//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_EVOLUTION_PARAMETERS
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_BODY_EVOLUTION_PARAMETERS<TV>::
RIGID_BODY_EVOLUTION_PARAMETERS()
    :simulate_rigid_bodies(false),write_rigid_bodies(true),rigid_body_ether_viscosity(0),max_rigid_body_rotation_per_time_step((T).1*(T)pi),
    max_rigid_body_linear_movement_fraction_per_time_step((T).1),minimum_rigid_body_time_step_fraction(0),maximum_rigid_body_time_step_fraction((T)1.1),clamp_rigid_body_velocities(false),
    max_rigid_body_linear_velocity(0),max_rigid_body_angular_velocity(0),rigid_cfl((T).5),rigid_minimum_dt(0),rigid_maximum_dt(FLT_MAX),
    correct_evolution_energy(false),residual_push_out_depth(0),correct_contact_energy(false),fix_velocities_in_position_update(false),threadid(1)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_BODY_EVOLUTION_PARAMETERS<TV>::
~RIGID_BODY_EVOLUTION_PARAMETERS()
{
}
//#####################################################################
template class RIGID_BODY_EVOLUTION_PARAMETERS<VECTOR<float,1> >;
template class RIGID_BODY_EVOLUTION_PARAMETERS<VECTOR<float,2> >;
template class RIGID_BODY_EVOLUTION_PARAMETERS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_BODY_EVOLUTION_PARAMETERS<VECTOR<double,1> >;
template class RIGID_BODY_EVOLUTION_PARAMETERS<VECTOR<double,2> >;
template class RIGID_BODY_EVOLUTION_PARAMETERS<VECTOR<double,3> >;
#endif
}

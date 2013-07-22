//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_EVOLUTION_PARAMETERS
//#####################################################################
#ifndef __RIGID_BODY_EVOLUTION_PARAMETERS__
#define __RIGID_BODY_EVOLUTION_PARAMETERS__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_EVOLUTION_PARAMETERS.h>
namespace PhysBAM{

template<class TV>
class RIGID_BODY_EVOLUTION_PARAMETERS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    bool simulate_rigid_bodies;
    bool write_rigid_bodies;
    T rigid_body_ether_viscosity;
    T max_rigid_body_rotation_per_time_step;
    T max_rigid_body_linear_movement_fraction_per_time_step;
    T minimum_rigid_body_time_step_fraction; // fraction of frame time
    T maximum_rigid_body_time_step_fraction;
    bool clamp_rigid_body_velocities;
    T max_rigid_body_linear_velocity,max_rigid_body_angular_velocity; // magnitudes of vectors
    T rigid_cfl;
    T rigid_minimum_dt;
    T rigid_maximum_dt;
    RIGID_GEOMETRY_EVOLUTION_PARAMETERS rigid_geometry_evolution_parameters;
    bool correct_evolution_energy;
    T residual_push_out_depth;
    bool correct_contact_energy;
    bool fix_velocities_in_position_update;

    int threadid;

    RIGID_BODY_EVOLUTION_PARAMETERS();
    virtual ~RIGID_BODY_EVOLUTION_PARAMETERS();
//#####################################################################
};
}
#endif

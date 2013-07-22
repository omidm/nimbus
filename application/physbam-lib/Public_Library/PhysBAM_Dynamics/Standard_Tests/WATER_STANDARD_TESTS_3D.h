//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class WATER_STANDARD_TESTS_3D
//#####################################################################
// Provides 5 tests that are consistent across all water codes:
//   1. Flat surface of water that should remain flat
//   2. Falling drop
//   3. 4 Sources spraying water at each of the walls
//   4. Sphere splashing into a pool of water
//   5. Glass filling
//   6. Wetting
//   7. drop moving upward with no gravity at constant velocity
//   8. Thin source spraying down with one-way SPH for removed negative particles
//   9. Two-way SPH with four falling drops of varying densities
//   10. Targeting densities for SPH particles.
//   11. plane moving upward under no gravity at constant velocity
//   12. falling drop constant viscosity
//   13. falling drop variable viscosity
//   15. falling drop analytic test downward velocity
//   16. drop expanding
//   20. sphere moving up from bottom out of water to test noise
//   21. flat surface of water
//   22. initial height of water noise on height between grid cells
// Also supports a variety of standard resolutions in powers of 2.
//#####################################################################
#ifndef __WATER_STANDARD_TESTS_3D__
#define __WATER_STANDARD_TESTS_3D__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/FLUID_COLLISION_BODY_INACCURATE_UNION.h>
namespace PhysBAM{

template<class T_GRID> class FLUIDS_PARAMETERS;
template<class TV> class SOLIDS_FLUIDS_EXAMPLE;
template<class TV> class RIGID_BODY_COLLECTION;

template<class T_GRID>
class WATER_STANDARD_TESTS_3D
{
    typedef typename T_GRID::SCALAR T;typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
public:
    SOLIDS_FLUIDS_EXAMPLE<TV>& example;
    FLUIDS_PARAMETERS<T_GRID>& fluids_parameters;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID> inaccurate_union;
    bool use_inaccurate_body_collisions;
    bool use_variable_density_for_sph,use_two_way_coupling_for_sph,convert_sph_particles_to_fluid,use_analytic_divergence,use_analytic_divergence_for_expansion_only,
         adjust_cell_weights_on_neumann_boundaries,enforce_density_near_interface;
    T flip_ratio,target_particles_per_unit_volume,neumann_boundary_slip_multiplier,ballistic_particles_as_percentage_of_target,particle_targeting_time;

    int test_number;
    GRID<TV> grid;
    ARRAY<CYLINDER<T> > sources;
    ARRAY<MATRIX<T,4> > world_to_source;
    ARRAY<TV> source_velocity;
    int sphere;
    INTERPOLATION_CURVE<T,TV> motion_curve;
    mutable RANDOM_NUMBERS<T> height_noise_random;

//#####################################################################
    WATER_STANDARD_TESTS_3D(SOLIDS_FLUIDS_EXAMPLE<TV>& example,FLUIDS_PARAMETERS<T_GRID>& fluids_parameters,RIGID_BODY_COLLECTION<TV>& rigid_body_collection);
    virtual ~WATER_STANDARD_TESTS_3D(){};
    void Initialize(const int test_number_input,const int resolution);
    void Initialize_Advection(const bool always_use_objects=false);
    virtual TV Initial_Velocity(const TV& X) const;
    virtual void Get_Variable_Viscosity(ARRAY<T,VECTOR<int,3> >& variable_viscosity,const T time) const;
    T Initial_Phi(const TV& X) const;
    T Initial_Phi_Object(const TV& X) const;
    void Initialize_Bodies();
    virtual void Update_Sources(const T time);
    virtual void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id);
    virtual bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id);
    virtual void Limit_Dt(T& dt,const T time);
    void Initialize_SPH_Particles();
    TV Analytic_Velocity(const T time,const TV& location) const;
//#####################################################################
};
}
#endif

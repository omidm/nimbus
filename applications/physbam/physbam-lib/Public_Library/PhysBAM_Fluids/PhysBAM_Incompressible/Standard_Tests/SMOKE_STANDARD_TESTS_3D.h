//#####################################################################
// Copyright 2006-2007, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SMOKE_STANDARD_TESTS_3D
//#####################################################################
// Provides 4 tests that are consistent across all smoke codes:
//   1. Plume
//   2. Plume past sphere
//   3. Explosion
//   4. Explosion with vortex particles
// Also supports a variety of resolutions.
//#####################################################################
#ifndef __SMOKE_STANDARD_TESTS_3D__
#define __SMOKE_STANDARD_TESTS_3D__

#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
namespace PhysBAM{

template<class TV> class BOX;
template<class TV> class GRID;
template<class TV> class SOLIDS_FLUIDS_EXAMPLE;
template<class TV> class RIGID_BODY_COLLECTION;
template<class T> class VORTEX_PARTICLE_EVOLUTION_3D;
template<class T_GRID> class FLUIDS_PARAMETERS_UNIFORM;
template<class T_GRID> class INCOMPRESSIBLE_FLUID_CONTAINER;

template<class T_GRID>
class SMOKE_STANDARD_TESTS_3D
{
    typedef typename T_GRID::SCALAR T;typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;
public:
    SOLIDS_FLUIDS_EXAMPLE<TV>& example;
    FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters;
    INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& incompressible_fluid_container;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    VORTEX_PARTICLE_EVOLUTION_3D<T>* vortex_particle_evolution;

    int test_number;
    GRID<TV> grid;
    BOX<TV> source;
    MATRIX<T,4> world_to_source;
    TV source_velocity;
    T rho;
    T explosion_divergence,explosion_end_time;
    T source_vorticity_magnitude,particle_vorticity_minimum_density,particle_radius;
    RANDOM_NUMBERS<T> random;

    SMOKE_STANDARD_TESTS_3D(SOLIDS_FLUIDS_EXAMPLE<TV>& example,FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters,INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& incompressible_fluid_container,
        RIGID_BODY_COLLECTION<TV>& rigid_body_collection);
    virtual ~SMOKE_STANDARD_TESTS_3D();

    bool Use_Source() const
    {return test_number==1 || test_number==2;}

//#####################################################################
    virtual TV Initial_Velocity(const TV& X) const;
    void Initialize(const int test_number_input,const int resolution);
    void Initialize_Bodies();
    void Get_Divergence(ARRAY<T,VECTOR<int,3> >& divergence,const T dt,const T time);
    void Get_Body_Force(ARRAY<T,FACE_INDEX<3> >& force,const T dt,const T time);
//#####################################################################
};
}
#endif

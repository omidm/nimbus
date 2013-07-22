//#####################################################################
// Copyright 2005, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class WATER_STANDARD_TESTS_MULTIPHASE
//#####################################################################
// Test descriptions:
//   11. Two drops colliding
//   12. Sphere splashing into a pool of two liquid phases
//   13. 4 phase splash
//   14. Rising air bubble in water
// Also supports a variety of standard resolutions in powers of 2.
//#####################################################################
#ifndef __WATER_STANDARD_TESTS_MULTIPHASE__
#define __WATER_STANDARD_TESTS_MULTIPHASE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID> class SOLIDS_FLUIDS_EXAMPLE_UNIFORM;
template<class T_GRID> class FLUIDS_PARAMETERS_UNIFORM;
template<class T_GRID> class INCOMPRESSIBLE_FLUID_CONTAINER;
template<class TV> class RIGID_BODY_COLLECTION;

template<class T_GRID,class T_WATER_STANDARD_TESTS>
class WATER_STANDARD_TESTS_MULTIPHASE:public T_WATER_STANDARD_TESTS
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_GRID::VECTOR_INT TV_INT;
public:
    typedef T_WATER_STANDARD_TESTS BASE;
    using BASE::world_to_source;using BASE::rigid_body_collection;using BASE::sources;using BASE::fluids_parameters;using BASE::grid;using BASE::example;using BASE::Initial_Phi;using BASE::sphere;
    using BASE::source_velocity;

    FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters_uniform;
    INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& incompressible_fluid_container;
    int test_number;
    ARRAY<int> source_region;
    bool use_open_wall;
    int air_region;

    T_LEVELSET* armadillo;
    
    WATER_STANDARD_TESTS_MULTIPHASE(SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>& example,FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters_input,INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& incompressible_fluid_container_input,
        RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input);
    ~WATER_STANDARD_TESTS_MULTIPHASE();

//#####################################################################
    void Initialize(const int test_number_input,const int resolution,const int restart_frame);
    virtual void Initialize_Advection(const bool always_use_objects=false);
    void Initialize_Bodies();
    static int Number_Of_Regions(int test_number);
    static int Non_Multiphase_Test_Number(int test_number);
    T Initial_Phi(const int region,const TV& X) const;
    TV Initial_Velocity(const TV& X) const;
    virtual void Update_Sources(const T time);
//#####################################################################
};
}
#endif

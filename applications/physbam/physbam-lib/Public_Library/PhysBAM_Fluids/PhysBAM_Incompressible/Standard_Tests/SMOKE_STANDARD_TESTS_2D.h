//#####################################################################
// Copyright 2006-2007, Jon Gretarsson, Avi Robinson-Mosher, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SMOKE_STANDARD_TESTS_2D
//#####################################################################
// Provides 3 tests that are consistent across all smoke codes:
//   1. Plume
//   2. Plume past circle
//   3. Explosion
//   4. Vortex in a box test
// Also supports a variety of resolutions.
//#####################################################################
#ifndef __SMOKE_STANDARD_TESTS_2D__
#define __SMOKE_STANDARD_TESTS_2D__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
namespace PhysBAM{

template<class TV> class SOLIDS_FLUIDS_EXAMPLE;
template<class TV> class RIGID_BODY_COLLECTION;
template<class T_GRID> class FLUIDS_PARAMETERS_UNIFORM;
template<class T_GRID> class INCOMPRESSIBLE_FLUID_CONTAINER;
template<class T_GRID> class PROJECTION_DYNAMICS_UNIFORM;
template<class T_GRID>
class SMOKE_STANDARD_TESTS_2D
{
    typedef typename T_GRID::SCALAR T;typedef VECTOR<T,2> TV;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
public:
    SOLIDS_FLUIDS_EXAMPLE<TV>& example;
    FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters;
    INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& incompressible_fluid_container;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;

    int test_number;
    T_GRID grid;
    BOX<TV> source;
    MATRIX<T,3> world_to_source;
    VECTOR<T,2> source_velocity;
    T rho;
    T explosion_divergence,explosion_end_time;
    int oriented_box;
    T rotation_angle;
    FRAME<TV> rotation_frame;
    int left_box,right_box;
    T_FACE_ARRAYS_SCALAR beta_face;
    T_FACE_ARRAYS_SCALAR divergence_face_weights;
    
    SMOKE_STANDARD_TESTS_2D(SOLIDS_FLUIDS_EXAMPLE<TV>& example,FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters,INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& incompressible_fluid_container,
        RIGID_BODY_COLLECTION<TV>& rigid_body_collection);
    virtual ~SMOKE_STANDARD_TESTS_2D();

//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
void Set_Dirichlet_Boundary_Conditions(const T time)
{
    //TODO This doesn't do anything.  To remove it, we'd need to change something both in Projects and in Public_Library
}
//#####################################################################
    virtual VECTOR<T,2> Initial_Velocity(const VECTOR<T,2>& X) const;
    void Initialize_Bodies();
    void Initialize(const int test_number_input,const int resolution,const T angle_fraction=0);
    void Get_Divergence(ARRAY<T,VECTOR<int,2> >& divergence,const T dt,const T time);
    void Get_Object_Velocities(PROJECTION_DYNAMICS_UNIFORM<T_GRID>& projection,const T dt,const T time);
//#####################################################################
};
}
#endif

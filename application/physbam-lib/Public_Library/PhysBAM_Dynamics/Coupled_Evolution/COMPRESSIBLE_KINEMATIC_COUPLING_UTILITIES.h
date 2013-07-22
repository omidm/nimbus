//#####################################################################
// Copyright 2009-2011, Mridul Aanjaneya, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES
//#####################################################################
#ifndef __COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES__
#define __COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/BOUNDARY_OBJECT.h>
namespace PhysBAM{

template<class T_GRID> class GRID_BASED_COLLISION_GEOMETRY;
template<class T_GRID> class GRID_BASED_COLLISION_GEOMETRY_UNIFORM;

template<class TV>
class COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES
{
    typedef typename TV::SCALAR T;typedef GRID<TV> T_GRID;typedef VECTOR<T,TV::dimension+2> TV_DIMENSION;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename REBIND<T_ARRAYS_SCALAR,bool>::TYPE T_ARRAYS_BOOL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY_SCALAR;
    typedef typename REBIND_LENGTH<T_BOUNDARY_SCALAR,TV::dimension+2>::TYPE T_BOUNDARY_DIMENSION_SCALAR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_SCALAR T_LINEAR_INTERPOLATION_SCALAR;
    typedef typename T_LINEAR_INTERPOLATION_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_LINEAR_INTERPOLATION_DIMENSION;
    typedef typename LEVELSET_POLICY<T_GRID>::EXTRAPOLATION_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_EXTRAPOLATION_SCALAR_DIMENSION;
    enum {d=TV::dimension+2};

  public:
    T_GRID& grid;
    T_ARRAYS_DIMENSION_SCALAR& U;
    T_ARRAYS_BOOL& psi;
    MPI_UNIFORM_GRID<T_GRID>* mpi_grid;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>* collision_bodies_affecting_fluid;
    BOUNDARY_OBJECT<T_GRID,TV_DIMENSION> *object_boundary;
    T_BOUNDARY_DIMENSION_SCALAR* boundary;
    
    bool use_fast_marching,use_higher_order_solid_extrapolation;
    int number_of_cells_to_extrapolate;
    T_ARRAYS_BOOL outside_fluid;
    T_ARRAYS_SCALAR phi_all_solids_negated;

    COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES(T_GRID& grid_input,T_ARRAYS_DIMENSION_SCALAR& U_input,T_ARRAYS_BOOL& psi_input,
                                              MPI_UNIFORM_GRID<T_GRID>* mpi_grid_input=0);
    ~COMPRESSIBLE_KINEMATIC_COUPLING_UTILITIES();

    void Set_Boundary(T_BOUNDARY_DIMENSION_SCALAR* boundary_input)
    {boundary=boundary_input;}
//#####################################################################
    void Initialize_Solid_Fluid_Coupling(GRID_BASED_COLLISION_GEOMETRY<T_GRID>* collision_bodies_affecting_fluid_input);
    void Update_Cut_Out_Grid();
    void Fill_Solid_Cells(const T dt,const T time);
    void Extrapolate_State_Into_Solids(T_ARRAYS_SCALAR& phi_all_solids_negated,const T dt,const T time,const int number_of_ghost_cells,const int number_of_cells_to_extrapolate);
    void Compute_Phi_Solids(const int number_of_ghost_cells);
    void Get_Neumann_Data(const TV& location,const T max_distance,TV& normal_direction,T& object_velocity_normal_component,TV& reflected_point) const;
};
}
#endif

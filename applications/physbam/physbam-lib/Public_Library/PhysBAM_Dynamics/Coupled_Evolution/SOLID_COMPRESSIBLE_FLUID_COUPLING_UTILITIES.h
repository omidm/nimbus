//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES
//#####################################################################
#ifndef __SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES__
#define __SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES__
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Basic_Geometry/POLYGON.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_COLLIDABLE_UNIFORM_FORWARD.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/CUT_CELL.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
namespace PhysBAM{
template<class TV> class GRID;
template<class T,int d> class VECTOR;
template<class T_GRID> class MPI_UNIFORM_GRID;
template<class T_GRID> class EULER_UNIFORM;
template<class T_GRID> class GRID_BASED_COLLISION_GEOMETRY;
template<class T_GRID> class GRID_BASED_COLLISION_GEOMETRY_UNIFORM;
template<class T_GRID> class EULER_FLUID_FORCES;

template<class TV>
class SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES
{
    typedef typename TV::SCALAR T;
    typedef GRID<TV> T_GRID;
    typedef VECTOR<T,TV::dimension+2> TV_DIMENSION;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_TV;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<CUT_CELLS<T,TV::dimension>*>::TYPE T_ARRAYS_CUT_CELLS;
    typedef typename REBIND<T_ARRAYS_SCALAR,bool>::TYPE T_ARRAYS_BOOL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<int>::TYPE T_FACE_ARRAYS_INT;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_FACE_ARRAYS_DIMENSION_SCALAR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_SCALAR T_LINEAR_INTERPOLATION_SCALAR;
    typedef typename T_LINEAR_INTERPOLATION_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_LINEAR_INTERPOLATION_DIMENSION;
    typedef typename LEVELSET_POLICY<T_GRID>::EXTRAPOLATION_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_EXTRAPOLATION_SCALAR_DIMENSION;
    typedef typename LEVELSET_POLICY<T_GRID>::EXTRAPOLATION_SCALAR T_EXTRAPOLATION_SCALAR;

public:
    EULER_UNIFORM<T_GRID>& euler;
    MPI_UNIFORM_GRID<T_GRID>* mpi_grid;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>* collision_bodies_affecting_fluid;

    T_ARRAYS_BOOL uncovered_cells;
    bool thinshell,use_fast_marching,use_higher_order_solid_extrapolation,fluid_affects_solid;
    int number_of_cells_to_extrapolate;
    TV_DIMENSION solid_state;

    T_ARRAYS_DIMENSION_SCALAR U_n;
    T_FACE_ARRAYS_BOOL solid_fluid_face_time_n;
    T_ARRAYS_BOOL cells_inside_fluid_time_n,outside_fluid;
    EULER_FLUID_FORCES<T_GRID>* euler_fluid_forces;
    T_FACE_ARRAYS_SCALAR pressure_at_faces;
    T_ARRAYS_SCALAR phi_all_solids_negated;

    T_ARRAYS_BOOL near_interface;
    T_ARRAYS_CUT_CELLS cut_cells_n,cut_cells_n_p_half,cut_cells_np1;
    T_ARRAYS_SCALAR cell_volumes_n,cell_volumes_n_p_half,cell_volumes_np1;

    T_ARRAYS_BOOL psi_n,psi_n_p_half,psi_np1,uncovered_cells_n_p_half;

    HASHTABLE<PAIR<TV_INT,int>,ARRAY<PAIR<TV_INT,int> > > visibility_n_p_half;
    HASHTABLE<PAIR<TV_INT,int>,ARRAY<PAIR<TV_INT,int> > > visibility_np1;

    T_FACE_ARRAYS_DIMENSION_SCALAR accumulated_flux;

    SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES(EULER_UNIFORM<T_GRID>& euler_input,MPI_UNIFORM_GRID<T_GRID>* mpi_grid_input=0);
    ~SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES();

//#####################################################################
    void Initialize_Solid_Fluid_Coupling(GRID_BASED_COLLISION_GEOMETRY<T_GRID>* collision_bodies_affecting_fluid_input);
    void Update_Cut_Out_Grid();
    void Fill_Uncovered_Cells();
    void Fill_Solid_Cells(bool fill_pressure_only=false);
    void Project_Fluid_Pressure_At_Neumann_Faces(const T_ARRAYS_SCALAR& p_ghost,T_FACE_ARRAYS_SCALAR& p_face) const;
    void Apply_Isobaric_Fix(const T dt,const T time);
    void Extract_Time_N_Data_For_Explicit_Fluid_Forces();
private:
    void Get_Neumann_Data(const TV& location,const T max_distance,TV& normal_direction,T& object_velocity_normal_component,TV& reflected_point) const;
    void Get_Neumann_Data(const TV& location,const T max_distance,TV& object_velocity,TV& reflected_point) const;
    void Extrapolate_State_Into_Solids(T_ARRAYS_SCALAR& phi_all_solids_negated,const int number_of_ghost_cells,const int number_of_cells_to_extrapolate);
    void Compute_Phi_Solids(const int number_of_ghost_cells);
//#####################################################################
public:
    void Snapshot_State(const T_ARRAYS_DIMENSION_SCALAR& U_ghost);
    void Initialize_Collision_Data();
    void Update_Np1_Collision_Data(const T dt);
    void Compute_Intermediate_Solid_Position_Data(const T dt);
    void Revert_Cells_Near_Interface(const int iteration_number);
    void Update_Cells_Near_Interface(const T dt,const int rk_order,const int rk_substep);
    void Compute_Post_Advected_Variables();
//#####################################################################
};
}
#endif

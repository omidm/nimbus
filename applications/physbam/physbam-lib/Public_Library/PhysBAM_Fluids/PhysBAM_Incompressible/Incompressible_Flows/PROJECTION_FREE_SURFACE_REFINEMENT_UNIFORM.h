//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM  
//#####################################################################
#ifndef __PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM__
#define __PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_REFINEMENT_UNIFORM.h>
#include <PhysBAM_Dynamics/Interpolation/FIRE_INTERPOLATION_POLICY.h>
namespace PhysBAM{

template<class T_GRID>
class PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM:public PROJECTION_REFINEMENT_UNIFORM<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_ARRAYS_SLIP T_FACE_ARRAYS_SLIP_SCALAR;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;typedef typename FIRE_INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP_FIRE_MULTIPHASE T_FACE_LOOKUP_FIRE_MULTIPHASE;
public:
    typedef PROJECTION_REFINEMENT_UNIFORM<T_GRID> BASE;
    using BASE::fine_mpi_grid;using BASE::collidable_solver;using BASE::elliptic_solver;using BASE::p;using BASE::coarse_mpi_grid;using BASE::solid_wall;using BASE::fine_psi_N;using BASE::poisson;
    using BASE::coarse_grid;using BASE::local_grid;using BASE::fast_local_projection;using BASE::fine_grid;using BASE::coarse_scale;using BASE::domain_boundary;using BASE::face_velocities_save;
    
    BOUNDARY_UNIFORM<GRID<TV>,T> *boundary,*phi_boundary;
public:
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T> phi_interpolation;
    PROJECTION_DYNAMICS_UNIFORM<T_GRID> levelset_projection;
    T_LEVELSET &levelset,coarse_levelset;
    T_ARRAYS_SCALAR phi_ghost,coarse_phi,local_phi;
    bool surface_solve;
    int buffer;
public:

    PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM(const T_GRID& mac_grid,T_LEVELSET& levelset_input,const int scale,const T alpha=1,const bool use_surface_solve=true,const bool flame_input=false,const bool multiphase=false,const bool use_variable_beta=false,const bool use_poisson=false);
    virtual ~PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM();

//#####################################################################
    virtual void Initialize_Grid(const T_GRID& mac_grid);
    virtual void Set_Coarse_Boundary_Conditions(T_FACE_ARRAYS_SCALAR& coarse_face_velocities);
    bool Set_Local_Boundary_Conditions(GRID<TV>& local_grid,PROJECTION_UNIFORM<GRID<TV> >& local_projection,TV_INT coarse_index);
    void Set_Local_Phi_From_Fine_Phi(GRID<TV>& local_mac_grid,ARRAY<T,TV_INT>& local_phi,const ARRAY<T,TV_INT>& fine_phi,TV_INT cell_index);
    virtual void Local_Projection_PCG(T_FACE_ARRAYS_SCALAR& fine_face_velocities,T_GRID& local_grid,T_FACE_ARRAYS_SCALAR& local_face_velocities,FAST_PROJECTION_DYNAMICS_UNIFORM<GRID<TV> >& local_projection,const T dt,const T time,TV_INT cell_index);
    bool Contains_Inside(TV_INT cell_index,const ARRAY<T,TV_INT>& levelset_phi,int buffer);
    bool Contains_Outside(TV_INT cell_index,const ARRAY<T,TV_INT>& levelset_phi,int buffer);
    void Set_Coarse_Phi_From_Fine_Phi(ARRAY<T,TV_INT>& coarse_phi,const ARRAY<T,TV_INT>& fine_phi);
    void Set_Levelset_Boundary_Conditions(const GRID<TV>& levelset_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& levelset_velocities,const ARRAY<T,TV_INT>& levelset_phi,const T time);
    void Map_Fine_To_Levelset_For_Constraints(T_FACE_ARRAYS_SCALAR& face_velocities);
    virtual void Map_Fine_To_Coarse(T_FACE_ARRAYS_SCALAR& coarse_face_velocities,const T_FACE_ARRAYS_SCALAR& face_velocities);
    virtual void Map_Coarse_To_Fine(const T_FACE_ARRAYS_SCALAR& coarse_face_velocities,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
//#####################################################################
};
}
#endif

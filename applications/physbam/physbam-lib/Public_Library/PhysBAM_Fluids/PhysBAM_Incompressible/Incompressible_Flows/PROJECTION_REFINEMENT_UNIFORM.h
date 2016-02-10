//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROJECTION_REFINEMENT_UNIFORM  
//#####################################################################
#ifndef __PROJECTION_REFINEMENT_UNIFORM__
#define __PROJECTION_REFINEMENT_UNIFORM__

#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/FAST_PROJECTION_DYNAMICS_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
#include <PhysBAM_Dynamics/Interpolation/FIRE_INTERPOLATION_POLICY.h>
namespace PhysBAM{

template<class T_GRID>
class PROJECTION_REFINEMENT_UNIFORM:public PROJECTION_DYNAMICS_UNIFORM<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_ARRAYS_SLIP T_FACE_ARRAYS_SLIP_SCALAR;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;typedef typename FIRE_INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP_FIRE_MULTIPHASE T_FACE_LOOKUP_FIRE_MULTIPHASE;
public:
    typedef PROJECTION_DYNAMICS_UNIFORM<T_GRID> BASE;
    using BASE::poisson;using BASE::elliptic_solver;using BASE::p;
    
    THREAD_QUEUE* thread_queue;
public:
    MPI_UNIFORM_GRID<T_GRID> *coarse_mpi_grid,*fine_mpi_grid;
    FAST_PROJECTION_DYNAMICS_UNIFORM<T_GRID> fast_local_projection;
    T_FACE_ARRAYS_SCALAR coarse_face_velocities,coarse_face_velocities_save,face_velocities_save,local_face_velocities;
    T_FACE_ARRAYS_SCALAR &beta_face;
    T_FACE_ARRAYS_BOOL fine_psi_N;
    GRID<TV> coarse_grid,fine_grid,local_grid;
    VECTOR<VECTOR<bool,2>,TV::dimension> domain_boundary,solid_wall;
    T scale_face_inverse,alpha;
public:
    int coarse_scale;

    PROJECTION_REFINEMENT_UNIFORM(const T_GRID& mac_grid,const int scale,const T alpha=1,const bool flame_input=false,const bool multiphase=false,const bool use_variable_beta=false,const bool use_poisson=false,THREAD_QUEUE* thread_queue_input=0);
    PROJECTION_REFINEMENT_UNIFORM(const T_GRID& mac_grid,T_LEVELSET& levelset_input,const int scale,const T alpha=1);
    virtual ~PROJECTION_REFINEMENT_UNIFORM();

//#####################################################################
    virtual void Initialize_Grid(const T_GRID& mac_grid);
    void Average_Velocities_From_Fine_To_Coarse(ARRAY<T,FACE_INDEX<TV::dimension> >& coarse_face_velocities,const ARRAY<T,FACE_INDEX<TV::dimension> >& fine_face_velocities);
    void Set_Beta_Face_For_Boundary_Conditions(T_FACE_ARRAYS_SCALAR& coarse_face_velocities);
    virtual void Set_Coarse_Boundary_Conditions(T_FACE_ARRAYS_SCALAR& coarse_face_velocities);
    virtual void Map_Fine_To_Coarse(T_FACE_ARRAYS_SCALAR& coarse_face_velocities,const T_FACE_ARRAYS_SCALAR& face_velocities);
    void Map_Coarse_To_Fine(const T_FACE_ARRAYS_SCALAR& coarse_face_velocities,T_FACE_ARRAYS_SCALAR& face_velocities);
    void Map_Fine_To_Local_Boundary_For_Cell(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,T_FACE_ARRAYS_SCALAR& fine_face_velocities,TV_INT coarse_cell_index);
    void Map_Fine_To_Local_Interior_For_Cell(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,T_FACE_ARRAYS_SCALAR& fine_face_velocities,TV_INT coarse_cell_index,bool zero_out);
    bool Map_Fine_To_Local_Boundaries_For_Cell(GRID<TV>& local_mac_grid,ARRAY<bool,FACE_INDEX<TV::dimension> >& local_psi_N,TV_INT cell_index);
    void Map_Local_To_Fine_Interior_For_Cell(GRID<TV>& local_mac_grid,ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities,T_FACE_ARRAYS_SCALAR& fine_face_velocities,TV_INT cell_index);
    virtual void Local_Projection_PCG(T_FACE_ARRAYS_SCALAR& fine_face_velocities,T_GRID& local_grid,T_FACE_ARRAYS_SCALAR& local_face_velocities,FAST_PROJECTION_DYNAMICS_UNIFORM<GRID<TV> >& local_projection,const T dt,const T time,TV_INT cell_index);
    void Threaded_Local_Projection_PCG(RANGE<TV_INT>& domain,T_FACE_ARRAYS_SCALAR& fine_face_velocities,const T dt,const T time);
    void Local_Projection_PCG(T_FACE_ARRAYS_SCALAR& fine_face_velocities,const T dt,const T time);
    virtual void Map_Coarse_To_Fine(const T_FACE_ARRAYS_SCALAR& coarse_face_velocities,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    virtual void Make_Divergence_Free(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
//#####################################################################
};
}
#endif

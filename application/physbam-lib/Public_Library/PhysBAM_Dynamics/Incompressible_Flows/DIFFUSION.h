//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIFFUSION_UNIFORM
//#####################################################################
#ifndef __DIFFUSION_UNIFORM__
#define __DIFFUSION_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_UNIFORM.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM_FORWARD.h>

namespace PhysBAM{

template<class T_GRID,class T2>
class DIFFUSION_UNIFORM
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;typedef typename REBIND<T_BOUNDARY,T2>::TYPE T_BOUNDARY_T2;
    typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;typedef typename MPI_GRID_POLICY<T_GRID>::MPI_GRID T_MPI_GRID;    
public:
    bool diffuse_weights,diffuse_errors;
    T_MPI_GRID* mpi_grid;
    int num_diffusion_iterations;
    int evenodd,evenodd_cell;
    int max_value,min_value,hard_max,hard_min;

    VECTOR<VECTOR<bool,2>,TV::dimension> solid_walls;
    VECTOR<VECTOR<bool,2>,TV::dimension> mpi_boundary;

//#####################################################################
    DIFFUSION_UNIFORM(T_MPI_GRID* mpi_grid_input);
    ~DIFFUSION_UNIFORM();
    bool Is_MPI_Boundary(const RANGE<TV_INT>& inside_domain,const FACE_INDEX<TV::dimension>& face);
    bool Is_MPI_Boundary(const RANGE<TV_INT>& inside_domain,const TV_INT& index);
    void Cell_Diffusion_Value_Helper(FACE_ITERATOR& iterator,T_ARRAYS_T2& Z,ARRAY<bool,TV_INT>* inside);
    void Cell_Diffusion_Sum_Helper(FACE_ITERATOR& iterator,ARRAY<T,TV_INT>& sum_jc_cell,T_ARRAYS_T2& Z,ARRAY<bool,TV_INT>* inside);
    void Cell_Diffusion_Error_Helper(FACE_ITERATOR& iterator,ARRAY<T,TV_INT>& sum_jc_cell,T_ARRAYS_T2& Z,ARRAY<bool,TV_INT>* inside);
    void Cell_Diffusion_Helper(FACE_ITERATOR& iterator,ARRAY<T,TV_INT>* sum_jc_cell,T_ARRAYS_T2& Z,ARRAY<bool,TV_INT>* inside);
    void Face_Diffusion_Sum_Helper(const GRID<TV>& grid,FACE_INDEX<TV::dimension>& first_face_index,FACE_INDEX<TV::dimension>& second_face_index,ARRAY<T,FACE_INDEX<TV::dimension> >& sum_jc,T_FACE_ARRAYS_SCALAR& Z,ARRAY<bool,FACE_INDEX<TV::dimension> >* inside);
    void Face_Diffusion_Helper(FACE_ITERATOR& iterator,int axis,ARRAY<T,FACE_INDEX<TV::dimension> >* sum_jc,T_FACE_ARRAYS_SCALAR& Z,ARRAY<bool,FACE_INDEX<TV::dimension> >* inside);
    void Face_Diffusion_Helper(CELL_ITERATOR& iterator,int axis,ARRAY<T,FACE_INDEX<TV::dimension> >* sum_jc,T_FACE_ARRAYS_SCALAR& Z,ARRAY<bool,FACE_INDEX<TV::dimension> >* inside);
    void Face_Diffusion(const T_GRID& grid,ARRAY<T,FACE_INDEX<TV::dimension> >* sum_jc,T_FACE_ARRAYS_SCALAR& Z,T_BOUNDARY& boundary,BOUNDARY_UNIFORM<T_GRID,T>* boundary_sum,ARRAY<bool,FACE_INDEX<TV::dimension> >* inside=0);
    void Cell_Diffusion(const T_GRID& grid,T_ARRAYS_T2& Z,T_BOUNDARY_T2& boundary,ARRAY<T,TV_INT>* sum_jc_cell=0,BOUNDARY_UNIFORM<T_GRID,T>* boundary_sum=0,ARRAY<bool,TV_INT>* inside=0);
//#####################################################################
};
}
#endif

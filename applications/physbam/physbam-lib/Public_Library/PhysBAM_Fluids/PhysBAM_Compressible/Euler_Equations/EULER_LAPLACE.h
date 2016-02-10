//#####################################################################
// Copyright 2007, Jon Gretarsson, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_LAPLACE
//#####################################################################
#ifndef __EULER_LAPLACE__
#define __EULER_LAPLACE__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UNIFORM_FORWARD.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
namespace PhysBAM{

template<class T_GRID> struct GRID_ARRAYS_POLICY;

template<class T_LAPLACE>
class EULER_LAPLACE:public T_LAPLACE
{
    typedef typename T_LAPLACE::GRID_T T_GRID;typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename REBIND<T_ARRAYS_SCALAR,int>::TYPE T_ARRAYS_INT;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::INDEX INDEX;
public:
    typedef T_LAPLACE BASE;
    using BASE::filled_region_touches_dirichlet;using BASE::Solve;using BASE::Find_Solution_Regions;using BASE::f;
    using BASE::grid;using BASE::solve_neumann_regions;using BASE::dt;using BASE::dt_is_set;using BASE::Initialize_Grid;

    const T_ARRAYS_SCALAR& one_over_rho_c_squared;
    T_ARRAYS_SCALAR* phi;

    EULER_LAPLACE(const T_GRID& grid_input,T_ARRAYS_SCALAR& u_input,const T_ARRAYS_SCALAR& one_over_rho_c_squared_input,T_ARRAYS_SCALAR* phi_input=0)
        :T_LAPLACE(grid_input,u_input,true,false,false),one_over_rho_c_squared(one_over_rho_c_squared_input),phi(phi_input)
    {}

    EULER_LAPLACE(const T_GRID& grid_input,T_ARRAYS_SCALAR& u_input,const T_ARRAYS_SCALAR& one_over_rho_c_squared_input,const bool flag)
        :T_LAPLACE(grid_input,u_input,flag,flag,false),one_over_rho_c_squared(one_over_rho_c_squared_input)
    {}

    void Solve(const T time,const bool solution_regions_already_computed=false)
    {if(!solution_regions_already_computed) Find_Solution_Regions();
    T_LAPLACE::Solve(time,true);}

    void Find_A(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,const ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,
            T_ARRAYS_INT& cell_index_to_matrix_index)
    {assert(dt_is_set);dt_is_set=false;

    T_LAPLACE::Find_A(domain,A_array,b_array,filled_region_cell_count,cell_index_to_matrix_index);
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){INDEX cell_index=iterator.Cell_Index();
        int color=filled_region_colors(cell_index);
        if(color!=-1 && (filled_region_touches_dirichlet(color)||solve_neumann_regions)){
            int matrix_index=cell_index_to_matrix_index(cell_index);
            SPARSE_MATRIX_FLAT_NXN<T>& A=A_array(filled_region_colors(cell_index));
            if(!phi || (phi && (*phi)(cell_index)>0)) A(matrix_index,matrix_index)-=one_over_rho_c_squared(cell_index)/(dt*dt);}}
    }

//#####################################################################
};
}
#endif

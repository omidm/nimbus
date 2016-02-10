//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/FAST_PROJECTION_DYNAMICS_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> FAST_PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
FAST_PROJECTION_DYNAMICS_UNIFORM(const int scale,const bool flame_input,const bool multiphase,const bool use_variable_beta,const bool use_poisson)
    :PROJECTION_DYNAMICS_UNIFORM<T_GRID>(GRID<TV>(TV_INT::All_Ones_Vector()*scale,RANGE<TV>(TV(),TV::All_Ones_Vector()),true),flame_input,multiphase,use_variable_beta,use_poisson)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> FAST_PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
FAST_PROJECTION_DYNAMICS_UNIFORM(const int scale,T_LEVELSET& levelset_input)
    :PROJECTION_DYNAMICS_UNIFORM<T_GRID>(GRID<TV>(TV_INT::All_Ones_Vector()*scale,RANGE<TV>(TV(),TV::All_Ones_Vector()),true),levelset_input)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> FAST_PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
FAST_PROJECTION_DYNAMICS_UNIFORM(FAST_PROJECTION_DYNAMICS_UNIFORM<T_GRID>& projection)
    :PROJECTION_DYNAMICS_UNIFORM<T_GRID>(projection.p_grid,projection.flame,projection.poisson!=NULL?projection.poisson->multiphase:false,projection.poisson!=NULL?projection.poisson->use_variable_beta:false,projection.poisson!=NULL)
{
    A=projection.A;
    b=projection.b;
    cell_index_to_matrix_index=projection.cell_index_to_matrix_index;
    matrix_index_to_cell_index=projection.matrix_index_to_cell_index;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> FAST_PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
~FAST_PROJECTION_DYNAMICS_UNIFORM()
{
}
//#####################################################################
// Function Initialize_Grid
//#####################################################################
template<class T_GRID> void FAST_PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Initialize_Grid(const T_GRID& mac_grid)
{
    BASE::Initialize_Grid(mac_grid);
    int scale=mac_grid.counts.x;
    elliptic_solver->Set_Neumann_Outer_Boundaries();
    int number_of_elements=TV::dimension==2?scale*scale:scale*scale*scale;
    cell_index_to_matrix_index.Resize(mac_grid.Domain_Indices());
    matrix_index_to_cell_index.Resize(number_of_elements);
    b.Resize(number_of_elements);
    ARRAY<int> row_counts;
    row_counts.Resize(number_of_elements);
    int count=0;
    for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        int matrix_index;
        count++;cell_index_to_matrix_index(cell_index)=matrix_index=count;
        matrix_index_to_cell_index(matrix_index)=cell_index;}
    for(int i=1;i<=row_counts.m;i++){
        int boundary=0;
        for(int j=1;j<=TV::dimension;j++) if(matrix_index_to_cell_index(i)(j)==1 || matrix_index_to_cell_index(i)(j)==mac_grid.Counts()(j)) boundary++;
        row_counts(i)=(2*TV::dimension+1)-boundary;}
    A.Set_Row_Lengths(row_counts);
    TV one_over_dx2=Inverse(mac_grid.dX*mac_grid.dX);
    T default_row_sum=-2*one_over_dx2.L1_Norm();
    TV_INT grid_counts=mac_grid.counts;
    for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        T row_sum=default_row_sum;
        int matrix_index=cell_index_to_matrix_index(cell_index);
        for(int axis=1;axis<=GRID<TV>::dimension;axis++){TV_INT offset;offset[axis]=1;
            if(elliptic_solver->psi_N.Component(axis)(cell_index)) row_sum+=one_over_dx2[axis];
            else A.Set_Element(matrix_index,cell_index_to_matrix_index(cell_index-offset),one_over_dx2[axis]);
            if(elliptic_solver->psi_N.Component(axis)(cell_index+offset)) row_sum+=one_over_dx2[axis];
            else A.Set_Element(matrix_index,cell_index_to_matrix_index(cell_index+offset),one_over_dx2[axis]);}
        A.Set_Element(matrix_index,matrix_index,row_sum);}
}
//#####################################################################
// Function Make_Divergence_Free
//#####################################################################
template<class T_GRID> void FAST_PROJECTION_DYNAMICS_UNIFORM<T_GRID>::
Make_Divergence_Free_Fast(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    Compute_Divergence(typename INTERPOLATION_POLICY<GRID<TV> >::FACE_LOOKUP(face_velocities),elliptic_solver);
    for(typename GRID<TV>::CELL_ITERATOR iterator(p_grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        int matrix_index=cell_index_to_matrix_index(cell_index);
        b(matrix_index)=elliptic_solver->f(cell_index);}
    elliptic_solver->Solve_Subregion(matrix_index_to_cell_index,A,b);
    Apply_Pressure(face_velocities,dt,time);
    A.Negate();
}
//#####################################################################
template class FAST_PROJECTION_DYNAMICS_UNIFORM<GRID<VECTOR<float,1> > >;
template class FAST_PROJECTION_DYNAMICS_UNIFORM<GRID<VECTOR<float,2> > >;
template class FAST_PROJECTION_DYNAMICS_UNIFORM<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FAST_PROJECTION_DYNAMICS_UNIFORM<GRID<VECTOR<double,1> > >;
template class FAST_PROJECTION_DYNAMICS_UNIFORM<GRID<VECTOR<double,2> > >;
template class FAST_PROJECTION_DYNAMICS_UNIFORM<GRID<VECTOR<double,3> > >;
#endif

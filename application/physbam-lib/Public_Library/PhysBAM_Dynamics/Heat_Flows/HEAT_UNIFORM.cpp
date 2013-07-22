//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h> // TODO(jontg): There doesn't appear to be any good reason for these dependencies to exist...
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Dynamics/Heat_Flows/HEAT_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> HEAT_UNIFORM<T_GRID>::
HEAT_UNIFORM(T_GRID& grid_input,T_ARRAYS_SCALAR& Q_input)
    :grid(grid_input),Q(Q_input),laplace(grid,Q)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> HEAT_UNIFORM<T_GRID>::
~HEAT_UNIFORM()
{}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T_GRID> void HEAT_UNIFORM<T_GRID>::
Euler_Step(const T dt,const T time)
{
    TV kappa_over_DX2=kappa*sqr(grid.one_over_dX);
    laplace.f.Fill(0);
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();if(!laplace.psi_D(cell))
        for(int axis=1;axis<=T_GRID::dimension;axis++){
            TV_INT face1=iterator.First_Face_Index(axis),face2=iterator.Second_Face_Index(axis),offset=TV_INT::Axis_Vector(axis);
            if(!laplace.psi_N.Component(axis)(face1)) laplace.f(cell)-=(Q(cell)-Q(cell-offset))*kappa_over_DX2[axis];
            if(!laplace.psi_N.Component(axis)(face2)) laplace.f(cell)+=(Q(cell+offset)-Q(cell))*kappa_over_DX2[axis];}}
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();Q(cell)+=dt*laplace.f(cell);}
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR HEAT_UNIFORM<T_GRID>::
CFL()
{
    T dt_overall=kappa*2*grid.one_over_dX.Magnitude_Squared();
    return 1/max(dt_overall,1/max_time_step);
}
//#####################################################################
// Function Backward_Euler_Calculate_Right_Hand_Side
//#####################################################################
// boundary conditions at time n
template<class T_GRID> void HEAT_UNIFORM<T_GRID>::
Backward_Euler_Calculate_Right_Hand_Side(const T dt,const T time)
{
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();laplace.f(cell)=Q(cell);}
}
//#####################################################################
// Function Backward_Euler_Step
//#####################################################################
// boundary conditions at time n+1
template<class T_GRID> void HEAT_UNIFORM<T_GRID>::
Backward_Euler_Step(const T dt,const T time)
{
    Implicit_Solve(dt*kappa,dt,time);
}
//#####################################################################
// Function Crank_Nicolson_Calculate_Right_Hand_Side
//#####################################################################
// boundary conditions at time n
template<class T_GRID> void HEAT_UNIFORM<T_GRID>::
Crank_Nicolson_Calculate_Right_Hand_Side(const T dt,const T time)
{
    TV half_dt_kappa_over_DX2=(T).5*dt*kappa*sqr(grid.one_over_dX);
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        laplace.f(cell)=Q(cell);
        if(!laplace.psi_D(cell)){
            for(int axis=1;axis<=T_GRID::dimension;axis++){
                TV_INT face1=iterator.First_Face_Index(axis),face2=iterator.Second_Face_Index(axis),offset=TV_INT::Axis_Vector(axis);
                if(!laplace.psi_N.Component(axis)(face1)) laplace.f(cell)-=(Q(cell)-Q(cell-offset))*half_dt_kappa_over_DX2[axis];
                if(!laplace.psi_N.Component(axis)(face2)) laplace.f(cell)+=(Q(cell+offset)-Q(cell))*half_dt_kappa_over_DX2[axis];}}}

    if(laplace.second_order_cut_cell_method)
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face=iterator.Face_Index(),cell1=iterator.First_Cell_Index(),cell2=iterator.Second_Cell_Index();
            if(!laplace.psi_N.Component(axis)(face) && LEVELSET_UTILITIES<T>::Interface(laplace.levelset->phi(cell1),laplace.levelset->phi(cell2))){
                T theta=LEVELSET_UTILITIES<T>::Theta(laplace.levelset->phi(cell1),laplace.levelset->phi(cell2));
                if(laplace.psi_D(cell1) && !laplace.psi_D(cell2)){
                    laplace.f(cell2)+=half_dt_kappa_over_DX2[axis]*(Q(cell2)-Q(cell1));laplace.f(cell2)-=half_dt_kappa_over_DX2[axis]*(Q(cell2)-laplace.u_interface.Component(axis)(face))/(1-theta);}
                else if(laplace.psi_D(cell2) && !laplace.psi_D(cell1)){
                    laplace.f(cell1)-=half_dt_kappa_over_DX2[axis]*(Q(cell2)-Q(cell1));laplace.f(cell1)+=half_dt_kappa_over_DX2[axis]*(laplace.u_interface.Component(axis)(face)-Q(cell1))/theta;}}}
}
//#####################################################################
// Function Crank_Nicolson_Step
//#####################################################################
// boundary conditions at time n+1
template<class T_GRID> void HEAT_UNIFORM<T_GRID>::
Crank_Nicolson_Step(const T dt,const T time)
{
    Implicit_Solve((T).5*dt*kappa,dt,time);
}
//#####################################################################
// Function Implicit_Solve
//#####################################################################
// duplicate of laplace.Solve() with some changes
template<class T_GRID> void HEAT_UNIFORM<T_GRID>::
Implicit_Solve(const T coefficient,const T dt,const T time)
{
    laplace.coefficient=coefficient;
    laplace.Solve(time);
}
//#####################################################################
template class HEAT_UNIFORM<GRID<VECTOR<float,1> > >;
template class HEAT_UNIFORM<GRID<VECTOR<float,2> > >;
template class HEAT_UNIFORM<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class HEAT_UNIFORM<GRID<VECTOR<double,1> > >;
template class HEAT_UNIFORM<GRID<VECTOR<double,2> > >;
template class HEAT_UNIFORM<GRID<VECTOR<double,3> > >;
#endif

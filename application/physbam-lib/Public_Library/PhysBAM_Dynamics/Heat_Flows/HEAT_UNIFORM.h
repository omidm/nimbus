//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HEAT_UNIFORM
//#####################################################################
//
// Solves T_t = k laplacian T.
// Input T with 1 ghost cell with values equal to the initial guess at Dirchlet boundary condition points.
// Set psi_N=true for each point that gets Neumann boundary conditions. Set psi_D=true for each point that gets Dirichlet boundary conditions.
//
// Euler_Step() - set boundary conditions at time n before calling this.
// Backward_Euler_Step() - set boundary conditions at time n, call Backward_Euler_Calculate_Right_Hand_Side(), then set boundary conditions at time n+1 before calling this.
// Crank_Nicolson() - set boundary conditions at time n, call Crank_Nicolson_Calculate_Right_Hand_Side(), then set boundary conditions at time n+1 before calling this.
//
//#####################################################################
#ifndef __HEAT_UNIFORM__
#define __HEAT_UNIFORM__

#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/LAPLACE_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Heat_Flows/HEAT.h>
#include <PhysBAM_Dynamics/Heat_Flows/HEAT_LAPLACE.h>
namespace PhysBAM{

template<class T_GRID>
class HEAT_UNIFORM:public HEAT<typename T_GRID::SCALAR>
{
private:
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;

    using HEAT<T>::max_time_step;using HEAT<T>::density;using HEAT<T>::specific_heat;using HEAT<T>::thermal_conductivity;using HEAT<T>::kappa;
public:

    T_GRID& grid;
    T_ARRAYS_SCALAR& Q;
    HEAT_LAPLACE<LAPLACE_COLLIDABLE_UNIFORM<T_GRID> > laplace;

    HEAT_UNIFORM(T_GRID& grid_input,T_ARRAYS_SCALAR& Q_input);
    ~HEAT_UNIFORM();

//#####################################################################
    // explicit update with CFL=1/2 and all Neumann outer boundaries, no smoothing where phi is negative if defined
    void Euler_Step(const T dt,const T time=0);
    T CFL();
    void Backward_Euler_Step(const T dt,const T time=0);
    void Backward_Euler_Calculate_Right_Hand_Side(const T dt,const T time=0);
    void Crank_Nicolson_Step(const T dt,const T time=0);
    void Crank_Nicolson_Calculate_Right_Hand_Side(const T dt,const T time=0);
private:
    void Implicit_Solve(const T coefficient,const T dt,const T time=0);
//#####################################################################
};
}
#endif

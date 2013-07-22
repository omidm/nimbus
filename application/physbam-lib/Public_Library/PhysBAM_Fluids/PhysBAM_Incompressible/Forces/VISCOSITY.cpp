//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/VISCOSITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/IMPLICIT_VISCOSITY_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> VISCOSITY<T_GRID>::
VISCOSITY(LAPLACE_UNIFORM<T_GRID>& elliptic_solver_input,const T_ARRAYS_SCALAR& variable_viscosity_input,const T density_input,const T viscosity_input,bool implicit_viscosity_input,
    bool use_explicit_part_of_implicit_viscosity_input,bool use_variable_viscosity_input,int maximum_implicit_viscosity_iterations_input,bool use_psi_R_input)
    :elliptic_solver(elliptic_solver_input),density(density_input),viscosity(viscosity_input),implicit_viscosity(implicit_viscosity_input),
    use_explicit_part_of_implicit_viscosity(use_explicit_part_of_implicit_viscosity_input),variable_viscosity(variable_viscosity_input),use_variable_viscosity(use_variable_viscosity_input),
    maximum_implicit_viscosity_iterations(maximum_implicit_viscosity_iterations_input),use_psi_R(use_psi_R_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> VISCOSITY<T_GRID>::
~VISCOSITY()
{
}
//#####################################################################
// Function Add_Explicit_Forces
//#####################################################################
template<class T_GRID> void VISCOSITY<T_GRID>::
Add_Explicit_Forces(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    if(dt && use_variable_viscosity && (!implicit_viscosity || use_explicit_part_of_implicit_viscosity)){
        if(!implicit_viscosity) PHYSBAM_NOT_IMPLEMENTED();
        IMPLICIT_VISCOSITY_UNIFORM<T_GRID>::Variable_Viscosity_Explicit_Part(density,variable_viscosity,grid,face_velocities,face_velocities_ghost,dt,time);}
}
//#####################################################################
// Function Add_Implicit_Forces_Before_Projection
//#####################################################################
template<class T_GRID> void VISCOSITY<T_GRID>::
Add_Implicit_Forces_Before_Projection(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    if(!dt || (!use_variable_viscosity && viscosity==0)) return;
    for(int axis=1;axis<=T_GRID::dimension;axis++){
        IMPLICIT_VISCOSITY_UNIFORM<T_GRID> implicit_viscosity(elliptic_solver,variable_viscosity,density,viscosity,elliptic_solver.mpi_grid,axis,use_variable_viscosity,use_psi_R);
        implicit_viscosity.Viscous_Update(grid,face_velocities,face_velocities_ghost,dt,time,maximum_implicit_viscosity_iterations);}
    if(elliptic_solver.mpi_grid) elliptic_solver.mpi_grid->Copy_Common_Face_Data(face_velocities);
}
//#####################################################################
// Function Initialize_Grids
//#####################################################################
template<class T_GRID> void VISCOSITY<T_GRID>::
Initialize_Grids(const T_GRID& grid)
{
    elliptic_solver.Initialize_Grid(grid);
}
//#####################################################################
template class VISCOSITY<GRID<VECTOR<float,1> > >;
template class VISCOSITY<GRID<VECTOR<float,2> > >;
template class VISCOSITY<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VISCOSITY<GRID<VECTOR<double,1> > >;
template class VISCOSITY<GRID<VECTOR<double,2> > >;
template class VISCOSITY<GRID<VECTOR<double,3> > >;
#endif

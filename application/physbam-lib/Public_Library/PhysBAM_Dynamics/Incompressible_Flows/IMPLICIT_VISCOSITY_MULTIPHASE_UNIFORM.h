//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM
//#####################################################################
#ifndef __IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM__
#define __IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/IMPLICIT_VISCOSITY_UNIFORM.h>
#include <string>
namespace PhysBAM{

template<class T_GRID> class POISSON_COLLIDABLE_UNIFORM;
template<class T_LAPLACE> class HEAT_LAPLACE;
template<class T_GRID> class PROJECTION_DYNAMICS_UNIFORM;

template<class T_GRID>
class IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM:public IMPLICIT_VISCOSITY_UNIFORM<T_GRID>
{
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::AVERAGING T_AVERAGING;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET_MULTIPLE T_LEVELSET_MULTIPLE;
    typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    typedef typename MPI_GRID_POLICY<T_GRID>::MPI_GRID T_MPI_GRID;
public:
    typedef IMPLICIT_VISCOSITY_UNIFORM<T_GRID> BASE;
    using BASE::face_grid;using BASE::heat_solver;using BASE::u;using BASE::axis;
    using BASE::mpi_grid;using BASE::use_variable_viscosity;

    PROJECTION_DYNAMICS_UNIFORM<T_GRID>& projection;
    ARRAY<T> densities;
    ARRAY<T> viscosities;

    IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM(PROJECTION_DYNAMICS_UNIFORM<T_GRID>& projection_input,const T_ARRAYS_SCALAR& variable_viscosity_input,const ARRAY<T>& densities_input,const ARRAY<T>& viscosities_input,T_MPI_GRID* mpi_grid_input,const int axis_input,bool use_variable_viscosity_input);
    virtual ~IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM();

//#####################################################################
private:
    void Allocate_Heat_Solver() PHYSBAM_OVERRIDE;
    void Setup_Viscosity(const T dt) PHYSBAM_OVERRIDE;
    void Setup_Boundary_Conditions(const T_FACE_ARRAYS_SCALAR& face_velocities) PHYSBAM_OVERRIDE;
    void Calculate_Velocity_Jump();
    void Debug_Write(const std::string& output_directory_input);
//#####################################################################
};
}
#endif

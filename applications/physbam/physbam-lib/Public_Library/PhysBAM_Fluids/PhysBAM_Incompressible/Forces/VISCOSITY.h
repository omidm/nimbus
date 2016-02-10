//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VISCOSITY
//#####################################################################
#ifndef __VISCOSITY__
#define __VISCOSITY__

#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/INCOMPRESSIBLE_FLUIDS_FORCES.h>
namespace PhysBAM{
template<class T_GRID> class LAPLACE_UNIFORM;

template<class T_GRID>
class VISCOSITY:public INCOMPRESSIBLE_FLUIDS_FORCES<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;
    typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    LAPLACE_UNIFORM<T_GRID>& elliptic_solver;
public:
    T density;
    T viscosity;
    bool implicit_viscosity;
    bool use_explicit_part_of_implicit_viscosity;
    const T_ARRAYS_SCALAR& variable_viscosity;
    bool use_variable_viscosity;
    int maximum_implicit_viscosity_iterations;
    bool use_psi_R;

    VISCOSITY(LAPLACE_UNIFORM<T_GRID>& elliptic_solver_input,const T_ARRAYS_SCALAR& variable_viscosity_input,const T density_input,const T viscosity_input,bool implicit_viscosity_input,
        bool use_explicit_part_of_implicit_viscosity_input,bool use_variable_viscosity_input,int maximum_implicit_viscosity_iterations_input,bool use_psi_R_input);
    virtual ~VISCOSITY();

//#####################################################################
    void Add_Explicit_Forces(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time) PHYSBAM_OVERRIDE;
    void Add_Implicit_Forces_Before_Projection(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Add_Implicit_Forces_Projection(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Initialize_Grids(const T_GRID& grid) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

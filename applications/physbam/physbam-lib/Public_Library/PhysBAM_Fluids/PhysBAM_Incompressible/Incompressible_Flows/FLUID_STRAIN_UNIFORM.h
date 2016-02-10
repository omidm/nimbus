//#####################################################################
// Copyright 2002-2006, Doug Enright, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andy Lutomirksi, Jonathan Su, Paul-James White.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_STRAIN_UNIFORM  
//#####################################################################
#ifndef __FLUID_STRAIN_UNIFORM__
#define __FLUID_STRAIN_UNIFORM__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Boundaries/BOUNDARY_FORWARD.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/FLUID_STRAIN.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_POLICY.h>
namespace PhysBAM{

template<class T> class EXTERNAL_STRAIN_ADJUSTMENT;
template<class T_GRID> class LEVELSET_MULTIPLE_UNIFORM;

template<class T_GRID>
class FLUID_STRAIN_UNIFORM:public FLUID_STRAIN<typename T_GRID::SCALAR>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX SYMMETRIC_MATRIX;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<SYMMETRIC_MATRIX>::TYPE T_ARRAYS_SYMMETRIC_MATRIX;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;typedef typename ADVECTION_POLICY<T_GRID>::ADVECTION_SEMI_LAGRANGIAN_SCALAR T_ADVECTION_SEMI_LAGRANGIAN_SCALAR;
    typedef typename REBIND<T_ADVECTION_SEMI_LAGRANGIAN_SCALAR,SYMMETRIC_MATRIX>::TYPE T_ADVECTION_SEMI_LAGRANGIAN_SYMMETRIC_MATRIX;
    typedef typename INCOMPRESSIBLE_POLICY<T_GRID>::PROJECTION T_PROJECTION;
public:
    using FLUID_STRAIN<T>::viscosity_index;using FLUID_STRAIN<T>::strainrate_time;using FLUID_STRAIN<T>::elastic_modulus;
    using FLUID_STRAIN<T>::plasticity_alpha;using FLUID_STRAIN<T>::plasticity_gamma;

    T_GRID grid;
    T_ARRAYS_SYMMETRIC_MATRIX e; // strain tensor
    BOUNDARY_UNIFORM<T_GRID,SYMMETRIC_MATRIX>* e_boundary;
    ADVECTION<T_GRID,SYMMETRIC_MATRIX>* e_advection;
    EXTERNAL_STRAIN_ADJUSTMENT<T>* external_strain_adjustment;
private:               
    BOUNDARY_UNIFORM<T_GRID,SYMMETRIC_MATRIX>& e_boundary_default;
    T_ADVECTION_SEMI_LAGRANGIAN_SYMMETRIC_MATRIX& e_advection_default;
    mutable bool cfl_called;
public:

    FLUID_STRAIN_UNIFORM(const T_GRID& grid_input);
    ~FLUID_STRAIN_UNIFORM();

    void Initialize_Grid(const T_GRID& grid_input)
    {assert(grid_input.Is_MAC_Grid());grid=grid_input;e.Resize(grid.Domain_Indices());}

    void Set_Custom_Boundary(BOUNDARY_UNIFORM<T_GRID,SYMMETRIC_MATRIX>& e_boundary_input)
    {e_boundary=&e_boundary_input;}

    void Set_Custom_Advection(ADVECTION<T_GRID,SYMMETRIC_MATRIX>& e_advection_input)
    {e_advection=&e_advection_input;}

    void Set_External_Strain_Adjustment(EXTERNAL_STRAIN_ADJUSTMENT<T>& external_strain_adjustment_input)
    {external_strain_adjustment=&external_strain_adjustment_input;}

//#####################################################################
    void Update_Strain_Equation(const T dt,const T time,const T density,T_FACE_ARRAYS_SCALAR& face_velocities,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,
        const T_ARRAYS_SCALAR& phi_ghost,const int number_of_ghost_cells);
    void Update_Strain_Equation_Multiphase(const T dt,const T time,const T density,T_FACE_ARRAYS_SCALAR& face_velocities,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,
        const LEVELSET_MULTIPLE_UNIFORM<T_GRID>& levelset,const int region,const int number_of_ghost_cells);
private:
    void Update_Strain_Equation_Helper_Cell_Centered(const T dt,const T time,const T density,const T heaviside_bandwidth,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,
        T_ARRAYS_VECTOR& V,const T_ARRAYS_SCALAR& phi_ghost,const int number_of_ghost_cells);
public:
    void Extrapolate_Strain_Across_Interface(T_ARRAYS_SCALAR& phi_ghost,const T band_width=3);
    T CFL(const T density) const;
//#####################################################################
};

// not implemented in one dimension
template<class T>
class FLUID_STRAIN_UNIFORM<GRID<VECTOR<T,1> > >:public FLUID_STRAIN<T>
{
    typedef VECTOR<T,1> TV;
public:
//#####################################################################
    typedef typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX SYMMETRIC_MATRIX;typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<SYMMETRIC_MATRIX>::TYPE T_ARRAYS_SYMMETRIC_MATRIX;

    T_ARRAYS_SYMMETRIC_MATRIX e; // strain tensor
    BOUNDARY_UNIFORM<GRID<TV>,SYMMETRIC_MATRIX>* e_boundary;

    FLUID_STRAIN_UNIFORM(const GRID<TV>& grid_input){PHYSBAM_NOT_IMPLEMENTED();}
    void Initialize_Grid(const GRID<TV>& grid_input){PHYSBAM_NOT_IMPLEMENTED();}
    void Set_Custom_Boundary(BOUNDARY_UNIFORM<GRID<TV>,SYMMETRIC_MATRIX>& e_boundary_input){e_boundary=&e_boundary_input;}
    void Update_Strain_Equation(const T dt,const T time,const T density,T_FACE_ARRAYS_SCALAR& face_velocities,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,
        const ARRAY<T,VECTOR<int,1> >& phi_ghost,const int number_of_ghost_cells){PHYSBAM_NOT_IMPLEMENTED();}
    void Update_Strain_Equation_Multiphase(const T dt,const T time,const T density,T_FACE_ARRAYS_SCALAR& face_velocities,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,
        const LEVELSET_MULTIPLE_UNIFORM<GRID<TV> >& levelset,const int region,const int number_of_ghost_cells){PHYSBAM_NOT_IMPLEMENTED();}
    void Extrapolate_Strain_Across_Interface(ARRAY<T,VECTOR<int,1> >& phi_ghost,const T band_width=3){PHYSBAM_NOT_IMPLEMENTED();}
    T CFL(const T density) const{PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
};

}
#endif

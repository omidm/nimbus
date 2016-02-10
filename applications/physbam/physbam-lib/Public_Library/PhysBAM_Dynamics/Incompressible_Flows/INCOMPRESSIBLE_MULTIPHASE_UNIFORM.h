//#####################################################################
// Copyright 2005-2007, Ron Fedkiw, Frank Losasso, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_MULTIPHASE_UNIFORM
//#####################################################################
#ifndef __INCOMPRESSIBLE_MULTIPHASE_UNIFORM__
#define __INCOMPRESSIBLE_MULTIPHASE_UNIFORM__

#include <PhysBAM_Fluids/PhysBAM_Incompressible/Boundaries/BOUNDARY_MAC_GRID_SOLID_WALL_SLIP.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/FLUID_STRAIN_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID>
class INCOMPRESSIBLE_MULTIPHASE_UNIFORM:public INCOMPRESSIBLE_UNIFORM<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;typedef typename T_GRID::VECTOR_INT T_VECTOR_INT;
    typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;typedef typename LEVELSET_POLICY<T_GRID>::EXTRAPOLATION_SCALAR T_EXTRAPOLATION_SCALAR;typedef typename INTERPOLATION_POLICY<T_GRID>::INTERPOLATION_SCALAR T_INTERPOLATION_SCALAR;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;typedef typename LEVELSET_POLICY<T_GRID>::FAST_LEVELSET_T T_FAST_LEVELSET;
    typedef typename T_ARRAYS_SCALAR::template REBIND<typename TV::SPIN>::TYPE T_ARRAYS_SPIN;
public:
    typedef INCOMPRESSIBLE_UNIFORM<T_GRID> BASE;
    using BASE::use_force;using BASE::viscosity;using BASE::use_variable_viscosity;using BASE::use_variable_vorticity_confinement;using BASE::dt_old;using BASE::gravity;
    using BASE::downward_direction;using BASE::vorticity_confinements;using BASE::nonzero_viscosity;using BASE::nonzero_surface_tension;using BASE::mpi_grid;
    using BASE::use_explicit_part_of_implicit_viscosity;using BASE::vorticity_confinement;using BASE::max_time_step;using BASE::advection;
    using BASE::Set_Custom_Advection;using BASE::GFM;using BASE::number_of_interface_cells;using BASE::viscosities;using BASE::surface_tensions;
    using BASE::projection;using BASE::grid;using BASE::boundary;using BASE::force;using BASE::variable_vorticity_confinement;using BASE::strain;using BASE::variable_viscosity;
    using BASE::maximum_implicit_viscosity_iterations;

    T_FACE_ARRAYS_SCALAR viscous_force;
    T_LEVELSET* levelset_for_dirichlet_regions;
    ARRAY<FLUID_STRAIN_UNIFORM<T_GRID>*> strains;

    INCOMPRESSIBLE_MULTIPHASE_UNIFORM(const T_GRID& grid_input,PROJECTION_DYNAMICS_UNIFORM<T_GRID>& projection_input);
    virtual ~INCOMPRESSIBLE_MULTIPHASE_UNIFORM();

    void Use_Strain(const ARRAY<bool>& use_multiphase_strain)
    {for(int i=1;i<=strains.m;i++)delete strains(i);
    strains.Resize(use_multiphase_strain.m);
    for(int i=1;i<=use_multiphase_strain.m;i++)if(use_multiphase_strain(i))strains(i)=new FLUID_STRAIN_UNIFORM<T_GRID>(grid);}

    // overrides version from BASE
    void Advance_One_Time_Step_Forces(const T dt,const T time,const bool implicit_viscosity=false,const T_ARRAYS_SCALAR* phi_ghost=0)
    {PHYSBAM_NOT_IMPLEMENTED();/*PHYSBAM_ASSERT(!phi_ghost);Advance_One_Time_Step_Forces(dt,time,implicit_viscosity,0,0);*/}

    // overrides version from BASE
    void Advance_One_Time_Step_Convection(const T dt,const T time,T_FACE_ARRAYS_SCALAR& face_velocities_to_advect)
    {PHYSBAM_NOT_IMPLEMENTED();/*Advance_One_Time_Step_Convection(dt,time,face_velocities_to_advect,0);*/}
    
//#####################################################################
    void Advance_One_Time_Step_Forces(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,const bool implicit_viscosity,const ARRAY<T_ARRAYS_SCALAR>* phi_ghost,
        const ARRAY<bool>* pseudo_dirichlet_regions,const int number_of_ghost_cells);
    void Advance_One_Time_Step_Convection(const T dt,const T time,T_FACE_ARRAYS_SCALAR& advecting_face_velocities,T_FACE_ARRAYS_SCALAR& face_velocities_to_advect,const ARRAY<bool>* pseudo_dirichlet_regions,const int number_of_ghost_cells);
    void Advance_One_Time_Step_Implicit_Part(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,const bool implicit_viscosity=false,BOUNDARY_UNIFORM<T_GRID,T>* projection_boundary=0,
        bool use_levelset_viscosity=false,BOUNDARY_CONDITIONS_CALLBACKS<TV>* bc_callbacks=0,bool print_viscosity_matrix=false) PHYSBAM_OVERRIDE;
    void Calculate_Pressure_Jump(const T dt,const T time);
    T CFL(T_FACE_ARRAYS_SCALAR& face_velocities,const bool inviscid=false,const bool viscous_only=false) const;
    void Set_Dirichlet_Boundary_Conditions(ARRAY<T_ARRAYS_SCALAR>& phis,const ARRAY<bool>& dirichlet_regions,const ARRAY<T>* pressures=0);
    void Add_Surface_Tension(T_LEVELSET& levelset,const T time);
    void Compute_Vorticity_Confinement_Force(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_ARRAYS_VECTOR& F) PHYSBAM_OVERRIDE;
protected:
    void Discretize_Explicit_Viscous_Terms(const T dt){PHYSBAM_NOT_IMPLEMENTED();}
    void Implicit_Viscous_Update(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
//#####################################################################
};
}
#endif

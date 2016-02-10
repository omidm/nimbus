#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Frank Losasso, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_DYADIC
//#####################################################################
//
// Solves V_t+V GRAD(V)+GRAD(p)=0.
//
//#####################################################################
#ifndef __INCOMPRESSIBLE_DYADIC__
#define __INCOMPRESSIBLE_DYADIC__

#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_POLICY_DYADIC.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_DYADIC.h>
namespace PhysBAM{

template<class T_GRID>
class INCOMPRESSIBLE_DYADIC:public INCOMPRESSIBLE<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::MAP_MESH MAP_MESH;typedef typename T_GRID::CELL CELL;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FLOOD_FILL FLOOD_FILL;typedef ARRAY<VECTOR<bool,T_GRID::dimension> > T_FACE_NEIGHBORS;
public:
    typedef INCOMPRESSIBLE<T_GRID> BASE;
    using BASE::use_force;using BASE::viscosity;using BASE::use_variable_viscosity;using BASE::use_variable_vorticity_confinement;using BASE::use_explicit_part_of_implicit_viscosity;
    using BASE::gravity;using BASE::downward_direction;using BASE::max_time_step;using BASE::vorticity_confinement;using BASE::maximum_implicit_viscosity_iterations;using BASE::advection;
    
    T_GRID& grid;
    BOUNDARY_DYADIC<T_GRID,T>* boundary; // TODO: this might change to a vector boundary class
    PROJECTION_DYADIC<T_GRID> projection;
    ARRAY<T> force;
    ARRAY<T> variable_vorticity_confinement;
    ARRAY<T> flame_speed;
private:
    BOUNDARY_DYADIC<T_GRID,T>& boundary_default;
public:

    INCOMPRESSIBLE_DYADIC(T_GRID& grid_input,const bool flame=false);
    virtual ~INCOMPRESSIBLE_DYADIC();

    void Initialize_Grids()
    {INCOMPRESSIBLE<T_GRID>::Initialize_Grids(grid);
    projection.Initialize_Grid();
    if(projection.flame)flame_speed.Resize(grid.number_of_cells);}

    void Initialize_Grids(const T_GRID& grid_input) PHYSBAM_OVERRIDE
    {PHYSBAM_ASSERT(&grid==&grid_input);Initialize_Grids();}

    void Set_Custom_Boundary(BOUNDARY_DYADIC<T_GRID,T>& boundary_input)
    {boundary=&boundary_input;}

    void Set_Body_Force(const bool use_force_input=true)
    {use_force=use_force_input;if(use_force) force.Resize(grid.number_of_faces);else force.Resize(0);}

    void Use_Variable_Vorticity_Confinement(const bool use_variable_vorticity_confinement_input=true)
    {use_variable_vorticity_confinement=use_variable_vorticity_confinement_input;
    if(use_variable_vorticity_confinement) variable_vorticity_confinement.Resize(grid.number_of_nodes);
    else variable_vorticity_confinement.Resize(0);}

//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    void Advance_One_Time_Step(const T dt,const T time,const bool implicit_viscosity=false,const bool crank_nicolson=false,const bool first_step=false,const ARRAY<T>* phi_ghost=0);
    void Advance_One_Time_Step_Explicit_Part(const T dt,const T time,const bool implicit_viscosity=false,const bool crank_nicolson=false,const bool first_step=false,const ARRAY<T>* phi_ghost=0);
    void Advance_One_Time_Step_Implicit_Part(const T dt,const T time,const bool implicit_viscosity=false,const bool crank_nicolson=false);
    T CFL(const bool inviscid=false,const bool viscous_only=false) const;
    void Calculate_Pressure_Jump(const T dt,const T time);
    void Set_Dirichlet_Boundary_Conditions(const ARRAY<T>& phi,const T pressure=0);
    void Calculate_Flame_Speed();
    void Extrapolate_Velocity_Across_Interface(ARRAY<T>& phi_ghost,const bool enforce_divergence_free=false,const T band_width=3,const T damping=0,const TV& air_speed=TV(),
        const T_FACE_NEIGHBORS* face_neighbors_visible=0);
    virtual void Compute_Vorticity_Confinement_Force(T_GRID& grid,const ARRAY<TV>& V_ghost,ARRAY<TV>& F)=0;
//#####################################################################
};
}
#endif
#endif

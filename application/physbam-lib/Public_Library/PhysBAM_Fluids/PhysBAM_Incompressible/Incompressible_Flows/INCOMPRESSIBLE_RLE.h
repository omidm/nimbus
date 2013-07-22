#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_RLE
//#####################################################################
#ifndef __INCOMPRESSIBLE_RLE__
#define __INCOMPRESSIBLE_RLE__

#include <PhysBAM_Tools/Grids_RLE_Boundaries/BOUNDARY_POLICY_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/LEVELSET_POLICY_RLE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_RLE.h>
namespace PhysBAM{

template<class T_GRID>
class INCOMPRESSIBLE_RLE:public INCOMPRESSIBLE<T_GRID>
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_Y_ITERATOR FACE_Y_ITERATOR;
    typedef typename TV::SPIN T_SPIN;typedef ARRAY<VECTOR<bool,T_GRID::dimension> > T_FACE_NEIGHBORS;
public:
    typedef INCOMPRESSIBLE<T_GRID> BASE;
    using BASE::use_force;using BASE::max_time_step;using BASE::viscosity;using BASE::use_variable_viscosity;using BASE::gravity;using BASE::downward_direction;using BASE::advection;
    using BASE::valid_mask;

    const T_GRID& grid;
    ARRAY<T> V;
    BOUNDARY_RLE<T_GRID,T>* boundary;
    PROJECTION_RLE<T_GRID> projection;
    ARRAY<T> force;
    ADVECTION_FACE_VELOCITY_DONOR_CELL_RLE<T>& donor_cell_advection;
    MPI_RLE_GRID<T_GRID>* mpi_grid;
private:
    BOUNDARY_RLE<T_GRID,T>& boundary_default;
public:

    INCOMPRESSIBLE_RLE(const T_GRID& grid_input);
    virtual ~INCOMPRESSIBLE_RLE();

    void Initialize_Grids()
    {INCOMPRESSIBLE<T_GRID>::Initialize_Grids(grid);
    V.Resize(grid.number_of_faces);projection.Initialize_Grid();}

    void Initialize_Grids(const T_GRID& grid_input)
    {PHYSBAM_ASSERT(&grid==&grid_input);Initialize_Grids();}

    void Set_Custom_Boundary(BOUNDARY_RLE<T_GRID,T>& boundary_input)
    {boundary=&boundary_input;}

    void Set_Body_Force(const bool use_force_input=true)
    {use_force=use_force_input;if(use_force) force.Resize(grid.number_of_faces);else force.Clean_Memory();}

//#####################################################################
    void Advance_One_Time_Step_Explicit_Part(const T dt,const T time,const int substeps,const bool implicit_viscosity,const ARRAY<T>* phi_ghost);
    void Advance_One_Time_Step_Implicit_Part(const T dt,const T time,const bool implicit_viscosity=false,const ARRAY<T>* phi_ghost=0);
    T CFL(const bool inviscid=false,const bool viscous_only=false) const;
    void Set_Dirichlet_Boundary_Conditions(const ARRAY<T>& phi,const T pressure=0);
    void Extrapolate_Velocity_Across_Interface(ARRAY<T>& phi_ghost,const bool enforce_divergence_free,const T band_width,const T_FACE_NEIGHBORS* face_neighbors_visible);
    void Transfer_Velocity(const T_GRID& new_grid);
    void Transfer_Velocity_Ghost(const T_GRID& new_grid,ARRAY<T>& V_ghost) const;
private:
    struct Add_Impulse{template<class T_FACE> static void Apply(const T_GRID& grid,ARRAY<T>& V,const TV& impulse);};
    struct Compute_Maximum_Cell_Velocity{template<class T_FACE> static void Apply(const T_GRID& grid,const ARRAY<T>& V,ARRAY<TV>& max_V_cell);};
    struct Find_Fixed_Faces{template<class T_FACE> static void Apply(const T_GRID& grid,const ARRAY<T>& phi_ghost,ARRAY<bool>& fixed);};
    struct Transfer_Horizontal_Velocity{template<class T_FACE> static void Apply(const T_GRID& old_grid,const T_GRID& new_grid,const ARRAY<T>& V,ARRAY<T>& new_V,
        const ARRAY<bool>& valid_mask,ARRAY<bool>& new_valid_mask);};
    static void Transfer_Vertical_Velocity(const T_GRID& old_grid,const T_GRID& new_grid,const ARRAY<T>& V,ARRAY<T>& new_V,
        const ARRAY<bool>& valid_mask,ARRAY<bool>& new_valid_mask);
//#####################################################################
};
}
#endif
#endif

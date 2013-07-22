//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE  
//#####################################################################
//
// Inherited by INCOMPRESSIBLE_2D, and INCOMPRESSIBLE_3D.
// Also inherited by INCOMPRESSIBLE_FLAME_1D, INCOMPRESSIBLE_FLAME_2D, and INCOMPRESSIBLE_FLAME_3D.
//
//#####################################################################
#ifndef __INCOMPRESSIBLE__
#define __INCOMPRESSIBLE__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Advection_Collidable/ADVECTION_COLLIDABLE_POLICY_DYADIC.h>
#include <PhysBAM_Geometry/Grids_RLE_Advection_Collidable/ADVECTION_COLLIDABLE_POLICY_RLE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_COLLIDABLE_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_POLICY.h>
#include <PhysBAM_Dynamics/Advection_Equations/FIRE_ADVECTION_POLICY.h>
namespace PhysBAM{

template<class T_GRID>
class INCOMPRESSIBLE:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;typedef typename ADVECTION_COLLIDABLE_POLICY<T_GRID>::ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_LOOKUP_COLLIDABLE T_FACE_LOOKUP_COLLIDABLE;typedef typename INCOMPRESSIBLE_POLICY<T_GRID>::PROJECTION T_PROJECTION;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET_MULTIPLE T_LEVELSET_MULTIPLE;typedef typename FIRE_ADVECTION_POLICY<T_GRID>::ADVECTION_SEMI_LAGRANGIAN_FIRE_MULTIPHASE T_ADVECTION_SEMI_LAGRANGIAN_FIRE_MULTIPHASE;
    typedef typename FIRE_ADVECTION_POLICY<T_GRID>::ADVECTION_SEMI_LAGRANGIAN_FIRE_MULTIPHASE_COLLIDABLE T_ADVECTION_SEMI_LAGRANGIAN_FIRE_MULTIPHASE_COLLIDABLE;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef typename FIRE_ADVECTION_POLICY<T_GRID>::ADVECTION_WRAPPER_SEMI_LAGRANGIAN_FIRE_MULTIPHASE T_ADVECTION_WRAPPER_SEMI_LAGRANGIAN_FIRE_MULTIPHASE;
    typedef typename FIRE_INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP_FIRE_MULTIPHASE_COLLIDABLE T_FACE_LOOKUP_FIRE_MULTIPHASE_COLLIDABLE;
    typedef typename FIRE_ADVECTION_POLICY<T_GRID>::NESTED_ADVECTION_WRAPPER_SEMI_LAGRANGIAN_FIRE_MULTIPHASE_COLLIDABLE T_NESTED_ADVECTION_WRAPPER_SEMI_LAGRANGIAN_FIRE_MULTIPHASE_COLLIDABLE;
    typedef typename FIRE_ADVECTION_POLICY<T_GRID>::ADVECTION_WRAPPER_SEMI_LAGRANGIAN_FIRE_MULTIPHASE_COLLIDABLE T_ADVECTION_WRAPPER_SEMI_LAGRANGIAN_FIRE_MULTIPHASE_COLLIDABLE;
    typedef typename ADVECTION_COLLIDABLE_POLICY<T_GRID>::ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_SLIP_FACE T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_SLIP_FACE;
    typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_LOOKUP_COLLIDABLE_SLIP T_FACE_LOOKUP_COLLIDABLE_SLIP;
public:
    ADVECTION<T_GRID,T>* advection;
protected:               
    T max_time_step;
    T gravity;
    TV downward_direction;
    T surface_tension;
    T viscosity;
    T vorticity_confinement;
    T dt_old; // saved for the central time method
    bool use_force_x,use_force_y,use_force_z; // TODO(jontg): deprecated?
    bool use_force;
    bool conserve_kinetic_energy;
public:
    bool use_variable_surface_tension;
    bool use_variable_viscosity;
    int maximum_implicit_viscosity_iterations;
protected:
    bool use_variable_vorticity_confinement;
    bool use_explicit_part_of_implicit_viscosity;

    // for multiphase
    bool GFM; // (true) for GFM, (false) for delta function smearing
    T number_of_interface_cells; // e.g. 3 interface cells gives half_width = 1.5*dx
    ARRAY<T> viscosities; // constant viscosity for each region
    ARRAY<T,VECTOR<int,2> > surface_tensions;
    ARRAY<T> vorticity_confinements;
    bool nonzero_viscosity,nonzero_surface_tension;

public:
    T_FACE_ARRAYS_BOOL valid_mask;

    T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE* nested_semi_lagrangian_collidable;
    ADVECTION_WRAPPER_COLLIDABLE_FACE<T_GRID,T,T_FACE_LOOKUP,T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE,T_FACE_LOOKUP_COLLIDABLE>* semi_lagrangian_collidable;
    T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_SLIP_FACE* nested_semi_lagrangian_collidable_slip;
    ADVECTION_WRAPPER_COLLIDABLE_FACE<T_GRID,T,T_FACE_LOOKUP,T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_SLIP_FACE,T_FACE_LOOKUP_COLLIDABLE_SLIP>* semi_lagrangian_collidable_slip;
    T_ADVECTION_SEMI_LAGRANGIAN_FIRE_MULTIPHASE* nested_semi_lagrangian_fire_multiphase;
    T_ADVECTION_WRAPPER_SEMI_LAGRANGIAN_FIRE_MULTIPHASE* semi_lagrangian_fire_multiphase;
    T_ADVECTION_SEMI_LAGRANGIAN_FIRE_MULTIPHASE_COLLIDABLE* nested_nested_semi_lagrangian_fire_multiphase_collidable;
    T_NESTED_ADVECTION_WRAPPER_SEMI_LAGRANGIAN_FIRE_MULTIPHASE_COLLIDABLE* nested_semi_lagrangian_fire_multiphase_collidable;
    T_ADVECTION_WRAPPER_SEMI_LAGRANGIAN_FIRE_MULTIPHASE_COLLIDABLE* semi_lagrangian_fire_multiphase_collidable;
public:

    INCOMPRESSIBLE();
    virtual ~INCOMPRESSIBLE();

    virtual void Initialize_Grids(const T_GRID& grid)
    {valid_mask.Resize(grid.Domain_Indices(3),true,true,true);}

    void Set_Custom_Advection(ADVECTION<T_GRID,T>& advection_input)
    {advection=&advection_input;}
        
    void Set_Max_Time_Step(const T max_time_step_input=1e8)
    {max_time_step=max_time_step_input;}

    void Set_Gravity(const T gravity_input=9.8)
    {gravity=gravity_input;downward_direction=T_GRID::dimension>1?-TV::Axis_Vector(2):TV();}

    void Set_Gravity(const T gravity_input,const TV& downward_direction_input)
    {gravity=gravity_input;downward_direction=downward_direction_input;}

    void Set_Surface_Tension(const T surface_tension_input=0)
    {surface_tension=surface_tension_input;use_variable_surface_tension=false;nonzero_surface_tension=(surface_tension!=0);}

    void Set_Viscosity(const T viscosity_input=.001137)
    {viscosity=viscosity_input;use_variable_viscosity=false;nonzero_viscosity=(viscosity!=0);}

    void Set_Vorticity_Confinement(const T vorticity_confinement_input=.3)
    {vorticity_confinement=vorticity_confinement_input;}

    void Set_Maximum_Implicit_Viscosity_Iterations(const int maximum_implicit_viscosity_iterations_input=0)
    {maximum_implicit_viscosity_iterations=maximum_implicit_viscosity_iterations_input;}

    void Use_Explicit_Part_Of_Implicit_Viscosity(const bool use_explicit_part_of_implicit_viscosity_input=true)
    {use_explicit_part_of_implicit_viscosity=use_explicit_part_of_implicit_viscosity_input;}

    // functions used only in multiphase
    void Use_GFM() 
    {GFM=true;number_of_interface_cells=0;}

    void Use_Delta_Function_Method(const T number_of_interface_cells_input=3) 
    {GFM=false;number_of_interface_cells=number_of_interface_cells_input;}

    void Set_Viscosity(const ARRAY<T>& viscosities_input)
    {viscosities=viscosities_input;
    nonzero_viscosity=false;for(int i=1;i<=viscosities.m;i++) if(viscosities(i)!=0){nonzero_viscosity=true;break;}}

    void Set_Surface_Tension(const ARRAY<T,VECTOR<int,2> >& surface_tensions_input)
    {surface_tensions=surface_tensions_input;
    nonzero_surface_tension=false;for(int j=1;j<=surface_tensions.counts.x;j++) if(surface_tensions(1,j)!=0){nonzero_surface_tension=true;break;}}

    void Set_Vorticity_Confinement(const ARRAY<T> vorticity_confinements_input)
    {vorticity_confinements=vorticity_confinements_input;}

//#####################################################################
    void Use_Semi_Lagrangian_Collidable_Advection(T_GRID_BASED_COLLISION_GEOMETRY& body_list);
    void Use_Semi_Lagrangian_Collidable_Slip_Advection(T_GRID_BASED_COLLISION_GEOMETRY& body_list);
    void Use_Semi_Lagrangian_Fire_Multiphase_Advection(const T_PROJECTION& projection,const T_LEVELSET_MULTIPLE& levelset_multiple_n_plus_one);
    void Use_Semi_Lagrangian_Fire_Multiphase_Collidable_Advection(T_GRID_BASED_COLLISION_GEOMETRY& body_list,const T_PROJECTION& projection,const T_LEVELSET_MULTIPLE& levelset_multiple_n_plus_one);
//#####################################################################
};
}
#endif

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2007, Geoffrey Irving, Frank Losasso, Tamar Shinar, Eftychios Sifakis, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_EXAMPLE_DYADIC
//#####################################################################
#ifndef __SOLIDS_FLUIDS_EXAMPLE_DYADIC__
#define __SOLIDS_FLUIDS_EXAMPLE_DYADIC__

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_OBJECT_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/SOLIDS_COLLISIONS_FORWARD.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/INCOMPRESSIBLE_COLLISIONS_FORWARD.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_POLICY.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_CALLBACKS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_DYADIC.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE.h>
namespace PhysBAM{

template<class T_GRID>
class SOLIDS_FLUIDS_EXAMPLE_DYADIC:public SOLIDS_FLUIDS_EXAMPLE<typename T_GRID::VECTOR_T>,public LEVELSET_CALLBACKS<T_GRID>,public FLUIDS_PARAMETERS_CALLBACKS<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::MAP_MESH MAP_MESH;typedef typename T_GRID::CELL CELL;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    typedef typename LEVELSET_ADVECTION_POLICY<T_GRID>::LEVELSET_ADVECTION T_LEVELSET_ADVECTION;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FLOOD_FILL FLOOD_FILL;
    typedef typename T_GRID::UNIFORM_GRID UNIFORM_GRID;
    typedef typename UNIFORM_GRID::CELL_ITERATOR UNIFORM_CELL_ITERATOR;
    typedef typename UNIFORM_GRID::NODE_ITERATOR UNIFORM_NODE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_UNIFORM_ARRAYS;
    typedef typename LEVELSET_POLICY<T_GRID>::PARTICLE_LEVELSET T_PARTICLE_LEVELSET;typedef typename TURBULENCE_POLICY<TV>::TURBULENCE T_TURBULENCE;
    typedef typename INCOMPRESSIBLE_POLICY<T_GRID>::INCOMPRESSIBLE T_INCOMPRESSIBLE;typedef typename LEVELSET_POLICY<T_GRID>::EXTRAPOLATION_SCALAR EXTRAPOLATION_SCALAR;
    typedef typename T_UNIFORM_ARRAYS::template REBIND<TV>::TYPE T_UNIFORM_ARRAYS_VECTOR;typedef typename MATRIX_POLICY<TV>::TRANSFORMATION_MATRIX TRANSFORMATION_MATRIX;
    typedef typename LEVELSET_POLICY<UNIFORM_GRID>::LEVELSET UNIFORM_LEVELSET;
    typedef typename LEVELSET_POLICY<T_GRID>::EXTRAPOLATION_VECTOR EXTRAPOLATION_VECTOR;
    typedef typename LEVELSET_POLICY<T_GRID>::PARTICLE_LEVELSET_EVOLUTION T_PARTICLE_LEVELSET_EVOLUTION;
    typedef typename TV::SPIN T_ANGULAR_VELOCITY;typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_HELPER T_LINEAR_INTERPOLATION_DYADIC_HELPER;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE<TV> BASE;
    using BASE::first_frame;using BASE::initial_time;using BASE::restart;using BASE::output_directory;using BASE::Write_Frame_Title;using BASE::stream_type;
    using BASE::Adjust_Phi_With_Sources;using FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::Get_Source_Reseed_Mask;
    using FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::Adjust_Density_And_Temperature_With_Sources;using BASE::solids_parameters;
    using FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::Get_Source_Velocities;using FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::Get_Object_Velocities; // silence -Woverloaded-virtual
    using LEVELSET_CALLBACKS<T_GRID>::Get_Levelset_Velocity; // silence -Woverloaded-virtual
protected:
    using BASE::minimum_collision_thickness;
public:

    FLUIDS_PARAMETERS_DYADIC<T_GRID> fluids_parameters;
    ARRAY<bool> collision_body_affected_by_fluid;
    
    SOLIDS_FLUIDS_EXAMPLE_DYADIC(const STREAM_TYPE stream_type,const typename FLUIDS_PARAMETERS<T_GRID>::TYPE type);
    virtual ~SOLIDS_FLUIDS_EXAMPLE_DYADIC();

    void Get_Levelset_Velocity(const T_GRID& grid,T_LEVELSET& levelset,ARRAY<T>& face_velocities_levelset,const T time) const PHYSBAM_OVERRIDE
    {fluids_parameters.Get_Levelset_Velocity(grid,levelset,face_velocities_levelset,time);}

    void Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time) PHYSBAM_OVERRIDE
    {fluids_parameters.Adjust_Particle_For_Domain_Boundaries(particles,index,V,particle_type,dt,time);}

    void Get_Body_Force(ARRAY<TV>& force,const T dt,const T time) // TODO: change this to match base class version below
    {if(fluids_parameters.fire || fluids_parameters.smoke) fluids_parameters.Get_Body_Force(force,dt,time);else PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}

    void Get_Body_Force(ARRAY<T>& force,const T dt,const T time) PHYSBAM_OVERRIDE
    {PHYSBAM_FATAL_ERROR("call other version of Get_Body_Force");}

    void Setup_Initial_Refinement() PHYSBAM_OVERRIDE
    {fluids_parameters.Setup_Initial_Refinement(initial_time);}

    virtual void Specify_Refinement(ARRAY<typename CELL::REFINE_ACTION>& refine_action,ARRAY<CELL*>* cells_to_specify_for,const T dt,const T time)
    {return fluids_parameters.Specify_Refinement(refine_action,cells_to_specify_for,dt,time);}

    virtual void Update_Fluid_Parameters(const T dt,const T time)
    {fluids_parameters.Update_Fluid_Parameters(dt,time);}

//#####################################################################
    void Add_Volumetric_Body_To_Fluid_Simulation(RIGID_BODY<TV>& rigid_body,const bool affects_fluid=true,const bool affected_by_fluid=true);
    void Add_Thin_Shell_To_Fluid_Simulation(RIGID_BODY<TV>& rigid_body,const bool affects_fluid=true,const bool affected_by_fluid=true);
    void Add_To_Fluid_Simulation(DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions,const bool affects_fluid=true,const bool affected_by_fluid=true);
    virtual void Initialize_Solid_Fluid_Coupling(); // called after grids are initialized
    void Set_Dirichlet_Boundary_Conditions(const T time) PHYSBAM_OVERRIDE;
    template<class GEOMETRY> void Get_Source_Velocities(const GEOMETRY& source,const TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity);
    template<class GEOMETRY> void Adjust_Phi_With_Source(const GEOMETRY& source,const TRANSFORMATION_MATRIX& world_to_source);
    template<class GEOMETRY> void Get_Source_Reseed_Mask(const GEOMETRY& source,const TRANSFORMATION_MATRIX& world_to_source,ARRAY<bool>*& cell_centered_mask,const bool reset_mask);
    template<class GEOMETRY> void Adjust_Density_And_Temperature_With_Sources(const GEOMETRY& source,const TRANSFORMATION_MATRIX& world_to_source,const T source_density,
        const T source_temperature);
    void Revalidate_Fluid_Scalars();
    void Revalidate_Phi_After_Modify_Levelset();
    void Revalidate_Fluid_Velocity();
    void Get_Object_Velocities(const T dt,const T time) PHYSBAM_OVERRIDE;
    void Get_Object_Velocities(PROJECTION_DYADIC<T_GRID>& projection,ARRAY<T> face_velocities,const T dt,const T time) PHYSBAM_OVERRIDE;
    void Initialize_Swept_Occupied_Blocks_For_Advection(const T dt,const T time,const bool include_removed_particle_velocities);
    void Read_Output_Files_Fluids(const int frame) PHYSBAM_OVERRIDE;
    void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
#endif

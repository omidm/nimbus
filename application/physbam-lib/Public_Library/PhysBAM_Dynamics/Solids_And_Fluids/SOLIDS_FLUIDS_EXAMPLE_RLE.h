#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_EXAMPLE_RLE
//#####################################################################
#ifndef __SOLIDS_FLUIDS_EXAMPLE_RLE__
#define __SOLIDS_FLUIDS_EXAMPLE_RLE__

#include <PhysBAM_Tools/Fourier_Transforms_Calculations/DEEP_WATER_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_OBJECT_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/SOLIDS_COLLISIONS_FORWARD.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/INCOMPRESSIBLE_COLLISIONS_FORWARD.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_RLE.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_RLE.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_RLE.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_CALLBACKS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE.h>
namespace PhysBAM{

template<class T_GRID>
class SOLIDS_FLUIDS_EXAMPLE_RLE:public SOLIDS_FLUIDS_EXAMPLE<typename T_GRID::VECTOR_T>,protected FLUIDS_PARAMETERS<T_GRID>,public LEVELSET_CALLBACKS<T_GRID>,
    public FLUIDS_PARAMETERS_CALLBACKS<T_GRID>
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::BLOCK_ITERATOR BLOCK_ITERATOR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_Y_ITERATOR FACE_Y_ITERATOR;
    typedef typename T_GRID::ARRAYS_HORIZONTAL T_ARRAYS_HORIZONTAL;typedef typename T_ARRAYS_HORIZONTAL::template REBIND<int>::TYPE T_ARRAYS_HORIZONTAL_INT;
    typedef typename T_GRID::VECTOR_HORIZONTAL TV_HORIZONTAL;
    typedef typename T_GRID::BOX_HORIZONTAL_INT T_BOX_HORIZONTAL_INT;typedef typename T_GRID::VECTOR_HORIZONTAL_INT TV_HORIZONTAL_INT;
    typedef typename T_GRID::HORIZONTAL_GRID::CELL_ITERATOR HORIZONTAL_CELL_ITERATOR;typedef typename MATRIX_POLICY<TV>::TRANSFORMATION_MATRIX T_TRANSFORMATION_MATRIX;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE<TV> BASE;
    using BASE::first_frame;using BASE::initial_time;using BASE::restart;using BASE::output_directory;using BASE::Write_Frame_Title;
    using BASE::minimum_collision_thickness;using BASE::solids_parameters;using BASE::stream_type;
    using FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::Get_Source_Velocities;using FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::Get_Source_Reseed_Mask;
    using FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::Get_Object_Velocities;using FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::Delete_Particles_Inside_Objects; // silence -Woverloaded-virtual

    FLUIDS_PARAMETERS<T_GRID>& fluids_parameters;

    INCOMPRESSIBLE_RLE<T_GRID> incompressible;
    T_ARRAYS_HORIZONTAL_INT ground_j;
    PARTICLE_LEVELSET_RLE<T_GRID> particle_levelset;
    LEVELSET_ADVECTION_RLE<T_GRID> levelset_advection;
    MPI_RLE_GRID<T_GRID>* mpi_grid;
    bool refine_all_water;
    bool use_phi_for_vertical_refinement;
    T vertical_refinement_depth;
    bool refine_near_some_objects;
    bool refine_near_negative_particles;
    bool refine_near_removed_negative_particles;
    bool enforce_refinement_slope;
    bool use_incompressible_cfl;
    bool clamp_long_velocities_to_short_velocities;
    DEEP_WATER_EVOLUTION<TV_HORIZONTAL> deep_water;
    bool use_deep_water;

    SOLIDS_FLUIDS_EXAMPLE_RLE(const STREAM_TYPE stream_type);
    virtual ~SOLIDS_FLUIDS_EXAMPLE_RLE();

//#####################################################################
    void Use_Fluid_Coupling_Defaults();
    void Use_Fluid_Slip_Coupling_Defaults();
    void Use_No_Fluid_Coupling_Defaults();
    virtual T Initial_Ground(const typename T_GRID::VECTOR_HORIZONTAL& X) const=0;
    virtual T Initial_Phi(const TV& X) const=0;
    virtual T Initial_Phi_Object(const TV& X) const=0;
    void Initialize_Ground();
    void Initialize_Grid();
    void Initialize_Phi() PHYSBAM_OVERRIDE;
    void Initialize_MPI();
    virtual void Initialize_Grids_Extra(){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Update_Fluid_Parameters(const T dt,const T time);
    void Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time) PHYSBAM_OVERRIDE;
private:
    struct Adjust_Particle_For_Domain_Boundaries_Helper{template<class T_FACE> static void Apply(const FLUIDS_PARAMETERS<T_GRID>& example,TV& X,TV& X_new,TV& V,
        const T max_collision_distance,const T dt);};
public:
    void Delete_Particles_Inside_Objects(const T time);
    template<class T_PARTICLES> void Delete_Particles_Inside_Objects(ARRAY<T_PARTICLES*>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time);
    virtual void Delete_Particles_Inside_Objects(const BLOCK_ITERATOR& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time);
    virtual void Get_Neumann_And_Dirichlet_Boundary_Conditions(const T dt,const T time);
    void Set_Domain_Boundary_Conditions(const T time);
private:
    struct Set_Domain_Boundary_Conditions_Helper{template<class T_FACE> static void Apply(SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>& example);};
public:
    void Set_Dirichlet_Boundary_Conditions(const T time) PHYSBAM_OVERRIDE;
    void Revalidate_Fluid_Scalars();
    void Revalidate_Phi_After_Modify_Levelset();
    void Revalidate_Fluid_Velocity();
    void Get_Object_Velocities(const T dt,const T time) PHYSBAM_OVERRIDE;
    virtual void Get_Ground_Velocities(const T dt,const T time);
protected:
    struct Get_Ground_Velocities_Horizontal{template<class T_FACE> static void Apply(SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>& example);};
    void Get_Ground_Velocities_Vertical();
    struct Clamp_Long_Velocities_To_Short_Velocities{template<class T_FACE> static void Apply(const T_GRID& grid,ARRAY<T>& V);};
public:
    void Initialize_Swept_Occupied_Blocks_For_Advection(const T dt,const T time,const bool include_removed_particle_velocities);
    virtual void Get_Cell_Should_Be_Long(ARRAY<bool>& cell_should_be_long,const T time) const;
private:
    struct Get_Cell_Should_Be_Long_For_Objects{template<class T_FACE> static void Apply(const SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>& example,ARRAY<bool>& cell_should_be_long);};
public:
    virtual void Modify_Grid_After_Rebuild(T_GRID& new_grid,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Transfer_Extra_State(const T_GRID& new_grid){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    void Add_Volumetric_Body_To_Fluid_Simulation(RIGID_BODY<TV>& rigid_body,const bool affects_fluid);
    void Add_Thin_Shell_To_Fluid_Simulation(RIGID_BODY<TV>& rigid_body,const bool affects_fluid);
    void Add_To_Fluid_Simulation(DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions,const bool affects_fluid);
    void Initialize_Solid_Fluid_Coupling();
    template<class GEOMETRY> void Get_Source_Velocities(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity);
private:
    template<class GEOMETRY> struct Get_Source_Velocities_Horizontal{template<class T_FACE> static void Apply(PROJECTION_RLE<T_GRID>& projection,const GEOMETRY& source,
        const T_TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity);};
    template<class GEOMETRY> void Get_Source_Velocities_Vertical(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity);
public:
    template<class GEOMETRY> void Adjust_Phi_With_Source(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const bool only_near_existing_water=false);
    template<class GEOMETRY> void Get_Source_Reseed_Mask(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,ARRAY<bool>*& cell_centered_mask,
        const bool reset_mask) const;
    template<class GEOMETRY> void Get_Source_Cell_Should_Be_Long(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,ARRAY<bool>& cell_should_be_long,
        const bool only_near_existing_water=false) const;
    template<class T_PARTICLES> int Total_Number_Of_Particles(const ARRAY<T_PARTICLES*>& particles) const;
    template<class T_PARTICLES> void Read_Particles(const T_PARTICLES& template_particles,ARRAY<T_PARTICLES*>& particles,const std::string& prefix,const int frame);
    template<class T_PARTICLES> void Read_Particles(const T_PARTICLES& template_particles,ARRAY<T_PARTICLES*>& particles,T_PARTICLES*& particles_in_long_cells,const std::string& prefix,
        const int frame);
    template<class T_PARTICLES> void Write_Particles(const ARRAY<T_PARTICLES*>& particles,const std::string& prefix,const int frame) const;
    template<class T_PARTICLES> void Write_Particles(const ARRAY<T_PARTICLES*>& particles,const T_PARTICLES* particles_in_long_cells,const std::string& prefix,const int frame) const;
    void Read_Output_Files_Fluids(const int frame) PHYSBAM_OVERRIDE;
    void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
#endif 

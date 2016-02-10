//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Eftychios Sifakis, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_EXAMPLE_UNIFORM
//#####################################################################
#ifndef __SOLIDS_FLUIDS_EXAMPLE_UNIFORM__
#define __SOLIDS_FLUIDS_EXAMPLE_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_OBJECT_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/SOLIDS_COLLISIONS_FORWARD.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_CALLBACKS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Grid_Based_Fields/INCOMPRESSIBLE_FLUID_CONTAINER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/INCOMPRESSIBLE_COLLISIONS_FORWARD.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_POLICY.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_CALLBACKS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SPH_CALLBACKS.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY;

template<class T_GRID>
class SOLIDS_FLUIDS_EXAMPLE_UNIFORM:public SOLIDS_FLUIDS_EXAMPLE<typename T_GRID::VECTOR_T>,public LEVELSET_CALLBACKS<T_GRID>,public SPH_CALLBACKS<T_GRID>,
    public FLUIDS_PARAMETERS_CALLBACKS<T_GRID>,public CONSERVATION_CALLBACKS<T_GRID,VECTOR<typename T_GRID::SCALAR,T_GRID::dimension+2> >
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef VECTOR<T,T_GRID::dimension+2> TV_DIMENSION;
    typedef typename T_GRID::VECTOR_INT T_VECTOR_INT;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;typedef typename T_ARRAYS_SCALAR::template REBIND<char>::TYPE T_ARRAYS_CHAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename REBIND<T_FACE_ARRAYS_SCALAR,bool>::TYPE T_FACE_ARRAYS_BOOL;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    typedef typename MATRIX_POLICY<TV>::TRANSFORMATION_MATRIX T_TRANSFORMATION_MATRIX;
    typedef typename TV::SPIN T_ANGULAR_VELOCITY;typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_LOOKUP_COLLIDABLE T_FACE_LOOKUP_COLLIDABLE;
    typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::AVERAGING T_AVERAGING;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_SCALAR T_LINEAR_INTERPOLATION_SCALAR;
    typedef typename LEVELSET_POLICY<T_GRID>::PARTICLE_LEVELSET_EVOLUTION T_PARTICLE_LEVELSET_EVOLUTION;
    typedef typename LEVELSET_POLICY<T_GRID>::FAST_LEVELSET_T T_FAST_LEVELSET;typedef typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX SYMMETRIC_MATRIX;
    typedef typename LEVELSET_ADVECTION_POLICY<T_GRID>::FAST_LEVELSET_ADVECTION_T T_FAST_LEVELSET_ADVECTION;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE<TV> BASE;
    using BASE::output_directory;using BASE::first_frame;using BASE::restart;using BASE::Write_Frame_Title;using BASE::solids_parameters;using BASE::stream_type;
    using BASE::solids_fluids_parameters;using BASE::solid_body_collection;using BASE::solids_evolution;
    using BASE::Adjust_Phi_With_Sources;using BASE::minimum_collision_thickness;using FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::Adjust_Density_And_Temperature_With_Sources;
    using FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::Get_Source_Reseed_Mask;
    using FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::Get_Source_Velocities;using FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::Get_Object_Velocities; // silence -Woverloaded-virtual

    INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID> incompressible_fluid_container;
    FLUIDS_PARAMETERS_UNIFORM<T_GRID> fluids_parameters;
    int resolution;

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM(const STREAM_TYPE stream_type,const int number_of_regions,const typename FLUIDS_PARAMETERS<T_GRID>::TYPE type,const int array_collection_type=0);
    ~SOLIDS_FLUIDS_EXAMPLE_UNIFORM();

    void Get_Levelset_Velocity(const T_GRID& grid,T_LEVELSET& levelset,T_FACE_ARRAYS_SCALAR& V_levelset,const T time) const PHYSBAM_OVERRIDE
    {if(fluids_parameters.analytic_test) Get_Analytic_Velocities(time);V_levelset=incompressible_fluid_container.face_velocities;}

    void Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time) PHYSBAM_OVERRIDE
    {fluids_parameters.Adjust_Particle_For_Domain_Boundaries(particles,index,V,particle_type,dt,time);}

    void Get_Body_Force(T_FACE_ARRAYS_SCALAR& force,const T dt,const T time) PHYSBAM_OVERRIDE
    {if(fluids_parameters.fire||fluids_parameters.smoke) fluids_parameters.Get_Body_Force(force,dt,time);else PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}

    virtual void Update_Fluid_Parameters(const T dt,const T time)
    {fluids_parameters.Update_Fluid_Parameters(dt,time);}

    virtual void Evolve_Density_And_Temperature(const T dt,const T time)
    {fluids_parameters.Evolve_Density_And_Temperature(dt,time);}

    virtual void Apply_Isobaric_Fix(const T dt,const T time)
    {fluids_parameters.Apply_Isobaric_Fix(dt,time);}

//#####################################################################
    void Add_Volumetric_Body_To_Fluid_Simulation(RIGID_BODY<TV>& rigid_body,bool add_collision=true,bool add_coupling=true);
    void Add_Thin_Shell_To_Fluid_Simulation(RIGID_BODY<TV>& rigid_body,bool add_collision=true,bool add_coupling=true);
    void Add_To_Fluid_Simulation(DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions,bool add_collision=true,bool add_coupling=true);
    virtual void Clamp_Dt_Adaptively(T_GRID& grid,const T_ARRAYS_DIMENSION_SCALAR& rhs,const T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_BOOL& psi,T_ARRAYS_SCALAR& rho_dt,T_ARRAYS_SCALAR& e_dt,const T& dt,T& clamp_rho,T& clamp_e);
    virtual void Initialize_MPI();
    virtual void Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization(); // called before grids are initialized
    virtual void Initialize_Solid_Fluid_Coupling_After_Grid_Initialization(); // called after grids are initialized
    virtual void Initialize_Compressible_Incompressible_Coupling();
    virtual void Set_Ghost_Density_And_Temperature_Inside_Flame_Core();
    void Set_Dirichlet_Boundary_Conditions(const T time) PHYSBAM_OVERRIDE;
    template<class GEOMETRY> void Get_Source_Velocities(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity);
    template<class GEOMETRY> void Get_Source_Velocities(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity,const T_FACE_ARRAYS_BOOL& invalid_mask);
    template<class GEOMETRY> void Adjust_Phi_With_Source(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source);
    template<class GEOMETRY> void Adjust_Phi_With_Source(const GEOMETRY& source,const int region,const T_TRANSFORMATION_MATRIX& world_to_source);
    template<class GEOMETRY> void Get_Source_Reseed_Mask(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,T_ARRAYS_BOOL*& cell_centered_mask,const bool reset_mask);
    template<class GEOMETRY> void Adjust_Density_And_Temperature_With_Sources(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const T source_density,
        const T source_temperature);
    void Revalidate_Fluid_Scalars();
    void Revalidate_Phi_After_Modify_Levelset();
    void Revalidate_Fluid_Velocity(T_FACE_ARRAYS_SCALAR& face_velocities);
    bool Get_Psi_D_Inside_Solids(T_ARRAYS_BOOL& psi_D) PHYSBAM_OVERRIDE; // return true if overwritten
    void Post_Velocity_Advection_Callback(const T dt,const T time){}
    void Get_Object_Velocities(LAPLACE_UNIFORM<T_GRID>* elliptic_solver,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time) PHYSBAM_OVERRIDE;
    void Get_Levelset_Velocity(const T_GRID& grid,LEVELSET_MULTIPLE_UNIFORM<T_GRID>& levelset_multiple,T_FACE_ARRAYS_SCALAR& V_levelset,const T time) const PHYSBAM_OVERRIDE;
    void Initialize_Swept_Occupied_Blocks_For_Advection(const T dt,const T time,T maximum_fluid_velocity,const bool include_removed_particle_velocities);
    void Read_Output_Files_Fluids(const int frame) PHYSBAM_OVERRIDE;
    virtual void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE;
    void Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time) PHYSBAM_OVERRIDE;
    void Log_Parameters() const PHYSBAM_OVERRIDE;
    void Register_Options() PHYSBAM_OVERRIDE;
    void Parse_Options() PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

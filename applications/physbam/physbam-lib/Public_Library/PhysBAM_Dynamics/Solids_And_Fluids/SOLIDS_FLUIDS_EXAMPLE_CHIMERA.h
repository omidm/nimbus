//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Eftychios Sifakis, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_EXAMPLE_CHIMERA
//#####################################################################
#ifndef __SOLIDS_FLUIDS_EXAMPLE_CHIMERA__
#define __SOLIDS_FLUIDS_EXAMPLE_CHIMERA__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_OBJECT_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/SOLIDS_COLLISIONS_FORWARD.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Grid_Based_Fields/INCOMPRESSIBLE_FLUID_CONTAINER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/INCOMPRESSIBLE_COLLISIONS_FORWARD.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_POLICY.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_CALLBACKS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SPH_CALLBACKS.h>
#include <PhysBAM_Dynamics/Grids_Chimera/Grids_Chimera_Boundaries/BOUNDARY_CHIMERA.h>
#include <PhysBAM_Dynamics/Grids_Chimera/Parallel_Computation/CHIMERA_GRID.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <PhysBAM_Dynamics/Boundaries/BOUNDARY_PHI_WATER.h>
#include <PhysBAM_Dynamics/Grids_Chimera_PDE_Linear/Parallel_Computation/LAPLACE_CHIMERA_MPI.h>
#include <PhysBAM_Geometry/Grids_ALE_Advection/ADVECTION_WRAPPER_ALE.h>
#include <PhysBAM_Geometry/Grids_ALE_Advection/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_ALE.h>
#include <PhysBAM_Dynamics/Grids_Chimera_Advection/ADVECTION_WRAPPER_MACCORMACK_CHIMERA.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY;

template<class T_GRID,class T2,class T_NESTED_ADVECTION,class T_NESTED_LOOKUP> class ADVECTION_WRAPPER_MACCORMACK_CHIMERA;

template<class T_GRID>
class SOLIDS_FLUIDS_EXAMPLE_CHIMERA:public SOLIDS_FLUIDS_EXAMPLE<typename T_GRID::VECTOR_T>,public LEVELSET_CALLBACKS<T_GRID>,public SPH_CALLBACKS<T_GRID>,
                                    public FLUIDS_PARAMETERS_CALLBACKS<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef VECTOR<T,T_GRID::dimension+2> TV_DIMENSION;typedef typename TV::SPIN T_SPIN;
    typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,T_GRID::dimension> TV_SIDES;
    typedef typename T_GRID::VECTOR_INT T_VECTOR_INT;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
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
    typedef ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T, AVERAGING_UNIFORM<GRID<TV>, FACE_LOOKUP_UNIFORM<GRID<TV> > >,LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T,FACE_LOOKUP_UNIFORM<GRID<TV> > > > T_ADVECTION_UNIFORM;
    typedef ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_ALE<GRID<TV>,T,FACE_LOOKUP_UNIFORM<T_GRID>,FACE_LOOKUP_COLLIDABLE_UNIFORM<T_GRID> > T_ADVECTION_COLLIDABLE_ALE;
    typedef ADVECTION_WRAPPER_ALE<GRID<TV>,T,T_ADVECTION_UNIFORM,FACE_LOOKUP_UNIFORM<T_GRID> > T_ADVECTION_WRAPPER_ALE;
    typedef ADVECTION_WRAPPER_MACCORMACK_CHIMERA<GRID<TV>,T,T_ADVECTION_UNIFORM,FACE_LOOKUP_UNIFORM<T_GRID> > T_ADVECTION_MACCORMACK_CHIMERA;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE<TV> BASE;
    using BASE::output_directory;using BASE::first_frame;using BASE::auto_restart;using BASE::restart;using BASE::restart_frame;using BASE::Write_Frame_Title;using BASE::solids_parameters;using BASE::stream_type;
    using BASE::solids_fluids_parameters;using BASE::solid_body_collection;using BASE::solids_evolution;
    using BASE::minimum_collision_thickness;using FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::Adjust_Density_And_Temperature_With_Sources;
    using FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::Get_Source_Reseed_Mask;
    using FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::Get_Source_Velocities;using FLUIDS_PARAMETERS_CALLBACKS<T_GRID>::Get_Object_Velocities; // silence -Woverloaded-virtual

    virtual bool Get_Source_Velocity(const TV& location,const TV& normal,T& normal_gradient,const T time,const T dt){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return false;}
    virtual bool Get_Dirichlet_Boundary_Condition(const TV& location,T& phi,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return false;}

    RIGID_GRID_COLLECTION<T_GRID> &rigid_grid_collection;
    ARRAY<INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>*> incompressible_fluid_containers;
    FLUIDS_PARAMETERS_UNIFORM<T_GRID> fluids_parameters;
    T cfl_alpha;//number of ghost cells from grid motion constraint
    T cfl_beta;//number of ghost cells from fluid velocity constraint
    CHIMERA_GRID<T_GRID>* chimera_grid;
    ARRAY<SOLID_BODY_COLLECTION<TV>*> object_space_solid_body_collection;
    BOUNDARY_CHIMERA<T_GRID,T>* boundary_chimera;
    BOUNDARY_PHI_WATER<T_GRID>* boundary_phi_water;
    VECTOR<int,2> particles_counter;
    int resolution;
    bool simulate_solids;
    bool kinematic_grid;
    bool run_advection_test;
    bool solve_poisson_equation;
    bool solve_heat_equation;
    bool solve_heat_equation_on_faces;
    bool solve_heat_equation_on_faces_coupled;
    bool use_flip_face_update;
    bool use_flip_face_update_consistent;
    bool use_flip_cell_update;
    bool use_flip_update_overlap;
    TV_INT test_numbers;
    bool enforce_second_order_maccormack_scalar;
    bool enforce_second_order_maccormack_vector;
    bool clamp_extrema_maccormack_scalar;
    bool clamp_extrema_maccormack_vector;
    bool fill_overlap_region_recursively;
    bool use_maccormack_advection_scalar;
    bool use_maccormack_advection_vector;
    bool use_collidable_advection;
    ARRAY<ADVECTION<T_GRID,T>*> advection;
    T_ADVECTION_MACCORMACK_CHIMERA* advection_maccormack_scalar;
    T_ADVECTION_MACCORMACK_CHIMERA* advection_maccormack_vector;
    ARRAY<T_ADVECTION_COLLIDABLE_ALE*> advection_collidable;
    T_ADVECTION_UNIFORM advection_uniform;
    bool revalidate_velocities_in_projection;
    bool define_slip_boundary;
    bool write_voronoi;

    //added here because it isn't used in real simulations only for multidirectional diffusion equation testing
    ARRAY<T_FACE_ARRAYS_SCALAR> face_temperatures;

    SOLIDS_FLUIDS_EXAMPLE_CHIMERA(const STREAM_TYPE stream_type,const int number_of_regions,const typename FLUIDS_PARAMETERS<T_GRID>::TYPE type,const int array_collection_type=0);
    ~SOLIDS_FLUIDS_EXAMPLE_CHIMERA();

    RIGID_GRID<T_GRID>& Add_Rigid_Grid()
    {
        INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>* new_container=new INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>(rigid_grid_collection,fluids_parameters.number_of_ghost_cells);
        incompressible_fluid_containers.Append(new_container);
        return new_container->rigid_grid;
    }

    void Add_Object_Space_Solid_Body_Collection()
    {
        if(!simulate_solids) return;
        for(int grid_index=1;grid_index<=rigid_grid_collection.particles.array_collection->Size();grid_index++)
            object_space_solid_body_collection.Append(new SOLID_BODY_COLLECTION<TV>(this,0));
    }

    void Update_Grid(const T time,const T dt)
    {            
        for(int i=1;i<=chimera_grid->global_grid.particles.array_collection->Size();i++) chimera_grid->global_grid.Rigid_Grid(i).Update_Frame_From_Twist(dt,time);
        if(kinematic_grid){
            Set_Grid_External_Positions(chimera_grid->global_grid.particles.X,chimera_grid->global_grid.particles.rotation,time+dt);
            Set_Grid_External_Velocities(chimera_grid->global_grid.particles.V,chimera_grid->global_grid.particles.angular,time+dt,time+dt);}
        chimera_grid->Synchronize_Grid_Position();// synchronize positions of global grids and local grids
        chimera_grid->Update_Oriented_Boxes();
    }

    void Update_Grid(ARRAY<TV>& X,ARRAY<TV>& V,ARRAY<ROTATION<TV> >& rotation,ARRAY<typename TV::SPIN>& angular,const T time,const T dt)
    {
        for(int i=1;i<=chimera_grid->global_grid.particles.array_collection->Size();i++){
            X(i)+=V(i)*dt;rotation(i)=ROTATION<TV>::From_Rotation_Vector(angular(i)*dt)*rotation(i);rotation(i).Normalize();}
        if(kinematic_grid){
            Set_Grid_External_Positions(X,rotation,time+dt);
            Set_Grid_External_Velocities(V,angular,time+dt,time+dt);}
    }

    void Update_Object_Space_Rigid_Body_Collection(bool use_previous_grid_frame)
    {
        if(!simulate_solids || !use_collidable_advection) return;
        RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
        for(int i=1;i<=rigid_grid_collection.particles.array_collection->Size();i++){
            RIGID_BODY_COLLECTION<TV>& object_space_rigid_body_collection=object_space_solid_body_collection(i)->rigid_body_collection;
            RIGID_GRID<T_GRID>& rigid_grid=rigid_grid_collection.Rigid_Grid(i);
            for(int b=1;b<=rigid_body_collection.rigid_geometry_collection.particles.array_collection->Size();b++){
                if(use_previous_grid_frame){
                    FRAME<TV> new_frame=rigid_grid.previous_state.frame.Inverse()*rigid_body_collection.Rigid_Body(b).Frame();
                    object_space_rigid_body_collection.Rigid_Body(b).X()=new_frame.t;
                    object_space_rigid_body_collection.Rigid_Body(b).Rotation()=new_frame.r;
                    object_space_rigid_body_collection.Rigid_Body(b).V()=rigid_grid.previous_state.Object_Space_Vector(rigid_body_collection.Rigid_Body(b).V());}
                else{
                    FRAME<TV> new_frame=rigid_grid.Frame().Inverse()*rigid_body_collection.Rigid_Body(b).Frame();
                    object_space_rigid_body_collection.Rigid_Body(b).X()=new_frame.t;
                    object_space_rigid_body_collection.Rigid_Body(b).Rotation()=new_frame.r;
                    object_space_rigid_body_collection.Rigid_Body(b).V()=rigid_grid.Object_Space_Vector(rigid_body_collection.Rigid_Body(b).V());}}}
    }

    void Update_Object_Space_Rigid_Body_Collection_To_New_Grid_Frame(int grid_index)
    {
        if(!simulate_solids || !use_collidable_advection) return;
        RIGID_BODY_COLLECTION<TV>& object_space_rigid_body_collection=object_space_solid_body_collection(grid_index)->rigid_body_collection;
        RIGID_GRID<T_GRID>& rigid_grid=rigid_grid_collection.Rigid_Grid(grid_index);
        for(int b=1;b<=object_space_rigid_body_collection.rigid_geometry_collection.particles.array_collection->Size();b++){
            FRAME<TV> new_frame=rigid_grid.Frame().Inverse()*rigid_grid.previous_state.frame*object_space_rigid_body_collection.Rigid_Body(b).Frame();
            object_space_rigid_body_collection.Rigid_Body(b).X()=new_frame.t;
            object_space_rigid_body_collection.Rigid_Body(b).Rotation()=new_frame.r;
            object_space_rigid_body_collection.Rigid_Body(b).V()=rigid_grid.Current_Object_Space_Vector(object_space_rigid_body_collection.Rigid_Body(b).V());}
    }

    void Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time) PHYSBAM_OVERRIDE
    {
    }

    void Extrapolate_Phi_Into_Objects(const T time)
    {
        if(!simulate_solids) return;
        for(int grid_index=1;grid_index<=rigid_grid_collection.particles.array_collection->Size();grid_index++){
            RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection=object_space_solid_body_collection(grid_index)->rigid_body_collection.rigid_geometry_collection;
            T_GRID grid=rigid_grid_collection.Rigid_Grid(grid_index).grid;
            T_PARTICLE_LEVELSET_EVOLUTION& particle_levelset_evolution=incompressible_fluid_containers(grid_index)->particle_levelset_evolution;
            for(int id=1;id<=rigid_geometry_collection.particles.array_collection->Size();id++){
                T_ARRAYS_SCALAR phi_object(grid.Domain_Indices(fluids_parameters.number_of_ghost_cells));
                for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
                    phi_object(iterator.Cell_Index())=-rigid_geometry_collection.Rigid_Geometry(id).Implicit_Geometry_Extended_Value(iterator.Location());
                EXTRAPOLATION_UNIFORM<T_GRID,T> extrapolate(grid,phi_object,particle_levelset_evolution.particle_levelset.levelset.phi,0);
                extrapolate.Set_Band_Width((T)fluids_parameters.number_of_ghost_cells);extrapolate.Extrapolate();}}
    }

    void Initialize_Boundary_Phi_Water(const T_FACE_ARRAYS_SCALAR* V_input,const TV_SIDES& constant_extrapolation=TV_SIDES::Constant_Vector(TV_BOOL2::Constant_Vector(true)))
    {boundary_phi_water=new BOUNDARY_PHI_WATER<T_GRID>(constant_extrapolation);boundary_phi_water->Set_Velocity_Pointer(*V_input);}
    
    void Get_Velocity_Pointer_For_Boundary_Phi_Water(const T_FACE_ARRAYS_SCALAR* V_input)
    {if(boundary_phi_water) boundary_phi_water->Set_Velocity_Pointer(*V_input);}

    void Get_Body_Force(T_FACE_ARRAYS_SCALAR& force,const T dt,const T time) PHYSBAM_OVERRIDE
    {if(fluids_parameters.fire||fluids_parameters.smoke) fluids_parameters.Get_Body_Force(force,dt,time);else PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}

    virtual void Update_Fluid_Parameters(const T dt,const T time)
    {fluids_parameters.Update_Fluid_Parameters(dt,time);}

    virtual void Evolve_Density_And_Temperature(const T dt,const T time)
    {fluids_parameters.Evolve_Density_And_Temperature(dt,time);}

    virtual void Apply_Isobaric_Fix(const T dt,const T time)
    {fluids_parameters.Apply_Isobaric_Fix(dt,time);}

    virtual void Initialize_Density()
    {}
    
    virtual void Set_Grid_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) {};
    virtual void Set_Grid_External_Velocities(ARRAY_VIEW<TV> V,ARRAY_VIEW<T_SPIN> angular_velocity,const T velocity_time,const T current_position_time) {};
    virtual void Initialize_Phi();
    virtual bool Adjust_Phi_With_Sources(const T time);
    virtual void Get_Levelset_Velocity(const T_GRID& grid,T_LEVELSET& levelset,T_FACE_ARRAYS_SCALAR& V_levelset,const T time) const PHYSBAM_OVERRIDE;
    virtual void Get_Levelset_Velocity(const T_GRID& grid,LEVELSET_MULTIPLE_UNIFORM<T_GRID>& levelset,T_FACE_ARRAYS_SCALAR& V_levelset,const T time) const PHYSBAM_OVERRIDE {}
    virtual T Get_Analytic_Laplacian(const TV& location){return (T)0;}
    virtual void Fill_Ghost_Regions_Analytic(ARRAY<T_ARRAYS_SCALAR>& density_ghost,ARRAY<T_FACE_ARRAYS_SCALAR>& face_velocities_ghost,const int number_of_ghost_cells,const T time,const T dt) {if(run_advection_test) PHYSBAM_FATAL_ERROR();};
//#####################################################################
    void Advect_Scalar_Field(ARRAY<T_ARRAYS_SCALAR*>& u_array,ARRAY<T_ARRAYS_SCALAR*>& u_ghost_array,ARRAY<T_FACE_ARRAYS_SCALAR*>& face_velocities_ghost,const T dt,const T time);
    void Advect_Vector_Field(ARRAY<T_FACE_ARRAYS_SCALAR*>& u_array,ARRAY<T_FACE_ARRAYS_SCALAR*>& u_ghost_array,ARRAY<T_FACE_ARRAYS_SCALAR*>& face_velocities_ghost,const T dt,const T time);
    void Extrapolate_Scalar_Field_Into_Objects(ARRAY<T_ARRAYS_SCALAR*>& u_array,const T dt,const T time);

    void Exchange_Split_Grid_Ghost_Cells(ARRAY<T_ARRAYS_SCALAR*>& u_array);

    void Fill_Ghost_Cells_Chimera(int number_of_ghost_cells,ARRAY<T_ARRAYS_SCALAR*>& u_array,ARRAY<T_ARRAYS_SCALAR*>& u_ghost_array,bool only_overwrite_nan=false);
    void Fill_Ghost_Cells_Face_Chimera(int number_of_ghost_cells,ARRAY<T_FACE_ARRAYS_SCALAR*>& u_array,ARRAY<T_FACE_ARRAYS_SCALAR*>& u_ghost_array);
    void Coupling_Overlap_Regions_Cell(int number_of_ghost_cells,ARRAY<T_ARRAYS_SCALAR*>& u_array,bool only_overwrite_nan=false);
    void Coupling_Overlap_Regions_Face(int number_of_ghost_cells,ARRAY<T_FACE_ARRAYS_SCALAR*>& u_array);

    void Add_Volumetric_Body_To_Fluid_Simulation(int grid_index,RIGID_BODY<TV>& rigid_body,bool add_collision=true,bool add_coupling=true);
    void Add_Thin_Shell_To_Fluid_Simulation(RIGID_BODY<TV>& rigid_body,bool add_collision=true,bool add_coupling=true);
    void Add_To_Fluid_Simulation(DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions,bool add_collision=true,bool add_coupling=true);
    virtual void Initialize_MPI();
    virtual void Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization(); // called before grids are initialized
    virtual void Initialize_Solid_Fluid_Coupling_After_Grid_Initialization(); // called after grids are initialized
    virtual void Initialize_Compressible_Incompressible_Coupling();
    virtual void Set_Ghost_Density_And_Temperature_Inside_Flame_Core();
    void Set_Dirichlet_Boundary_Conditions(const T time) PHYSBAM_OVERRIDE;
    void Get_Source_Velocities_Chimera(const T time,const T dt);
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
    void Post_Velocity_Advection_Callback(const T dt,const T time){}
    void Get_Object_Velocities(LAPLACE_UNIFORM<T_GRID>* elliptic_solver,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time) PHYSBAM_OVERRIDE;
    void Initialize_Swept_Occupied_Blocks_For_Advection(const T dt,const T time,const bool include_removed_particle_velocities);
    void Read_Output_Files_Fluids(const int frame) PHYSBAM_OVERRIDE;
    virtual void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE;
    void Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time) PHYSBAM_OVERRIDE;
    void Log_Parameters() const PHYSBAM_OVERRIDE;
    void Register_Options() PHYSBAM_OVERRIDE;
    void Parse_Options() PHYSBAM_OVERRIDE;

    virtual T Get_Density(const TV& location)
    {return fluids_parameters.density;}
//#####################################################################
};
}
#endif

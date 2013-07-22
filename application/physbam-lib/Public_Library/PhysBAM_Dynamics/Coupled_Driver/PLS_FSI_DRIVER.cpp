//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLS_FSI_DRIVER
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Tools/Utilities/INTERRUPTS.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Level_Sets/EXTRAPOLATION_HIGHER_ORDER.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SURFACE_TENSION_FORCE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/FRACTURE_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Boundaries/BOUNDARY_MAC_GRID_SOLID_WALL_SLIP.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/DETONATION_SHOCK_DYNAMICS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Boundaries/BOUNDARY_PHI_WATER.h>
#include <PhysBAM_Dynamics/Coupled_Driver/PLS_FSI_DRIVER.h>
#include <PhysBAM_Dynamics/Coupled_Driver/PLS_FSI_EXAMPLE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_UNIFORM_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PLS_FSI_DRIVER<TV>::
PLS_FSI_DRIVER(PLS_FSI_EXAMPLE<TV>& example_input)
    :BASE(example_input),example(example_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PLS_FSI_DRIVER<TV>::
~PLS_FSI_DRIVER()
{}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Preprocess_Frame(const int frame)
{
    if(example.substeps_delay_frame==frame){
        example.Set_Write_Substeps_Level(example.substeps_delay_level);
        output_number=frame-1;}
    example.Preprocess_Frame(frame);
}
//#####################################################################
// Function Execute_Main_Program
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Execute_Main_Program()
{
    {LOG::SCOPE scope("INITIALIZING","Initializing");
    Initialize();
    example.Post_Initialization();
    example.Log_Parameters();
    if(!example.restart) Write_Output_Files(example.first_frame);}
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Function Simulate_To_Frame
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Simulate_To_Frame(const int frame_input)
{
    while(current_frame<frame_input){
        LOG::SCOPE scope("FRAME","Frame %d",current_frame+1);
        Preprocess_Frame(current_frame+1);
        Advance_To_Target_Time(example.Time_At_Frame(current_frame+1));
        Postprocess_Frame(++current_frame);
        if(example.write_output_files && example.write_substeps_level==-1) Write_Output_Files(current_frame);
        else if(example.write_substeps_level!=-1) Write_Substep(STRING_UTILITIES::string_sprintf("END Frame %d",current_frame),0,example.write_substeps_level);
        std::stringstream ss;
        ss<<"TIME = "<<time<<std::endl;
        LOG::filecout(ss.str());}
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Initialize()
{
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>& collision_bodies_affecting_fluid=*example.fluids_parameters.collision_bodies_affecting_fluid;

    if(example.auto_restart){
        std::string last_frame_file=example.output_directory+"/common/last_frame";
        int last_frame;FILE_UTILITIES::Read_From_Text_File(last_frame_file,last_frame);
        example.restart=true;example.restart_frame=last_frame;
        std::stringstream ss;
        ss<<"Auto Restart from frame "<<last_frame<<" (from file "<<last_frame_file<<")"<<std::endl;
        LOG::filecout(ss.str());}
    if(example.restart){current_frame=example.restart_frame;Read_Time(current_frame);}
    else current_frame=example.first_frame;
    output_number=current_frame;
    time=example.Time_At_Frame(current_frame);
    example.fluids_parameters.callbacks=&example;

    *example.fluids_parameters.grid=example.fluids_parameters.grid->Get_MAC_Grid();
    example.fluids_parameters.p_grid=*example.fluids_parameters.grid;
    example.Initialize_Fluids_Grids();
    GRID<TV>& grid=*example.fluids_parameters.grid;

    example.fluids_parameters.use_poisson=true;
    SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* coupled_evolution=new SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>(example.solids_parameters,example.solid_body_collection,
        example.fluids_parameters,example.solids_fluids_parameters,example.incompressible_fluid_container);
    delete example.solids_evolution;
    example.solids_evolution=coupled_evolution;
    example.fluids_parameters.projection=coupled_evolution;
    SOLIDS_EVOLUTION<TV>& solids_evolution=*coupled_evolution;

    solids_evolution.Set_Solids_Evolution_Callbacks(example);
    example.Initialize_Bodies();

    example.fluids_parameters.particle_levelset_evolution=new typename LEVELSET_POLICY<GRID<TV> >::PARTICLE_LEVELSET_EVOLUTION(*example.fluids_parameters.grid,example.fluids_parameters.number_of_ghost_cells);
    example.fluids_parameters.projection=coupled_evolution;
    example.fluids_parameters.incompressible=new INCOMPRESSIBLE_UNIFORM<GRID<TV> >(*example.fluids_parameters.grid,*example.fluids_parameters.projection);
    example.fluids_parameters.projection=0;
    example.fluids_parameters.phi_boundary=&example.fluids_parameters.phi_boundary_water; // override default
    example.fluids_parameters.phi_boundary_water.Set_Velocity_Pointer(example.incompressible_fluid_container.face_velocities);
    example.fluids_parameters.boundary_mac_slip.Set_Phi(example.fluids_parameters.particle_levelset_evolution->phi);

    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<TV> >* particle_levelset_evolution=example.fluids_parameters.particle_levelset_evolution;
    INCOMPRESSIBLE_UNIFORM<GRID<TV> >* incompressible=example.fluids_parameters.incompressible;

    example.Parse_Late_Options();

    Initialize_Fluids_Grids(); // this needs to be here because the arrays have to be resized for multiphase

    example.fluids_parameters.collision_bodies_affecting_fluid->Initialize_Grids();
    coupled_evolution->Setup_Boundary_Condition_Collection();
    solids_evolution.time=time;

    if(example.restart){
        LOG::SCOPE scope("reading solids data");
        example.Read_Output_Files_Solids(example.restart_frame);
        solids_evolution.time=time=example.Time_At_Frame(example.restart_frame);}

    solids_evolution.Initialize_Deformable_Objects(example.frame_rate,example.restart);

    solids_evolution.Initialize_Rigid_Bodies(example.frame_rate,example.restart);

    // time
    particle_levelset_evolution->Set_Time(time);
    particle_levelset_evolution->Set_CFL_Number(example.fluids_parameters.cfl);

    // sets up the proper wall states
    VECTOR<VECTOR<bool,2>,TV::dimension> domain_open_boundaries=VECTOR_UTILITIES::Complement(example.fluids_parameters.domain_walls);
    if(example.fluids_parameters.phi_boundary) example.fluids_parameters.phi_boundary->Set_Constant_Extrapolation(domain_open_boundaries);
    if(example.fluids_parameters.fluid_boundary) example.fluids_parameters.fluid_boundary->Set_Constant_Extrapolation(domain_open_boundaries);

    example.Initialize_Advection();
    Initialize_Fluids_Grids(); // initialize the valid masks

    // initialize levelset
    particle_levelset_evolution->Set_Number_Particles_Per_Cell(example.fluids_parameters.number_particles_per_cell);
    particle_levelset_evolution->Set_Levelset_Callbacks(example);
    particle_levelset_evolution->Initialize_FMM_Initialization_Iterative_Solver(example.fluids_parameters.refine_fmm_initialization_with_iterative_solver);

    particle_levelset_evolution->particle_levelset.levelset.Set_Custom_Boundary(*example.fluids_parameters.phi_boundary);
    particle_levelset_evolution->Bias_Towards_Negative_Particles(example.fluids_parameters.bias_towards_negative_particles);
    if(example.fluids_parameters.use_removed_positive_particles) particle_levelset_evolution->Particle_Levelset(1).Use_Removed_Positive_Particles();
    if(example.fluids_parameters.use_removed_negative_particles) particle_levelset_evolution->Particle_Levelset(1).Use_Removed_Negative_Particles();
    if(example.fluids_parameters.store_particle_ids){
        particle_levelset_evolution->Particle_Levelset(1).Store_Unique_Particle_Id();}
    particle_levelset_evolution->Use_Particle_Levelset(example.fluids_parameters.use_particle_levelset);

    // solid fluid coupling
    particle_levelset_evolution->particle_levelset.levelset.Set_Collision_Body_List(collision_bodies_affecting_fluid);
    particle_levelset_evolution->particle_levelset.levelset.Set_Face_Velocities_Valid_Mask(&incompressible->valid_mask);
    particle_levelset_evolution->particle_levelset.Set_Collision_Distance_Factors(example.fluids_parameters.min_collision_distance_factor,
        example.fluids_parameters.max_collision_distance_factor);

    // incompressible flow
    incompressible->Set_Custom_Boundary(*example.fluids_parameters.fluid_boundary);

    // set up the initial state
    if(example.restart){
        example.Read_Output_Files_Fluids(current_frame);
        Initialize_Fluids_Grids();
        collision_bodies_affecting_fluid.Rasterize_Objects();
        collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.Minimum_Edge_Length(),5);

        collision_bodies_affecting_fluid.Compute_Grid_Visibility(); // compute grid visibility (for advection later)
        particle_levelset_evolution->Set_Seed(2606);
        if(!example.fluids_parameters.write_particles) particle_levelset_evolution->Seed_Particles(time);
        particle_levelset_evolution->Delete_Particles_Outside_Grid();
        if(example.fluids_parameters.delete_fluid_inside_objects) Delete_Particles_Inside_Objects(time);
        example.Update_Fluid_Parameters((T)1./example.frame_rate,time);}
    else{
        Initialize_Fluids_Grids();
        collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false);
        collision_bodies_affecting_fluid.Rasterize_Objects();
        collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.Minimum_Edge_Length(),5);
        example.Initialize_Phi();
        example.Adjust_Phi_With_Sources(time);
        particle_levelset_evolution->Make_Signed_Distance();
        particle_levelset_evolution->Fill_Levelset_Ghost_Cells(time);

        collision_bodies_affecting_fluid.Compute_Grid_Visibility(); // compute grid visibility (for averaging face velocities to nodes below)
        particle_levelset_evolution->Set_Seed(2606);
        particle_levelset_evolution->Seed_Particles(time);
        particle_levelset_evolution->Delete_Particles_Outside_Grid();
        if(example.fluids_parameters.delete_fluid_inside_objects) Delete_Particles_Inside_Objects(time);

        example.Initialize_Velocities();
        example.Update_Fluid_Parameters((T)1./example.frame_rate,time);

        int extrapolation_ghost_cells=2*example.fluids_parameters.number_of_ghost_cells+2;
        T extrapolation_bandwidth=(T)(extrapolation_ghost_cells-1);
        T_ARRAYS_SCALAR exchanged_phi_ghost(grid.Domain_Indices(extrapolation_ghost_cells));
        particle_levelset_evolution->particle_levelset.levelset.boundary->Fill_Ghost_Cells(grid,particle_levelset_evolution->phi,exchanged_phi_ghost,0,time,extrapolation_ghost_cells);
//        Extrapolate_Velocity_Across_Interface(example.incompressible_fluid_container.face_velocities,particle_levelset_evolution->Particle_Levelset(1).levelset,extrapolation_bandwidth);
        if(!example.two_phase) 
            incompressible->Extrapolate_Velocity_Across_Interface(example.incompressible_fluid_container.face_velocities,exchanged_phi_ghost,
                example.fluids_parameters.enforce_divergence_free_extrapolation,extrapolation_bandwidth,0,TV(),&collision_bodies_affecting_fluid.face_neighbors_visible);
    }
}
//#####################################################################
// Function Initialize_Fluids_Grids
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Initialize_Fluids_Grids()
{
    T_GRID& grid=*example.fluids_parameters.grid;
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters=example.fluids_parameters;
    *example.fluids_parameters.grid=example.fluids_parameters.grid->Get_MAC_Grid();
    example.fluids_parameters.p_grid=*example.fluids_parameters.grid;
    example.Initialize_Fluids_Grids();
    example.incompressible_fluid_container.Initialize_Grids();
    fluids_parameters.particle_levelset_evolution->Initialize_Domain(fluids_parameters.p_grid);
    fluids_parameters.particle_levelset_evolution->particle_levelset.Set_Band_Width((T)2*fluids_parameters.particle_half_bandwidth);
    fluids_parameters.incompressible->valid_mask.Resize(grid.Domain_Indices(example.fluids_parameters.number_of_ghost_cells),true,true,true);
    fluids_parameters.incompressible->grid=grid.Get_MAC_Grid();
//    p.Resize(p_grid.Domain_Indices(1));p_save_for_projection.Resize(p_grid.Domain_Indices(1));face_velocities_save_for_projection.Resize(p_grid);
    dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(solids_evolution).Initialize_Grid_Arrays();
}
//#####################################################################
// Function First_Order_Time_Step
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
First_Order_Time_Step(int substep,T dt)
{
    FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters=example.fluids_parameters;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>& collision_bodies_affecting_fluid=*fluids_parameters.collision_bodies_affecting_fluid;
    EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities=*example.solid_body_collection.example_forces_and_velocities;
    SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>& slip=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*example.solids_evolution);
    T_GRID& grid=*fluids_parameters.grid;
    INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& incompressible_fluid_container=example.incompressible_fluid_container;
    INCOMPRESSIBLE_UNIFORM<T_GRID>* incompressible=fluids_parameters.incompressible;

    example.solid_body_collection.Print_Energy(time,1);
    collision_bodies_affecting_fluid.Rasterize_Objects(); // non-swept
    collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.Minimum_Edge_Length(),5);  // static occupied blocks
    // swept occupied blocks
    example.Initialize_Swept_Occupied_Blocks_For_Advection(dt,time,example.incompressible_fluid_container.face_velocities);

    if(example.use_pls_evolution_for_structure) Advance_Particles_With_PLS(dt);

    Write_Substep("start step",substep,1);
    Advect_Fluid(dt,substep);
    std::stringstream ss;
    ss<<"Maximum face velocity (after advect) = ("<<incompressible_fluid_container.face_velocities.Maxabs().Magnitude()<<": "<<incompressible_fluid_container.face_velocities.Maxabs()<<std::endl;
    LOG::filecout(ss.str());
    Write_Substep("advect",substep,1);

    example_forces_and_velocities.Update_Time_Varying_Material_Properties(time+dt);
    example.solid_body_collection.Update_Position_Based_State(time+dt,true);
//    example.solid_body_collection.deformable_body_collection.template Find_Force<SURFACE_TENSION_FORCE<VECTOR<T,2> >*>()->Dump_Curvatures();
    slip.two_phase=example.two_phase;
    slip.Solve(incompressible_fluid_container.face_velocities,dt,time,time+dt,false,false);
    Write_Substep("pressure solve",substep,1);

    slip.Print_Maximum_Velocities(time);

//        if(incompressible) T_ARRAYS_SCALAR::Copy(incompressible->projection.p,incompressible->projection.p_save_for_projection); // save the good pressure for later
    if(!example.use_pls_evolution_for_structure) slip.Euler_Step_Position(dt,time+dt);

    Write_Substep("euler step position",substep,1);

    LOG::Time("extrapolating velocity across interface");
//    Extrapolate_Velocity_Across_Interface(example.incompressible_fluid_container.face_velocities,particle_levelset_evolution->Particle_Levelset(1).levelset,extrapolation_bandwidth);
    if(!example.two_phase) Extrapolate_Velocity_Across_Interface(time,dt);
        //incompressible->Extrapolate_Velocity_Across_Interface(example.incompressible_fluid_container.face_velocities,exchanged_phi_ghost,
        //    fluids_parameters.enforce_divergence_free_extrapolation,extrapolation_bandwidth,0,TV(),&collision_bodies_affecting_fluid.face_neighbors_visible);
    Write_Substep("extrapolate about interface",substep,1);
    incompressible->boundary->Apply_Boundary_Condition_Face(incompressible->grid,example.incompressible_fluid_container.face_velocities,time+dt);
    std::stringstream ss1;
    ss1<<"Maximum face velocity = ("<<incompressible_fluid_container.face_velocities.Maxabs().Magnitude()<<": "<<incompressible_fluid_container.face_velocities.Maxabs()<<std::endl;
    LOG::filecout(ss1.str());

    Write_Substep("end step",substep,1);
    example.solids_evolution->time+=dt;
    time+=dt;
}
//#####################################################################
// Function Extrapolate_Velocity_Across_Interface
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Extrapolate_Velocity_Across_Interface(T time,T dt)
{
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>* particle_levelset_evolution=example.fluids_parameters.particle_levelset_evolution;
    T_GRID& grid=*example.fluids_parameters.grid;
    SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>& slip=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*example.solids_evolution);
    ARRAY<bool,FACE_INDEX<TV::m> >& valid_faces=slip.solved_faces;
    int extrapolation_ghost_cells=2*example.fluids_parameters.number_of_ghost_cells+2;
    T_ARRAYS_SCALAR phi_ghost(grid.Domain_Indices(extrapolation_ghost_cells));
    particle_levelset_evolution->particle_levelset.levelset.boundary->Fill_Ghost_Cells(grid,particle_levelset_evolution->phi,phi_ghost,0,time+dt,extrapolation_ghost_cells);
    int band_width=extrapolation_ghost_cells-1;
    T delta=band_width*grid.dX.Max();
    for(int axis=1;axis<=TV::m;axis++){
        GRID<TV> face_grid=grid.Get_Face_Grid(axis);
        ARRAY<T,TV_INT> phi_face(face_grid.Domain_Indices(),false);
        ARRAYS_ND_BASE<TV>& face_velocity=example.incompressible_fluid_container.face_velocities.Component(axis);
        ARRAYS_ND_BASE<VECTOR<bool,TV::m> >& fixed_face=valid_faces.Component(axis);
        for(FACE_ITERATOR iterator(grid,0,T_GRID::WHOLE_REGION,0,axis);iterator.Valid();iterator.Next()){
            TV_INT index=iterator.Face_Index();
            T phi1=phi_ghost(iterator.First_Cell_Index()),phi2=phi_ghost(iterator.Second_Cell_Index());
            phi_face(index)=(T).5*(phi1+phi2);
            if(phi_face(index) >= delta && !fixed_face(index)) face_velocity(index)=(T)0;}

        EXTRAPOLATION_UNIFORM<GRID<TV>,T> extrapolate(face_grid,phi_face,face_velocity,extrapolation_ghost_cells);
        extrapolate.Set_Band_Width((T)band_width);
        extrapolate.Set_Custom_Seed_Done(&fixed_face);
        extrapolate.Extrapolate();}
}
//#####################################################################
// Function Advance_To_Target_Time
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Advance_Particles_With_PLS(T dt)
{
    ARRAY_VIEW<TV> X=example.solid_body_collection.deformable_body_collection.particles.X;
    ARRAY<T,FACE_INDEX<TV::m> >& face_velocities=example.incompressible_fluid_container.face_velocities;
    LINEAR_INTERPOLATION_MAC<TV,T> interp(*example.fluids_parameters.grid);
    RUNGEKUTTA<ARRAY_VIEW<TV> >* rk=RUNGEKUTTA<ARRAY_VIEW<TV> >::Create(X,example.fluids_parameters.particle_levelset_evolution->runge_kutta_order_particles,dt,0);
    for(int k=1;k<=example.fluids_parameters.particle_levelset_evolution->runge_kutta_order_particles;k++){
        for(int p=1;p<=X.m;p++) X(p)+=dt*interp.Clamped_To_Array(face_velocities,X(p));
        rk->Main();}
    delete rk;
}
//#####################################################################
// Function Advance_To_Target_Time
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters=example.fluids_parameters;
    SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callbacks=example.solids_evolution->solids_evolution_callbacks;
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>* particle_levelset_evolution=fluids_parameters.particle_levelset_evolution;
    EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities=*example.solid_body_collection.example_forces_and_velocities;

    bool done=false;for(int substep=1;!done;substep++){
        LOG::SCOPE scope("SUBSTEP","substep %d",substep);
        solids_evolution_callbacks->Preprocess_Solids_Substep(time,substep);
        particle_levelset_evolution->Set_Number_Particles_Per_Cell(fluids_parameters.number_particles_per_cell);
        T dt=Compute_Dt(time,target_time,done);
        example.Preprocess_Substep(dt,time);
        example.Update_Fluid_Parameters(dt,time);

        First_Order_Time_Step(substep,dt);

        example_forces_and_velocities.Advance_One_Time_Step_End_Callback(dt,time);
        solids_evolution_callbacks->Postprocess_Solids_Substep(example.solids_evolution->time,substep);
        example.Postprocess_Substep(dt,time);
        Write_Substep(STRING_UTILITIES::string_sprintf("END Substep %d",substep),substep,0);}
}
//#####################################################################
// Function Advect_Fluid
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Advect_Fluid(const T dt,const int substep)
{
    FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters=example.fluids_parameters;
    INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& incompressible_fluid_container=example.incompressible_fluid_container;
    T_GRID& grid=*fluids_parameters.grid;
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>* particle_levelset_evolution=fluids_parameters.particle_levelset_evolution;
    INCOMPRESSIBLE_UNIFORM<T_GRID>* incompressible=fluids_parameters.incompressible;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>& collision_bodies_affecting_fluid=*fluids_parameters.collision_bodies_affecting_fluid;

    LOG::SCOPE scalar_scope("SCALAR SCOPE");

    // important to compute ghost velocity values for particles near domain boundaries
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities_ghost;
    ARRAY<T,FACE_INDEX<TV::m> >& face_velocities=incompressible_fluid_container.face_velocities;
    face_velocities_ghost.Resize(incompressible->grid,example.fluids_parameters.number_of_ghost_cells,false);
    incompressible->boundary->Fill_Ghost_Cells_Face(grid,face_velocities,face_velocities_ghost,time+dt,example.fluids_parameters.number_of_ghost_cells);

    fluids_parameters.phi_boundary_water.Use_Extrapolation_Mode(false);
    example.Adjust_Phi_With_Objects(time);
    LOG::Time("advecting levelset");
    particle_levelset_evolution->Advance_Levelset(dt);
    fluids_parameters.phi_boundary_water.Use_Extrapolation_Mode(true);
    example.Extrapolate_Phi_Into_Objects(time+dt);
    Write_Substep("after levelset advection",0,1);
    LOG::Time("advecting particles");
    particle_levelset_evolution->Advance_Particles(face_velocities_ghost,dt,false);
    Write_Substep("after particle advection",0,1);

    example.Scalar_Advection_Callback(dt,time);

    LOG::Time("updating removed particle velocities");
    example.Modify_Removed_Particles_Before_Advection(dt,time);
    particle_levelset_evolution->Fill_Levelset_Ghost_Cells(time);
    PARTICLE_LEVELSET_UNIFORM<T_GRID>& pls=particle_levelset_evolution->Particle_Levelset(1);
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,TV> interpolation;
    if(pls.use_removed_positive_particles) for(NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) if(pls.removed_positive_particles(iterator.Node_Index())){
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_positive_particles(iterator.Node_Index());
        for(int p=1;p<=particles.array_collection->Size();p++){
            TV X=particles.X(p),V=interpolation.Clamped_To_Array_Face(grid,face_velocities_ghost,X);
            if(-pls.levelset.Phi(X)>1.5*particles.radius(p)) V-=fluids_parameters.removed_positive_particle_buoyancy_constant*fluids_parameters.gravity_direction; // buoyancy
            particles.V(p)=V;}}
    if(pls.use_removed_negative_particles) for(NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) if(pls.removed_negative_particles(iterator.Node_Index())){
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_negative_particles(iterator.Node_Index());
        for(int p=1;p<=particles.array_collection->Size();p++) particles.V(p)+=dt*fluids_parameters.gravity*fluids_parameters.gravity_direction; // ballistic
        if(fluids_parameters.use_body_force) for(int p=1;p<=particles.array_collection->Size();p++)
            particles.V(p)+=dt*interpolation.Clamped_To_Array_Face(grid,incompressible->force,particles.X(p));} // external forces

    LOG::Time("updating velocity (explicit part)");
    Write_Substep("before advection",substep,1);

    int extrapolation_ghost_cells=2*example.fluids_parameters.number_of_ghost_cells+2;
    T extrapolation_bandwidth=(T)(extrapolation_ghost_cells-1);
    T_ARRAYS_SCALAR exchanged_phi_ghost(grid.Domain_Indices(extrapolation_ghost_cells));
    particle_levelset_evolution->particle_levelset.levelset.boundary->Fill_Ghost_Cells(grid,particle_levelset_evolution->phi,exchanged_phi_ghost,0,time+dt,extrapolation_ghost_cells);
    if(example.convection_order>1){
        T rk_time=time;
        RUNGEKUTTA<ARRAY<T,FACE_INDEX<TV::dimension> > > rk(face_velocities);
        rk.Set_Order(example.convection_order);
        rk.Set_Time(time);
        rk.Start(dt);
        for(int i=1;i<=rk.order;i++){
//            Extrapolate_Velocity_Across_Interface(face_velocities,particle_levelset_evolution->Particle_Levelset(1).levelset,extrapolation_bandwidth);
            if(!example.two_phase)
                incompressible->Extrapolate_Velocity_Across_Interface(example.incompressible_fluid_container.face_velocities,exchanged_phi_ghost,
                    fluids_parameters.enforce_divergence_free_extrapolation,extrapolation_bandwidth,0,TV(),&collision_bodies_affecting_fluid.face_neighbors_visible);
            incompressible->Advance_One_Time_Step_Convection(dt,rk_time,face_velocities,face_velocities,example.fluids_parameters.number_of_ghost_cells);
            rk_time=rk.Main();}}
    else incompressible->Advance_One_Time_Step_Convection(dt,time,face_velocities_ghost,face_velocities,example.fluids_parameters.number_of_ghost_cells);
    Write_Substep("after advection",substep,1);

    if(fluids_parameters.use_body_force) fluids_parameters.callbacks->Get_Body_Force(fluids_parameters.incompressible->force,dt,time);

    LOG::Time("updating velocity (explicit part without convection)");
    Write_Substep("before viscosity",substep,1);
    if(SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(example.solids_evolution))
        coupled_evolution->Apply_Viscosity(face_velocities,dt,time);
    Write_Substep("after viscosity",substep,1);
    incompressible->Advance_One_Time_Step_Forces(face_velocities,dt,time,fluids_parameters.implicit_viscosity,&particle_levelset_evolution->phi,example.fluids_parameters.number_of_ghost_cells);
    Write_Substep("after integrate non advection forces",substep,1);

    LOG::Time("effective velocity acceleration structures");
    // revalidate scalars and velocity in body's new position
    Write_Substep("before scalar revalidation",0,1);
    collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false); // NON-swept acceleration structures
    collision_bodies_affecting_fluid.Rasterize_Objects(); // non-swept
    collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.Minimum_Edge_Length(),5);  // static occupied blocks
    collision_bodies_affecting_fluid.Compute_Grid_Visibility(); // used in fast marching and extrapolation too... NOTE: this requires that objects_in_cell be current!
    example.Revalidate_Fluid_Scalars(); // uses visibility

    example.Revalidate_Fluid_Velocity(incompressible_fluid_container.face_velocities); // uses visibility
    Write_Substep("after scalar revalidation",0,1);

    LOG::Time("modifying levelset");
    particle_levelset_evolution->Fill_Levelset_Ghost_Cells(time+dt);
    Write_Substep("after filling level set ghost cells and exchanging overlap particles",0,1);
    particle_levelset_evolution->Modify_Levelset_And_Particles(&face_velocities_ghost);
    Write_Substep("after modify levelset and particles",0,1);
    example.Revalidate_Phi_After_Modify_Levelset(); // second revalidation -- uses visibility too
    Write_Substep("after revalidate phi",0,1);

    LOG::Time("adding sources");
    if(example.Adjust_Phi_With_Sources(time+dt)) particle_levelset_evolution->Make_Signed_Distance();
    particle_levelset_evolution->Fill_Levelset_Ghost_Cells(time+dt);
    LOG::Time("getting sources");
    T_ARRAYS_BOOL* source_mask=0;example.Get_Source_Reseed_Mask(source_mask,time+dt);
    if(source_mask){LOG::Time("reseeding sources");particle_levelset_evolution->Reseed_Particles(time+dt,0,source_mask);delete source_mask;}
    Write_Substep("after adding sources",0,1);

    LOG::Time("deleting particles"); // needs to be after adding sources, since it does a reseed
    particle_levelset_evolution->Delete_Particles_Outside_Grid();
    if(fluids_parameters.delete_fluid_inside_objects) Delete_Particles_Inside_Objects(time+dt);
    LOG::Time("deleting particles in local maxima");
    particle_levelset_evolution->particle_levelset.Delete_Particles_In_Local_Maximum_Phi_Cells(1);
    LOG::Time("deleting particles far from interface");
    particle_levelset_evolution->Particle_Levelset(1).Delete_Particles_Far_From_Interface(); // uses visibility
    Write_Substep("after delete particles far from interface",0,1);

    LOG::Time("re-incorporating removed particles");
    // TODO: if your particles fall entirely within the grid this shouldn't need ghost cells, but it may need them for MPI
    particle_levelset_evolution->Particle_Levelset(1).Identify_And_Remove_Escaped_Particles(face_velocities_ghost,1.5,time+dt);
    if(particle_levelset_evolution->Particle_Levelset(1).use_removed_positive_particles || particle_levelset_evolution->Particle_Levelset(1).use_removed_negative_particles)
        particle_levelset_evolution->Particle_Levelset(1).Reincorporate_Removed_Particles(1,fluids_parameters.removed_particle_mass_scaling,
            fluids_parameters.reincorporate_removed_particle_velocity?&example.incompressible_fluid_container.face_velocities:0,true);

    // update ghost phi values
    particle_levelset_evolution->Fill_Levelset_Ghost_Cells(time+dt);

    Write_Substep("after scalar update",substep,1);
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Postprocess_Frame(const int frame)
{
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>* particle_levelset_evolution=example.fluids_parameters.particle_levelset_evolution;

    example.Postprocess_Phi(time);

    if(particle_levelset_evolution->use_particle_levelset && (frame-example.first_frame)%example.fluids_parameters.reseeding_frame_rate==0){
        LOG::Time("Reseeding... ");
        particle_levelset_evolution->Reseed_Particles(time);
        particle_levelset_evolution->Delete_Particles_Outside_Grid();
        if(example.fluids_parameters.delete_fluid_inside_objects) Delete_Particles_Inside_Objects(time);
        LOG::Stop_Time();}

    example.Postprocess_Frame(frame);
    example.solids_evolution->Postprocess_Frame(frame);
}
//#####################################################################
// Function Compute_Dt
//#####################################################################
template<class TV> typename TV::SCALAR PLS_FSI_DRIVER<TV>::
Compute_Dt(const T time,const T target_time,bool& done)
{
    SOLIDS_PARAMETERS<TV>& solids_parameters=example.solids_parameters;
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters=example.fluids_parameters;
    SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callbacks=solids_evolution.solids_evolution_callbacks;

    T fluids_dt=FLT_MAX;
    {T dt_levelset=fluids_parameters.particle_levelset_evolution->CFL(true,false);fluids_dt=min(fluids_dt,dt_levelset);}
    std::stringstream ss;
    ss<<"before example cfl clamping:  "<<fluids_dt<<std::endl;
    LOG::filecout(ss.str());
    example.Limit_Dt(fluids_dt,time);

    // solids dt
    T solids_dt=FLT_MAX;
    if(solids_evolution.Use_CFL()) solids_dt=min(solids_dt,example.solid_body_collection.CFL(solids_parameters.verbose_dt));
    solids_dt=min(solids_dt,solids_evolution_callbacks->Constraints_CFL());
    solids_dt=min(solids_dt,example.solid_body_collection.rigid_body_collection.CFL_Rigid(solids_parameters.rigid_body_evolution_parameters,solids_parameters.verbose_dt));
    solids_evolution_callbacks->Limit_Solids_Dt(solids_dt,time);

    if(example.fixed_dt){fluids_dt=example.fixed_dt;solids_dt=example.fixed_dt;}
    T dt=min(fluids_dt,solids_dt);
    std::stringstream ss1;
    ss1<<"fluids_dt = "<<fluids_dt<<", solids_dt = "<<solids_dt<<" dt="<<dt<<std::endl;
    LOG::filecout(ss1.str());
    done=false;
    PLS_FSI_EXAMPLE<TV>::Clamp_Time_Step_With_Target_Time(time,target_time,dt,done,solids_parameters.min_dt);
    return dt;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Write_Output_Files(const int frame)
{
    LOG::SCOPE scope("writing output files");
    FILE_UTILITIES::Create_Directory(example.output_directory);
    FILE_UTILITIES::Create_Directory(example.output_directory+STRING_UTILITIES::string_sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
    Write_First_Frame(frame);

    example.fluids_parameters.phi_boundary_water.Use_Extrapolation_Mode(false);
    example.Write_Output_Files(frame);
    example.fluids_parameters.phi_boundary_water.Use_Extrapolation_Mode(true);

    T_GRID& grid=*example.fluids_parameters.grid;
    int number_of_positive_particles=0,number_of_negative_particles=0,number_of_removed_positive_particles=0,number_of_removed_negative_particles=0;
    PARTICLE_LEVELSET_UNIFORM<T_GRID>* pls=0;
    pls=&example.fluids_parameters.particle_levelset_evolution->particle_levelset;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) if(pls->positive_particles(iterator.Cell_Index()))
        number_of_positive_particles+=pls->positive_particles(iterator.Cell_Index())->array_collection->Size();
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) if(pls->negative_particles(iterator.Cell_Index()))
        number_of_negative_particles+=pls->negative_particles(iterator.Cell_Index())->array_collection->Size();
    std::stringstream ss;
    ss<<number_of_positive_particles<<" positive and "<<number_of_negative_particles<<" negative particles "<<std::endl;
    if(pls->use_removed_positive_particles)
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) if(pls->removed_positive_particles(iterator.Cell_Index()))
            number_of_removed_positive_particles+=pls->removed_positive_particles(iterator.Cell_Index())->array_collection->Size();
    if(pls->use_removed_negative_particles)
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) if(pls->removed_negative_particles(iterator.Cell_Index()))
            number_of_removed_negative_particles+=pls->removed_negative_particles(iterator.Cell_Index())->array_collection->Size();
    ss<<number_of_removed_positive_particles<<" positive and "<<number_of_removed_negative_particles<<" negative removed particles "<<std::endl;
    LOG::filecout(ss.str());

    Write_Time(frame);
    Write_Last_Frame(frame);
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Delete_Particles_Inside_Objects(const T time)
{
    PARTICLE_LEVELSET_UNIFORM<T_GRID>* particle_levelset=&example.fluids_parameters.particle_levelset_evolution->particle_levelset;
    Delete_Particles_Inside_Objects<PARTICLE_LEVELSET_PARTICLES<TV> >(particle_levelset->positive_particles,PARTICLE_LEVELSET_POSITIVE,time);
    Delete_Particles_Inside_Objects<PARTICLE_LEVELSET_PARTICLES<TV> >(particle_levelset->negative_particles,PARTICLE_LEVELSET_NEGATIVE,time);
    if(particle_levelset->use_removed_positive_particles)
        Delete_Particles_Inside_Objects<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >(particle_levelset->removed_positive_particles,PARTICLE_LEVELSET_REMOVED_POSITIVE,time);
    if(particle_levelset->use_removed_negative_particles)
        Delete_Particles_Inside_Objects<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >(particle_levelset->removed_negative_particles,PARTICLE_LEVELSET_REMOVED_NEGATIVE,time);
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
template<class TV> template<class T_PARTICLES> void PLS_FSI_DRIVER<TV>::
Delete_Particles_Inside_Objects(ARRAY<T_PARTICLES*,TV_INT>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time)
{
    for(NODE_ITERATOR iterator(*example.fluids_parameters.grid);iterator.Valid();iterator.Next()){TV_INT block_index=iterator.Node_Index();if(particles(block_index)){
        BLOCK_UNIFORM<T_GRID> block(*example.fluids_parameters.grid,block_index);
        COLLISION_GEOMETRY_ID body_id;int aggregate_id;
        T_PARTICLES& block_particles=*particles(block_index);
        if(example.fluids_parameters.collision_bodies_affecting_fluid->Occupied_Block(block)){
            for(int k=block_particles.array_collection->Size();k>=1;k--)
                if(example.fluids_parameters.collision_bodies_affecting_fluid->Inside_Any_Simplex_Of_Any_Body(block_particles.X(k),body_id,aggregate_id))
                    block_particles.array_collection->Delete_Element(k);}
        example.fluids_parameters.callbacks->Delete_Particles_Inside_Objects(block_particles,particle_type,time);}}
}
//#####################################################################
// Function Extrapolate_Velocity_Across_Interface
//#####################################################################
template<class TV> void PLS_FSI_DRIVER<TV>::
Extrapolate_Velocity_Across_Interface(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T_FAST_LEVELSET& phi,const T band_width)
{
    T_GRID& grid=*example.fluids_parameters.grid;
    EXTRAPOLATION_HIGHER_ORDER<TV,T>::Extrapolate_Face(grid,phi,(int)ceil(band_width)+1,face_velocities,20,example.fluids_parameters.number_of_ghost_cells,band_width);
}
//#####################################################################
template class PLS_FSI_DRIVER<VECTOR<float,1> >;
template class PLS_FSI_DRIVER<VECTOR<float,2> >;
template class PLS_FSI_DRIVER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PLS_FSI_DRIVER<VECTOR<double,1> >;
template class PLS_FSI_DRIVER<VECTOR<double,2> >;
template class PLS_FSI_DRIVER<VECTOR<double,3> >;
#endif

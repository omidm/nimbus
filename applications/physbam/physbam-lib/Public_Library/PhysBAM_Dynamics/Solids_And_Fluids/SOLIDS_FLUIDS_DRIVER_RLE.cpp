#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_DRIVER_RLE
//#####################################################################
#include <PhysBAM_Tools/Grids_RLE_Boundaries/BOUNDARY_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/AVERAGING_RLE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_RLE_GRID.h>
#include <PhysBAM_Geometry/Grids_RLE_Collisions/GRID_BASED_COLLISION_GEOMETRY_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_RLE.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_RLE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> SOLIDS_FLUIDS_DRIVER_RLE<T_GRID>::
SOLIDS_FLUIDS_DRIVER_RLE(SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>& example_input)
    :SOLIDS_FLUIDS_DRIVER<TV>(example_input),example(example_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> SOLIDS_FLUIDS_DRIVER_RLE<T_GRID>::
~SOLIDS_FLUIDS_DRIVER_RLE()
{}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_RLE<T_GRID>::
Initialize()
{
    FLUIDS_PARAMETERS<T_GRID>& fluids_parameters=example.fluids_parameters;
    PARTICLE_LEVELSET_RLE<T_GRID>& particle_levelset=example.particle_levelset;
    INCOMPRESSIBLE_RLE<T_GRID>& incompressible=example.incompressible;
    DEEP_WATER_EVOLUTION<TV_HORIZONTAL>& deep_water=example.deep_water;
    GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>& collision_bodies_affecting_fluid=*fluids_parameters.collision_bodies_affecting_fluid;
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_parameters.solids_evolution;
    T_GRID& grid=*fluids_parameters.grid;

    SOLIDS_FLUIDS_DRIVER<TV>::Initialize();
    fluids_parameters.Set_Fluids_Parameters_Callbacks(example);
    solids_evolution.Set_Solids_Evolution_Callbacks(example);

    // initialize MPI
    if(example.mpi_grid) example.Initialize_MPI();

    solids_evolution.time=time;
    solids_evolution.Initialize_Deformable_Objects(example.frame_rate,example.restart);

    if(example.restart){
        example.Read_Output_Files_Solids(example.restart_frame);
        solids_evolution.time=time=example.initial_time+(example.restart_frame-example.first_frame)/example.frame_rate;}

    solids_evolution.Initialize_Rigid_Bodies(example.frame_rate,example.restart);

    // sets up the proper wall states
    fluids_parameters.Initialize_Domain_Boundary_Conditions();

    // initialize grid and levelset
    particle_levelset.Set_Number_Particles_Per_Cell(fluids_parameters.number_particles_per_cell);
    particle_levelset.levelset.Set_Levelset_Callbacks(example);
    particle_levelset.levelset.Initialize_FMM_Initialization_Iterative_Solver(fluids_parameters.refine_fmm_initialization_with_iterative_solver);
    if(fluids_parameters.phi_boundary) particle_levelset.levelset.Set_Custom_Boundary(*fluids_parameters.phi_boundary);
    particle_levelset.Bias_Towards_Negative_Particles(fluids_parameters.bias_towards_negative_particles);
    if(fluids_parameters.use_removed_positive_particles) particle_levelset.Use_Removed_Positive_Particles();
    if(fluids_parameters.use_removed_negative_particles) particle_levelset.Use_Removed_Negative_Particles();
    if(fluids_parameters.store_particle_ids) particle_levelset.Store_Unique_Particle_Id();
    particle_levelset.cfl_number=fluids_parameters.cfl;

    // solid fluid coupling
    particle_levelset.levelset.Set_Collision_Body_List(collision_bodies_affecting_fluid);
    particle_levelset.levelset.Set_Face_Velocities_Valid_Mask(&incompressible.valid_mask);
    particle_levelset.Set_Collision_Distance_Factors(fluids_parameters.min_collision_distance_factor,fluids_parameters.max_collision_distance_factor);

    // incompressible flow
    incompressible.Set_Custom_Boundary(*fluids_parameters.fluid_boundary);
    incompressible.projection.laplace.Set_Relative_Tolerance(fluids_parameters.incompressible_tolerance);
    incompressible.projection.laplace.pcg.Set_Maximum_Iterations(fluids_parameters.incompressible_iterations);
    incompressible.projection.laplace.pcg.Show_Results();
    incompressible.projection.laplace.pcg.Use_Incomplete_Cholesky(); // 3d definitely needs this, but 2d might work with modified incomplete cholesky

    // deep water
    if(example.use_deep_water){
        deep_water.grid=grid.horizontal_grid.Get_MAC_Grid();
        deep_water.Use_Surface_Pressure();
        deep_water.Initialize();}

    // set up the initial state
    if(example.restart){
        example.Read_Output_Files_Fluids(current_frame);
        example.Initialize_Ground();
        Initialize_Grids();
        collision_bodies_affecting_fluid.Rasterize_Objects();
        collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.Minimum_Edge_Length(),5);
        if(!fluids_parameters.write_particles){
            particle_levelset.Seed_Particles(time);
            example.particle_levelset.Delete_Particles_Outside_Grid();
            if(fluids_parameters.delete_fluid_inside_objects) example.Delete_Particles_Inside_Objects(time);}}
    else{
        example.Initialize_Ground();
        example.Initialize_Grid();
        Initialize_Grids();
        example.Initialize_Phi();
        particle_levelset.levelset.Fast_Marching_Method(time);
        collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false);
        collision_bodies_affecting_fluid.Rasterize_Objects();
        Rebuild_Grid(time,0);
        collision_bodies_affecting_fluid.Rasterize_Objects();
        collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.Minimum_Edge_Length(),5);
        example.Adjust_Phi_With_Sources(time);
        particle_levelset.levelset.Fast_Marching_Method(time);
        collision_bodies_affecting_fluid.Compute_Grid_Visibility(); // compute grid visibility (for averaging face velocities to nodes below)
        particle_levelset.random.Set_Seed(2606);
        particle_levelset.Seed_Particles(time);
        example.particle_levelset.Delete_Particles_Outside_Grid();
        if(fluids_parameters.delete_fluid_inside_objects) example.Delete_Particles_Inside_Objects(time);
        example.Initialize_Velocities();
        example.Update_Fluid_Parameters((T)1/example.frame_rate,time);
        particle_levelset.levelset.boundary->Fill_Ghost_Cells_Cell(grid,particle_levelset.phi,particle_levelset.phi,time);
        incompressible.Extrapolate_Velocity_Across_Interface(particle_levelset.phi,fluids_parameters.enforce_divergence_free_extrapolation,3,
            &collision_bodies_affecting_fluid.face_neighbors_visible);}

    if(fluids_parameters.monitor_mass){
        fluids_parameters.mass=particle_levelset.levelset.Approximate_Negative_Size();
        LOG::cout<<"Initial Size="<<fluids_parameters.mass<<std::endl;}

    collision_bodies_affecting_fluid.Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,time);      
}
//#####################################################################
// Function Initialize_Grids
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_RLE<T_GRID>::
Initialize_Grids()
{
    FLUIDS_PARAMETERS<T_GRID>& fluids_parameters=example.fluids_parameters;
    example.particle_levelset.Initialize_Particle_Levelset_Grid_Values();
    if(example.levelset_advection.semi_lagrangian_collidable || example.levelset_advection.semi_lagrangian_collidable_slip)
        example.particle_levelset.levelset.Initialize_Valid_Masks(example.particle_levelset.grid);
    example.particle_levelset.Set_Band_Width((T)2*fluids_parameters.particle_half_bandwidth);
    example.incompressible.Initialize_Grids();
    fluids_parameters.collision_bodies_affecting_fluid->Initialize_Grids();
    example.Initialize_Grids_Extra();
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_RLE<T_GRID>::
Postprocess_Frame(const int frame)
{
    FLUIDS_PARAMETERS<T_GRID>& fluids_parameters=example.fluids_parameters;
    PARTICLE_LEVELSET_RLE<T_GRID>& particle_levelset=example.particle_levelset;
    example.Postprocess_Phi(time); // actually changes phi !!!

    if((frame-example.first_frame)%fluids_parameters.reseeding_frame_rate == 0){
        LOG::Time("Reseeding... ");
        particle_levelset.Reseed_Particles(time);
        example.particle_levelset.Delete_Particles_Outside_Grid();
        if(fluids_parameters.delete_fluid_inside_objects) example.Delete_Particles_Inside_Objects(time);
        LOG::Stop_Time();}

    if(fluids_parameters.monitor_mass){
        T mass_new=particle_levelset.levelset.Approximate_Negative_Size();
        LOG::cout<<"Size = "<<mass_new<<" - change = "<<mass_new-fluids_parameters.mass<<std::endl;
        fluids_parameters.mass=mass_new;}

    SOLIDS_FLUIDS_DRIVER<TV>::Postprocess_Frame(frame);
}
//#####################################################################
// Function Advance_To_Target_Time
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_RLE<T_GRID>::
Advance_To_Target_Time(const T target_time)
{
    typedef typename T_GRID::VECTOR_HORIZONTAL_INT TV_HORIZONTAL_INT;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::HORIZONTAL_CELL_ITERATOR HORIZONTAL_CELL_ITERATOR;

    FLUIDS_PARAMETERS<T_GRID>& fluids_parameters=example.fluids_parameters;
    PARTICLE_LEVELSET_RLE<T_GRID>& particle_levelset=example.particle_levelset;
    INCOMPRESSIBLE_RLE<T_GRID>& incompressible=example.incompressible;
    GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>& collision_bodies_affecting_fluid=*fluids_parameters.collision_bodies_affecting_fluid;
//    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_parameters.solids_evolution;
    T_GRID& grid=*fluids_parameters.grid;

    bool done=false;for(int substep=1;!done;substep++){
        LOG::SCOPE scope("SUBSTEP","substep %d",substep);

        T dt=next_dt;done=next_done;
        if(example.abort_when_dt_below && dt<example.abort_when_dt_below) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("ABORTING BECAUSE dt < %g",example.abort_when_dt_below));

        particle_levelset.Set_Number_Particles_Per_Cell(fluids_parameters.number_particles_per_cell);
        example.Update_Fluid_Parameters(dt,time);
        //fluids_parameters.phi_boundary_water.Use_Extrapolation_Mode(false); // TODO: uncomment this once rle phi boundary water is written

        {LOG::SCOPE scope("SCALAR UPDATE","scalar update");

        LOG::Time("rasterize objects");
        if(substep==1) collision_bodies_affecting_fluid.Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);
        collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(true,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE,
            COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);
        if(substep==1) collision_bodies_affecting_fluid.Rasterize_Objects();
        if(substep==1) collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.Minimum_Edge_Length(),5);  // static occupied blocks
        LOG::Time("initializing swept occupied blocks");
        example.Initialize_Swept_Occupied_Blocks_For_Advection(dt,time,true); // computes swept occupied blocks

        example.Scalar_Advection_Callback(dt,time);

        if(example.use_deep_water){LOG::Time("deep water evolution");
            for(HORIZONTAL_CELL_ITERATOR iterator(grid.horizontal_grid);iterator.Valid();iterator.Next()){TV_HORIZONTAL_INT column=iterator.Cell_Index();
                CELL_ITERATOR cell(grid,column.Insert(example.ground_j(column),2));
                example.deep_water.surface_pressure(column)=incompressible.projection.p(cell.Cell());}
            example.deep_water.Advance_Height(dt);}

        // important to compute ghost velocity values for particles near domain boundaries
        ARRAY<T> V_ghost(grid.number_of_faces,false);
        incompressible.boundary->Fill_Ghost_Cells_Face(grid,incompressible.V,V_ghost,time);

        LOG::Time("updating removed particle velocities");
        example.Modify_Removed_Particles_Before_Advection(dt,time);
        particle_levelset.levelset.boundary->Fill_Ghost_Cells(grid,particle_levelset.phi,particle_levelset.phi,dt,time);
        if(particle_levelset.use_removed_positive_particles){
            LINEAR_INTERPOLATION_RLE<T_GRID,TV> interpolation;
            for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();if(particle_levelset.removed_positive_particles(b)){
                PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*particle_levelset.removed_positive_particles(b);
                for(int p=1;p<=particles.array_collection->Size();p++){
                    TV X=particles.X(p),V=interpolation.From_Block_Face(block,V_ghost,X);
                    if(-particle_levelset.levelset.Phi(block,X) > 1.5*particles.radius(p)) V.y+=fluids_parameters.removed_positive_particle_buoyancy_constant; // buoyancy
                    particles.V(p)=V;}}}
            if(particle_levelset.use_removed_positive_particles_in_long_cells) PHYSBAM_NOT_IMPLEMENTED();}
        if(particle_levelset.use_removed_negative_particles){
            TV gravity_impulse=dt*fluids_parameters.gravity*fluids_parameters.gravity_direction;
            if(fluids_parameters.use_body_force) PHYSBAM_NOT_IMPLEMENTED();
            for(int b=1;b<=grid.number_of_blocks;b++)if(particle_levelset.removed_negative_particles(b)){
                PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*particle_levelset.removed_negative_particles(b);
                for(int p=1;p<=particles.array_collection->Size();p++) particles.V(p)+=gravity_impulse;} // ballistic
            if(particle_levelset.use_removed_negative_particles_in_long_cells)
                for(int p=1;p<=particle_levelset.removed_negative_particles_in_long_cells->array_collection->Size();p++)particle_levelset.removed_negative_particles_in_long_cells->V(p)+=gravity_impulse;}

        LOG::Time("advecting levelset");
        example.levelset_advection.Euler_Step(incompressible.V,dt,time);
        Write_Substep("after levelset advection",0,1);
        LOG::Time("advecting particles");
        particle_levelset.Euler_Step_Particles(V_ghost,dt,time,true,true);
        Write_Substep("after particle advection",0,1);

        {LOG::SCOPE scope("EXPLICIT","updating explicit part");
        incompressible.Advance_One_Time_Step_Explicit_Part(dt,time,1,fluids_parameters.implicit_viscosity,0);
        time+=dt;
        Write_Substep("after explicit part",substep,0);}

        // rebuild grid
        Rebuild_Grid(time+dt,&V_ghost);

        // revalidate scalars and velocity in body's new position
        collision_bodies_affecting_fluid.Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE);
        collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false); // NON-swept acceleration structures
        collision_bodies_affecting_fluid.Rasterize_Objects();
        collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.Minimum_Edge_Length(),5);  // static occupied blocks
        collision_bodies_affecting_fluid.Compute_Grid_Visibility(); // used in fast marching and extrapolation too... NOTE: this requires that objects_in_cell be current!
        example.Revalidate_Fluid_Scalars(); // uses visibility
        example.Revalidate_Fluid_Velocity(); // uses visibility
        Write_Substep("after scalar revalidation",0,1);

        example.Extrapolate_Phi_Into_Objects(time);

        LOG::Time("modifying levelset");
        particle_levelset.Modify_Levelset_Using_Escaped_Particles();
        Write_Substep("after first modify",0,3);
        LOG::Time("reinitializing levelset");
        particle_levelset.levelset.Fast_Marching_Method(time);
        Write_Substep("after fast march",0,3);
        LOG::Time("modifying levelset");
        particle_levelset.Modify_Levelset_Using_Escaped_Particles();
        Write_Substep("after second modify",0,3);
        LOG::Time("compacting particles");
        particle_levelset.Compact_Particles_Into_Single_Particle_Bin();
        LOG::Time("adjusting particle radii");
        particle_levelset.Adjust_Particle_Radii();
        example.Revalidate_Phi_After_Modify_Levelset(); // second revalidation -- uses visibility too
        Write_Substep("after revalidate phi",0,1);

        LOG::Time("adding sources");
        if(example.Adjust_Phi_With_Sources(time)) particle_levelset.levelset.Fast_Marching_Method(time);
        particle_levelset.levelset.boundary->Fill_Ghost_Cells(grid,particle_levelset.phi,particle_levelset.phi,dt,time);
        LOG::Time("getting sources");
        ARRAY<bool>* source_mask=0;example.Get_Source_Reseed_Mask(source_mask,time);
        if(source_mask){LOG::Time("reseeding sources");particle_levelset.Reseed_Particles(time,source_mask);delete source_mask;}
        Write_Substep("after adding sources",0,1);

        LOG::Time("deleting particles"); // needs to be after adding sources, since it does a reseed
        particle_levelset.Delete_Particles_Outside_Grid();
        if(fluids_parameters.delete_fluid_inside_objects) example.Delete_Particles_Inside_Objects(time);
        LOG::Time("deleting particles in local maxima");
        particle_levelset.Delete_Particles_In_Local_Maximum_Phi_Cells();
        LOG::Time("deleting particles far from interface");
        particle_levelset.Delete_Particles_Far_From_Interface(); // uses visibility

        LOG::Time("reincorporating removed particles");
        example.Modify_Removed_Particles_Before_Reincorporation(dt,time);
        particle_levelset.Identify_And_Remove_Escaped_Particles(V_ghost,1.5);
        particle_levelset.Reincorporate_Removed_Particles(1);
        example.Modify_Removed_Particles_After_Reincorporation(dt,time);

        // update ghost phi values
        particle_levelset.levelset.boundary->Fill_Ghost_Cells(grid,particle_levelset.phi,particle_levelset.phi,dt,time);

        LOG::cout<<"substep = "<<substep<<", dt = "<<dt<<std::endl;
        Write_Substep("after scalar update",substep,0);}

        // finish fluid update
        time-=dt; // rewinding to time n
        Compute_Next_Dt(time+dt,done,next_dt,next_done);

        // update bodies (using next_dt)
        collision_bodies_affecting_fluid.Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,time+dt);
        PHYSBAM_FATAL_ERROR("SOLIDS_EVOLUTION::Advance_To_Target_Time no longer exists.");
        //solids_evolution.Advance_To_Target_Time(time+dt+next_dt);
        collision_bodies_affecting_fluid.Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE,time+dt+next_dt);
        collision_bodies_affecting_fluid.Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);
        collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false); // NON-swept acceleration structures

        // resize u_interface here in case we want nonzero Dirichlet boundary conditions
        if(fluids_parameters.second_order_cut_cell_method) incompressible.projection.laplace.Set_Up_Second_Order_Cut_Cell_Method();

        LOG::Time("getting boundary conditions");
        example.Get_Neumann_And_Dirichlet_Boundary_Conditions(dt,time+dt);

        incompressible.projection.p*=dt; // RESCALE PRESSURE FOR A BETTER INITIAL GUESS!
        if(fluids_parameters.surface_tension || fluids_parameters.variable_surface_tension) PHYSBAM_NOT_IMPLEMENTED();
        LOG::Time("solving for the pressure and viscosity");
        Write_Substep("before laplace solve",substep,0);
        incompressible.Advance_One_Time_Step_Implicit_Part(dt,time,fluids_parameters.implicit_viscosity,&particle_levelset.phi);
        Write_Substep("after laplace solve",substep,0);
        incompressible.projection.p*=1/dt; // scale pressure back to get a real pressure
        incompressible.boundary->Apply_Boundary_Condition_Face(grid,incompressible.V,time+dt);
        if(fluids_parameters.second_order_cut_cell_method) incompressible.projection.laplace.Set_Up_Second_Order_Cut_Cell_Method(false);
        time+=dt;

        LOG::Time("extrapolating velocity across interface");
        incompressible.Extrapolate_Velocity_Across_Interface(particle_levelset.phi,fluids_parameters.enforce_divergence_free_extrapolation,3,
            &collision_bodies_affecting_fluid.face_neighbors_visible);

        Write_Substep("end of iteration",substep,0);}
}
//#####################################################################
// Function Compute_Next_Dt
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_RLE<T_GRID>::
Compute_Next_Dt(const T next_time,const bool done_this_frame,T& next_dt,bool& next_done)
{
    PHYSBAM_FATAL_ERROR("DRIVERS NOW WORK OFF OF DT, NOT NEXT DT");
    FLUIDS_PARAMETERS<T_GRID>& fluids_parameters=example.fluids_parameters;
    if(!example.use_incompressible_cfl) next_dt=fluids_parameters.cfl*example.particle_levelset.levelset.CFL(example.incompressible.V);
    else next_dt=fluids_parameters.cfl*example.incompressible.CFL();
    example.Limit_Dt(next_dt,next_time);
    if(example.mpi_grid) example.mpi_grid->Synchronize_Dt(next_dt);
    if(example.abort_when_dt_below && next_dt<example.abort_when_dt_below) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("ABORTING BECAUSE dt < %g",example.abort_when_dt_below));
    next_done=false;
    if(done_this_frame) SOLIDS_FLUIDS_EXAMPLE<TV>::Clamp_Time_Step_With_Target_Time(next_time,example.Time_At_Frame(current_frame+2),next_dt,next_done);
    else SOLIDS_FLUIDS_EXAMPLE<TV>::Clamp_Time_Step_With_Target_Time(next_time,example.Time_At_Frame(current_frame+1),next_dt,next_done);
    LOG::cout<<"next_dt = "<<next_dt<<", next_done = "<<next_done<<std::endl;
}
//#####################################################################
// Function Rebuild_Grid
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_RLE<T_GRID>::
Rebuild_Grid(const T time,ARRAY<T>* V_ghost)
{
    FLUIDS_PARAMETERS<T_GRID>& fluids_parameters=example.fluids_parameters;
    PARTICLE_LEVELSET_RLE<T_GRID>& particle_levelset=example.particle_levelset;
    INCOMPRESSIBLE_RLE<T_GRID>& incompressible=example.incompressible;
    T_GRID& grid=*fluids_parameters.grid;
    // build new grid
    LOG::Time("rebuilding grid");
    ARRAY<bool> cell_should_be_long(grid.number_of_cells);
    example.Get_Cell_Should_Be_Long(cell_should_be_long,time);
    T_GRID new_grid;
    T_GRID::Rebuild(grid,new_grid,cell_should_be_long,(int)ceil(fluids_parameters.cfl),&example.ground_j);
    example.Modify_Grid_After_Rebuild(new_grid,time);
    LOG::Stop_Time();
    LOG::cout<<"new grid: cells = "<<new_grid.number_of_cells<<", blocks = "<<new_grid.number_of_blocks<<", faces = "<<new_grid.number_of_faces<<std::endl;
    // transfer data to new grid
    LOG::Time("transferring data to new grid");
    particle_levelset.levelset.Transfer_Phi(new_grid);
    particle_levelset.Transfer_Particles(new_grid);
    incompressible.projection.Transfer_Pressure(new_grid);
    incompressible.Transfer_Velocity(new_grid);
    if(V_ghost) incompressible.Transfer_Velocity_Ghost(new_grid,*V_ghost);
    example.Transfer_Extra_State(new_grid);
    // replace old grid with new grid and adjust vertical space
    LOG::Time("finishing grid transfer");
    T_GRID::Transfer(new_grid,grid);
    grid.Adjust_Vertical_Space();
    Initialize_Grids();
    LOG::Stop_Time();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_RLE<T_GRID>::
Write_Output_Files(const int frame)
{
    LOG::SCOPE scope("WRITING OUTPUT FILES","writing output files");
    FILE_UTILITIES::Create_Directory(example.output_directory);
    FILE_UTILITIES::Create_Directory(example.output_directory+STRING_UTILITIES::string_sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
    Write_First_Frame(frame);
    example.Write_Output_Files(frame);
    Write_Time(frame);
    Write_Last_Frame(frame);
}
//#####################################################################
template class SOLIDS_FLUIDS_DRIVER_RLE<RLE_GRID_2D<float> >;
template class SOLIDS_FLUIDS_DRIVER_RLE<RLE_GRID_3D<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SOLIDS_FLUIDS_DRIVER_RLE<RLE_GRID_2D<double> >;
template class SOLIDS_FLUIDS_DRIVER_RLE<RLE_GRID_3D<double> >;
#endif
#endif 

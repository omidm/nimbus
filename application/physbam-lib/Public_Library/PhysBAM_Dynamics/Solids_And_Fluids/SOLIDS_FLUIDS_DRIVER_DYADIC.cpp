#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2007, Frank Losasso, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_DRIVER_DYADIC
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_OCTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_QUADTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_CELL_HELPER.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_CELL_HELPER.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Collisions/GRID_BASED_COLLISION_GEOMETRY_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_DYADIC.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Boundaries/BOUNDARY_PHI_WATER.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_DYADIC.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_DYADIC.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> SOLIDS_FLUIDS_DRIVER_DYADIC<T_GRID>::
SOLIDS_FLUIDS_DRIVER_DYADIC(SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>& example_input)
    :SOLIDS_FLUIDS_DRIVER<TV>(example_input),example(example_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> SOLIDS_FLUIDS_DRIVER_DYADIC<T_GRID>::
~SOLIDS_FLUIDS_DRIVER_DYADIC()
{}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_DYADIC<T_GRID>::
Initialize()
{
    T_GRID& grid=*example.fluids_parameters.grid;
    T_INCOMPRESSIBLE* incompressible=&example.fluids_parameters.incompressible;
    PARTICLE_LEVELSET_EVOLUTION_DYADIC<T_GRID>* particle_levelset_evolution=&example.fluids_parameters.particle_levelset_evolution;
    GRID_BASED_COLLISION_GEOMETRY_DYADIC<T_GRID>& collision_bodies_affecting_fluid=*example.fluids_parameters.collision_bodies_affecting_fluid;
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_parameters.solids_evolution;

    // fake what we are simulating with number of regions (matches the uniform code better). 0=smoke, 1=water, 2=fire
    int number_of_regions=0;
    if(example.fluids_parameters.water) number_of_regions=1;
    else if(example.fluids_parameters.fire) number_of_regions=2;

    SOLIDS_FLUIDS_DRIVER<TV>::Initialize();
    example.fluids_parameters.Set_Fluids_Parameters_Callbacks(example);
    solids_evolution.Set_Solids_Evolution_Callbacks(example);
    Initialize_Fluids_Grids(); // this needs to be here because the arrays have to be resized for multiphase

    example.Initialize_Bodies();
    solids_evolution.time=time;
    solids_evolution.Initialize_Deformable_Objects(example.frame_rate,example.restart);
    example.Parse_Late_Options();

    if(T_GRID::dimension==3 && !example.solids_parameters.fracture && !example.use_melting) // fracture and melting initialize collisions in Initialize_Bodies
        example.solids_parameters.Initialize_Triangle_Collisions();

    if(example.restart){
        LOG::SCOPE scope("Reading solids data","Reading solids data");
        example.Read_Output_Files_Solids(example.restart_frame);
        solids_evolution.time=time=example.Time_At_Frame(example.restart_frame);}

    solids_evolution.Initialize_Rigid_Bodies(example.frame_rate,example.restart);

    // NOTHING MORE TO DO FOR SOLIDS SIMS
    if(!example.fluids_parameters.smoke && !example.fluids_parameters.fire && !example.fluids_parameters.water) return;

    // fire smoke and water stuff...
    if(example.fluids_parameters.fire) assert(example.fluids_parameters.use_density && example.fluids_parameters.use_temperature); // both assumed to be defined
    if(example.fluids_parameters.solid_affects_fluid && example.fluids_parameters.fluid_affects_solid){
        assert(!example.fluids_parameters.second_order_cut_cell_method);assert(!example.fluids_parameters.surface_tension&&!example.fluids_parameters.variable_surface_tension);
        assert(!example.fluids_parameters.viscosity && !example.fluids_parameters.variable_viscosity);assert(example.fluids_parameters.smoke);}

    // time
    if(number_of_regions>=1){
        particle_levelset_evolution->Set_Time(time);
        particle_levelset_evolution->Set_CFL_Number(example.fluids_parameters.cfl);}

    // sets up the proper wall states
    example.fluids_parameters.Initialize_Domain_Boundary_Conditions();

    example.Initialize_Advection();

    // initialize levelset
    particle_levelset_evolution->Use_Particle_Levelset(example.fluids_parameters.use_particle_levelset);
    if(example.fluids_parameters.fire) particle_levelset_evolution->Use_Particle_Levelset(false);
    if(example.fluids_parameters.water || example.fluids_parameters.fire){
        particle_levelset_evolution->particle_levelset.Set_Number_Particles_Per_Cell(example.fluids_parameters.number_particles_per_cell);
        particle_levelset_evolution->particle_levelset.levelset.Set_Levelset_Callbacks(example);
        particle_levelset_evolution->particle_levelset.levelset.Initialize_FMM_Initialization_Iterative_Solver(
            example.fluids_parameters.refine_fmm_initialization_with_iterative_solver);
        if(example.fluids_parameters.phi_boundary) particle_levelset_evolution->particle_levelset.levelset.Set_Custom_Boundary(*example.fluids_parameters.phi_boundary);
        particle_levelset_evolution->particle_levelset.Bias_Towards_Negative_Particles(example.fluids_parameters.bias_towards_negative_particles);}
    if(example.fluids_parameters.water){
        if(example.fluids_parameters.use_removed_positive_particles) particle_levelset_evolution->particle_levelset.Use_Removed_Positive_Particles();
        if(example.fluids_parameters.use_removed_negative_particles) particle_levelset_evolution->particle_levelset.Use_Removed_Negative_Particles();
        if(example.fluids_parameters.store_particle_ids) particle_levelset_evolution->particle_levelset.Store_Unique_Particle_Id();}

    // solid fluid coupling
    if(example.fluids_parameters.water || example.fluids_parameters.fire){
        particle_levelset_evolution->particle_levelset.levelset.Set_Collision_Body_List(collision_bodies_affecting_fluid);
        particle_levelset_evolution->particle_levelset.levelset.Set_Face_Velocities_Valid_Mask(&incompressible->valid_mask);}
    if(example.fluids_parameters.water){
        particle_levelset_evolution->particle_levelset.Set_Collision_Distance_Factors(example.fluids_parameters.min_collision_distance_factor,
            example.fluids_parameters.max_collision_distance_factor);}

    // incompressible flow
    incompressible->Set_Custom_Boundary(*example.fluids_parameters.fluid_boundary);
    incompressible->projection.elliptic_solver->Set_Relative_Tolerance(example.fluids_parameters.incompressible_tolerance);
    incompressible->projection.elliptic_solver->pcg.Set_Maximum_Iterations(example.fluids_parameters.incompressible_iterations);
    incompressible->projection.elliptic_solver->pcg.Show_Results();
    incompressible->projection.elliptic_solver->pcg.Use_Incomplete_Cholesky();
    if(number_of_regions==1 && example.fluids_parameters.second_order_cut_cell_method)
        incompressible->projection.laplace->Use_External_Level_Set(particle_levelset_evolution->particle_levelset.levelset);
    if(example.fluids_parameters.fire){
        incompressible->projection.elliptic_solver->Use_Internal_Level_Set();
        incompressible->projection.poisson->Set_Jump();}

    // set up Kolmogorov's spectrum
    if(example.fluids_parameters.kolmogorov) example.fluids_parameters.Initialize_Turbulence(time,example.frame_rate);

    // set up the initial state
    if(example.restart){
        example.Read_Output_Files_Fluids(current_frame);
        Initialize_Fluids_Grids();
        grid.Node_Iterator_Data();
        grid.Face_Iterator_Data();
        collision_bodies_affecting_fluid.Rasterize_Objects();
        collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.Minimum_Edge_Length(),5);
        if(example.fluids_parameters.water){
            particle_levelset_evolution->particle_levelset.random.Set_Seed(2606);
            if(!example.fluids_parameters.write_particles) particle_levelset_evolution->particle_levelset.Seed_Particles(time);
            particle_levelset_evolution->particle_levelset.Delete_Particles_Outside_Grid();
            if(example.fluids_parameters.delete_fluid_inside_objects) example.fluids_parameters.Delete_Particles_Inside_Objects(time);}
        if(example.fluids_parameters.fire) particle_levelset_evolution->particle_levelset.levelset.Compute_Normals(time);}
    else{
        grid.Initialize(grid.uniform_grid,example.fluids_parameters.maximum_tree_depth,3,true,true);
        example.Setup_Initial_Refinement();
        Write_Substep("After Initial Refinement",0);

        LEVELSET_DYADIC<T_GRID>& levelset=particle_levelset_evolution->particle_levelset.levelset;
        ARRAY<T> phi_ghost(grid.number_of_cells);
        levelset.boundary->Fill_Ghost_Cells_Cell(grid,levelset.phi,phi_ghost,time);
        levelset.phi=phi_ghost;
        Write_Substep("After Filling Ghost Cells",0);

        Initialize_Fluids_Grids();
        grid.Node_Iterator_Data();
        grid.Face_Iterator_Data();
        collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false);
        collision_bodies_affecting_fluid.Rasterize_Objects();
        collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.Minimum_Edge_Length(),5);
        if(example.fluids_parameters.water || example.fluids_parameters.fire){
            example.Initialize_Phi();
            collision_bodies_affecting_fluid.Initialize_Grids();
            particle_levelset_evolution->Make_Signed_Distance();
            Write_Substep("After Fast March",0);
            Update_Tree_Topology((T)1./example.frame_rate,time);
            Write_Substep("After Updating Tree Topology",0);}
        if(example.fluids_parameters.water || example.fluids_parameters.fire){
            example.Adjust_Phi_With_Sources(time);
            Write_Substep("After Adjusting For Phi and Extrapolate Phi Into Objects",0);
            particle_levelset_evolution->Make_Signed_Distance();}
        collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false); // needed before particle seeding
        collision_bodies_affecting_fluid.Compute_Grid_Visibility(); // compute grid visibility (for averaging face velocities to nodes below)
        if(example.fluids_parameters.water && particle_levelset_evolution->use_particle_levelset){
            particle_levelset_evolution->particle_levelset.random.Set_Seed(2606);
            particle_levelset_evolution->particle_levelset.Seed_Particles(time);
            particle_levelset_evolution->particle_levelset.Delete_Particles_Outside_Grid();
            if(example.fluids_parameters.delete_fluid_inside_objects) example.fluids_parameters.Delete_Particles_Inside_Objects(time);}
        if(example.fluids_parameters.use_density || example.fluids_parameters.use_temperature){
            example.fluids_parameters.Initialize_Density_And_Temperature(time);
            example.Adjust_Density_And_Temperature_With_Sources(time);}
        if(example.fluids_parameters.fire) Set_Ghost_Density_And_Temperature_Inside_Flame_Core();
        example.Initialize_Velocities();
        example.Update_Fluid_Parameters((T)1./example.frame_rate,time);
        // TODO: cannot do this at initial time because we do not have two states to use to compute pseudo velocity
        //example.fluids_parameters.Get_Neumann_And_Dirichlet_Boundary_Conditions((T)1./example.frame_rate,time); // fictitious dt
        if(example.fluids_parameters.fire){
            incompressible->projection.Update_Phi_And_Move_Velocity_Discontinuity(particle_levelset_evolution->particle_levelset.levelset,time,true);
            particle_levelset_evolution->particle_levelset.levelset.Compute_Normals();
            incompressible->Calculate_Flame_Speed();}
        if(example.fluids_parameters.water){
            ARRAY<T> phi_ghost(grid.number_of_cells);
            particle_levelset_evolution->particle_levelset.levelset.boundary->Fill_Ghost_Cells_Cell(grid,
                particle_levelset_evolution->phi,phi_ghost,time); // TODO: why do we fill ghost cells here when we also fill ghost cells in extrapolation???
            incompressible->Extrapolate_Velocity_Across_Interface(phi_ghost,example.fluids_parameters.enforce_divergence_free_extrapolation,3,0,TV(),
                &collision_bodies_affecting_fluid.face_neighbors_visible);}}

    if(example.fluids_parameters.water && example.fluids_parameters.monitor_mass){
        example.fluids_parameters.mass=particle_levelset_evolution->particle_levelset.levelset.Approximate_Negative_Material();
        LOG::cout<<"Initial Material = "<<example.fluids_parameters.mass<<std::endl;}

    // don't need to advance bodies here if we're skipping fluid simulation
    if(example.fluids_parameters.simulate && (example.fluids_parameters.smoke || example.fluids_parameters.water || example.fluids_parameters.fire))
        collision_bodies_affecting_fluid.Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,time);
}
//#####################################################################
// Function Initialize_Fluids_Grids
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_DYADIC<T_GRID>::
Initialize_Fluids_Grids()
{
    LOG::SCOPE scope("INITIALIZE FLUIDS GRIDS","initialize fluids grids");
    example.fluids_parameters.Initialize_Grids();
    if(example.fluids_parameters.water || example.fluids_parameters.fire) example.fluids_parameters.particle_levelset_evolution.Initialize_Domain();
    if(example.fluids_parameters.water) example.fluids_parameters.particle_levelset_evolution.particle_levelset.Set_Band_Width((T)2*example.fluids_parameters.particle_half_bandwidth);
    if(example.fluids_parameters.use_density) example.fluids_parameters.density_container.Initialize_Array();
    if(example.fluids_parameters.use_temperature) example.fluids_parameters.temperature_container.Initialize_Array();
    example.fluids_parameters.incompressible.Initialize_Grids();
    if(example.fluids_parameters.solid_affects_fluid && example.fluids_parameters.fluid_affects_solid) example.Initialize_Solid_Fluid_Coupling();
    example.fluids_parameters.grid->Tree_Topology_Changed();
    example.Topology_Changed();
}
//#####################################################################
// Function Set_Ghost_Density_And_Temperature_Inside_Flame_Core
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_DYADIC<T_GRID>::
Set_Ghost_Density_And_Temperature_Inside_Flame_Core()
{
    for(int i=1;i<=example.fluids_parameters.grid->number_of_nodes;i++) if(example.fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i)<=0){
        example.fluids_parameters.density_container.density(i)=example.fluids_parameters.density;
        example.fluids_parameters.temperature_container.temperature(i)=example.fluids_parameters.temperature_products;}
}
//#####################################################################
// Function Advance_To_Target_Time
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_DYADIC<T_GRID>::
Advance_To_Target_Time(const T target_time)
{
//    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_parameters.solids_evolution;

    // fake what we are simulating with number of regions (matches the uniform code better). 0=smoke, 1=water, 2=fire
    int number_of_regions=0;
    if(example.fluids_parameters.water) number_of_regions=1;
    else if(example.fluids_parameters.fire) number_of_regions=2;

    if(!example.fluids_parameters.simulate || !(example.fluids_parameters.smoke || example.fluids_parameters.water || example.fluids_parameters.fire)){
        T dt=target_time-time;
        PHYSBAM_FATAL_ERROR("SOLIDS_EVOLUTION::Advance_To_Target_Time no longer exists.");
        //solids_evolution.Advance_To_Target_Time(target_time);
        time=target_time;
        if(example.use_melting){
            if(number_of_regions>=1) example.Update_Melting_Substep_Parameters(dt,time);
            example.Melting_Substep(dt,time);}
        return;}

    T_GRID& grid=*example.fluids_parameters.grid;
    PARTICLE_LEVELSET_EVOLUTION_DYADIC<T_GRID>* particle_levelset_evolution=&example.fluids_parameters.particle_levelset_evolution;
    T_INCOMPRESSIBLE* incompressible=&example.fluids_parameters.incompressible;
    GRID_BASED_COLLISION_GEOMETRY_DYADIC<T_GRID>& collision_bodies_affecting_fluid=*example.fluids_parameters.collision_bodies_affecting_fluid;

    bool done=false;for(int substep=1;!done;substep++){
        LOG::SCOPE scope("SUBSTEP","substep %d",substep);

        T dt=next_dt;done=next_done;
        if(example.abort_when_dt_below && dt<example.abort_when_dt_below) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("ABORTING BECAUSE dt < %g",example.abort_when_dt_below));

        // TODO: Figure out where update tology should be
        Write_Substep("after scalar update",substep);
        Update_Tree_Topology(dt,time);
        Write_Substep("after update_tree_topology",substep);

        if(number_of_regions) particle_levelset_evolution->particle_levelset.Set_Number_Particles_Per_Cell(example.fluids_parameters.number_particles_per_cell);
        example.Update_Fluid_Parameters(dt,time);
        if(example.fluids_parameters.use_body_force) example.fluids_parameters.callbacks->Get_Body_Force(incompressible->force,dt,time);
        if(example.use_melting && number_of_regions==1) example.Update_Melting_Substep_Parameters(dt,time);
        if(number_of_regions==1) example.fluids_parameters.phi_boundary_water.Use_Extrapolation_Mode(false);

        {LOG::SCOPE scope("SCALAR UPDATE","scalar update");

        LOG::Time("rasterize objects");
        if(substep==1) collision_bodies_affecting_fluid.Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);
        collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(true,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);
        if(substep==1) collision_bodies_affecting_fluid.Rasterize_Objects(); // non-swept
        if(substep==1) collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.Minimum_Edge_Length(),5);  // static occupied blocks
        LOG::Time("initializing swept occupied blocks");
        example.Initialize_Swept_Occupied_Blocks_For_Advection(dt,time,true); // swept occupied blocks

        example.Scalar_Advection_Callback(dt,time);

        // important to compute ghost velocity values for particles near domain boundaries
        ARRAY<T> face_velocities_ghost;
        if(number_of_regions==1){
            face_velocities_ghost.Resize(incompressible->grid.number_of_faces,false);
            incompressible->boundary->Fill_Ghost_Cells_Face(grid,incompressible->projection.face_velocities,face_velocities_ghost,time);}

        LOG::Time("updating removed particle velocities");
        example.Modify_Removed_Particles_Before_Advection(dt,time);
        if(number_of_regions==1) example.fluids_parameters.phi_boundary_water.Use_Extrapolation_Mode(true);
        particle_levelset_evolution->Fill_Levelset_Ghost_Cells(time);
        if(number_of_regions){
            PARTICLE_LEVELSET_DYADIC<T_GRID>& pls=particle_levelset_evolution->particle_levelset;
            ARRAY<T_CELL*>& cell_pointer_from_index=grid.Cell_Pointer_From_Index();
            LINEAR_INTERPOLATION_DYADIC<T_GRID,TV> interpolation;
            if(pls.use_removed_positive_particles){
                ARRAY<TV> V_node_ghost(grid.number_of_nodes,false);
                LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID>::Interpolate_From_Faces_To_Nodes(grid,face_velocities_ghost,V_node_ghost);
                ARRAY<T> phi_node(grid.number_of_nodes,false);
                LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID>::Interpolate_From_Cells_To_Nodes(grid,particle_levelset_evolution->phi,phi_node);
                for(int cell_index=1;cell_index<=pls.removed_positive_particles.m;cell_index++) if(pls.removed_positive_particles(cell_index)){
                    T_CELL* cell=cell_pointer_from_index(cell_index);
                    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_positive_particles(cell_index);
                    for(int p=1;p<=particles.array_collection->Size();p++){
                        TV X=particles.X(p),V=interpolation.From_Close_Cell_Face(grid,cell,face_velocities_ghost,&V_node_ghost,X);
                        if(-pls.levelset.Phi(cell,X,phi_node)>1.5*particles.radius(p))
                            V.y+=example.fluids_parameters.removed_positive_particle_buoyancy_constant; // add buoyancy
                        particles.V(p)=V;}}}
            if(pls.use_removed_negative_particles) for(int cell_index=1;cell_index<=particle_levelset_evolution->particle_levelset.removed_negative_particles.m;cell_index++)
                if(pls.removed_negative_particles(cell_index)){
                    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_negative_particles(cell_index);
                    for(int p=1;p<=particles.array_collection->Size();p++) particles.V(p)+=dt*example.fluids_parameters.gravity*example.fluids_parameters.gravity_direction; // ballistic
                    // TODO: fix the body forces
                    /*if(example.fluids_parameters.use_body_force){
                      T_CELL* cell=cell_pointer_from_index(cell_index);
                      for(int p=1;p<=particles.array_collection->Size();p++)
                      particles.V(p)+=dt*interpolation.Interpolate_Nodes(cell,incompressible->force,particles.X(p));}*/}} // from external forces

        if(number_of_regions){
            LOG::Time("advecting levelset");
            if(example.fluids_parameters.fire) particle_levelset_evolution->particle_levelset.levelset.Compute_Normals(time);
            particle_levelset_evolution->Advance_Levelset(dt);
            Write_Substep("after levelset advection",0,1);
            particle_levelset_evolution->particle_levelset.Euler_Step_Particles(face_velocities_ghost,dt,time,true,true);
            Write_Substep("after particle advection",0,1);}
        if(example.fluids_parameters.use_density || example.fluids_parameters.use_temperature){LOG::Time("advecting density and temperature");
            example.fluids_parameters.Evolve_Density_And_Temperature(dt,time);}
        LOG::Time("updating velocity (explicit part)"); // TODO: fix vorticity confinement for thin shells
        incompressible->Advance_One_Time_Step_Explicit_Part(dt,time,example.fluids_parameters.implicit_viscosity,false,substep==1,&particle_levelset_evolution->phi);
        example.fluids_parameters.Blend_In_External_Velocity(dt,time);
        if(!number_of_regions) time+=dt;else time=particle_levelset_evolution->time;
        Write_Substep("after explicit part",substep,0);

        // revalidate scalars and velocity in body's new position
        collision_bodies_affecting_fluid.Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE);
        collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false); // NON-swept acceleration structures
        collision_bodies_affecting_fluid.Rasterize_Objects(); // non-swept
        collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.Minimum_Edge_Length(),5);  // static occupied blocks
        collision_bodies_affecting_fluid.Compute_Grid_Visibility(); // used in fast marching and extrapolation too... NOTE: this requires that objects_in_cell be current!
        example.Revalidate_Fluid_Scalars(); // uses visibility
        example.Revalidate_Fluid_Velocity(); // uses visibility
        Write_Substep("after scalar revalidation",0,1);

        if(number_of_regions){
            /* TODO: fix due to multiphase
            if(example.fluids_parameters.use_reacting_flow){
                example.Set_Ghost_Density_And_Temperature_Inside_Flame_Core();}
            */
           example.Extrapolate_Phi_Into_Objects(time);

            LOG::Time("modifying levelset");
            particle_levelset_evolution->Fill_Levelset_Ghost_Cells(time);
            particle_levelset_evolution->Modify_Levelset_And_Particles();
            Write_Substep("after modify levelset and particles",0,1);
            example.Revalidate_Phi_After_Modify_Levelset(); // second revalidation -- uses visibility too
            Write_Substep("after revalidate phi",0,1);

            LOG::Time("adding sources");
            if(example.Adjust_Phi_With_Sources(time)) particle_levelset_evolution->Make_Signed_Distance();
            particle_levelset_evolution->Fill_Levelset_Ghost_Cells(time);
            LOG::Time("getting sources");
            ARRAY<bool>* source_mask=0;example.Get_Source_Reseed_Mask(source_mask,time);
            if(source_mask){LOG::Time("reseeding sources");particle_levelset_evolution->particle_levelset.Reseed_Particles(time,source_mask);delete source_mask;}
            Write_Substep("after adding sources",0,1);

            LOG::Time("deleting particles"); // needs to be after adding sources, since it does a reseed
            particle_levelset_evolution->particle_levelset.Delete_Particles_Outside_Grid();
            if(example.fluids_parameters.delete_fluid_inside_objects) example.fluids_parameters.Delete_Particles_Inside_Objects(time);
            LOG::Time("deleting particles in local maxima");
            /* TODO: Not implemented in dyadic:
            if(number_of_regions==1) particle_levelset_evolution->particle_levelset.Delete_Particles_In_Local_Maximum_Phi_Cells(1);
            else for(int i=1;i<=number_of_regions;i++) if(example.fluids_parameters.dirichlet_regions(i))
                particle_levelset_evolution->Particle_Levelset(i).Delete_Particles_In_Local_Maximum_Phi_Cells(-1);
            */
            LOG::Time("deleting particles far from interface");
            for(int i=1;i<=number_of_regions;i++) particle_levelset_evolution->particle_levelset.Delete_Particles_Far_From_Interface(); // uses visibility

            LOG::Time("re-incorporating removed particles");
            example.Modify_Removed_Particles_Before_Reincorporation(dt,time);
            if(number_of_regions==1){
                particle_levelset_evolution->particle_levelset.Identify_And_Remove_Escaped_Particles(face_velocities_ghost,1.5,true);
                if(particle_levelset_evolution->particle_levelset.use_removed_positive_particles || particle_levelset_evolution->particle_levelset.use_removed_negative_particles)
                    particle_levelset_evolution->particle_levelset.Reincorporate_Removed_Particles(1);}
            example.Modify_Removed_Particles_After_Reincorporation(dt,time);}

        // update ghost phi values
        if(number_of_regions>=1) particle_levelset_evolution->Fill_Levelset_Ghost_Cells(time);

        LOG::cout<<"substep = "<<substep<<", dt = "<<dt<<std::endl;
        Write_Substep("after scalar update",substep,0);}

        // finish fluid update
        time-=dt; // rewinding to time n

        /* TODO: Not implemented in dyadic
        if(example.fluids_parameters.use_non_zero_divergence){LOG::Time("getting Divergence");
            example.fluids_parameters.callbacks->Get_Divergence(example.fluids_parameters.incompressible->projection.divergence,dt,time);}
        */

        Compute_Next_Dt(time+dt,done,next_dt,next_done);

        // update bodies (using next_dt)
        if(example.use_melting){if(number_of_regions){example.Melting_Substep(dt,time);example.Modify_Fluid_For_Melting(dt,time);}else example.Melting_Substep(dt,time);}

        collision_bodies_affecting_fluid.Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,time+dt);
        PHYSBAM_FATAL_ERROR("SOLIDS_EVOLUTION::Advance_To_Target_Time no longer exists.");
        //solids_evolution.Advance_To_Target_Time(time+dt+next_dt); // TODO: allow example to override method of evolving bodies
        if(example.use_melting && !number_of_regions) example.Modify_Fluid_For_Melting(dt,time);
        collision_bodies_affecting_fluid.Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE,time+dt+next_dt);
        collision_bodies_affecting_fluid.Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);
        collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false); // NON-swept acceleration structures
        if(example.use_melting && number_of_regions){ // re-rasterize changed bodies
            collision_bodies_affecting_fluid.Rasterize_Objects(); // non-swept
            collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.Minimum_Edge_Length(),5);}  // static occupied blocks
        Write_Substep("after poisson solve, body update",substep);

        LOG::Time("getting Neumann and Dirichlet boundary conditions");
        example.fluids_parameters.Get_Neumann_And_Dirichlet_Boundary_Conditions(dt,time+dt);
        if(!example.fluids_parameters.fire && !example.fluids_parameters.surface_tension && !example.fluids_parameters.variable_surface_tension) 
            incompressible->projection.p*=dt; // RESCALE PRESSURE FOR A BETTER INITIAL GUESS!
        if(number_of_regions==1){
            if(example.fluids_parameters.second_order_cut_cell_method) incompressible->projection.laplace->Set_Up_Second_Order_Cut_Cell_Method();
            if(example.fluids_parameters.surface_tension || example.fluids_parameters.variable_surface_tension) PHYSBAM_NOT_IMPLEMENTED();}
        LOG::Time("solving for the pressure and viscosity");
        Write_Substep("before laplace solve",substep,0);
        incompressible->Advance_One_Time_Step_Implicit_Part(dt,time,example.fluids_parameters.implicit_viscosity);
        Write_Substep("after laplace solve",substep,0);
        if(!example.fluids_parameters.fire && !example.fluids_parameters.surface_tension && !example.fluids_parameters.variable_surface_tension)
            incompressible->projection.p*=(1/dt); // scale pressure back to get a real pressure
        incompressible->boundary->Apply_Boundary_Condition(incompressible->grid,incompressible->projection.face_velocities,time+dt);
        if(example.fluids_parameters.second_order_cut_cell_method) incompressible->projection.elliptic_solver->Set_Up_Second_Order_Cut_Cell_Method(false);
        time+=dt;

        LOG::Time("extrapolating velocity across interface");
        if(number_of_regions==1){
            incompressible->Extrapolate_Velocity_Across_Interface(particle_levelset_evolution->phi,example.fluids_parameters.enforce_divergence_free_extrapolation,3,0,TV(),
                &collision_bodies_affecting_fluid.face_neighbors_visible);
            if(example.fluids_parameters.use_strain) PHYSBAM_NOT_IMPLEMENTED();}

        if(example.fluids_parameters.move_grid){example.fluids_parameters.Move_Grid(time);Initialize_Fluids_Grids();}

        Write_Substep("end of iteration",substep,0);
    }
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_DYADIC<T_GRID>::
Postprocess_Frame(const int frame)
{
    PARTICLE_LEVELSET_EVOLUTION_DYADIC<T_GRID>& particle_levelset_evolution=example.fluids_parameters.particle_levelset_evolution;
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_parameters.solids_evolution;

    if(example.fluids_parameters.water){
        example.Postprocess_Phi(time); // actually changes phi !!!
        if((frame-example.first_frame)%example.fluids_parameters.reseeding_frame_rate==0){
            LOG::Time("Reseeding... ");
            particle_levelset_evolution.particle_levelset.Reseed_Particles(time);
            particle_levelset_evolution.particle_levelset.Delete_Particles_Outside_Grid();
            if(example.fluids_parameters.delete_fluid_inside_objects) example.fluids_parameters.Delete_Particles_Inside_Objects(time);
            LOG::Stop_Time();}}

    if(example.fluids_parameters.water && example.fluids_parameters.monitor_mass){
        T mass_new=particle_levelset_evolution.particle_levelset.levelset.Approximate_Negative_Material();
        LOG::cout<<"Material = "<<mass_new<<" - change = "<<mass_new-example.fluids_parameters.mass<<std::endl;
        example.fluids_parameters.mass=mass_new;}

    SOLIDS_FLUIDS_DRIVER<TV>::Postprocess_Frame(frame);
    if(example.solids_parameters.rigid_body_collision_parameters.rigid_collisions_print_interpenetration_statistics)
        solids_evolution.rigid_body_collisions->Print_Interpenetration_Statistics();
}
//#####################################################################
// Function Compute_Next_Dt
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_DYADIC<T_GRID>::
Compute_Next_Dt(const T next_time,const bool done_this_frame,T& next_dt,bool& next_done)
{
    PHYSBAM_FATAL_ERROR("DRIVERS NOW WORK OFF OF DT, NOT NEXT DT");
    assert(example.fluids_parameters.smoke || example.fluids_parameters.fire || example.fluids_parameters.water);
    if(example.fluids_parameters.fire) example.fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.Compute_Normals(next_time);
    next_dt=example.fluids_parameters.water?FLT_MAX:example.fluids_parameters.cfl*example.fluids_parameters.incompressible.CFL();
    if(example.fluids_parameters.water || example.fluids_parameters.fire){
        T next_dt_levelset=example.fluids_parameters.particle_levelset_evolution.CFL();
        next_dt=min(next_dt,next_dt_levelset);}
    example.Limit_Dt(next_dt,time);
    next_done=false;
    if(done_this_frame) SOLIDS_FLUIDS_EXAMPLE<TV>::Clamp_Time_Step_With_Target_Time(next_time,example.Time_At_Frame(current_frame+2),next_dt,next_done);
    else SOLIDS_FLUIDS_EXAMPLE<TV>::Clamp_Time_Step_With_Target_Time(next_time,example.Time_At_Frame(current_frame+1),next_dt,next_done);
}
//#####################################################################
// Function Find_Cells_To_Refine
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_DYADIC<T_GRID>::
Find_Ghost_Cells_To_Refine(void* data,const T_CELL* cell1,const T_CELL* cell2,int axis) // nonconst int to bypass windows compiler error
{
    ARRAY<T_CELL*>* cells_to_refine=(ARRAY<T_CELL*>*)data;
    if(cell1->Depth_Of_This_Cell()<cell2->Depth_Of_This_Cell()) cells_to_refine->Append((T_CELL*)cell1);
    else if(cell2->Depth_Of_This_Cell()<cell1->Depth_Of_This_Cell()) cells_to_refine->Append((T_CELL*)cell2);
}
//#####################################################################
// Function Coarsen_Tree
//#####################################################################
template<class T_GRID> bool SOLIDS_FLUIDS_DRIVER_DYADIC<T_GRID>::
Coarsen_Tree(T_CELL* cell,const ARRAY<typename T_CELL::REFINE_ACTION>& refine_action)
{
    if(!cell->Has_Children())return refine_action(cell->Cell())==T_CELL::COARSEN;
    bool can_coarsen=true;for(int i=0;i<T_GRID::number_of_children_per_cell;i++) if(!Coarsen_Tree(cell->Child(i),refine_action))can_coarsen=false;
    if(can_coarsen){
        PROJECTION_DYADIC<T_GRID>& projection=example.fluids_parameters.incompressible.projection;
        if(example.fluids_parameters.fire){
            ARRAY<T>& phi_for_flame=projection.elliptic_solver->levelset->phi;
            T phi=0;for(int i=0;i<T_GRID::number_of_children_per_cell;i++) phi+=phi_for_flame(cell->Child(i)->Cell());
            phi_for_flame(cell->Cell())=phi*T_GRID::one_over_number_of_children_per_cell;
            CELL_HELPER::Interpolate_Face_Values_From_Direct_Children(*cell,projection.phi_face_for_flame,FACE_LOOKUP_DYADIC<T_GRID>(projection.phi_face_for_flame));
            FACE_LOOKUP_FIRE_DYADIC<T_GRID> face_lookup_fire(projection.face_velocities,projection,&phi_for_flame,&projection.phi_face_for_flame);
            CELL_HELPER::Interpolate_Face_Values_From_Direct_Children(*cell,projection.face_velocities,face_lookup_fire);}
        else CELL_HELPER::Interpolate_Face_Values_From_Direct_Children(*cell,projection.face_velocities,FACE_LOOKUP_DYADIC<T_GRID>(projection.face_velocities));
        if(example.fluids_parameters.water){
            T phi=0;for(int i=0;i<T_GRID::number_of_children_per_cell;i++) phi+=example.fluids_parameters.particle_levelset_evolution.phi(cell->Child(i)->Cell());
            example.fluids_parameters.particle_levelset_evolution.phi(cell->Cell())=phi*T_GRID::one_over_number_of_children_per_cell;}
        cell->Delete_Children();return refine_action(cell->Cell())==T_CELL::COARSEN;}
    return false;
}
//#####################################################################
// Function Coarsen_Ghost_Cells
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_DYADIC<T_GRID>::
Coarsen_Ghost_Cells(T_CELL* cell)
{
    if(cell->Has_Children()){for(int i=0;i<T_GRID::number_of_children_per_cell;i++) Coarsen_Ghost_Cells(cell->Child(i));cell->Delete_Children();}
    if(example.fluids_parameters.water) Delete_All_Particles_In_Cell(cell);
}
//#####################################################################
// Function Delete_All_Particles_In_Cell
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_DYADIC<T_GRID>::
Delete_All_Particles_In_Cell(T_CELL* cell)
{
    T_PARTICLE_LEVELSET& particle_levelset=example.fluids_parameters.particle_levelset_evolution.particle_levelset;
    if(particle_levelset.positive_particles(cell->Cell())){
        particle_levelset.particle_pool.Free_Particle(particle_levelset.positive_particles(cell->Cell()));
        particle_levelset.positive_particles(cell->Cell())=0;}
    if(particle_levelset.negative_particles(cell->Cell())){
        particle_levelset.particle_pool.Free_Particle(particle_levelset.negative_particles(cell->Cell()));
        particle_levelset.negative_particles(cell->Cell())=0;}
    if(particle_levelset.use_removed_positive_particles)
        if(particle_levelset.removed_positive_particles(cell->Cell())){
            delete particle_levelset.removed_positive_particles(cell->Cell());
            particle_levelset.removed_positive_particles(cell->Cell())=0;}
    if(particle_levelset.use_removed_negative_particles)
        if(particle_levelset.removed_negative_particles(cell->Cell())){
            delete particle_levelset.removed_negative_particles(cell->Cell());
            particle_levelset.removed_negative_particles(cell->Cell())=0;}
}
//#####################################################################
// Function Interpolate_To_Direct_Children
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_DYADIC<T_GRID>::
Interpolate_To_Direct_Children(T_CELL* cell)
{
    T_PARTICLE_LEVELSET& particle_levelset=example.fluids_parameters.particle_levelset_evolution.particle_levelset;
    T_INCOMPRESSIBLE& incompressible=example.fluids_parameters.incompressible;
    if(example.fluids_parameters.use_density) cell->Interpolate_Node_Values_To_Direct_Children(example.fluids_parameters.density_container.density);
    if(example.fluids_parameters.use_temperature) cell->Interpolate_Node_Values_To_Direct_Children(example.fluids_parameters.temperature_container.temperature);
    cell->Interpolate_Face_Values_To_Direct_Children(incompressible.projection.face_velocities);
    if(example.fluids_parameters.water || example.fluids_parameters.fire){
        for(int i=0;i<T_GRID::number_of_children_per_cell;i++)
            particle_levelset.levelset.phi(cell->Child(i)->Cell())=particle_levelset.levelset.phi(cell->Cell());
        if(example.fluids_parameters.fire){
            T_LEVELSET& levelset_fire=*incompressible.projection.elliptic_solver->levelset;
            cell->Interpolate_Face_Values_To_Direct_Children(incompressible.projection.phi_face_for_flame);
            for(int i=0;i<T_GRID::number_of_children_per_cell;i++)
                levelset_fire.phi(cell->Child(i)->Cell())=levelset_fire.phi(cell->Cell());}}
}
//#####################################################################
// Function Update_Tree_Topology
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_DYADIC<T_GRID>::
Update_Tree_Topology(const T dt,const T time)
{
    FLUIDS_PARAMETERS_DYADIC<T_GRID>& fluids_parameters=example.fluids_parameters;
    T_GRID& grid=*fluids_parameters.grid;
    T_INCOMPRESSIBLE& incompressible=fluids_parameters.incompressible;
    T_PARTICLE_LEVELSET& particle_levelset=fluids_parameters.particle_levelset_evolution.particle_levelset;
    PROJECTION_DYADIC<T_GRID>& projection=incompressible.projection;
    LOG::SCOPE scope("UPDATE TREE TOPOLOGY","update tree topology");
    if(fluids_parameters.move_grid){fluids_parameters.Move_Grid(time);Initialize_Fluids_Grids();}

    int old_number_of_cells=grid.number_of_cells,old_number_of_nodes=grid.number_of_nodes,old_number_of_faces=grid.number_of_faces;

    // coarsen the tree
    LOG::Time("coarsening tree");
    ARRAY<typename T_CELL::REFINE_ACTION> refine_action(grid.number_of_cells);example.Specify_Refinement(refine_action,0,dt,time);
    for(UNIFORM_CELL_ITERATOR iterator(grid.uniform_grid,0);iterator.Valid();iterator.Next()) if(grid.cells(iterator.Cell_Index())->Has_Children())
        Coarsen_Tree(grid.cells(iterator.Cell_Index()),refine_action);
    // delete all the children in the innermost band of ghost cells (they will be regenerated below)
    for(UNIFORM_CELL_ITERATOR iterator(grid.uniform_grid,1,UNIFORM_GRID::GHOST_REGION);iterator.Valid();iterator.Next())Coarsen_Ghost_Cells(grid.cells(iterator.Cell_Index()));

    LOG::Time("compacting array indices");
    ARRAY<int> cell_mapping_array,node_mapping_array,face_mapping_array;
    grid.Compact_Array_Indices(&cell_mapping_array,&node_mapping_array,&face_mapping_array);

    LOG::Time("compacting arrays");
    {ARRAY<T> temp(max(old_number_of_cells,old_number_of_nodes,old_number_of_faces),false);
    ARRAY<T>::Compact_Array_Using_Compaction_Array(projection.face_velocities,face_mapping_array,&temp);
    ARRAY<T>::Compact_Array_Using_Compaction_Array(projection.p,cell_mapping_array,&temp);
    if(fluids_parameters.fire || fluids_parameters.water){
        ARRAY<T>::Compact_Array_Using_Compaction_Array(fluids_parameters.particle_levelset_evolution.phi,cell_mapping_array,&temp);}
    if(fluids_parameters.fire){
        ARRAY<T>::Compact_Array_Using_Compaction_Array(projection.elliptic_solver->levelset->phi,cell_mapping_array,&temp);
        ARRAY<T>::Compact_Array_Using_Compaction_Array(projection.phi_face_for_flame,face_mapping_array,&temp);}
    if(fluids_parameters.use_density) ARRAY<T>::Compact_Array_Using_Compaction_Array(fluids_parameters.density_container.density,node_mapping_array,&temp);
    if(fluids_parameters.use_temperature) ARRAY<T>::Compact_Array_Using_Compaction_Array(fluids_parameters.temperature_container.temperature,node_mapping_array,&temp);}
    if(fluids_parameters.water){
        {ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*> temp(old_number_of_cells,false);
        ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>::Compact_Array_Using_Compaction_Array(particle_levelset.positive_particles,cell_mapping_array,&temp);
        ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>::Compact_Array_Using_Compaction_Array(particle_levelset.negative_particles,cell_mapping_array,&temp);}
        {ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*> temp(old_number_of_cells,false);
        if(particle_levelset.use_removed_positive_particles)
            ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>::Compact_Array_Using_Compaction_Array(particle_levelset.removed_positive_particles,cell_mapping_array,&temp);
        if(particle_levelset.use_removed_negative_particles)
        ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>::Compact_Array_Using_Compaction_Array(particle_levelset.removed_negative_particles,cell_mapping_array,&temp);}}
    {ARRAY<typename T_CELL::REFINE_ACTION> temp(old_number_of_cells,false);
    ARRAY<typename T_CELL::REFINE_ACTION>::Compact_Array_Using_Compaction_Array(refine_action,cell_mapping_array,&temp);}
    Initialize_Fluids_Grids();

    // refine the tree
    LOG::Time("refining interior cells");
    ARRAY<T_CELL*> new_cells;new_cells.Preallocate(grid.number_of_cells/2);
    while(true){
        bool refined_any_cells=false;int number_of_cells=grid.number_of_cells;ARRAY<T_CELL*>& cell_pointer_from_index=grid.Cell_Pointer_From_Index();
        for(int i=1;i<=number_of_cells;i++) if(refine_action(i)==T_CELL::REFINE){
            T_CELL* cell=cell_pointer_from_index(i);
            if(grid.Domain().Lazy_Inside(cell->Center()) && cell->Depth_Of_This_Cell()<grid.maximum_depth && !cell->Has_Children()){
                refined_any_cells=true;cell->Create_Children(grid.number_of_cells,&new_cells,grid.number_of_nodes,0,grid.number_of_faces,0,&grid);}
            else refine_action(i)=T_CELL::NONE;} // so that we know which cells to interpolate from
        if(!refined_any_cells) break;
        Initialize_Fluids_Grids();grid.Cell_Pointer_From_Index(); // refreshes the cell_pointer_from_index_array
        for(int i=1;i<=number_of_cells;i++) if(refine_action(i)==T_CELL::REFINE)Interpolate_To_Direct_Children(cell_pointer_from_index(i));
        refine_action.Resize(0);refine_action.Resize(grid.number_of_cells);example.Specify_Refinement(refine_action,&new_cells,dt,time);new_cells.Remove_All();}

    LOG::Time("refining ghost cells");
    while(true){
        ARRAY<T_CELL*> cells_to_refine;
        MAP_MESH::Map_Left_Boundary_Faces(grid.uniform_grid,grid.cells,&cells_to_refine,Find_Ghost_Cells_To_Refine);
        MAP_MESH::Map_Right_Boundary_Faces(grid.uniform_grid,grid.cells,&cells_to_refine,Find_Ghost_Cells_To_Refine);
        for(int i=2;i<=T_GRID::number_of_faces_per_cell;i++)
            MAP_MESH::Map_Domain_Side_Faces(grid.uniform_grid,grid.cells,i,grid.number_of_ghost_cells,&cells_to_refine,Find_Ghost_Cells_To_Refine);
        if(cells_to_refine.m==0) break;
        for(int i=1;i<=cells_to_refine.m;i++) if(!cells_to_refine(i)->Has_Children()){
            cells_to_refine(i)->Create_Children(grid.number_of_cells,0,grid.number_of_nodes,0,grid.number_of_faces,0,&grid);}}

    LOG::cout<<"Cells: "<<old_number_of_cells<<"->"<<grid.number_of_cells<<", Nodes: "<<old_number_of_nodes<<"->"<<grid.number_of_nodes<<", Faces: "<<old_number_of_faces<<"->"<<grid.number_of_faces<<std::endl;

    Initialize_Fluids_Grids();

    LOG::Time("updating iterator data");
    grid.Node_Iterator_Data();grid.Face_Iterator_Data();grid.Fully_Refined_Block();
    example.fluids_parameters.collision_bodies_affecting_fluid->Initialize_Grids();

    // delete the particles that are no longer in valid blocks
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){
        if(!grid.fully_refined_block(iterator.Cell_Index())){
            Delete_All_Particles_In_Cell(iterator.Cell_Pointer());}}
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_DYADIC<T_GRID>::
Write_Output_Files(const int frame)
{
    LOG::Time("writing output files");
    FILE_UTILITIES::Create_Directory(example.output_directory);
    FILE_UTILITIES::Create_Directory(example.output_directory+STRING_UTILITIES::string_sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
    Write_First_Frame(frame);
    if(example.fluids_parameters.water) example.fluids_parameters.phi_boundary_water.Use_Extrapolation_Mode(false);
    example.Write_Output_Files(frame);
    if(example.fluids_parameters.water) example.fluids_parameters.phi_boundary_water.Use_Extrapolation_Mode(true);

    if(example.fluids_parameters.water){
        T_PARTICLE_LEVELSET& pls=example.fluids_parameters.particle_levelset_evolution.particle_levelset;
        int number_of_positive_particles=0;for(int i=1;i<=pls.positive_particles.m;i++) if(pls.positive_particles(i))
            number_of_positive_particles+=pls.positive_particles(i)->array_collection->Size();
        int number_of_negative_particles=0;for(int i=1;i<=pls.negative_particles.m;i++) if(pls.negative_particles(i))
            number_of_negative_particles+=pls.negative_particles(i)->array_collection->Size();
        LOG::cout<<number_of_positive_particles<<" positive and "<<number_of_negative_particles<<" negative particles "<<std::endl;

        int number_of_removed_positive_particles=0,number_of_removed_negative_particles=0;
        if(pls.use_removed_positive_particles)
            for(int i=1;i<=pls.removed_positive_particles.m;i++) if(pls.removed_positive_particles(i))
            number_of_removed_positive_particles+=pls.removed_positive_particles(i)->array_collection->Size();
        if(pls.use_removed_negative_particles)
            for(int i=1;i<=pls.removed_negative_particles.m;i++) if(pls.removed_negative_particles(i))
            number_of_removed_negative_particles+=pls.removed_negative_particles(i)->array_collection->Size();
            LOG::cout<<number_of_removed_positive_particles<<" positive and "<<number_of_removed_negative_particles<<" negative removed particles "<<std::endl;}
    Write_Time(frame);
    Write_Last_Frame(frame);
    LOG::Stop_Time();
}
//#####################################################################
template class SOLIDS_FLUIDS_DRIVER_DYADIC<QUADTREE_GRID<float> >;
template class SOLIDS_FLUIDS_DRIVER_DYADIC<OCTREE_GRID<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SOLIDS_FLUIDS_DRIVER_DYADIC<QUADTREE_GRID<double> >;
template class SOLIDS_FLUIDS_DRIVER_DYADIC<OCTREE_GRID<double> >;
#endif
#endif

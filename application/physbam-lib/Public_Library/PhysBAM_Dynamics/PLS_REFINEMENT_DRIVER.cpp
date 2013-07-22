//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Tools/Vectors/VECTOR_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Dynamics/Boundaries/BOUNDARY_PHI_WATER.h>
#include <PhysBAM_Dynamics/PLS_REFINEMENT_DRIVER.h>
#include <PhysBAM_Dynamics/PLS_REFINEMENT_EXAMPLE.h>
using namespace PhysBAM;
namespace{
template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
{
    ((PLS_REFINEMENT_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
}
};
//#####################################################################
// Initialize
//#####################################################################
template<class TV> PLS_REFINEMENT_DRIVER<TV>::
PLS_REFINEMENT_DRIVER(PLS_REFINEMENT_EXAMPLE<TV>& example)
    :example(example),kinematic_evolution(example.rigid_geometry_collection,true)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> PLS_REFINEMENT_DRIVER<TV>::
~PLS_REFINEMENT_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void PLS_REFINEMENT_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void PLS_REFINEMENT_DRIVER<TV>::
Initialize()
{
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);

    // setup time
    if(example.restart) current_frame=example.restart;else current_frame=example.first_frame;
    output_number=current_frame;
    time=example.Time_At_Frame(current_frame);
 
    // initialize collision objects
    example.Initialize_Bodies();
    example.collision_bodies_affecting_fluid.Add_Bodies(example.rigid_geometry_collection);
    kinematic_evolution.Get_Current_Kinematic_Keyframes(1/example.frame_rate,time);
    kinematic_evolution.Set_External_Positions(example.rigid_geometry_collection.particles.X,example.rigid_geometry_collection.particles.rotation,time);
    kinematic_evolution.Set_External_Velocities(example.rigid_geometry_collection.particles.V,example.rigid_geometry_collection.particles.angular_velocity,time,time);

    example.phi_boundary_water.Set_Velocity_Pointer(example.fine_face_velocities);

    {
        example.particle_levelset_evolution.Initialize_Domain(example.fine_mac_grid);
        example.Set_Band_Width(6);
        example.incompressible.Initialize_Grids(example.fine_mac_grid);
        example.projection.Initialize_Grid(example.coarse_mac_grid);
        example.collision_bodies_affecting_fluid.Initialize_Grids();
    }
    example.coarse_face_velocities.Resize(example.coarse_mac_grid);
    example.fine_face_velocities.Resize(example.fine_mac_grid);
    example.coarse_phi.Resize(example.coarse_mac_grid.Domain_Indices(1));

    example.particle_levelset_evolution.Set_Time(time);
    example.particle_levelset_evolution.Set_CFL_Number(example.cfl);

    if(example.coarse_mpi_grid) example.coarse_mpi_grid->Initialize(example.domain_boundary);
    if(example.coarse_mpi_grid) example.coarse_mpi_grid->Initialize(example.non_mpi_boundary);
    example.incompressible.mpi_grid=example.fine_mpi_grid;
    example.projection.elliptic_solver->mpi_grid=example.coarse_mpi_grid;
    example.particle_levelset_evolution.Particle_Levelset(1).mpi_grid=example.fine_mpi_grid;
    if(example.fine_mpi_grid){
        example.boundary=new BOUNDARY_MPI<GRID<TV> >(example.fine_mpi_grid,example.boundary_scalar);
        example.boundary_coarse=new BOUNDARY_MPI<GRID<TV> >(example.coarse_mpi_grid,example.boundary_scalar);
        example.phi_boundary=new BOUNDARY_MPI<GRID<TV> >(example.fine_mpi_grid,example.phi_boundary_water);
        example.particle_levelset_evolution.Particle_Levelset(1).last_unique_particle_id=example.fine_mpi_grid->rank*30000000;}
    else{
        example.boundary=&example.boundary_scalar;
        example.phi_boundary=&example.phi_boundary_water;}
    example.Initialize_MPI();

    VECTOR<VECTOR<bool,2>,TV::dimension> domain_open_boundaries=VECTOR_UTILITIES::Complement(example.domain_boundary);
    example.phi_boundary->Set_Constant_Extrapolation(domain_open_boundaries);
    example.boundary->Set_Constant_Extrapolation(domain_open_boundaries);
    if(example.use_collidable_advection){
        example.particle_levelset_evolution.Levelset_Advection(1).Use_Semi_Lagrangian_Collidable_Advection(example.collision_bodies_affecting_fluid,(T)1e-5,example.incompressible.valid_mask);
        example.particle_levelset_evolution.Levelset_Advection(1).Set_Custom_Advection(*new T_ADVECTION_SEMI_LAGRANGIAN_SCALAR);}
    else example.particle_levelset_evolution.Levelset_Advection(1).Set_Custom_Advection(example.advection_scalar);
    example.incompressible.Set_Custom_Advection(example.advection_scalar);

    {
        example.particle_levelset_evolution.Initialize_Domain(example.fine_mac_grid);
        example.Set_Band_Width(6);
        example.incompressible.Initialize_Grids(example.fine_mac_grid);
        example.projection.Initialize_Grid(example.coarse_mac_grid);
        example.collision_bodies_affecting_fluid.Initialize_Grids();
    }

    example.particle_levelset_evolution.Set_Number_Particles_Per_Cell(16);
    example.particle_levelset_evolution.Set_Levelset_Callbacks(example);
    example.particle_levelset_evolution.Initialize_FMM_Initialization_Iterative_Solver(true);

    example.particle_levelset_evolution.particle_levelset.levelset.Set_Custom_Boundary(*example.phi_boundary);
    example.particle_levelset_evolution.Bias_Towards_Negative_Particles(false);
    example.particle_levelset_evolution.Particle_Levelset(1).Use_Removed_Positive_Particles();
    example.particle_levelset_evolution.Particle_Levelset(1).Use_Removed_Negative_Particles();
    example.particle_levelset_evolution.Particle_Levelset(1).Store_Unique_Particle_Id();
    example.particle_levelset_evolution.Use_Particle_Levelset(true);
    example.particle_levelset_evolution.particle_levelset.levelset.Set_Collision_Body_List(example.collision_bodies_affecting_fluid);
    example.particle_levelset_evolution.particle_levelset.levelset.Set_Face_Velocities_Valid_Mask(&example.incompressible.valid_mask);
    example.particle_levelset_evolution.particle_levelset.Set_Collision_Distance_Factors((T).1,1);

    example.incompressible.Set_Custom_Boundary(*example.boundary);
    example.incompressible.projection.elliptic_solver->Set_Relative_Tolerance((T)1e-8);
    example.incompressible.projection.elliptic_solver->pcg.Set_Maximum_Iterations(1000);
    example.incompressible.projection.elliptic_solver->pcg.evolution_solver_type=krylov_solver_cg;
    example.incompressible.projection.elliptic_solver->pcg.cg_restart_iterations=40;
    example.incompressible.projection.elliptic_solver->pcg.Show_Results();
    example.incompressible.projection.collidable_solver->Use_External_Level_Set(*new T_LEVELSET(example.coarse_mac_grid,example.coarse_phi));

    {
        example.particle_levelset_evolution.Initialize_Domain(example.fine_mac_grid);
        example.Set_Band_Width(6);
        example.incompressible.Initialize_Grids(example.fine_mac_grid);
        example.projection.Initialize_Grid(example.coarse_mac_grid);
        example.collision_bodies_affecting_fluid.Initialize_Grids();
    }

    if(example.restart){
        example.Read_Output_Files(example.restart);
        example.collision_bodies_affecting_fluid.Rasterize_Objects();
        example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*example.fine_mac_grid.Minimum_Edge_Length(),5);} // compute grid visibility (for advection later)
    else{
        example.collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false);
        example.collision_bodies_affecting_fluid.Rasterize_Objects();
        example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*example.fine_mac_grid.Minimum_Edge_Length(),5);
        example.Initialize_Phi();
        example.Adjust_Phi_With_Sources(time);
        example.particle_levelset_evolution.Make_Signed_Distance();
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time);}

    example.collision_bodies_affecting_fluid.Compute_Grid_Visibility();
    example.particle_levelset_evolution.Set_Seed(2606);
    if((!example.write_debug_data && example.restart) || !example.restart) example.particle_levelset_evolution.Reseed_Particles(time);
    example.particle_levelset_evolution.Delete_Particles_Outside_Grid();
    
    //add forces
    example.incompressible.Set_Gravity(example.gravity);
    example.incompressible.Set_Body_Force(true);
    example.incompressible.projection.Use_Non_Zero_Divergence(false);
    example.incompressible.projection.elliptic_solver->Solve_Neumann_Regions(true);
    example.incompressible.projection.elliptic_solver->solve_single_cell_neumann_regions=false;
    example.incompressible.Use_Explicit_Part_Of_Implicit_Viscosity(false);
    example.incompressible.Set_Maximum_Implicit_Viscosity_Iterations(40);
    example.incompressible.Use_Variable_Vorticity_Confinement(false);
    example.incompressible.Set_Surface_Tension(0);
    example.incompressible.Set_Variable_Surface_Tension(false);
    example.incompressible.Set_Viscosity(0);
    example.incompressible.Set_Variable_Viscosity(false);
    example.incompressible.projection.Set_Density(1e3);

    example.Extrapolate_Velocity_Across_Interface(3,time);

    if(example.restart) example.Set_Coarse_Phi_From_Fine_Phi(example.coarse_phi,example.particle_levelset_evolution.phi);
    example.Set_Boundary_Conditions(time); // get so CFL is correct
    if(!example.restart) Write_Output_Files(example.first_frame);
    if(example.split_dir!=""){Write_Output_Files(current_frame);exit(0);}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after init",0,1);
}
//#####################################################################
// Advance_To_Target_Time
//#####################################################################
template<class TV> void PLS_REFINEMENT_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;for(int substep=1;!done;substep++){
        LOG::SCOPE scope("SUBSTEP","substep %d",substep);
        // choose time step
        example.particle_levelset_evolution.Set_Number_Particles_Per_Cell(16);
        T dt=example.cfl*example.incompressible.CFL(example.fine_face_velocities);dt=min(dt,example.particle_levelset_evolution.CFL(true,false));
        if(example.fine_mpi_grid) example.fine_mpi_grid->Synchronize_Dt(dt);
        if(time+dt>=target_time){dt=target_time-time;done=true;}
        else if(time+2*dt>=target_time){dt=(T).5*(target_time-time);}
        std::stringstream ss;ss<<"dt is "<<dt<<std::endl;LOG::filecout(ss.str());

        if(example.use_collidable_advection) example.collision_bodies_affecting_fluid.Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,time);

        // kinematic_update
        kinematic_evolution.Set_External_Positions(example.rigid_geometry_collection.particles.X,example.rigid_geometry_collection.particles.rotation,time);
        kinematic_evolution.Set_External_Velocities(example.rigid_geometry_collection.particles.V,example.rigid_geometry_collection.particles.angular_velocity,time,time);
        for(int i=1;i<=example.rigid_geometry_collection.kinematic_rigid_geometry.m;i++){
            RIGID_GEOMETRY<TV>& rigid_geometry=example.rigid_geometry_collection.Rigid_Geometry(i);            
            rigid_geometry.X()+=dt*rigid_geometry.V();
            rigid_geometry.Rotation()=ROTATION<TV>::From_Rotation_Vector(dt*rigid_geometry.Angular_Velocity())*rigid_geometry.Rotation();rigid_geometry.Rotation().Normalize();}
        kinematic_evolution.Set_External_Positions(example.rigid_geometry_collection.particles.X,example.rigid_geometry_collection.particles.rotation,time+dt);
 
        LOG::Time("Collision setup 1");
        PARTICLE_LEVELSET_UNIFORM<GRID<TV> >& pls=example.particle_levelset_evolution.Particle_Levelset(1);
        if(example.use_collidable_advection){
            example.collision_bodies_affecting_fluid.Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE,time+dt);
            example.collision_bodies_affecting_fluid.Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);
            example.collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(true,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE,
                COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);
            example.collision_bodies_affecting_fluid.Rasterize_Objects(); // non-swept
            example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*example.fine_mac_grid.Minimum_Edge_Length(),5);  // static occupied blocks
            T maximum_particle_speed=0,max_particle_collision_distance=pls.max_collision_distance_factor*example.fine_mac_grid.dX.Max();
            if(pls.use_removed_negative_particles) for(typename GRID<TV>::CELL_ITERATOR iterator(pls.levelset.grid);iterator.Valid();iterator.Next()){
                PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles=pls.removed_negative_particles(iterator.Cell_Index());
                if(particles) maximum_particle_speed=max(maximum_particle_speed,ARRAYS_COMPUTATIONS::Maximum_Magnitude(particles->V));}
            if(pls.use_removed_positive_particles) for(typename GRID<TV>::CELL_ITERATOR iterator(pls.levelset.grid);iterator.Valid();iterator.Next()){
                PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles=pls.removed_positive_particles(iterator.Cell_Index());
                if(particles) maximum_particle_speed=max(maximum_particle_speed,ARRAYS_COMPUTATIONS::Maximum_Magnitude(particles->V));}
            example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(true,dt*maximum_particle_speed+2*max_particle_collision_distance+(T).5*example.fine_mac_grid.dX.Max(),10);}
        else{
            T maximum_fluid_speed=example.fine_face_velocities.Maxabs().Max();
            T max_particle_collision_distance=example.particle_levelset_evolution.Particle_Levelset(1).max_collision_distance_factor*example.fine_mac_grid.dX.Max();
            example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(true,dt*maximum_fluid_speed+2*max_particle_collision_distance+(T).5*example.fine_mac_grid.dX.Max(),10);}

        LOG::Time("Fill ghost");
        PHYSBAM_DEBUG_WRITE_SUBSTEP("before advection",0,1);
        ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities_ghost;face_velocities_ghost.Resize(example.incompressible.grid,example.number_of_ghost_cells,false);
        example.boundary->Fill_Ghost_Cells_Face(example.fine_mac_grid,example.fine_face_velocities,face_velocities_ghost,time+dt,example.number_of_ghost_cells);

        LOG::Time("Advect Levelset");
        example.phi_boundary_water.Use_Extrapolation_Mode(false);
        example.Adjust_Phi_With_Objects(time);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("before phi",0,1);
        example.particle_levelset_evolution.Advance_Levelset(dt);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after phi",0,1);
        example.phi_boundary_water.Use_Extrapolation_Mode(true);
        example.Extrapolate_Phi_Into_Objects(time+dt);
        LOG::Time("Step Particles");
        example.particle_levelset_evolution.Particle_Levelset(1).Euler_Step_Particles(face_velocities_ghost,dt,time,true,true,false,false);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after advection",0,1);
        
        LOG::Time("Fill ghost");
        example.phi_boundary_water.Use_Extrapolation_Mode(true);
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time);
        LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,TV> interpolation;
        TV gravity_direction=-TV::Axis_Vector(2);
        LOG::Time("Update Particles");
        if(pls.use_removed_positive_particles) for(typename GRID<TV>::NODE_ITERATOR iterator(example.fine_mac_grid);iterator.Valid();iterator.Next()) if(pls.removed_positive_particles(iterator.Node_Index())){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_positive_particles(iterator.Node_Index());
            for(int p=1;p<=particles.array_collection->Size();p++){
                TV X=particles.X(p),V=interpolation.Clamped_To_Array_Face(example.fine_mac_grid,face_velocities_ghost,X);
                if(-pls.levelset.Phi(X)>1.5*particles.radius(p)) V-=gravity_direction*(T).3; // buoyancy
                particles.V(p)=V;}}
        if(pls.use_removed_negative_particles) for(typename GRID<TV>::NODE_ITERATOR iterator(example.fine_mac_grid);iterator.Valid();iterator.Next()) if(pls.removed_negative_particles(iterator.Node_Index())){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_negative_particles(iterator.Node_Index());
            for(int p=1;p<=particles.array_collection->Size();p++) particles.V(p)+=dt*example.gravity*gravity_direction;
            for(int p=1;p<=particles.array_collection->Size();p++) particles.V(p)+=dt*interpolation.Clamped_To_Array_Face(example.fine_mac_grid,example.incompressible.force,particles.X(p));} // ballistic
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after particles",0,1);

        LOG::Time("Convection");
        example.incompressible.Advance_One_Time_Step_Convection(dt,time,example.fine_face_velocities,example.fine_face_velocities,example.number_of_ghost_cells);
        example.boundary->Apply_Boundary_Condition_Face(example.fine_mac_grid,example.fine_face_velocities,time+dt);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after convection",0,1);
        LOG::Time("Forces");
        example.incompressible.Advance_One_Time_Step_Forces(example.fine_face_velocities,dt,time,false,0,example.number_of_ghost_cells);
        example.boundary->Fill_Ghost_Cells_Face(example.fine_mac_grid,example.fine_face_velocities,face_velocities_ghost,time+dt,example.number_of_ghost_cells);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after forces",0,1);

        LOG::Time("Update Collidable 2");
        if(example.use_collidable_advection){
            example.collision_bodies_affecting_fluid.Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE);
            example.collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false); // NON-swept acceleration structures
            example.collision_bodies_affecting_fluid.Rasterize_Objects(); // non-swept
            example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*example.fine_mac_grid.Minimum_Edge_Length(),5);  // static occupied blocks
            example.collision_bodies_affecting_fluid.Compute_Grid_Visibility(); // used in fast marching and extrapolation too... NOTE: this requires that objects_in_cell be current!
            if(example.particle_levelset_evolution.Levelset_Advection(1).nested_semi_lagrangian_collidable)  
                example.particle_levelset_evolution.Levelset_Advection(1).nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(example.fine_mac_grid,(T)1e-5,example.particle_levelset_evolution.Levelset(1).phi);
            if(example.incompressible.nested_semi_lagrangian_collidable) 
                example.incompressible.nested_semi_lagrangian_collidable->Average_To_Invalidated_Face(example.fine_mac_grid,example.fine_face_velocities);
        }

        LOG::Time("Fill ghost");
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time+dt);
        LOG::Time("Exchange Particles");
        pls.Exchange_Overlap_Particles();
        LOG::Time("Modifying Levelset");
        example.particle_levelset_evolution.Modify_Levelset_And_Particles(&face_velocities_ghost);
        example.particle_levelset_evolution.Make_Signed_Distance();
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after modify particles",0,1);

        LOG::Time("adding sources");
        example.Adjust_Phi_With_Sources(time+dt);
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time+dt);
        LOG::Time("getting sources");

        LOG::Time("deleting particles"); // needs to be after adding sources, since it does a reseed
        example.particle_levelset_evolution.Delete_Particles_Outside_Grid();
        LOG::Time("deleting particles in local maxima");
        example.particle_levelset_evolution.particle_levelset.Delete_Particles_In_Local_Maximum_Phi_Cells(1);
        LOG::Time("deleting particles far from interface");
        example.particle_levelset_evolution.Particle_Levelset(1).Delete_Particles_Far_From_Interface(); // uses visibility
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete",0,1);

        LOG::Time("re-incorporating removed particles");
        example.particle_levelset_evolution.Particle_Levelset(1).Identify_And_Remove_Escaped_Particles(face_velocities_ghost,1.5,time+dt);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after remove",0,1);
        if(example.particle_levelset_evolution.Particle_Levelset(1).use_removed_positive_particles || example.particle_levelset_evolution.Particle_Levelset(1).use_removed_negative_particles)
            example.particle_levelset_evolution.Particle_Levelset(1).Reincorporate_Removed_Particles(1,1,0,true);

        LOG::Time("Fill ghost");
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time+dt);
        
        //Project
        kinematic_evolution.Set_External_Velocities(example.rigid_geometry_collection.particles.V,example.rigid_geometry_collection.particles.angular_velocity,time+dt,time+dt);
        
        LOG::Time("getting Neumann and Dirichlet boundary conditions");
        ARRAY<T,TV_INT> phi_ghost(example.fine_mac_grid.Domain_Indices(example.number_of_ghost_cells));
        example.phi_boundary->Fill_Ghost_Cells(example.fine_mac_grid,example.particle_levelset_evolution.phi,phi_ghost,dt,time,example.number_of_ghost_cells);
        example.Set_Coarse_Phi_From_Fine_Phi(example.coarse_phi,phi_ghost);
        example.Set_Boundary_Conditions(time+dt);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after boundary",0,1);
        
        LOG::Time("Coarse Projection");
        example.Preprocess_Projection(dt,time);
        //example.incompressible.Set_Dirichlet_Boundary_Conditions(&example.coarse_phi,0);
        example.projection.p*=dt; // rescale pressure for guess
        example.incompressible.projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method();
        PHYSBAM_DEBUG_WRITE_SUBSTEP("before projection",0,1);
        example.incompressible.Advance_One_Time_Step_Implicit_Part(example.coarse_face_velocities,dt,time,false,example.boundary_coarse);
        example.projection.p*=(1/dt); // unscale pressure
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after projection",0,1);
        if(example.boundary_coarse) example.boundary_coarse->Apply_Boundary_Condition_Face(example.coarse_mac_grid,example.coarse_face_velocities,time+dt);
        LOG::Time("Fine Projection");
        example.incompressible.projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method(false);
        example.Postprocess_Projection(dt,time);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after postprocess",0,1);

        LOG::Time("Boundary");
        example.boundary->Apply_Boundary_Condition_Face(example.fine_mac_grid,example.fine_face_velocities,time+dt);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after boundary",0,1);

        LOG::Time("Extrapolation");
        example.Extrapolate_Velocity_Across_Interface(3,time);

        time+=dt;}
}
//#####################################################################
// Simulate_To_Frame
//#####################################################################
template<class TV> void PLS_REFINEMENT_DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    while(current_frame<frame){
        LOG::SCOPE scope("FRAME","Frame %d",current_frame+1);
        kinematic_evolution.Get_Current_Kinematic_Keyframes(example.Time_At_Frame(current_frame+1)-time,time);
        Advance_To_Target_Time(example.Time_At_Frame(current_frame+1));
        if((current_frame-example.first_frame)%20==0){
            example.particle_levelset_evolution.Reseed_Particles(time);
            example.particle_levelset_evolution.Delete_Particles_Outside_Grid();}
        Write_Output_Files(++output_number);
        current_frame++;}
} 
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void PLS_REFINEMENT_DRIVER<TV>::
Write_Substep(const std::string& title,const int substep,const int level)
{
    if(level<=example.write_substeps_level){
        example.frame_title=title;
        std::stringstream ss;ss<<"Writing substep ["<<example.frame_title<<"]: output_number="<<output_number+1<<", time="<<time<<", frame="<<current_frame<<", substep="<<substep<<std::endl;LOG::filecout(ss.str());
        Write_Output_Files(++output_number);example.frame_title="";}
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void PLS_REFINEMENT_DRIVER<TV>::
Write_Output_Files(const int frame)
{
    FILE_UTILITIES::Create_Directory(example.output_directory);
    FILE_UTILITIES::Create_Directory(example.output_directory+STRING_UTILITIES::string_sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+STRING_UTILITIES::string_sprintf("/%d/frame_title",frame),example.frame_title);
    if(frame==example.first_frame) 
        FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/first_frame",frame,"\n");
    example.Write_Output_Files(frame);
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/last_frame",frame,"\n");
}
//#####################################################################
template class PLS_REFINEMENT_DRIVER<VECTOR<float,2> >;
template class PLS_REFINEMENT_DRIVER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PLS_REFINEMENT_DRIVER<VECTOR<double,2> >;
template class PLS_REFINEMENT_DRIVER<VECTOR<double,3> >;
#endif

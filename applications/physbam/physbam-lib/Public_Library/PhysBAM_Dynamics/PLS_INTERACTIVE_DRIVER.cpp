//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Tools/Vectors/VECTOR_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM_FORWARD.h>
#include <PhysBAM_Dynamics/Boundaries/BOUNDARY_PHI_WATER.h>
#include <PhysBAM_Dynamics/PLS_INTERACTIVE_DRIVER.h>
#include <PhysBAM_Dynamics/PLS_EXAMPLE.h>
#ifdef USE_PTHREADS
    #include <pthread.h>
#endif

using namespace PhysBAM;
namespace{
template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
{
    ((PLS_INTERACTIVE_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
}
};
//#####################################################################
// Initialize
//#####################################################################
template<class TV> PLS_INTERACTIVE_DRIVER<TV>::
PLS_INTERACTIVE_DRIVER(PLS_EXAMPLE<TV>& example)
    :example(example)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> PLS_INTERACTIVE_DRIVER<TV>::
~PLS_INTERACTIVE_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void PLS_INTERACTIVE_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void PLS_INTERACTIVE_DRIVER<TV>::
Initialize()
{
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);

    // setup time
    if(example.restart) current_frame=example.restart;else current_frame=example.first_frame;
    output_number=current_frame;
    time=example.Time_At_Frame(current_frame);

    example.phi_boundary_water.Set_Velocity_Pointer(example.face_velocities);

    {
        example.particle_levelset_evolution.Initialize_Domain(example.mac_grid);
        example.particle_levelset_evolution.particle_levelset.Set_Band_Width(6);
        example.incompressible.Initialize_Grids(example.mac_grid);
        example.projection.Initialize_Grid(example.mac_grid);
        example.collision_bodies_affecting_fluid.Initialize_Grids();
    }
    example.face_velocities.Resize(example.mac_grid);
    if(example.incompressible.conserve_kinetic_energy) example.incompressible.Conserve_Kinetic_Energy();

    example.particle_levelset_evolution.Set_Time(time);
    example.particle_levelset_evolution.Set_CFL_Number((T).9);

    if(example.mpi_grid) example.mpi_grid->Initialize(example.domain_boundary);
    ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>* advection_conservative=dynamic_cast<ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>*>(example.incompressible.advection);
    if(advection_conservative){
        VECTOR<VECTOR<bool,2>,TV::dimension> domain_boundary_local;
        for(int i=1;i<=TV::dimension;i++){domain_boundary_local(i)(1)=true;domain_boundary_local(i)(2)=true;}        
        if(example.mpi_grid) example.mpi_grid->Initialize(domain_boundary_local);
        if(example.mpi_grid) example.mpi_grid->Initialize(advection_conservative->solid_walls);
        for(int i=1;i<=TV::dimension;i++){advection_conservative->mpi_boundary(i)(1)=!domain_boundary_local(i)(1);advection_conservative->mpi_boundary(i)(2)=!domain_boundary_local(i)(2);}}
    example.incompressible.mpi_grid=example.mpi_grid;
    example.projection.elliptic_solver->mpi_grid=example.mpi_grid;
    example.particle_levelset_evolution.Particle_Levelset(1).mpi_grid=example.mpi_grid;
    if(example.mpi_grid){
        example.boundary=new BOUNDARY_MPI<GRID<TV> >(example.mpi_grid,example.boundary_scalar);
        example.phi_boundary=new BOUNDARY_MPI<GRID<TV> >(example.mpi_grid,example.phi_boundary_water);
        example.particle_levelset_evolution.Particle_Levelset(1).last_unique_particle_id=example.mpi_grid->rank*30000000;}
    else{
        example.boundary=&example.boundary_scalar;
        example.phi_boundary=&example.phi_boundary_water;}

    if(PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<GRID<TV> > *refine=dynamic_cast<PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<GRID<TV> >*>(&example.projection)){
        refine->boundary=example.boundary;
        refine->phi_boundary=example.phi_boundary;}

    VECTOR<VECTOR<bool,2>,TV::dimension> domain_open_boundaries=VECTOR_UTILITIES::Complement(example.domain_boundary);
    example.phi_boundary->Set_Constant_Extrapolation(domain_open_boundaries);
    example.boundary->Set_Constant_Extrapolation(domain_open_boundaries);
    example.particle_levelset_evolution.Levelset_Advection(1).Set_Custom_Advection(example.advection_scalar);
    //example.incompressible.Set_Custom_Advection(example.advection_scalar);

    {
        example.particle_levelset_evolution.Initialize_Domain(example.mac_grid);
        example.particle_levelset_evolution.particle_levelset.Set_Band_Width(6);
        example.incompressible.Initialize_Grids(example.mac_grid);
        example.projection.Initialize_Grid(example.mac_grid);
        example.collision_bodies_affecting_fluid.Initialize_Grids();
    }
    
    example.particle_levelset_evolution.Particle_Levelset(1).number_of_ghost_cells=example.number_of_ghost_cells;
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
    example.incompressible.projection.elliptic_solver->pcg.Set_Maximum_Iterations(40);
    example.incompressible.projection.elliptic_solver->pcg.evolution_solver_type=krylov_solver_cg;
    example.incompressible.projection.elliptic_solver->pcg.cg_restart_iterations=0;
    example.incompressible.projection.elliptic_solver->pcg.Show_Results();
    example.incompressible.projection.collidable_solver->Use_External_Level_Set(example.particle_levelset_evolution.particle_levelset.levelset);

    {
        example.particle_levelset_evolution.Initialize_Domain(example.mac_grid);
        example.particle_levelset_evolution.particle_levelset.Set_Band_Width(6);
        example.incompressible.Initialize_Grids(example.mac_grid);
        example.projection.Initialize_Grid(example.mac_grid);
        example.collision_bodies_affecting_fluid.Initialize_Grids();
    }

    if(example.restart){
        example.Read_Output_Files(example.restart);
        example.collision_bodies_affecting_fluid.Rasterize_Objects();
        example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*example.mac_grid.Minimum_Edge_Length(),5);} // compute grid visibility (for advection later)
    else{
        example.collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false);
        example.collision_bodies_affecting_fluid.Rasterize_Objects();
        example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*example.mac_grid.Minimum_Edge_Length(),5);
        example.Initialize_Phi();
        example.Adjust_Phi_With_Sources(time);
        example.particle_levelset_evolution.Make_Signed_Distance();
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time);}

    example.collision_bodies_affecting_fluid.Compute_Grid_Visibility();
    example.particle_levelset_evolution.Set_Seed(2606);
    if(!example.restart) example.particle_levelset_evolution.Seed_Particles(time);
    example.particle_levelset_evolution.Delete_Particles_Outside_Grid();
    
    //add forces
    example.incompressible.Set_Gravity();
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

    int extrapolation_cells=2*example.number_of_ghost_cells+2;
    ARRAY<T,TV_INT> exchanged_phi_ghost(example.mac_grid.Domain_Indices(extrapolation_cells));
    example.particle_levelset_evolution.particle_levelset.levelset.boundary->Fill_Ghost_Cells(example.mac_grid,example.particle_levelset_evolution.phi,exchanged_phi_ghost,0,time,extrapolation_cells);
    example.incompressible.Extrapolate_Velocity_Across_Interface(example.face_velocities,exchanged_phi_ghost,false,(T)example.number_of_ghost_cells,0,TV());

    example.Set_Boundary_Conditions(time); // get so CFL is correct
    if(!example.restart && example.write_output_files) Write_Output_Files(example.first_frame);
#if USE_UNIX_TIME
    gettimeofday(&since_last_frame,NULL);
#endif
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after init",0,1);
}
//#####################################################################
// Advance_To_Target_Time
//#####################################################################
template<class TV> void PLS_INTERACTIVE_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;for(int substep=1;!done;substep++){
        T dt=FLT_MAX;
        if(example.cfl>0){example.incompressible.CFL(example.face_velocities);dt=min(dt,example.particle_levelset_evolution.CFL(false,false));dt=example.cfl*dt;}
        if(example.mpi_grid) example.mpi_grid->Synchronize_Dt(dt);
        if(time+dt>=target_time){dt=target_time-time;done=true;}
        else if(time+2*dt>=target_time){dt=(T).5*(target_time-time);}
        LOG::cout<<"dt is "<<dt<<std::endl;
        ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>* advection_conservative=dynamic_cast<ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>*>(example.incompressible.advection);
        example.number_of_ghost_cells=max(3,1+example.incompressible.Real_CFL(example.face_velocities,false,false,dt));
        if(example.mpi_grid) example.mpi_grid->Synchronize_Ghost_Cells(example.number_of_ghost_cells);
        if(advection_conservative){
            advection_conservative->number_of_ghost_cells=example.number_of_ghost_cells;}

        T maximum_fluid_speed=example.face_velocities.Maxabs().Max();
        T max_particle_collision_distance=example.particle_levelset_evolution.Particle_Levelset(1).max_collision_distance_factor*example.mac_grid.dX.Max();
        example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(true,dt*maximum_fluid_speed+2*max_particle_collision_distance+(T).5*example.mac_grid.dX.Max(),10);

PHYSBAM_DEBUG_WRITE_SUBSTEP("before advection",0,1);
        T_FACE_ARRAYS_SCALAR face_velocities_ghost;face_velocities_ghost.Resize(example.incompressible.grid,example.number_of_ghost_cells,false);
        example.incompressible.boundary->Fill_Ghost_Cells_Face(example.mac_grid,example.face_velocities,face_velocities_ghost,time+dt,example.number_of_ghost_cells);

        T_ARRAYS_SCALAR phi_back(example.mac_grid.Domain_Indices(example.number_of_ghost_cells));
        example.phi_boundary->Fill_Ghost_Cells(example.mac_grid,example.particle_levelset_evolution.Particle_Levelset(1).levelset.phi,phi_back,dt,time,example.number_of_ghost_cells);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("before phi",0,1);
        example.phi_boundary_water.Use_Extrapolation_Mode(false);
        example.particle_levelset_evolution.Advance_Levelset(dt);
        ARRAY<T,TV_INT> phi_ghost(example.mac_grid.Domain_Indices(example.number_of_ghost_cells));
        //example.phi_boundary_water.Fill_Ghost_Cells(example.mac_grid,example.particle_levelset_evolution.phi,phi_ghost,dt,time,example.number_of_ghost_cells);
        //example.advection_scalar.Update_Advection_Equation_Cell(example.mac_grid,example.particle_levelset_evolution.phi,phi_ghost,example.face_velocities,example.phi_boundary_water,dt,time);
        example.phi_boundary_water.Use_Extrapolation_Mode(true);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after phi",0,1);
        example.particle_levelset_evolution.Particle_Levelset(1).Euler_Step_Particles(face_velocities_ghost,dt,time,true,true,false,false);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after advection",0,1);
        
        example.phi_boundary_water.Use_Extrapolation_Mode(true);
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time);
        PARTICLE_LEVELSET_UNIFORM<GRID<TV> >& pls=example.particle_levelset_evolution.Particle_Levelset(1);
        LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,TV> interpolation;
        if(pls.use_removed_positive_particles) for(typename GRID<TV>::NODE_ITERATOR iterator(example.mac_grid);iterator.Valid();iterator.Next()) if(pls.removed_positive_particles(iterator.Node_Index())){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_positive_particles(iterator.Node_Index());
            for(int p=1;p<=particles.array_collection->Size();p++){
                TV X=particles.X(p),V=interpolation.Clamped_To_Array_Face(example.mac_grid,face_velocities_ghost,X);
                if(-pls.levelset.Phi(X)>1.5*particles.radius(p)) V-=-TV::Axis_Vector(2)*(T).3; // buoyancy
                particles.V(p)=V;}}
        if(pls.use_removed_negative_particles) for(typename GRID<TV>::NODE_ITERATOR iterator(example.mac_grid);iterator.Valid();iterator.Next()) if(pls.removed_negative_particles(iterator.Node_Index())){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_negative_particles(iterator.Node_Index());
            for(int p=1;p<=particles.array_collection->Size();p++) particles.V(p)+=-TV::Axis_Vector(2)*dt*(T)9.8; // ballistic
            for(int p=1;p<=particles.array_collection->Size();p++) particles.V(p)+=dt*interpolation.Clamped_To_Array_Face(example.mac_grid,example.incompressible.force,particles.X(p));} // external forces
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after particles",0,1);

        example.particle_levelset_evolution.Make_Signed_Distance();
        example.phi_boundary->Fill_Ghost_Cells(example.mac_grid,pls.levelset.phi,phi_ghost,dt,time,example.number_of_ghost_cells);
        if(advection_conservative){PHYSBAM_FATAL_ERROR();}
            //advection_conservative->pls=&pls;
            //T_FACE_ARRAYS_SCALAR face_velocities_ghost(example.mac_grid,2*example.number_of_ghost_cells+1,false);
            //example.incompressible.boundary->Fill_Ghost_Cells_Face(example.mac_grid,example.face_velocities,face_velocities_ghost,time,2*example.number_of_ghost_cells+1);
            //advection_conservative->Update_Advection_Equation_Face_Lookup(example.mac_grid,phi_back,phi_ghost,example.face_velocities,face_velocities_ghost,face_velocities_ghost,*example.incompressible.boundary,dt,time,0,0,0,0);}
        else example.incompressible.Advance_One_Time_Step_Convection(dt,time,example.face_velocities,example.face_velocities,example.number_of_ghost_cells);
        //example.incompressible.Advance_One_Time_Step_Convection(dt,time,example.face_velocities,example.face_velocities,example.number_of_ghost_cells);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after convection",0,1);
        example.Advect_Particles(dt,time);
        example.incompressible.Add_Energy_With_Vorticity(example.face_velocities,example.domain_boundary,dt,time,example.number_of_ghost_cells,&pls.levelset);
        example.incompressible.Advance_One_Time_Step_Forces(example.face_velocities,dt,time,false,0,example.number_of_ghost_cells);
        example.boundary->Fill_Ghost_Cells_Face(example.mac_grid,example.face_velocities,face_velocities_ghost,time+dt,example.number_of_ghost_cells);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after forces",0,1);

        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time+dt);
        example.particle_levelset_evolution.Particle_Levelset(1).Exchange_Overlap_Particles();
        example.particle_levelset_evolution.Modify_Levelset_And_Particles(&face_velocities_ghost);
        example.particle_levelset_evolution.Make_Signed_Distance();
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after particles",0,1);

        example.Adjust_Phi_With_Sources(time+dt);
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time+dt);

        example.particle_levelset_evolution.Delete_Particles_Outside_Grid();
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete 1",0,1);
        example.particle_levelset_evolution.particle_levelset.Delete_Particles_In_Local_Maximum_Phi_Cells(1);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete 2",0,1);
        example.particle_levelset_evolution.Particle_Levelset(1).Delete_Particles_Far_From_Interface(); // uses visibility
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after delete 3",0,1);

        example.particle_levelset_evolution.Particle_Levelset(1).Identify_And_Remove_Escaped_Particles(face_velocities_ghost,1.5,time+dt);
        if(advection_conservative){
            example.particle_levelset_evolution.Particle_Levelset(1).Fix_Momentum_With_Escaped_Particles(face_velocities_ghost,advection_conservative->momentum_lost,1.5,1,time+dt);}
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after remove",0,1);
        if(example.particle_levelset_evolution.Particle_Levelset(1).use_removed_positive_particles || example.particle_levelset_evolution.Particle_Levelset(1).use_removed_negative_particles)
            example.particle_levelset_evolution.Particle_Levelset(1).Reincorporate_Removed_Particles(1,1,0,true);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after reincorporate",0,1);

        // update ghost phi values
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time+dt);

        PHYSBAM_DEBUG_WRITE_SUBSTEP("before boundary",0,1);
        example.Set_Boundary_Conditions(time);
        example.incompressible.Set_Dirichlet_Boundary_Conditions(&example.particle_levelset_evolution.phi,0);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after boundary",0,1);

        example.projection.p*=dt;
        example.projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method();
        example.incompressible.Advance_One_Time_Step_Implicit_Part(example.face_velocities,dt,time,true);
        example.projection.p*=(1/dt);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after solve",0,1);

        example.incompressible.boundary->Apply_Boundary_Condition_Face(example.incompressible.grid,example.face_velocities,time+dt);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after boundary",0,1);
        example.projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method(false);

        int band_width=example.number_of_ghost_cells+1;
        T_ARRAYS_SCALAR exchanged_phi_ghost(example.mac_grid.Domain_Indices(2*band_width+2));
        example.particle_levelset_evolution.particle_levelset.levelset.boundary->Fill_Ghost_Cells(example.mac_grid,example.particle_levelset_evolution.phi,exchanged_phi_ghost,0,time+dt,2*band_width+2);
        example.incompressible.Extrapolate_Velocity_Across_Interface(example.face_velocities,exchanged_phi_ghost,false,(T)band_width,0,TV());
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after extrapolate",0,1);

        time+=dt;}
}
//#####################################################################
// Function Update_Positions_From_Interactive
//#####################################################################
template<class TV> void PLS_INTERACTIVE_DRIVER<TV>::Update_Positions_From_Interactive()
{ 
    /*if(visualizer->current_selection && visualizer->current_selection->type==OPENGL_SELECTION::COMPONENT_RIGID_BODIES_3D){
        OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T> *real_selection=dynamic_cast<OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T>*>(visualizer->current_selection);
        example.rigid_geometry_collection.particles.X(real_selection->body_id)[1]=visualizer->opengl_world.object_location(1);
        example.rigid_geometry_collection.particles.X(real_selection->body_id)[2]=visualizer->opengl_world.object_location(2);
        example.rigid_geometry_collection.particles.X(real_selection->body_id)[3]=visualizer->opengl_world.object_location(3);
        if(rigid_component) rigid_component->Reinitialize(true,false);}*/
}
//#####################################################################
// Simulate_To_Frame
//#####################################################################
template<class TV> void Sync_Scalar(THREADED_UNIFORM_GRID<TV>& threaded_grid,ARRAY<typename TV::SCALAR,VECTOR<int,TV::dimension> >& phi,OPENGL_COMPONENT_LEVELSET_3D<typename TV::SCALAR>& component)
{threaded_grid.Sync_Scalar(phi,component.opengl_levelset_multiview->levelset->phi);}
template<class TV> void Sync_Scalar(THREADED_UNIFORM_GRID<TV>& threaded_grid,ARRAY<typename TV::SCALAR,VECTOR<int,TV::dimension> >& phi,OPENGL_COMPONENT_LEVELSET_2D<typename TV::SCALAR>& component)
{PHYSBAM_FATAL_ERROR();}
template<class TV> void Sync_Scalar(THREADED_UNIFORM_GRID<TV>& threaded_grid,ARRAY<typename TV::SCALAR,VECTOR<int,TV::dimension> >& phi,OPENGL_COMPONENT_LEVELSET_1D<typename TV::SCALAR>& component)
{PHYSBAM_FATAL_ERROR();}
template<class TV> void PLS_INTERACTIVE_DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    while(frame<0 || current_frame<frame){
        Advance_To_Target_Time(example.Time_At_Frame(current_frame+1));
        if((current_frame-example.first_frame)%1==0){
            example.particle_levelset_evolution.Reseed_Particles(time);
            example.particle_levelset_evolution.Delete_Particles_Outside_Grid();}
        if(example.mpi_grid && example.mpi_grid->threaded_grid && visualizer){
            Sync_Scalar(*example.mpi_grid->threaded_grid,example.particle_levelset_evolution.particle_levelset.levelset.phi,*levelset_component);}
#ifdef USE_PTHREADS
        if(example.mpi_grid && example.mpi_grid->threaded_grid) pthread_barrier_wait(example.mpi_grid->threaded_grid->barr);
#endif
        if(((example.mpi_grid && example.mpi_grid->threaded_grid && example.mpi_grid->threaded_grid->tid==1) || !example.mpi_grid || !example.mpi_grid->threaded_grid) && visualizer){
#if USE_UNIX_TIME
            gettimeofday(&now,NULL);
            while((now.tv_sec + now.tv_usec/1000000.0)-(since_last_frame.tv_sec+since_last_frame.tv_usec/1000000.0)<(double)1/24) gettimeofday(&now,NULL);
            since_last_frame.tv_sec=now.tv_sec;
            since_last_frame.tv_usec=now.tv_usec;
#endif
            visualizer->Set_Frame(current_frame);}
        if(example.write_output_files) Write_Output_Files(++output_number);
        current_frame++;}
} 
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void PLS_INTERACTIVE_DRIVER<TV>::
Write_Substep(const std::string& title,const int substep,const int level)
{
    if(level<=example.write_substeps_level){
        example.frame_title=title;
        LOG::cout<<"Writing substep ["<<example.frame_title<<"]: output_number="<<output_number+1<<", time="<<time<<", frame="<<current_frame<<", substep="<<substep<<std::endl;
        Write_Output_Files(++output_number);example.frame_title="";}
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void PLS_INTERACTIVE_DRIVER<TV>::
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
template class PLS_INTERACTIVE_DRIVER<VECTOR<float,2> >;
template class PLS_INTERACTIVE_DRIVER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PLS_INTERACTIVE_DRIVER<VECTOR<double,2> >;
template class PLS_INTERACTIVE_DRIVER<VECTOR<double,3> >;
#endif

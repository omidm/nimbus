//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Jon Gretarsson, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Sergey Levine, Nick Rasmussen, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_DRIVER_CHIMERA
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Utilities/INTERRUPTS.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/FRACTURE_EVOLUTION.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Fluids/Coupled_Evolution/COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/DETONATION_SHOCK_DYNAMICS.h>
#include <PhysBAM_Dynamics/Boundaries/BOUNDARY_PHI_WATER.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA.h>
#include <PhysBAM_Dynamics/Incompressible_Flows/INCOMPRESSIBLE_MULTIPHASE_UNIFORM.h>
#include <PhysBAM_Dynamics/Incompressible_Flows/SPH_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/VOF_ADVECTION.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_UNIFORM_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_CHIMERA.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_CHIMERA.h>
using namespace PhysBAM;
//#####################################################################
// Macros for chimera grid and data access
//#####################################################################
#define FACE_VELOCITIES(I) example.incompressible_fluid_containers(I)->face_velocities
#define PSI_N(I) example.incompressible_fluid_containers(I)->psi_N
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> SOLIDS_FLUIDS_DRIVER_CHIMERA<T_GRID>::
SOLIDS_FLUIDS_DRIVER_CHIMERA(SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>& example)
    :SOLIDS_FLUIDS_DRIVER<TV>(example),example(example),last_dt(0),restart_dt(0),reset_with_restart(false),evolution(new SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>(example))
{
    evolution->driver=this;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> SOLIDS_FLUIDS_DRIVER_CHIMERA<T_GRID>::
~SOLIDS_FLUIDS_DRIVER_CHIMERA()
{
    //delete evolution;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_CHIMERA<T_GRID>::
Initialize()
{
    SOLIDS_FLUIDS_DRIVER<TV>::Initialize();
    
    example.chimera_grid->Get_Pointer_Of_Incompressible_Fluid_Containers(&(example.incompressible_fluid_containers));

    example.fluids_parameters.Set_Fluids_Parameters_Callbacks(example);
    
    SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
    solids_evolution.Set_Solids_Evolution_Callbacks(example);
    example.Initialize_Bodies();
    solids_evolution.Initialize_Rigid_Bodies(example.frame_rate,example.restart);
    example.Update_Object_Space_Rigid_Body_Collection(false);

    //initialize valid masks
    int number_of_ghost_cells=example.fluids_parameters.number_of_ghost_cells;
    if(example.use_collidable_advection) for(int grid_index=1;grid_index<=example.rigid_grid_collection.particles.array_collection->Size();grid_index++){
        example.incompressible_fluid_containers(grid_index)->face_velocities_valid_mask.Resize(example.rigid_grid_collection.Rigid_Grid(grid_index).grid.Domain_Indices(number_of_ghost_cells));
        if(example.fluids_parameters.smoke){
            example.incompressible_fluid_containers(grid_index)->density_container.valid_mask_current.Resize(example.rigid_grid_collection.Rigid_Grid(grid_index).grid.Domain_Indices(number_of_ghost_cells));
            example.incompressible_fluid_containers(grid_index)->density_container.valid_mask_next.Resize(example.rigid_grid_collection.Rigid_Grid(grid_index).grid.Domain_Indices(number_of_ghost_cells));}
        if(example.fluids_parameters.water)
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.levelset.Initialize_Valid_Masks(example.rigid_grid_collection.Rigid_Grid(grid_index).grid);}

    example.Parse_Late_Options();

    if(example.restart){
        example.Read_Output_Files_Solids(example.restart_frame);
        example.Read_Output_Files_Fluids(example.restart_frame);}
    
    if(example.use_collidable_advection) for(int grid_index=1;grid_index<=example.rigid_grid_collection.particles.array_collection->Size();grid_index++){
        example.incompressible_fluid_containers(grid_index)->collision_bodies_affecting_fluid.Initialize_Grids();
        example.incompressible_fluid_containers(grid_index)->collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false);
        example.incompressible_fluid_containers(grid_index)->collision_bodies_affecting_fluid.Rasterize_Objects();
        example.incompressible_fluid_containers(grid_index)->collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*example.rigid_grid_collection.Rigid_Grid(grid_index).grid.Minimum_Edge_Length(),5);
        example.incompressible_fluid_containers(grid_index)->collision_bodies_affecting_fluid.Compute_Grid_Visibility();}

    if(example.solve_heat_equation || example.solve_heat_equation_on_faces || example.solve_heat_equation_on_faces_coupled){
        if(!example.restart)
            example.Initialize_Velocities();}
    if(example.fluids_parameters.smoke){
        if(!example.restart){
            example.Initialize_Density();
            example.Initialize_Velocities();}}
    else if(example.fluids_parameters.water){
        example.particles_counter.Fill(0);
        for(int grid_index=1;grid_index<=example.rigid_grid_collection.particles.array_collection->Size();grid_index++){
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Set_Time(time);
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Set_CFL_Number(example.fluids_parameters.cfl);
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Set_Number_Particles_Per_Cell(example.fluids_parameters.number_particles_per_cell);
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.Set_Band_Width((T)2*example.fluids_parameters.particle_half_bandwidth);
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Set_Levelset_Callbacks(example);
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Initialize_FMM_Initialization_Iterative_Solver(example.fluids_parameters.refine_fmm_initialization_with_iterative_solver);
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Bias_Towards_Negative_Particles(example.fluids_parameters.bias_towards_negative_particles);
            if(example.fluids_parameters.use_removed_positive_particles) example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Particle_Levelset(grid_index).Use_Removed_Positive_Particles();
            if(example.fluids_parameters.use_removed_negative_particles) example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Particle_Levelset(grid_index).Use_Removed_Negative_Particles();
            if(example.fluids_parameters.store_particle_ids)
                example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Particle_Levelset(grid_index).Store_Unique_Particle_Id();
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Use_Particle_Levelset(example.fluids_parameters.use_particle_levelset);}

        if(!example.restart){
            example.Initialize_Phi();
            example.Adjust_Phi_With_Sources(time);}
  
        for(int grid_index=1;grid_index<=example.rigid_grid_collection.particles.array_collection->Size();grid_index++){
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Initialize_Domain(example.incompressible_fluid_containers(grid_index)->rigid_grid.grid);
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.levelset_advection.Set_Custom_Advection(*(example.advection(grid_index)));
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.levelset.Set_Custom_Boundary(*(evolution->boundary)); //yuey
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.levelset.Set_Collision_Body_List(example.incompressible_fluid_containers(grid_index)->collision_bodies_affecting_fluid);
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.levelset.Set_Face_Velocities_Valid_Mask(&example.incompressible_fluid_containers(grid_index)->face_velocities_valid_mask);
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.Set_Collision_Distance_Factors(example.fluids_parameters.min_collision_distance_factor,example.fluids_parameters.max_collision_distance_factor);}
        ARRAY<T_ARRAYS_SCALAR*> scalar_array(example.rigid_grid_collection.particles.array_collection->Size());
        for(int grid_index=1;grid_index<=example.rigid_grid_collection.particles.array_collection->Size();grid_index++){
            scalar_array(grid_index)=&example.incompressible_fluid_containers(grid_index)->density_container.density;}

        example.Fill_Ghost_Cells_Chimera(2*example.fluids_parameters.number_of_ghost_cells+1,scalar_array,scalar_array);

        for(int grid_index=1;grid_index<=example.rigid_grid_collection.particles.array_collection->Size();grid_index++){
            if(!example.restart){
                //example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Use_Reinitialization();
                example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Use_Fast_Marching_Method();
                example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Make_Signed_Distance();}
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Set_Seed(2606);
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Seed_Particles(time);
            example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.Delete_Particles_Outside_Grid();
            if(example.fluids_parameters.delete_fluid_inside_objects) example.fluids_parameters.Delete_Particles_Inside_Objects(time);}
        
        example.Fill_Ghost_Cells_Chimera(example.fluids_parameters.number_of_ghost_cells,scalar_array,scalar_array);
    
        if(!example.restart) example.Initialize_Velocities();}

    example.Update_Grid(time,0);
}
//#####################################################################
// Function Advance_To_Target_Time
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_CHIMERA<T_GRID>::
Advance_To_Target_Time(const T target_time)
{
    //T dt_full_advance=target_time-time;

    LOG::SCOPE scope("Advance_To_Target_Time");

    bool done=false;for(int substep=1;!done;substep++){
        LOG::SCOPE scope("SUBSTEP",STRING_UTILITIES::string_sprintf("substep %d",substep));
        example.Get_Source_Velocities_Chimera(time,0);
        T dt=Compute_Dt(time,target_time,done);
        example.Preprocess_Substep(dt,time);
        if(example.run_advection_test) example.Get_Source_Velocities_Chimera(time,dt);//only makes a difference in Enright advection tests
        evolution->Advance_One_Time_Step(time,dt);
        last_dt=restart_dt?restart_dt:dt;
        example.Postprocess_Substep(last_dt,time);
        time+=last_dt;restart_dt=0;

        Write_Substep(STRING_UTILITIES::string_sprintf("END Substep %d",substep),substep,0);}
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_CHIMERA<T_GRID>::
Postprocess_Frame(const int frame)
{
    example.Postprocess_Frame(frame);
    if(example.fluids_parameters.water){
        for(int grid_index=1;grid_index<=example.rigid_grid_collection.particles.array_collection->Size();grid_index++){
            PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>& particle_levelset_evolution=example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution;
            //SOLIDS_EVOLUTION<TV>& solids_evolution=*example.solids_evolution;
            //int number_of_regions=example.fluids_parameters.number_of_regions;
            
            //if(number_of_regions>=1){
            //  example.Postprocess_Phi(time); // actually changes phi !!!
            
            //if(example.fluids_parameters.mass_conservation && (frame-example.first_frame)%example.fluids_parameters.reinitialize_geometry_frame_rate==0){
            //  LOG::Time("Reinitializing VOF geometry...");
            //  particle_levelset_evolution->Reinitialize_Geometry(number_of_regions);}
            LOG::Time("Reseeding... ");
            if(particle_levelset_evolution.use_particle_levelset && ((frame-example.first_frame)%(example.fluids_parameters.reseeding_frame_rate)==0)){
                particle_levelset_evolution.Reseed_Particles(time);
                particle_levelset_evolution.Delete_Particles_Outside_Grid();
                if(example.fluids_parameters.delete_fluid_inside_objects) example.fluids_parameters.Delete_Particles_Inside_Objects(time);}
            LOG::Stop_Time();}}
    
    //if(example.fluids_parameters.monitor_mass){
    //  if(number_of_regions==1){
    //      T mass_new=particle_levelset_evolution->levelset_advection.Approximate_Negative_Material();
    //      LOG::cout<<"Material = "<<mass_new<<" - change = "<<mass_new-example.fluids_parameters.mass<<std::endl;
    //      example.fluids_parameters.mass=mass_new;}
    //  else if(number_of_regions>=2){
    //      for(int i=1;i<=example.fluids_parameters.number_of_regions;i++){
    //          T mass_new=particle_levelset_evolution_multiple->Levelset_Advection(i).Approximate_Negative_Material();
    //          LOG::cout<<"Region "<<i<<": Material = "<<mass_new<<", change = "<<mass_new-example.fluids_parameters.masses(i)<<std::endl;
    //          example.fluids_parameters.masses(i)=mass_new;}}}
    
    //SOLIDS_FLUIDS_DRIVER<TV>::Postprocess_Frame(frame);
    //solids_evolution.Postprocess_Frame(frame);
    
    //if(solids_evolution.solids_parameters.rigid_body_collision_parameters.rigid_collisions_print_interpenetration_statistics)
    //  solids_evolution.rigid_body_collisions->Print_Interpenetration_Statistics();
}
//#####################################################################
// Function Compute_Dt
//#####################################################################
template<class T_GRID> typename T_GRID::VECTOR_T::SCALAR SOLIDS_FLUIDS_DRIVER_CHIMERA<T_GRID>::
Compute_Dt(const T time,const T target_time,bool& done)
{
    //COMPUTE DT HERE
    T dt,min_dt=0; //NEED TO SET min_dt

    //if(example.fixed_dt)
    //  dt=example.fixed_dt;
    //else
    dt=evolution->Calculate_Maximum_Allowable_dt(time);
    LOG::cout << "maximum allowable dt " << dt << std::endl;
    if(example.abort_when_dt_below && dt<example.abort_when_dt_below) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("dt too small (%g < %g)",dt,example.abort_when_dt_below));
    done=false;
    SOLIDS_FLUIDS_EXAMPLE<TV>::Clamp_Time_Step_With_Target_Time(time,target_time,dt,done,min_dt);
    LOG::cout << "time " << time << " dt " << dt << std::endl;
    return dt;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_DRIVER_CHIMERA<T_GRID>::
Write_Output_Files(const int frame)
{
    LOG::SCOPE scope("writing output files");
    FILE_UTILITIES::Create_Directory(example.output_directory);
    FILE_UTILITIES::Create_Directory(example.output_directory+STRING_UTILITIES::string_sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
    Write_First_Frame(frame);

    example.Write_Output_Files(frame);

    //WRITE DATA HERE

    Write_Time(frame);
    Write_Last_Frame(frame);
}
//#####################################################################
template class SOLIDS_FLUIDS_DRIVER_CHIMERA<GRID<VECTOR<float,1> > >;
template class SOLIDS_FLUIDS_DRIVER_CHIMERA<GRID<VECTOR<float,2> > >;
template class SOLIDS_FLUIDS_DRIVER_CHIMERA<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SOLIDS_FLUIDS_DRIVER_CHIMERA<GRID<VECTOR<double,1> > >;
template class SOLIDS_FLUIDS_DRIVER_CHIMERA<GRID<VECTOR<double,2> > >;
template class SOLIDS_FLUIDS_DRIVER_CHIMERA<GRID<VECTOR<double,3> > >;
#endif

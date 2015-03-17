//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/DOMAIN_ITERATOR_THREADED.h>
#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_THREADED.h>
#include <PhysBAM_Tools/Vectors/VECTOR_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM.h>
#include <PhysBAM_Dynamics/Boundaries/BOUNDARY_PHI_WATER.h>
#include "application/physbam-app/translator_test/data_face_array.h"
#include "application/physbam-app/translator_test/data_particle_array.h"
#include "application/physbam-app/translator_test/face_array_test.h"
#include "application/physbam-app/translator_test/particle_test.h"
#include "application/physbam-app/translator_test/water_driver.h"
#include "application/physbam-app/translator_test/water_example.h"
#include "shared/fast_log.hh"
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"
#include "worker/data.h"
using namespace PhysBAM;
namespace{
template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
{
    ((WATER_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
}
};
//#####################################################################
// Initialize
//#####################################################################
template<class TV> WATER_DRIVER<TV>::
WATER_DRIVER(WATER_EXAMPLE<TV>& example)
    :example(example),kinematic_evolution(example.rigid_geometry_collection,true),thread_queue(example.thread_queue)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> WATER_DRIVER<TV>::
~WATER_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Initialize()
{
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);

    // setup time
    if(example.restart) current_frame=example.restart;else current_frame=example.first_frame;
    output_number=current_frame;
    time=example.Time_At_Frame(current_frame);
    
    // initialize collision objects
    kinematic_evolution.Get_Current_Kinematic_Keyframes(1/example.frame_rate,time);
    kinematic_evolution.Set_External_Positions(example.rigid_geometry_collection.particles.X,example.rigid_geometry_collection.particles.rotation,time);
    kinematic_evolution.Set_External_Velocities(example.rigid_geometry_collection.particles.V,example.rigid_geometry_collection.particles.angular_velocity,time,time);

    example.phi_boundary_water.Set_Velocity_Pointer(example.face_velocities);

    {
        example.particle_levelset_evolution.Initialize_Domain(example.mac_grid);
        example.particle_levelset_evolution.particle_levelset.Set_Band_Width(6);
        example.incompressible.Initialize_Grids(example.mac_grid);
        example.projection.Initialize_Grid(example.mac_grid);
        example.collision_bodies_affecting_fluid.Initialize_Grids();
    }
    example.face_velocities.Resize(example.mac_grid);

    example.particle_levelset_evolution.Set_Time(time);
    example.particle_levelset_evolution.Set_CFL_Number((T).9);

    if(example.mpi_grid) example.mpi_grid->Initialize(example.domain_boundary);
    example.incompressible.mpi_grid=example.mpi_grid;
    example.projection.elliptic_solver->mpi_grid=example.mpi_grid;
    example.particle_levelset_evolution.particle_levelset.mpi_grid=example.mpi_grid;
    if(example.mpi_grid){
        example.boundary=new BOUNDARY_MPI<GRID<TV> >(example.mpi_grid,example.boundary_scalar);
        example.phi_boundary=new BOUNDARY_MPI<GRID<TV> >(example.mpi_grid,example.phi_boundary_water);
        example.particle_levelset_evolution.particle_levelset.last_unique_particle_id=example.mpi_grid->rank*30000000;}
    else{
        example.boundary=&example.boundary_scalar;
        example.phi_boundary=&example.phi_boundary_water;}

    if(PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<GRID<TV> > *refine=dynamic_cast<PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<GRID<TV> >*>(&example.projection)){
        refine->boundary=example.boundary;
        refine->phi_boundary=example.phi_boundary;}
    example.rigid_geometry_collection.Update_Kinematic_Particles();

    VECTOR<VECTOR<bool,2>,TV::dimension> domain_open_boundaries=VECTOR_UTILITIES::Complement(example.domain_boundary);
    example.phi_boundary->Set_Constant_Extrapolation(domain_open_boundaries);
    example.boundary->Set_Constant_Extrapolation(domain_open_boundaries);
    if(thread_queue){
        ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<GRID<TV>,T>* threaded_advection_scalar=new ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<GRID<TV>,T>(thread_queue);
        example.particle_levelset_evolution.Levelset_Advection(1).Set_Custom_Advection(*threaded_advection_scalar);
        example.incompressible.Set_Custom_Advection(*threaded_advection_scalar);
        example.particle_levelset_evolution.particle_levelset.Set_Thread_Queue(thread_queue);
        example.particle_levelset_evolution.particle_levelset.levelset.thread_queue=thread_queue;
        if(PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<GRID<TV> >* refinement=dynamic_cast<PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<GRID<TV> >*>(&example.projection))
            refinement->thread_queue=thread_queue;}
    else{
        example.particle_levelset_evolution.Levelset_Advection(1).Set_Custom_Advection(example.advection_scalar);
        example.incompressible.Set_Custom_Advection(example.advection_scalar);}

    example.particle_levelset_evolution.Set_Number_Particles_Per_Cell(16);
    example.particle_levelset_evolution.Set_Levelset_Callbacks(example);
    example.particle_levelset_evolution.Initialize_FMM_Initialization_Iterative_Solver(true);

    example.particle_levelset_evolution.particle_levelset.levelset.Set_Custom_Boundary(*example.phi_boundary);
    example.particle_levelset_evolution.Bias_Towards_Negative_Particles(false);
    example.particle_levelset_evolution.particle_levelset.Use_Removed_Positive_Particles();
    example.particle_levelset_evolution.particle_levelset.Use_Removed_Negative_Particles();
    example.particle_levelset_evolution.particle_levelset.Store_Unique_Particle_Id();
    example.particle_levelset_evolution.Use_Particle_Levelset(true);
    example.particle_levelset_evolution.particle_levelset.levelset.Set_Collision_Body_List(example.collision_bodies_affecting_fluid);
    example.particle_levelset_evolution.particle_levelset.levelset.Set_Face_Velocities_Valid_Mask(&example.incompressible.valid_mask);
    example.particle_levelset_evolution.particle_levelset.Set_Collision_Distance_Factors(.1,1);

    example.incompressible.Set_Custom_Boundary(*example.boundary);
    example.incompressible.projection.elliptic_solver->Set_Relative_Tolerance(1e-8);
    example.incompressible.projection.elliptic_solver->pcg.Set_Maximum_Iterations(40);
    example.incompressible.projection.elliptic_solver->pcg.evolution_solver_type=krylov_solver_cg;
    example.incompressible.projection.elliptic_solver->pcg.cg_restart_iterations=0;
    example.incompressible.projection.elliptic_solver->pcg.Show_Results();
    example.incompressible.projection.collidable_solver->Use_External_Level_Set(example.particle_levelset_evolution.particle_levelset.levelset);

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

    ARRAY<T,TV_INT> exchanged_phi_ghost(example.mac_grid.Domain_Indices(8));
    example.particle_levelset_evolution.particle_levelset.levelset.boundary->Fill_Ghost_Cells(example.mac_grid,example.particle_levelset_evolution.phi,exchanged_phi_ghost,0,time,8);
    example.incompressible.Extrapolate_Velocity_Across_Interface(example.face_velocities,exchanged_phi_ghost,false,3,0,TV());

    example.Set_Boundary_Conditions(time); // get so CFL is correct
    if(!example.restart) {
        printf("\nWrite initialized data ...\n");
        Write_Output_Files(example.first_frame);
    }

    // extra variables
    example.face_velocities_ghost_flag.Resize(example.incompressible.grid,example.number_of_ghost_cells,false);
    example.face_velocities_ghost.Resize(example.incompressible.grid,example.number_of_ghost_cells,false);
    example.incompressible.boundary->Fill_Ghost_Cells_Face(example.mac_grid,example.face_velocities,example.face_velocities_ghost,time,example.number_of_ghost_cells);
}
//#####################################################################
// Test Particles
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
TestParticles()
{
    printf("\nTesting particles ...\n\n");

    // Define/ create data
    int start[] = {-2, 1, 4, scale-2, scale+1};
    int delta[] = {3, 3, scale-6, 3, 3};
    nimbus::DataArray particle_array(5*5*5-1);
    nimbus::DataArray particle_out;
    nimbus::DataArray particle_in;
    size_t i = 0;
    for (size_t x = 0; x < 5; ++x) {
        for (size_t y = 0; y < 5; ++y) {
            for (size_t z = 0; z < 5; ++z) {
                if (! (x == y && y == z) ) {
                    particle_array[i] = new test::DataParticleArray("pos-particles");
                    particle_array[i]->set_region(
                            nimbus::GeometricRegion(start[x], start[y], start[z],
                                                    delta[x], delta[y], delta[z]));
                    particle_array[i]->Create();
                    if (x == 0 || x == 4 || y == 0 || y == 4 || z == 0 || z == 4) {
                        particle_out.push_back(particle_array[i]);
                    } else {
                        particle_in.push_back(particle_array[i]);
                    }
                    i++;
                }
            }
        }
    }

    // Particle test class
    test::ParticleTest particle_test;
    particle_test.shift = nimbus::Coord(0, 0, 0);
    particle_test.loc_region = nimbus::GeometricRegion(1, 1, 1, scale, scale, scale);
    particle_test.enl_region = nimbus::GeometricRegion(-2, -2, -2, scale+6, scale+6, scale+6);

    for (size_t t = 0; t < 10; ++t) {
        particle_test.WriteParticles(particle_test.enl_region, particle_out,
                                     &example.particle_levelset_evolution.particle_levelset,
                                     true);
        particle_test.DeleteOutsideParticles(&example.particle_levelset_evolution.particle_levelset, true);
        particle_test.ReadParticles(particle_test.enl_region, particle_out,
                                    &example.particle_levelset_evolution.particle_levelset,
                                    true);
    }

    // TODO(chinmayee): add test for delete in array of regions
}
//#####################################################################
// Test FaceArray
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
TestFaceArray()
{
    printf("\nTesting face arrays ...\n\n");

    // Define/ create data
    int start[] = {-2, 1, 4, scale-2, scale+1};
    int delta[] = {3, 3, scale-6, 3, 3};
    nimbus::DataArray face_array_all(5*5*5-1);
    nimbus::DataArray face_array_out;
    nimbus::DataArray face_array_in;
    size_t i = 0;
    for (size_t x = 0; x < 5; ++x) {
        for (size_t y = 0; y < 5; ++y) {
            for (size_t z = 0; z < 5; ++z) {
                if (! (x == y && y == z) ) {
                    face_array_all[i] = new test::DataFaceArray<float>("facearray");
                    face_array_all[i]->set_region(
                            nimbus::GeometricRegion(start[x], start[y], start[z],
                                                    delta[x], delta[y], delta[z]));
                    face_array_all[i]->Create();
                    if (x == 0 || x == 4 || y == 0 || y == 4 || z == 0 || z == 4) {
                        face_array_out.push_back(face_array_all[i]);
                    } else {
                        face_array_in.push_back(face_array_all[i]);
                    }
                    i++;
                }
            }
        }
    }

    // face_array test class
    test::FaceArrayTest<float> face_array_test;
    face_array_test.shift = nimbus::Coord(0, 0, 0);
    face_array_test.loc_region = nimbus::GeometricRegion(1, 1, 1, scale, scale, scale);
    face_array_test.enl_region = nimbus::GeometricRegion(-2, -2, -2, scale+6, scale+6, scale+6);

//    for (size_t t = 0; t < 10; ++t) {
//        face_array_test.WriteFaceArray(face_array_test.enl_region, face_array_out,
//                                       &example.face_velocities_ghost);
//        face_array_test.ReadFaceArray(face_array_test.enl_region, face_array_out,
//                                      &example.face_velocities_ghost);
//    }


    for (size_t t = 0; t < 10; ++t) {
        face_array_test.WriteFaceArray(face_array_test.enl_region, face_array_in,
                                       &example.face_velocities_ghost);
        example.face_velocities_ghost_flag.Fill(0);
        face_array_test.ReadFaceArray(face_array_test.enl_region, face_array_in,
                                      &example.face_velocities_ghost, example.face_velocities_ghost_flag);
    }
}
//#####################################################################
// Test
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Test()
{
    Initialize();

    nimbus::timer::InitializeKeys();
    nimbus::timer::InitializeTimers();

    printf("\nTesting translator ...\n");

    TestParticles();
    TestFaceArray();

    std::string file_name = "times.txt";
    FILE* temp = fopen(file_name.c_str(), "w");
    nimbus::timer::PrintTimerSummary(temp);
    fclose(temp);

    printf("\n... Test translator complete\n\n");
} 
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
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
template<class TV> void WATER_DRIVER<TV>::
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
template class WATER_DRIVER<VECTOR<float,3> >;

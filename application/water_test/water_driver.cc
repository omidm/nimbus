/*
 * Copyright 2013 Stanford University.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the
 *   distribution.
 *
 * - Neither the name of the copyright holders nor the names of
 *   its contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
 * THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "./water_app.h"
#include "./water_driver.h"

using namespace PhysBAM;    // NOLINT

template <class TV> WaterDriver<TV> ::
WaterDriver()
{
    //TODO: Initialize the example here
    //TODO: initialize all data and corresponding pointers

    mac_grid->initialize(TV_INT(), RANGE<TV>::Unit_Box(), true);

    stream_type = new STREAM_TYPE(float());

    // setup time
    initial_time = 0;
    first_frame = 0;
    last_frame = 100;
    frame_rate = 24;
    current_frame = 0;
    output_number = 0;
    time = Time_At_Frame(current_frame);

    // other parameters
    write_substeps_level = -1;
    write_output_files_flag = true;
    number_of_ghost_cells = 3;
    cfl = 0.9;
    mpi_grid.data = NULL;
}

template <class TV> WaterDriver<TV> ::
~WaterDriver()
{
}

template <class TV> void Write_Substep_Helper(
        void *writer,
        const std::string &title,
        int substep,
        int level)
{
    ((WaterDriver<TV>*)writer)->Write_Substep(title, substep, level);
}

template <class TV> void WaterDriver<TV>::
initialize(bool distributed)
{
    // initialize mpi grid and file names
    if (distributed)
    {
        mpi_grid.data = new MPI_UNIFORM_GRID<GRID<TV> >(
                *mac_grid.data, 
                number_of_ghost_cells);
    }

    //TODO: eventually change the way output is saved
    output_directory = "output";
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this, &Write_Substep_Helper<TV>);
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(write_substeps_level);
    if (distributed)
    {
        if (mpi_grid->Number_Of_Processors() > 1)
            output_directory += STRING_UTILITIES::string_sprintf("/%d",
                    (mpi_grid.data->rank+1));
    }
    FILE_UTILITIES::Create_Directory(output_directory + "/common");
    LOG::Instance()->Copy_Log_To_File(output_directory +
            "/common/log.txt", false);

    // initialize collision objects
    sim_data.kinematic_evolution.Get_Current_Kinematic_Keyframes(
            1/frame_rate,
            time);
    sim_data.kinematic_evolution.Set_External_Positions(
            sim_data.rigid_geometry_collection.particles.X,
            sim_data.rigid_geometry_collection.particles.rotation,
            time);
    sim_data.kinematic_evolution.Set_External_Velocities(
            sim_data.rigid_geometry_collection.particles.V,
            sim_data.rigid_geometry_collection.particles.angular_velocity,
            time,
            time);

    sim_data.phi_boundary_water.Set_Velocity_Pointer(*face_velocities.data);

    // initialize mac grid
    sim_data.particle_levelset_evolution.Initialize_Domain(
            *sim_data.mac_grid.data);
    sim_data.particle_levelset_evolution.particle_levelset.Set_Band_Width(6);
    sim_data.incompressible.Initialize_Grids(*sim_data.mac_grid.data);
    sim_data.projection.Initialize_Grid(*sim_data.mac_grid.data);
    sim_data.collision_bodies_affecting_fluid.Initialize_Grids();
    sim_data.face_velocities.Resize(*sim_data.mac_grid.data);

    sim_data.particle_levelset_evolution.Set_Time(time);
    sim_data.particle_levelset_evolution.Set_CFL_Number((T).9);

    if (distributed)
        sim_data.mpi_grid->Initialize(sim_data.domain_boundary);
    sim_data.incompressible.mpi_grid = mpi_grid.data;
    sim_data.projection.elliptic_solver->mpi_grid = mpi_grid.data;
    sim_data.particle_levelset_evolution.particle_levelset.mpi_grid =
        mpi_grid.data;
    if (distributed)
    {
        sim_data.boundary = new BOUNDARY_MPI<GRID<TV> >(
                mpi_grid.data,
                sim_data.boundary_scalar);
        sim_data.phi_boundary = new BOUNDARY_MPI<GRID<TV> >(
                mpi_grid.data,
                sim_data.phi_boundary_water);
        sim_data.particle_levelset_evolution.particle_levelset.
            last_unique_particle_id = mpi_grid.data->rank*30000000;
    }
    else
    {
        sim_data.boundary = &sim_data.boundary_scalar;
        sim_data.phi_boundary = &sim_data.phi_boundary_water;
    }

    if (PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<GRID<TV> > *refine
            = dynamic_cast
            <PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<GRID<TV> >*>
            (&sim_data.projection))
    {
        refine->boundary = sim_data.boundary;
        refine->phi_boundary = sim_data.phi_boundary;
    }
    sim_data.rigid_geometry_collection.Update_Kinematic_Particles();

    VECTOR<VECTOR<bool,2>,TV::dimension> domain_open_boundaries =
        VECTOR_UTILITIES::Complement(sim_data.domain_boundary);
    sim_data.phi_boundary->Set_Constant_Extrapolation(domain_open_boundaries);
    sim_data.boundary->Set_Constant_Extrapolation(domain_open_boundaries);

    sim_data.particle_levelset_evolution.Levelset_Advection(1).
        Set_Custom_Advection(sim_data.advection_scalar);
    sim_data.incompressible.Set_Custom_Advection(sim_data.advection_scalar);

    sim_data.particle_levelset_evolution.Set_Number_Particles_Per_Cell(16);
    sim_data.particle_levelset_evolution.Set_Levelset_Callbacks(this);
    sim_data.particle_levelset_evolution.
        Initialize_FMM_Initialization_Iterative_Solver(true);

    sim_data.particle_levelset_evolution.particle_levelset.levelset.
        Set_Custom_Boundary(*sim_data.phi_boundary);
    sim_data.particle_levelset_evolution.
        Bias_Towards_Negative_Particles(false);
    sim_data.particle_levelset_evolution.particle_levelset.
        Use_Removed_Positive_Particles();
    sim_data.particle_levelset_evolution.particle_levelset.
        Use_Removed_Negative_Particles();
    sim_data.particle_levelset_evolution.particle_levelset.
        Store_Unique_Particle_Id();
    sim_data.particle_levelset_evolution.Use_Particle_Levelset(true);
    sim_data.particle_levelset_evolution.particle_levelset.levelset.
        Set_Collision_Body_List(sim_data.collision_bodies_affecting_fluid);
    sim_data.particle_levelset_evolution.particle_levelset.levelset.
        Set_Face_Velocities_Valid_Mask(&sim_data.incompressible.valid_mask);
    sim_data.particle_levelset_evolution.particle_levelset.
        Set_Collision_Distance_Factors(.1,1);

    sim_data.incompressible.Set_Custom_Boundary(*sim_data.boundary);
    sim_data.incompressible.projection.elliptic_solver->
        Set_Relative_Tolerance(1e-8);
    sim_data.incompressible.projection.elliptic_solver->
        pcg.Set_Maximum_Iterations(40);
    sim_data.incompressible.projection.elliptic_solver->
        pcg.evolution_solver_type = krylov_solver_cg;
    sim_data.incompressible.projection.elliptic_solver->
        pcg.cg_restart_iterations=0;
    sim_data.incompressible.projection.elliptic_solver->
        pcg.Show_Results();
    sim_data.incompressible.projection.collidable_solver->
        Use_External_Level_Set(sim_data.particle_levelset_evolution.
                particle_levelset.levelset);

    sim_data.collision_bodies_affecting_fluid.
        Update_Intersection_Acceleration_Structures(false);
    sim_data.collision_bodies_affecting_fluid.Rasterize_Objects();
    sim_data.collision_bodies_affecting_fluid.
        Compute_Occupied_Blocks(false,
            (T)2*sim_data.mac_grid.Minimum_Edge_Length(), 5);
    // initialize_phi();
    sim_data.Adjust_Phi_With_Sources(time);
    sim_data.particle_levelset_evolution.Make_Signed_Distance();
    sim_data.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time);

    sim_data.collision_bodies_affecting_fluid.Compute_Grid_Visibility();
    sim_data.particle_levelset_evolution.Set_Seed(2606);
    sim_data.particle_levelset_evolution.Seed_Particles(time);
    sim_data.particle_levelset_evolution.Delete_Particles_Outside_Grid();

    //add forces
    sim_data.incompressible.Set_Gravity();
    sim_data.incompressible.Set_Body_Force(true);
    sim_data.incompressible.projection.Use_Non_Zero_Divergence(false);
    sim_data.incompressible.projection.elliptic_solver->
        Solve_Neumann_Regions(true);
    sim_data.incompressible.projection.elliptic_solver->
        solve_single_cell_neumann_regions = false;
    sim_data.incompressible.Use_Explicit_Part_Of_Implicit_Viscosity(false);
    sim_data.incompressible.Set_Maximum_Implicit_Viscosity_Iterations(40);
    sim_data.incompressible.Use_Variable_Vorticity_Confinement(false);
    sim_data.incompressible.Set_Surface_Tension(0);
    sim_data.incompressible.Set_Variable_Surface_Tension(false);
    sim_data.incompressible.Set_Viscosity(0);
    sim_data.incompressible.Set_Variable_Viscosity(false);
    sim_data.incompressible.projection.Set_Density(1e3);

    ARRAY<T,TV_INT> exchanged_phi_ghost(mac_grid.data->Domain_Indices(8));
    sim_data.particle_levelset_evolution.particle_levelset.levelset.
        boundary->Fill_Ghost_Cells(
                sim_data.mac_grid,
                sim_data.particle_levelset_evolution.phi,
                exchanged_phi_ghost,
                0, time, 8);
    sim_data.incompressible.Extrapolate_Velocity_Across_Interface(
            sim_data.face_velocities,
            exchanged_phi_ghost,
            false, 3, 0, TV() );

    // get so CFL is correct
    sim_data.Set_Boundary_Conditions(time);

    Write_Output_Files(sim_data.first_frame);
}

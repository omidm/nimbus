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

/*
 * Helper functions in water_driver.
 *
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#include <iostream>
#include "assert.h"
#include "shared/nimbus.h"
#include "./water_data_types.h"

using namespace PhysBAM;
using nimbus::Data;

// TODO(someone): right now, create may/ may not allocate required amount of
// memory. Will have to dig deep into PhysBAM to allocate required amount of
// memory at the beginning itself.

template <class TV> FaceArray<TV>::
FaceArray(int size)
{
    id_debug = face_array_id;
    this->size_ = size;
    grid = NULL;
    data = NULL;
}

template <class TV> void FaceArray<TV>::
create()
{
    std::cout << "Creating FaceArray\n";

    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    grid = new GRID<TV> (TV_INT::All_Ones_Vector()*size_,
            RANGE<TV>::Unit_Box(), true);
    assert(grid);

    data = new  ARRAY<T,FACE_INDEX<TV::dimension> >();
    assert(data);
    data->Resize(*grid);
}

template <class TV> Data* FaceArray<TV>::
clone()
{
    std::cout << "Cloning facearray\n";
    return new FaceArray<TV>(size_);
}

    template <class TV>
int FaceArray<TV> :: get_debug_info()
{
    return id_debug;
}

template <class TV> FaceArrayGhost<TV>::
FaceArrayGhost(int size)
{
    id_debug = face_array_ghost_id;
    this->size_ = size;
    grid = NULL;
    data = NULL;
}

template <class TV> void FaceArrayGhost<TV>::
create()
{
    std::cout << "Creating FaceArrayGhost\n";

    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    grid = new GRID<TV> (TV_INT::All_Ones_Vector()*size_,
            RANGE<TV>::Unit_Box(), true);
    assert(grid);

    data = new typename GRID_ARRAYS_POLICY<GRID<TV> >::
        FACE_ARRAYS();
    assert(data);
}

template <class TV> Data* FaceArrayGhost<TV>::
clone()
{
    std::cout << "Cloning facearrayghost\n";
    return new FaceArrayGhost<TV>(size_);
}

    template <class TV>
int FaceArrayGhost<TV> :: get_debug_info()
{
    return id_debug;
}

template <class TV, class T> NonAdvData<TV, T>::
NonAdvData(int size)
{
    id_debug = non_adv_id;

    this->size_ = size;

    number_of_ghost_cells = 3;
    time = (T)0;
    current_frame = 0;

    grid = NULL;

    boundary_scalar = NULL;
    boundary = NULL;
    phi_boundary = NULL;
    phi_boundary_water = NULL;
    domain_boundary = NULL;

    sources = NULL;

    particle_levelset_evolution = NULL;
    advection_scalar = NULL;

    collision_bodies_affecting_fluid = NULL;

    projection = NULL;
    incompressible = NULL;
}

template <class TV, class T> void NonAdvData<TV, T>::
create()
{
    std::cout << "Creating NonAdvData\n";

    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    grid = new GRID<TV> (TV_INT::All_Ones_Vector()*size_,
            RANGE<TV>::Unit_Box(), true);
    assert(grid);

    boundary_scalar = new BOUNDARY_UNIFORM<GRID<TV>, T>();
    phi_boundary_water = new
        typename GEOMETRY_BOUNDARY_POLICY<GRID<TV> >::
        BOUNDARY_PHI_WATER();
    domain_boundary = new VECTOR<VECTOR<bool,2>,TV::dimension>();

    sources = new ARRAY<IMPLICIT_OBJECT<TV> *>();

    particle_levelset_evolution = new
        PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<TV> >
        (*grid, number_of_ghost_cells);
    advection_scalar = new ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T>();

    collision_bodies_affecting_fluid = new
        typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<TV> >::
        GRID_BASED_COLLISION_GEOMETRY(*grid);

    projection = new PROJECTION_DYNAMICS_UNIFORM< GRID<TV> >
        (*grid, false, false, false, false, NULL);
    incompressible = new INCOMPRESSIBLE_UNIFORM<GRID<TV> >(*grid, *projection);
}

template <class TV, class T> Data* NonAdvData<TV, T>::
clone()
{
    std::cout << "Cloning nonadvdata\n";
    return new NonAdvData<TV, T>(size_);
}

    template <class TV, class T>
int NonAdvData<TV, T> :: get_debug_info()
{
    return id_debug;
}

template <class TV, class T> bool NonAdvData<TV, T>::
    initialize
(WaterDriver<TV> *driver, FaceArray<TV> *face_velocities, const int frame)
{
    std::cout << "Initializaing non advection data ...\n";

    typedef typename TV::template REBIND<int>::TYPE TV_INT;

    current_frame = frame;
    time = driver->Time_At_Frame(frame);

    driver->current_frame = current_frame;
    driver->output_number = current_frame;
    driver->time = time;

    for(int i=1;i<=TV::dimension;i++)
    {
        (*domain_boundary)(i)(1)=true;
        (*domain_boundary)(i)(2)=true;
    }
    (*domain_boundary)(2)(2)=false;

    phi_boundary_water->Set_Velocity_Pointer(*face_velocities->data);

    VECTOR<VECTOR<bool,2>,TV::dimension> domain_open_boundaries = 
        VECTOR_UTILITIES::Complement(*domain_boundary);

    phi_boundary = phi_boundary_water;
    phi_boundary->Set_Constant_Extrapolation(domain_open_boundaries);

    boundary = boundary_scalar;
    boundary->Set_Constant_Extrapolation(domain_open_boundaries);

    std::cout << "Moving to incompressible ...\n";

    incompressible->Initialize_Grids(*grid);
    incompressible->Set_Custom_Advection(*advection_scalar);
    incompressible->Set_Custom_Advection(*advection_scalar);
    incompressible->Set_Custom_Boundary(*boundary);
    incompressible->projection.elliptic_solver->Set_Relative_Tolerance(1e-8);
    incompressible->projection.elliptic_solver->pcg.Set_Maximum_Iterations(40);
    incompressible->projection.elliptic_solver->pcg.evolution_solver_type =
        krylov_solver_cg;
    incompressible->projection.elliptic_solver->pcg.cg_restart_iterations=0;
    incompressible->projection.elliptic_solver->pcg.Show_Results();
    incompressible->projection.collidable_solver->Use_External_Level_Set
        (particle_levelset_evolution->particle_levelset.levelset);
    //add forces
    incompressible->Set_Gravity();
    incompressible->Set_Body_Force(true);
    incompressible->projection.Use_Non_Zero_Divergence(false);
    incompressible->projection.elliptic_solver->Solve_Neumann_Regions(true);
    incompressible->projection.elliptic_solver->solve_single_cell_neumann_regions=false;
    incompressible->Use_Explicit_Part_Of_Implicit_Viscosity(false);
    incompressible->Set_Maximum_Implicit_Viscosity_Iterations(40);
    incompressible->Use_Variable_Vorticity_Confinement(false);
    incompressible->Set_Surface_Tension(0);
    incompressible->Set_Variable_Surface_Tension(false);
    incompressible->Set_Viscosity(0);
    incompressible->Set_Variable_Viscosity(false);
    incompressible->projection.Set_Density(1e3);

    std::cout << "Moving to collision bodies affecting fluid ...\n";

    collision_bodies_affecting_fluid->Initialize_Grids();
    collision_bodies_affecting_fluid->Update_Intersection_Acceleration_Structures(false);
    collision_bodies_affecting_fluid->Rasterize_Objects();
    collision_bodies_affecting_fluid->Compute_Occupied_Blocks
        (false, (T)2*grid->Minimum_Edge_Length(), 5);
    collision_bodies_affecting_fluid->Compute_Grid_Visibility();

    std::cout << "Moving to particle levelset evolution ...\n";

    particle_levelset_evolution->Initialize_Domain(*grid);
    particle_levelset_evolution->particle_levelset.Set_Band_Width(6);
    particle_levelset_evolution->Set_Time(time);
    particle_levelset_evolution->Set_CFL_Number((T).9);
    particle_levelset_evolution->Levelset_Advection(1).
        Set_Custom_Advection(*advection_scalar);
    particle_levelset_evolution->Set_Number_Particles_Per_Cell(16);
    particle_levelset_evolution->Set_Levelset_Callbacks(*driver);
    particle_levelset_evolution->Initialize_FMM_Initialization_Iterative_Solver(true);
    particle_levelset_evolution->particle_levelset.levelset.
        Set_Custom_Boundary(*phi_boundary);
    particle_levelset_evolution->Bias_Towards_Negative_Particles(false);
    particle_levelset_evolution->particle_levelset.Use_Removed_Positive_Particles();
    particle_levelset_evolution->particle_levelset.Use_Removed_Negative_Particles();
    particle_levelset_evolution->particle_levelset.Store_Unique_Particle_Id();
    particle_levelset_evolution->Use_Particle_Levelset(true);
    particle_levelset_evolution->particle_levelset.levelset.
        Set_Collision_Body_List(*collision_bodies_affecting_fluid);
    particle_levelset_evolution->particle_levelset.levelset.
        Set_Face_Velocities_Valid_Mask(&incompressible->valid_mask);
    particle_levelset_evolution->particle_levelset.Set_Collision_Distance_Factors(.1,1);

    Initialize_Phi();

    std::cout << "Initialized phi ...\n";

    std::cout << "Just before adjust phi with sources ...\n";
    //Adjust_Phi_With_Sources(time);
    std::cout << "After adjust phi with sources ...\n";

    std::cout << "1\n";
    particle_levelset_evolution->Make_Signed_Distance();
    std::cout << "2\n";
    particle_levelset_evolution->Set_Seed(2606);
    std::cout << "3\n";
    particle_levelset_evolution->Seed_Particles(time);
    std::cout << "4\n";
    particle_levelset_evolution->Delete_Particles_Outside_Grid();
    std::cout << "5\n";

    std::cout << "Extrapolate etc ...\n";

    ARRAY<T,TV_INT> exchanged_phi_ghost(grid->Domain_Indices(8));
    particle_levelset_evolution->particle_levelset.levelset.boundary->
        Fill_Ghost_Cells(*grid, particle_levelset_evolution->phi,
                exchanged_phi_ghost, 0, time, 8);
    incompressible->Extrapolate_Velocity_Across_Interface
        (*face_velocities->data, exchanged_phi_ghost, false, 3, 0, TV());

    projection->Initialize_Grid(*grid);

    std::cout << "Moving to incomplete implementation ...\n";

    Set_Boundary_Conditions(driver, time, face_velocities); // get so CFL is correct

    driver->Write_Output_Files(driver->first_frame);

    std::cout << "Successfully initialized non advection data\n";
    return false;
}

template <class TV, class T> void NonAdvData<TV, T>::
BeforeAdvection
(WaterDriver<TV> *driver, FaceArray<TV> *face_velocities)
{
    LOG::Time("Compute Occupied Blocks");
    T maximum_fluid_speed = face_velocities->data->Maxabs().Max();
    T max_particle_collision_distance = particle_levelset_evolution->
        particle_levelset.max_collision_distance_factor * grid->dX.Max();
    collision_bodies_affecting_fluid->Compute_Occupied_Blocks(true, dt *
            maximum_fluid_speed + 2 * max_particle_collision_distance + (T).5 *
            grid->dX.Max(), 10);

    T_FACE_ARRAYS_SCALAR face_velocities_ghost;
    face_velocities_ghost.Resize(incompressible->grid, number_of_ghost_cells, false);
    incompressible->boundary->Fill_Ghost_Cells_Face(incompressible->grid,
            *face_velocities->data, face_velocities_ghost, time + dt,
            number_of_ghost_cells);

    //Advect Phi 3.6% (Parallelized)
    LOG::Time("Advect Phi");
    phi_boundary_water->Use_Extrapolation_Mode(false);
    particle_levelset_evolution->Advance_Levelset(dt);
    phi_boundary_water->Use_Extrapolation_Mode(true);

    //Advect Particles 12.1% (Parallelized)
    LOG::Time("Step Particles");
    particle_levelset_evolution->particle_levelset.Euler_Step_Particles
        (face_velocities_ghost, dt, time, true, true, false, false);

    //Advect removed particles (Parallelized)
    LOG::Time("Advect Removed Particles");
    RANGE<TV_INT> domain(grid->Domain_Indices());
    domain.max_corner += TV_INT::All_Ones_Vector();

    incompressible->boundary->Fill_Ghost_Cells_Face(*grid,
            *face_velocities->data, face_velocities_ghost, time + dt,
            number_of_ghost_cells);
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,TV> interpolation;
    PARTICLE_LEVELSET_UNIFORM<GRID<TV> > &pls =
        particle_levelset_evolution->particle_levelset;
    if (pls.use_removed_positive_particles)
        for(typename GRID<TV>::NODE_ITERATOR iterator(*grid, domain);
                iterator.Valid();iterator.Next())
            if (pls.removed_positive_particles(iterator.Node_Index()))
            {
                PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> &particles = 
                    *pls.removed_positive_particles(iterator.Node_Index());
                for (int p=1; p<=particles.array_collection->Size(); p++)
                {
                    TV X = particles.X(p),
                       V = interpolation.Clamped_To_Array_Face
                           (*grid, face_velocities_ghost, X);
                    if (-pls.levelset.Phi(X) > 1.5*particles.radius(p))
                        V-=-TV::Axis_Vector(2)*.3; // buoyancy
                    particles.V(p) = V;
                }
            }
    if (pls.use_removed_negative_particles)
        for(typename GRID<TV>::NODE_ITERATOR iterator(*grid,domain);
                iterator.Valid();iterator.Next())
            if(pls.removed_negative_particles(iterator.Node_Index()))
            {
                PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> &particles =
                    *pls.removed_negative_particles(iterator.Node_Index());
                for (int p=1; p <= particles.array_collection->Size(); p++)
                    particles.V(p) += -TV::Axis_Vector(2)*dt*9.8; // ballistic
                for(int p=1; p<=particles.array_collection->Size(); p++)
                    particles.V(p) += dt * interpolation.Clamped_To_Array_Face
                        (*grid, incompressible->force, particles.X(p));
            } // external forces
}

template <class TV> void FaceArray<TV>::
Advection (WaterDriver<TV> *driver,
        NonAdvData<TV, T> *sim_data)
{
    T_FACE_ARRAYS_SCALAR face_velocities_ghost;
    face_velocities_ghost.Resize
        (sim_data->incompressible->grid, sim_data->number_of_ghost_cells, false);
    sim_data->incompressible->boundary->Fill_Ghost_Cells_Face(
            sim_data->incompressible->grid, *data,
            face_velocities_ghost, sim_data->time + sim_data->dt,
            sim_data->number_of_ghost_cells);

    sim_data->incompressible->advection->Update_Advection_Equation_Face(
            *grid, *data, face_velocities_ghost, face_velocities_ghost,
            *sim_data->incompressible->boundary, sim_data->dt, sim_data->time);
}

template <class TV, class T> void NonAdvData<TV, T>::
AfterAdvection
(WaterDriver<TV> *driver, FaceArray<TV> *face_velocities)
{
    T_FACE_ARRAYS_SCALAR face_velocities_ghost;
    face_velocities_ghost.Resize(incompressible->grid, number_of_ghost_cells, false);
    incompressible->boundary->Fill_Ghost_Cells_Face(incompressible->grid,
            *face_velocities->data, face_velocities_ghost, time + dt,
            number_of_ghost_cells);

    //Add Forces 0%
    LOG::Time("Forces");
    incompressible->Advance_One_Time_Step_Forces
        (*face_velocities->data, dt, time, true, 0, number_of_ghost_cells);

    //Modify Levelset with Particles 15% (Parallelizedish)
    LOG::Time("Modify Levelset");
    particle_levelset_evolution->particle_levelset.Exchange_Overlap_Particles();
    particle_levelset_evolution->Modify_Levelset_And_Particles(&face_velocities_ghost);

    //Adjust Phi 0%
    LOG::Time("Adjust Phi");
    Adjust_Phi_With_Sources(time + dt);

    //Delete Particles 12.5 (Parallelized)
    LOG::Time("Delete Particles");
    particle_levelset_evolution->Delete_Particles_Outside_Grid();                                                        //0.1%
    particle_levelset_evolution->particle_levelset.
        Delete_Particles_In_Local_Maximum_Phi_Cells(1);                           //4.9%
    particle_levelset_evolution->particle_levelset.
        Delete_Particles_Far_From_Interface(); // uses visibility                 //7.6%
    particle_levelset_evolution->particle_levelset.
        Identify_And_Remove_Escaped_Particles
        (face_velocities_ghost, 1.5, time + dt); //2.4%

    //Reincorporate Particles 0% (Parallelized)
    LOG::Time("Reincorporate Particles");
    if (particle_levelset_evolution->particle_levelset.
            use_removed_positive_particles ||
            particle_levelset_evolution->particle_levelset.
            use_removed_negative_particles)
        particle_levelset_evolution->particle_levelset.
            Reincorporate_Removed_Particles(1, 1, 0, true);

    //Project 7% (Parallelizedish)
    LOG::SCOPE *scope=0;
    scope=new LOG::SCOPE("Project");
    Set_Boundary_Conditions(driver, time, face_velocities);
    incompressible->Set_Dirichlet_Boundary_Conditions
        (&particle_levelset_evolution->phi, 0);
    projection->p *= dt;
    projection->collidable_solver->Set_Up_Second_Order_Cut_Cell_Method();
    incompressible->Advance_One_Time_Step_Implicit_Part
        (*face_velocities->data, dt, time, true);
    projection->p *= (1/dt);
    incompressible->boundary->Apply_Boundary_Condition_Face
        (incompressible->grid, *face_velocities->data, time + dt);
    projection->collidable_solver->Set_Up_Second_Order_Cut_Cell_Method(false);
    delete scope;

    //Extrapolate Velocity 7%
    LOG::Time("Extrapolate Velocity");
    T_ARRAYS_SCALAR exchanged_phi_ghost(grid->Domain_Indices(8));
    particle_levelset_evolution->particle_levelset.levelset.boundary->
        Fill_Ghost_Cells(*grid, particle_levelset_evolution->phi,
                exchanged_phi_ghost, 0, time + dt, 8);
    incompressible->Extrapolate_Velocity_Across_Interface
        (*face_velocities->data, exchanged_phi_ghost, false, 3, 0, TV());
}

template <class TV, class T> void NonAdvData<TV, T>::
Add_Source()
{
    TV point1, point2;
    BOX<TV> source;
    point1 = TV::All_Ones_Vector()*(T).5;point1(1) = .95;
    point1(2)=.6;point2=TV::All_Ones_Vector()*(T).65;
    point2(1)=1;point2(2)=.75;
    source.min_corner = point1;
    source.max_corner = point2;
    sources->Append(new ANALYTIC_IMPLICIT_OBJECT<BOX<TV> >(source));
}

#ifndef TEMPLATE_USE
#define TEMPLATE_USE
typedef VECTOR<float, 2> TVF2;
typedef float TF;
#endif  // TEMPLATE_USE

template class FaceArray<TVF2>;
template class FaceArrayGhost<TVF2>;
template class NonAdvData<TVF2, TF>;

//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "stdio.h"
#include "string.h"

#include "applications/physbam/smoke/app_utils.h"
#include "applications/physbam/smoke/data_names.h"
#include "applications/physbam/smoke/parameters.h"
#include "applications/physbam/smoke/physbam_include.h"
#include "applications/physbam/smoke/smoke_driver.h"
#include "applications/physbam/smoke/smoke_example.h"
#include "src/shared/dbg.h"
#include "src/shared/geometric_region.h"
#include "src/shared/nimbus.h"

using namespace PhysBAM;

template<class TV> void SMOKE_DRIVER<TV>::InitializeFirstDistributed(
    const nimbus::Job *job,
    const nimbus::DataArray &da)
{
  typedef application::DataConfig DataConfig;
  DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);
  output_number=current_frame;
  time=example.Time_At_Frame(current_frame);

  // allocates arrays for velocity and density
  {
    InitializeProjectionHelper(
       example.data_config,
       example.mac_grid,
       &example.projection);
    
    if (example.data_config.GetFlag(DataConfig::VELOCITY)) {
      LOG::Time("Velocity memory allocated.\n");
      example.face_velocities.Resize(example.mac_grid);
    }
    if (example.data_config.GetFlag(DataConfig::VELOCITY_GHOST)) {
      LOG::Time("Ghost Velocity memory allocated.\n");
      example.face_velocities_ghost.Resize(example.mac_grid,
          example.number_of_ghost_cells, false);
    }
    if (example.data_config.GetFlag(DataConfig::DENSITY)) {
      LOG::Time("Density memory allocated.\n");
      example.density.Resize(example.mac_grid.Domain_Indices(0));
    }
    if (example.data_config.GetFlag(DataConfig::DENSITY_GHOST)) {
      LOG::Time("Ghost Density memory allocated.\n");
      example.density_ghost.Resize(example.mac_grid.Domain_Indices(3));
    }
  }

  {
    // policies, etc.
    example.Initialize_Fields();

    if (example.nimbus_thread_queue) {
      example.projection.elliptic_solver->thread_queue=example.nimbus_thread_queue;
    }

    // setup laplace
    example.projection.elliptic_solver->Set_Relative_Tolerance(1e-9);
    example.projection.elliptic_solver->pcg.Set_Maximum_Iterations(1000);
    example.projection.elliptic_solver->pcg.evolution_solver_type = krylov_solver_cg;
    example.projection.elliptic_solver->pcg.cg_restart_iterations = 40;

    // domain boundaries
    {
      for (int i = 1; i <= TV::dimension; i++) {
	example.domain_boundary(i)(1) = true;
	example.domain_boundary(i)(2) = true;
      }

      if (example.nimbus_thread_queue) {
	example.boundary=new BOUNDARY_THREADED<GRID<TV> >(*example.nimbus_thread_queue,example.boundary_scalar);
      } else {
	example.boundary = &example.boundary_scalar;
      }

      VECTOR<VECTOR<bool, 2>, TV::dimension> constant_extrapolation;
      constant_extrapolation.Fill(VECTOR<bool, 2>::Constant_Vector(true));
      example.boundary->Set_Constant_Extrapolation(constant_extrapolation);
    }


    LAPLACE_UNIFORM<T_GRID>* laplace_solver =
      dynamic_cast<LAPLACE_UNIFORM<T_GRID>* >(
          example.projection.elliptic_solver);
    example.laplace_solver_wrapper.BindLaplaceAndInitialize(laplace_solver);

    example.Set_Boundary_Conditions(time);
  }
  // write, save
  // Write_Output_Files(example.first_frame);
  example.Save_To_Nimbus_No_Cache(job, da, current_frame);
}

template<class TV> void SMOKE_DRIVER<TV>::Initialize(
    const nimbus::Job *job,
    const nimbus::DataArray &da)
{
  typedef application::DataConfig DataConfig;
  DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);
  output_number=current_frame;

  // allocates arrays for velocity and density
  {
    InitializeProjectionHelper(
       example.data_config,
       example.mac_grid,
       &example.projection);

    if (example.data_config.GetFlag(DataConfig::VELOCITY)) {
      LOG::Time("Velocity memory allocated.\n");
      example.face_velocities.Resize(example.mac_grid);
    }
    if (example.data_config.GetFlag(DataConfig::VELOCITY_GHOST)) {
      LOG::Time("Ghost Velocity memory allocated.\n");
      example.face_velocities_ghost.Resize(example.mac_grid,
          example.number_of_ghost_cells, false);
    }
    if (example.data_config.GetFlag(DataConfig::DENSITY)) {
      LOG::Time("Density memory allocated.\n");
      example.density.Resize(example.mac_grid.Domain_Indices(0));
    }
    if (example.data_config.GetFlag(DataConfig::DENSITY_GHOST)) {
      LOG::Time("Ghost Density memory allocated.\n");
      example.density_ghost.Resize(example.mac_grid.Domain_Indices(3));
    }
  }
  {
    
    // policies, etc.
    if (example.nimbus_thread_queue) {
      example.projection.elliptic_solver->thread_queue=example.nimbus_thread_queue;
    }

    // setup laplace
    example.projection.elliptic_solver->Set_Relative_Tolerance(1e-9);
    example.projection.elliptic_solver->pcg.Set_Maximum_Iterations(1000);
    example.projection.elliptic_solver->pcg.evolution_solver_type = krylov_solver_cg;
    example.projection.elliptic_solver->pcg.cg_restart_iterations = 40;

    // load
    example.Load_From_Nimbus(job, da, current_frame);

    //domain boundaries
    {                                                                                                                      
      for (int i = 1; i <= TV::dimension; i++) {
	example.domain_boundary(i)(1) = true;
	example.domain_boundary(i)(2) = true;
      }

      if (example.nimbus_thread_queue) {
        example.boundary=new BOUNDARY_THREADED<GRID<TV> >(*example.nimbus_thread_queue,example.boundary_scalar);
      } else {
        example.boundary = &example.boundary_scalar;
      }

      VECTOR<VECTOR<bool, 2>, TV::dimension> constant_extrapolation;
      constant_extrapolation.Fill(VECTOR<bool, 2>::Constant_Vector(true));
      example.boundary->Set_Constant_Extrapolation(constant_extrapolation);
    }

    LAPLACE_UNIFORM<T_GRID>* laplace_solver =
      dynamic_cast<LAPLACE_UNIFORM<T_GRID>* >(
					      example.projection.elliptic_solver);
    example.laplace_solver_wrapper.BindLaplaceAndInitialize(laplace_solver);
  }
}

template<class TV> void SMOKE_DRIVER<TV>::InitializeUseCache(
    const nimbus::Job *job,
    const nimbus::DataArray &da)
{
  typedef application::DataConfig DataConfig;
  DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);
  output_number=current_frame;

  {
    InitializeProjectionHelper(
        example.data_config,
	example.mac_grid,
	&example.projection);
  }

  {
    // policies, etc.
    if (example.nimbus_thread_queue) {
      example.projection.elliptic_solver->thread_queue=example.nimbus_thread_queue;
    }

    //setup laplace
    example.projection.elliptic_solver->Set_Relative_Tolerance(1e-9);
    example.projection.elliptic_solver->pcg.Set_Maximum_Iterations(1000);
    example.projection.elliptic_solver->pcg.evolution_solver_type = krylov_solver_cg;
    example.projection.elliptic_solver->pcg.cg_restart_iterations = 40;

    // load
    example.Load_From_Nimbus(job, da, current_frame);

    // domain boundaries   
    {
      for (int i = 1; i <= TV::dimension; i++) {
	example.domain_boundary(i)(1) = true;
	example.domain_boundary(i)(2) = true;
      }

      if (example.nimbus_thread_queue) {
        example.boundary=new BOUNDARY_THREADED<GRID<TV> >(*example.nimbus_thread_queue,example.boundary_scalar);
      } else {
        example.boundary = &example.boundary_scalar;
      }

      VECTOR<VECTOR<bool,2>, TV::dimension> constant_extrapolation;
      constant_extrapolation.Fill(VECTOR<bool, 2>::Constant_Vector(true));
      example.boundary->Set_Constant_Extrapolation(constant_extrapolation);
    }

    LAPLACE_UNIFORM<T_GRID>* laplace_solver = 
      dynamic_cast<LAPLACE_UNIFORM<T_GRID>* >(
	  example.projection.elliptic_solver);
    example.laplace_solver_wrapper.BindLaplaceAndInitialize(laplace_solver);
  }
}

template<class TV> bool SMOKE_DRIVER<TV>::InitializeProjectionHelper(
    const application::DataConfig& data_config,
    const GRID<TV>& grid_input,
    PROJECTION_UNIFORM<GRID<TV> >* projection) {
  typedef application::DataConfig DataConfig;
  // assert(grid_input.Is_MAC_GRID());
  projection->p_grid = grid_input;
  // Laplace solver is used.
  assert(projection->poisson == NULL);
  assert(projection->laplace != NULL);
  assert(grid_input.DX() == TV() || grid_input.Is_MAC_Grid());

  // projection->laplace->Initialize_Grid(grid_input);
  {
    LAPLACE_UNIFORM<GRID<TV> >* laplace =
      dynamic_cast<LAPLACE_UNIFORM<GRID<TV> >*>(projection->laplace);
    laplace->grid = grid_input;
    // T_ARRAYS_SCALAR.
    if (data_config.GetFlag(DataConfig::DIVERGENCE)) {
      laplace->f.Resize(grid_input.Domain_Indices(1));
    }
    // T_FACE_ARRAYS_BOOL.
    if (data_config.GetFlag(DataConfig::PSI_N)) {
      laplace->psi_N.Resize(grid_input, 1);
    }
    // T_ARRAYS_BOOL.
    if (data_config.GetFlag(DataConfig::PSI_D)) {
      laplace->psi_D.Resize(grid_input.Domain_Indices(1));
    }
    // T_ARRAYS_INT.
    if (data_config.GetFlag(DataConfig::REGION_COLORS)) {
      laplace->filled_region_colors.Resize(
					   grid_input.Domain_Indices(1));
    }
    // Assume uniform region coloring.
    laplace->number_of_regions = 1;
    laplace->filled_region_touches_dirichlet.Resize(1);
    laplace->filled_region_touches_dirichlet(1) = true;
  }

  //Flag use_non_zero_divergence is expected to be false.
  assert(!projection->use_non_zero_divergence);
  projection->divergence.Clean_Memory();

  //T_ARRAYS_SCALAR.
  if (data_config.GetFlag(DataConfig::PRESSURE)) {
    projection->p.Resize(grid_input.Domain_Indices(1));
  }

  //T_ARRAYS_SCALAR.
  if (data_config.GetFlag(DataConfig::PRESSURE_SAVE)) {
    projection->p_save_for_projection.Resize(grid_input.Domain_Indices(1));
  }

  //T_FACE_ARRAYS_SCALAR.
  if (data_config.GetFlag(DataConfig::VELOCITY_SAVE)) {
    projection->face_velocities_save_for_projection.Resize(grid_input);
  }

  //dsd is not considered.
  // assert(projection->dsd == NULL);
  return true;
}

template class SMOKE_DRIVER<VECTOR<float,3> >;

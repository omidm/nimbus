//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "stdio.h"
#include "string.h"

#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_THREADED.h>

#include "application/smoke/app_utils.h"
#include "application/smoke/data_names.h"
#include "application/smoke/job_names.h"
#include "application/smoke/parameters.h"
#include "application/smoke/projection/laplace_solver_wrapper.h"
#include "application/smoke/projection/projection_helper.h"
#include "application/smoke/smoke_driver.h"
#include "application/smoke/smoke_example.h"
#include "data/physbam/translator_physbam.h"
#include "shared/dbg_modes.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "shared/nimbus.h"

using namespace PhysBAM;
namespace{
    template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
    {
        ((SMOKE_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
    }
};
//#####################################################################
// Initialize
//#####################################################################
template<class TV> SMOKE_DRIVER<TV>::
    SMOKE_DRIVER(SMOKE_EXAMPLE<TV>& example)
:example(example)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> SMOKE_DRIVER<TV>::
~SMOKE_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}

template<class TV> void SMOKE_DRIVER<TV>::
WriteOutputSplitImpl(const nimbus::Job *job,
			  const nimbus::DataArray &da,
			  const bool set_boundary_conditions,
			  const T dt, 
			  const int rank) {
  if (application::kUseGlobalWrite) {
    Write_Output_Files(++output_number);
  } else {
    Write_Output_Files(++output_number, rank);
  }
  // Save State
  // example.Save_To_Nimbus(job, da, current_frame+1);
}


template<class TV> bool SMOKE_DRIVER<TV>::
ScalarAdvanceImpl(const nimbus::Job *job,
		  const nimbus::DataArray &da,
		  const T dt) {
  LOG::Time("Scalar Advance");
  example.Get_Scalar_Field_Sources(time);
  // Array<T, TV_INT> density_ghost(example.mac_grid.Domain_Indices(3));
  example.boundary->Fill_Ghost_Cells(example.mac_grid, example.density, example.density_ghost,
				     dt, time, 3);
  example.advection_scalar.Update_Advection_Equation_Cell(example.mac_grid, 
    example.density, example.density_ghost, example.face_velocities, *example.boundary, 
    dt, time);
  // Save State.
  // example.Save_To_Nimbus(job, da, current_frame + 1);
  return true;
}

template<class TV> bool SMOKE_DRIVER<TV>::
ConvectImpl(const nimbus::Job *job,
	    const nimbus::DataArray &da,
	    const T dt) {
  LOG::Time("Convect");
  // Array<T, FACE_INDEX<TV::dimension> > face_velocities_ghost(example.mac_grid, 3, false);
  example.boundary->Fill_Ghost_Cells_Face(example.mac_grid, example.face_velocities,
					  example.face_velocities_ghost, time, 3);
  example.advection_scalar.Update_Advection_Equation_Face(example.mac_grid, 
    example.face_velocities, example.face_velocities_ghost, example.face_velocities_ghost,
    *example.boundary, dt, time);
  // Save State.
  // example.Save_To_Nimbus(job, da, current_frame + 1);
  return true;
}

template<class TV> bool SMOKE_DRIVER<TV>::
ProjectionCalculateBoundaryConditionPartOneImpl (
						 const nimbus::Job *job,
						 const nimbus::DataArray &da,
						 T dt) {
  // Sets boundary conditions.                                                                                                   
  // Local.                                                                                                                      
  // Read velocity and pressure. Write velocity, pressure, psi_D, and psi_N.                                                     
  example.Set_Boundary_Conditions(time + dt);

  return true;
}

template<class TV> bool SMOKE_DRIVER<TV>::
ProjectionCalculateBoundaryConditionPartTwoImpl (
						 const nimbus::Job *job,
						 const nimbus::DataArray &da,
						 T dt) {
  PROJECTION_UNIFORM<GRID<TV> >& projection = example.projection;

  LAPLACE_UNIFORM<T_GRID>& laplace_solver =
    *dynamic_cast<LAPLACE_UNIFORM<T_GRID>* >(
					    projection.elliptic_solver);

  // Scales pressure.                                                                                                            
  // Read/Write pressure.                                                                                                        
  projection.p *= dt;

  typedef typename INTERPOLATION_POLICY<GRID<TV> >::FACE_LOOKUP
    T_FACE_LOOKUP;

  // Computes divergence.                                                                                                        
  // Local.                                                                                                                      
  // Read velocity. Write divergence(solver->f).                                                                                 
  projection.Compute_Divergence(
				T_FACE_LOOKUP(example.face_velocities),
				&laplace_solver);

  // Coloring.                                                                                                                   
  // Local.                                                                                                                      
  // Read psi_D, psi_N.                                                                                                          
  // Write filled_region_colors.                                                                                                 
  FillUniformRegionColor(
			 laplace_solver.grid,
			 laplace_solver.psi_D, laplace_solver.psi_N,
			 false, &laplace_solver.filled_region_colors);

  return true;
}

// TAG_PROJECTION                                                                                                               
template<class TV> bool SMOKE_DRIVER<TV>::
ProjectionConstructMatrixImpl (
			       const nimbus::Job *job,
			       const nimbus::DataArray &da,
			       T dt) {
  // Read psi_N, psi_D, filled_region_colors, divergence, pressure.                                                             
  // Write A, b, x, tolerance, indexing.                                                                                        
  example.laplace_solver_wrapper.PrepareProjectionInput();
  return true;
}

template<class TV> bool SMOKE_DRIVER<TV>::
ProjectionWrapupImpl(
		     const nimbus::Job *job,
		     const nimbus::DataArray &da,
		     T dt) {
  // Applies pressure.                                                                                                          
  // Local.                                                                                                                     
  // Read pressure(u/p), levelset, psi_D, psi_N, u_interface, velocity.                                                         
  // Write velocity.                                                                                                            
  example.projection.Apply_Pressure(example.face_velocities, dt, time); //Last step in Make_Divergence_Free routine

  // Scales pressure.                                                                                                           
  // Read/Write pressure.                                                                                                       
  example.projection.p *= (1/dt);
  return true;
}

/*
template<class TV> bool SMOKE_DRIVER<TV>::
ProjectImpl(const nimbus::Job *job,
	    const nimbus::DataArray &da,
	    T dt) {
  LOG::Time("Project");
  example.Set_Boundary_Conditions(time+dt);
  example.projection.p* = dt;
  example.boundary->Apply_Boundary_Condition_Face(example.mac_grid, 
						  example.face_velocities, time + dt);
  example.projection.Make_Divergence_Free(example.face_velocities, dt, time);
  example.projection.p* = (1/dt);
  // Save State.
  // example.Save_To_Nimubs(job, da, current_frame + 1);
  return true;
}
*/

//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void SMOKE_DRIVER<TV>::
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
template<class TV> void SMOKE_DRIVER<TV>::
Write_Output_Files(const int frame, int rank)
{
    if (rank != -1) {
      std::string rank_name = "";
      std::stringstream temp_ss;
      temp_ss << "split_output/" << rank;
      rank_name = temp_ss.str();
      FILE_UTILITIES::Create_Directory("split_output/");
      FILE_UTILITIES::Create_Directory(rank_name);
      FILE_UTILITIES::Create_Directory(rank_name+STRING_UTILITIES::string_sprintf("/%d",frame));
      FILE_UTILITIES::Create_Directory(rank_name+"/common");
      FILE_UTILITIES::Write_To_Text_File(rank_name+STRING_UTILITIES::string_sprintf("/%d/frame_title",frame),example.frame_title);
      if(frame==example.first_frame)
        FILE_UTILITIES::Write_To_Text_File(rank_name+"/common/first_frame",frame,"\n");
      example.Write_Output_Files(frame, rank);
      FILE_UTILITIES::Write_To_Text_File(rank_name+"/common/last_frame",frame,"\n");
    } else {
      FILE_UTILITIES::Create_Directory(example.output_directory);
      FILE_UTILITIES::Create_Directory(example.output_directory+STRING_UTILITIES::string_sprintf("/%d",frame));
      FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
      FILE_UTILITIES::Write_To_Text_File(example.output_directory+STRING_UTILITIES::string_sprintf("/%d/frame_title",frame),example.frame_title);
      if(frame==example.first_frame)
        FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/first_frame",frame,"\n");
      example.Write_Output_Files(frame);
      FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/last_frame",frame,"\n");
    }
}
//#####################################################################
template class SMOKE_DRIVER<VECTOR<float,3> >;

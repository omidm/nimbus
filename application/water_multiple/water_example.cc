//#####################################################################
// Copyright 2009, Michael Lentine, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <set>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_1D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_3D.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Uniform/READ_WRITE_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/FLUID_GRAVITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/INCOMPRESSIBILITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM.h>
#include <PhysBAM_Dynamics/Geometry/GENERAL_GEOMETRY_FORWARD.h>

#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/cache_prototypes.h"
#include "application/water_multiple/cache_options.h"
#include "application/water_multiple/data_include.h"
#include "application/water_multiple/options.h"
#include "application/water_multiple/parameters.h"
#include "application/water_multiple/reg_def.h"
#include "application/water_multiple/water_example.h"
#include "data/physbam/translator_physbam_old.h"
#include "data/scalar_data.h"
#include "shared/nimbus.h"
#include "worker/physical_data_instance.h"

// TODO(quhang) In three places where nimbus_thread_queue is introduced.

using namespace PhysBAM;
//#####################################################################
// WATER_EXAMPLE
//#####################################################################
template<class TV> WATER_EXAMPLE<TV>::
WATER_EXAMPLE(const STREAM_TYPE stream_type_input,
              bool use_threading,
              int thread_quota) :
    nimbus_thread_queue(use_threading ?
                        new nimbus::NimbusThreadQueue(thread_quota) : NULL),
    stream_type(stream_type_input),
    initial_time(0),
    first_frame(0),
    last_frame(application::kLastFrame),
    frame_rate(24),
    write_substeps_level(-1),
    write_output_files(true),
    output_directory(application::kOutputDir),
    number_of_ghost_cells(application::kGhostNum),
    cfl(1),
    mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),//incompressible_fluid_collection(mac_grid),
    projection(*new PROJECTION_DYNAMICS_UNIFORM<GRID<TV> >(mac_grid,false,false,false,false,nimbus_thread_queue)),
    particle_levelset_evolution(*new  PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<TV> >(mac_grid,number_of_ghost_cells)),
    incompressible(mac_grid,projection,nimbus_thread_queue),
    boundary(0),
    collision_bodies_affecting_fluid(mac_grid)
{
    use_cache   = false;
    cache_fv    = NULL;
    cache_fvg   = NULL;
    cache_psi_n = NULL;
    cache_phi3  = NULL;
    cache_phi7  = NULL;
    cache_phi8  = NULL;
    cache_psi_d = NULL;
    cache_ple   = NULL;
    cache_pressure = NULL;
    cache_colors = NULL;
    cache_divergence = NULL;
    cache_matrix_a = NULL;
    create_destroy_ple = true;
    Initialize_Particles();
    Initialize_Read_Write_General_Structures();
}
//#####################################################################
// WATER_EXAMPLE
//#####################################################################
template<class TV> WATER_EXAMPLE<TV>::
WATER_EXAMPLE(const STREAM_TYPE stream_type_input,
              application::AppCacheObjects *cache,
              bool use_threading,
              int thread_quota) :
    nimbus_thread_queue(use_threading ?
                        new nimbus::NimbusThreadQueue(thread_quota) : NULL),
    stream_type(stream_type_input),
    initial_time(0),
    first_frame(0),
    last_frame(application::kLastFrame),
    frame_rate(24),
    write_substeps_level(-1),
    write_output_files(true),
    output_directory(application::kOutputDir),
    number_of_ghost_cells(application::kGhostNum),
    cfl(1),
    mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),//incompressible_fluid_collection(mac_grid),
    projection(*new PROJECTION_DYNAMICS_UNIFORM<GRID<TV> >(mac_grid,false,false,false,false,nimbus_thread_queue)),
    particle_levelset_evolution(*new  PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<TV> >(mac_grid,number_of_ghost_cells)),
    incompressible(mac_grid,projection,nimbus_thread_queue),
    boundary(0),
    collision_bodies_affecting_fluid(mac_grid)
{
    use_cache   = false;
    cache_fv    = cache->fv;
    cache_fvg   = cache->fvg;
    cache_psi_n = cache->psi_n;
    cache_phi3  = cache->phi3;
    cache_phi7  = cache->phi7;
    cache_phi8  = cache->phi8;
    cache_psi_d = cache->psi_d;
    cache_ple   = cache->ple;
    cache_pressure = cache->pressure;
    cache_colors = cache->color;
    cache_divergence = cache->divergence;
    cache_matrix_a = cache->matrix_a;
    create_destroy_ple = true;
    Initialize_Particles();
    Initialize_Read_Write_General_Structures();
}
//#####################################################################
// WATER_EXAMPLE
//#####################################################################
template<class TV> WATER_EXAMPLE<TV>::
WATER_EXAMPLE(const STREAM_TYPE stream_type_input,
              application::AppCacheObjects *cache,
              PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<TV> > *ple,
              bool use_threading,
              int thread_quota) :
    nimbus_thread_queue(use_threading ?
                        new nimbus::NimbusThreadQueue(thread_quota) : NULL),
    stream_type(stream_type_input),
    initial_time(0),
    first_frame(0),
    last_frame(application::kLastFrame),
    frame_rate(24),
    write_substeps_level(-1),
    write_output_files(true),
    output_directory(application::kOutputDir),
    number_of_ghost_cells(application::kGhostNum),
    cfl(1),
    mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),//incompressible_fluid_collection(mac_grid),
    projection(*new PROJECTION_DYNAMICS_UNIFORM<GRID<TV> >(mac_grid,false,false,false,false,nimbus_thread_queue)),
    particle_levelset_evolution(*ple),
    incompressible(mac_grid,projection,nimbus_thread_queue),
    boundary(0),
    collision_bodies_affecting_fluid(mac_grid)
{
    use_cache   = false;
    cache_fv    = cache->fv;
    cache_fvg   = cache->fvg;
    cache_psi_n = cache->psi_n;
    cache_phi3  = cache->phi3;
    cache_phi7  = cache->phi7;
    cache_phi8  = cache->phi8;
    cache_psi_d = cache->psi_d;
    cache_ple   = cache->ple;
    cache_pressure = cache->pressure;
    cache_colors = cache->color;
    cache_divergence = cache->divergence;
    cache_matrix_a = cache->matrix_a;
    create_destroy_ple = false;
    Initialize_Particles();
    Initialize_Read_Write_General_Structures();
}
//#####################################################################
// ~WATER_EXAMPLE
//#####################################################################
template<class TV> WATER_EXAMPLE<TV>::
~WATER_EXAMPLE()
{
    delete &projection;
    if (nimbus_thread_queue) {
      delete nimbus_thread_queue;
    }
}

// Initializes the initial levelset function.
template<class TV> void WATER_EXAMPLE<TV>::
Initialize_Phi() {
  ARRAY<T,TV_INT>& phi = particle_levelset_evolution.phi;
  for (typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);
       iterator.Valid();
       iterator.Next()) {
    const TV& X = iterator.Location();
    phi(iterator.Cell_Index()) = X.y - (T)0.35;
  }
}
//#####################################################################
// Initialize_Phi
//#####################################################################
// NOT REQUIRED IN NIMBUS
template<class TV> void WATER_EXAMPLE<TV>::
Initialize_Grid(TV_INT counts,RANGE<TV> domain)
{
    mac_grid.Initialize(counts,domain,true);
}

// Sets the boundary conditions before projection. It might read levelset and
// velocity near the boundary. For velocity, it only reads to check if an index
// is valid. It will try to write to psi_D and psi_N of the whole region, and
// write to pressure and velocity in the boundary region.
// From what I read. --quhang
template<class TV> void WATER_EXAMPLE<TV>::
Set_Boundary_Conditions(const T time) {
  projection.elliptic_solver->psi_D.Fill(false);
  projection.elliptic_solver->psi_N.Fill(false);
  for (int axis = 1; axis <= TV::dimension; axis++) {
    for (int axis_side = 1; axis_side <= 2; axis_side++) {
      int side = 2 * (axis-1) + axis_side;
      TV_INT interior_cell_offset =
          axis_side==1 ? TV_INT() : -TV_INT::Axis_Vector(axis);
      TV_INT exterior_cell_offset =
          axis_side==1 ? -TV_INT::Axis_Vector(axis) : TV_INT();
      TV_INT boundary_face_offset =
          axis_side==1 ? TV_INT::Axis_Vector(axis) : -TV_INT::Axis_Vector(axis);
      if (domain_boundary(axis)(axis_side)) {
        for (typename GRID<TV>::FACE_ITERATOR iterator(
                mac_grid, 1, GRID<TV>::BOUNDARY_REGION, side);
            iterator.Valid();
            iterator.Next()) {
          TV_INT face = iterator.Face_Index() + boundary_face_offset;
          if (particle_levelset_evolution.phi(face + interior_cell_offset)
              <= 0) {
            if (face_velocities.Component(axis).Valid_Index(face)){
              projection.elliptic_solver->psi_N.Component(axis)(face) = true;
              face_velocities.Component(axis)(face) = 0;
            }
          } else {
            TV_INT cell = face + exterior_cell_offset;
            projection.elliptic_solver->psi_D(cell) = true;
            projection.p(cell) = 0;
          }
        }
      } else {
        for (typename GRID<TV>::FACE_ITERATOR iterator(
                mac_grid, 1, GRID<TV>::BOUNDARY_REGION, side);
            iterator.Valid();
            iterator.Next()) {
          TV_INT cell = iterator.Face_Index() + interior_cell_offset;
          projection.elliptic_solver->psi_D(cell) = true;
          projection.p(cell) = 0;
        }
      }
    }
  }
  for (typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);
       iterator.Valid();
       iterator.Next()) {
    for (int i = 1; i <= sources.m; i++) {
      if (time <= 3 && sources(i)->Lazy_Inside(iterator.Location())) {
        projection.elliptic_solver->psi_N(iterator.Full_Index()) = true;
        if ((TV::dimension==2 && iterator.Axis()==1) ||
            (TV::dimension==3 && iterator.Axis()==3)) {
          face_velocities(iterator.Full_Index()) = -1;
        } else {
          face_velocities(iterator.Full_Index()) = 0;
        }
      }
    }
  }
}

// Enforces the boundary condition of levelset.
template<class TV> void WATER_EXAMPLE<TV>::
Adjust_Phi_With_Sources(const T time)
{
  if (time > 3) return;
  for (typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);
       iterator.Valid();
       iterator.Next()) {
    TV_INT index = iterator.Cell_Index();
    for (int i = 1; i <= sources.m; i++)
      particle_levelset_evolution.phi(index) = min(
          particle_levelset_evolution.phi(index),
          sources(i)->Extended_Phi(iterator.Location()));
  }
}

// Enforces the boundary condition of particles.
template<class TV> void WATER_EXAMPLE<TV>::
Adjust_Particle_For_Domain_Boundaries(
    PARTICLE_LEVELSET_PARTICLES<TV>& particles,
    const int index,
    TV& V,
    const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,
    const T dt,
    const T time) {
  if (particle_type == PARTICLE_LEVELSET_POSITIVE ||
      particle_type == PARTICLE_LEVELSET_REMOVED_POSITIVE)
    return;

  TV& X = particles.X(index);
  TV X_new = X + dt*V;
  T max_collision_distance =
      particle_levelset_evolution.particle_levelset.Particle_Collision_Distance(
          particles.quantized_collision_distance(index));
  T min_collision_distance =
      particle_levelset_evolution.particle_levelset.
          min_collision_distance_factor *
      max_collision_distance;
  TV min_corner = mac_grid.domain.Minimum_Corner();
  TV max_corner = mac_grid.domain.Maximum_Corner();
  for (int axis = 1; axis <= GRID<TV>::dimension; axis++) {
    if (domain_boundary[axis][1] &&
        X_new[axis] < min_corner[axis] + max_collision_distance) {
      T collision_distance = X[axis] - min_corner[axis];
      if (collision_distance > max_collision_distance)
        collision_distance = X_new[axis] - min_corner[axis];
      collision_distance = max(min_collision_distance, collision_distance);
      X_new[axis] += max((T)0, min_corner[axis]-X_new[axis]+collision_distance);
      V[axis] = max((T)0, V[axis]);
      X = X_new - dt*V;
    }
    if (domain_boundary[axis][2] &&
        X_new[axis] > max_corner[axis] - max_collision_distance) {
      T collision_distance = max_corner[axis] - X[axis];
      if (collision_distance > max_collision_distance)
        collision_distance = max_corner[axis] - X_new[axis];
      collision_distance = max(min_collision_distance, collision_distance);
      X_new[axis] -= max((T)0, X_new[axis]-max_corner[axis]+collision_distance);
      V[axis] = min((T)0,V[axis]);
      X = X_new-dt*V;
    }
  }
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Write_Output_Files(const int frame, int rank)
{
    if(!write_output_files) return;
    if (rank != -1) {
      std::string rank_name = "";
      std::stringstream temp_ss;
      temp_ss << "split_output/" << rank;
      rank_name = temp_ss.str();
      std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
      FILE_UTILITIES::Write_To_File(stream_type,rank_name+"/"+f+"/mac_velocities",face_velocities);
      if (rank != -1) {
        FILE_UTILITIES::Write_To_File(stream_type,rank_name+"/common/global_grid",
            GRID<TV>(kScale, kScale, kScale, 0, 1, 0, 1, 0, 1, true));
      }
      FILE_UTILITIES::Write_To_File(stream_type,rank_name+"/common/grid",mac_grid);
      FILE_UTILITIES::Write_To_File(stream_type,rank_name+"/"+f+"/pressure",incompressible.projection.p);
      FILE_UTILITIES::Write_To_File(stream_type,rank_name+"/"+f+"/psi_N",projection.elliptic_solver->psi_N);
      FILE_UTILITIES::Write_To_File(stream_type,rank_name+"/"+f+"/psi_D",projection.elliptic_solver->psi_D);
      T_PARTICLE_LEVELSET& particle_levelset=particle_levelset_evolution.particle_levelset;
      FILE_UTILITIES::Write_To_File(stream_type,rank_name+"/"+f+"/levelset",particle_levelset.levelset);
      FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf(
              "%s/%d/%s",rank_name.c_str(),frame,"positive_particles"),particle_levelset.positive_particles);
      FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf(
              "%s/%d/%s",rank_name.c_str(),frame,"negative_particles"),particle_levelset.negative_particles);
      FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf(
              "%s/%d/%s",rank_name.c_str(),frame,"removed_positive_particles"),particle_levelset.removed_positive_particles);
      FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf(
              "%s/%d/%s",rank_name.c_str(),frame,"removed_negative_particles"),particle_levelset.removed_negative_particles);
      FILE_UTILITIES::Write_To_Text_File(rank_name+"/"+f+"/last_unique_particle_id",particle_levelset.last_unique_particle_id);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
#else
      PHYSBAM_FATAL_ERROR("Cannot read doubles");
#endif
    } else {
      std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
      FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities);
      FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",mac_grid);
      FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/pressure",incompressible.projection.p);
      FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_N",projection.elliptic_solver->psi_N);
      FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_D",projection.elliptic_solver->psi_D);
      T_PARTICLE_LEVELSET& particle_levelset=particle_levelset_evolution.particle_levelset;
      FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/levelset",particle_levelset.levelset);
      FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf(
              "%s/%d/%s",output_directory.c_str(),frame,"positive_particles"),particle_levelset.positive_particles);
      FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf(
              "%s/%d/%s",output_directory.c_str(),frame,"negative_particles"),particle_levelset.negative_particles);
      FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf(
              "%s/%d/%s",output_directory.c_str(),frame,"removed_positive_particles"),particle_levelset.removed_positive_particles);
      FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf(
              "%s/%d/%s",output_directory.c_str(),frame,"removed_negative_particles"),particle_levelset.removed_negative_particles);
      FILE_UTILITIES::Write_To_Text_File(output_directory+"/"+f+"/last_unique_particle_id",particle_levelset.last_unique_particle_id);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
#else
      PHYSBAM_FATAL_ERROR("Cannot read doubles");
#endif
    }
}
//#####################################################################
// Read_Output_Files
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    T_PARTICLE_LEVELSET& particle_levelset=particle_levelset_evolution.particle_levelset;
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/"+f+"/levelset",particle_levelset.levelset);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"positive_particles"),particle_levelset.positive_particles);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"negative_particles"),particle_levelset.negative_particles);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_positive_particles"),particle_levelset.removed_positive_particles);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_negative_particles"),particle_levelset.removed_negative_particles);
    FILE_UTILITIES::Read_From_Text_File(output_directory+"/"+f+"/last_unique_particle_id",particle_levelset.last_unique_particle_id);
    std::string filename;
    filename=output_directory+"/"+f+"/pressure";
    if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading pressure "<<filename<<std::endl;LOG::filecout(ss.str());FILE_UTILITIES::Read_From_File(stream_type,filename,incompressible.projection.p);}
    filename=output_directory+"/"+f+"/mac_velocities";
    if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading mac_velocities "<<filename<<std::endl;LOG::filecout(ss.str());FILE_UTILITIES::Read_From_File(stream_type,filename,face_velocities);}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
#else
        PHYSBAM_FATAL_ERROR("Cannot read doubles");
#endif
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Save_To_Nimbus_No_Cache(const nimbus::Job *job, const nimbus::DataArray &da, const int frame)
{
    nimbus::int_dimension_t array_shift[3] = {
        local_region.x() - 1, local_region.y() - 1, local_region.z() - 1};
    nimbus::PdiVector pdv;

    GeometricRegion array_reg_central(local_region);
    GeometricRegion array_reg_outer(array_reg_central.NewEnlarged(application::kGhostNum));
    GeometricRegion array_reg_thin_outer(array_reg_central.NewEnlarged(1));
    GeometricRegion array_reg_outer_7(array_reg_central.NewEnlarged(7));
    GeometricRegion array_reg_outer_8(array_reg_central.NewEnlarged(8));

    GeometricRegion enlarge(1-application::kGhostNum,
                            1-application::kGhostNum,
                            1-application::kGhostNum,
                            local_region.dx()+2*application::kGhostNum,
                            local_region.dy()+2*application::kGhostNum,
                            local_region.dz()+2*application::kGhostNum);

    // mac velocities ghost
    const std::string fvstring = std::string(APP_FACE_VEL);
    if (application::GetTranslatorData(job, fvstring, da, &pdv, application::WRITE_ACCESS)
        && data_config.GetFlag(DataConfig::VELOCITY)) {
      translator.WriteFaceArrayFloat(
          &array_reg_central, array_shift, &pdv, &face_velocities);
    }
    application::DestroyTranslatorObjects(&pdv);

    // mac velocities ghost
    const std::string fvgstring = std::string(APP_FACE_VEL_GHOST);
    if (application::GetTranslatorData(job, fvgstring, da, &pdv, application::WRITE_ACCESS)
        && data_config.GetFlag(DataConfig::VELOCITY_GHOST)) {
      translator.WriteFaceArrayFloat(
          &array_reg_outer, array_shift, &pdv, &face_velocities_ghost);
    }
    application::DestroyTranslatorObjects(&pdv);

    // particle leveset quantities
    T_PARTICLE_LEVELSET& particle_levelset = particle_levelset_evolution.particle_levelset;

    // levelset
    const std::string lsstring = std::string(APP_PHI);
    if (application::GetTranslatorData(job, lsstring, da, &pdv, application::WRITE_ACCESS)
        && data_config.GetFlag(DataConfig::LEVELSET)) {
      translator.WriteScalarArrayFloat(
          &array_reg_outer,
          array_shift,
          &pdv,
          &particle_levelset.levelset.phi);
      std::cout << "OMID: write 3.\n";
    }
    application::DestroyTranslatorObjects(&pdv);
    if (application::GetTranslatorData(job, lsstring, da, &pdv, application::WRITE_ACCESS)
        && data_config.GetFlag(DataConfig::LEVELSET_WRITE)) {
      translator.WriteScalarArrayFloat(
          &array_reg_outer,
          array_shift,
          &pdv,
          &particle_levelset.levelset.phi);
      std::cout << "OMID: write 3.\n";
    }
    application::DestroyTranslatorObjects(&pdv);
    if (application::GetTranslatorData(job, lsstring, da, &pdv, application::WRITE_ACCESS)
        && data_config.GetFlag(DataConfig::LEVELSET_BW_SEVEN_WRITE)) {
      translator.WriteScalarArrayFloat(
          &array_reg_outer_7,
          array_shift,
          &pdv,
          &phi_ghost_bandwidth_seven);
      std::cout << "OMID: write 7.\n";
    }
    application::DestroyTranslatorObjects(&pdv);
    if (application::GetTranslatorData(job, lsstring, da, &pdv, application::WRITE_ACCESS)
        && data_config.GetFlag(DataConfig::LEVELSET_BW_EIGHT_WRITE)) {
      translator.WriteScalarArrayFloat(
          &array_reg_outer_8,
          array_shift,
          &pdv,
          &phi_ghost_bandwidth_eight);
      std::cout << "OMID: write 8.\n";
    }
    application::DestroyTranslatorObjects(&pdv);

    // positive particles
    const std::string ppstring = std::string(APP_POS_PARTICLES);
    if (application::GetTranslatorData(job, ppstring, da, &pdv, application::WRITE_ACCESS)
        && data_config.GetFlag(DataConfig::POSITIVE_PARTICLE)) {
      translator.WriteParticles(
          &enlarge, array_shift,
          &pdv, particle_levelset, kScale, true);
    }
    application::DestroyTranslatorObjects(&pdv);

    // negative particles
    const std::string npstring = std::string(APP_NEG_PARTICLES);
    if (application::GetTranslatorData(job, npstring, da, &pdv, application::WRITE_ACCESS)
        && data_config.GetFlag(DataConfig::NEGATIVE_PARTICLE)) {
      translator.WriteParticles(
          &enlarge, array_shift,
          &pdv, particle_levelset, kScale, false);
    }
    application::DestroyTranslatorObjects(&pdv);

    // Removed positive particles.
    const std::string prpstring = std::string(APP_POS_REM_PARTICLES);
    if (application::GetTranslatorData(job, prpstring, da, &pdv, application::WRITE_ACCESS)
        && data_config.GetFlag(DataConfig::REMOVED_POSITIVE_PARTICLE)) {
      translator.WriteRemovedParticles(
          &enlarge, array_shift,
          &pdv, particle_levelset, kScale, true);
    }
    application::DestroyTranslatorObjects(&pdv);

    // Removed negative particles.
    const std::string nrpstring = std::string(APP_NEG_REM_PARTICLES);
    if (application::GetTranslatorData(job, nrpstring, da, &pdv, application::WRITE_ACCESS)
        && data_config.GetFlag(DataConfig::REMOVED_POSITIVE_PARTICLE)) {
      translator.WriteRemovedParticles(
          &enlarge, array_shift,
          &pdv, particle_levelset, kScale, false);
    }
    application::DestroyTranslatorObjects(&pdv);

    // last unique particle id
    const std::string lupistring = std::string(APP_LAST_UNIQUE_PARTICLE_ID);
    if (Data *d = application::GetFirstData(lupistring, da)) {
        nimbus::ScalarData<int> *sd = static_cast<nimbus::ScalarData<int> * >(d);
        sd->set_scalar(particle_levelset.last_unique_particle_id);
    }

    Data* data_temp = application::GetTheOnlyData(
        job, std::string(APP_DT), da, application::WRITE_ACCESS);
    if (data_temp) {
      nimbus::ScalarData<float>* data_real =
          dynamic_cast<nimbus::ScalarData<float>*>(data_temp);
      data_real->set_scalar(dt_buffer);
      dbg(APP_LOG, "[Data Saving]%s: %0.9f\n", APP_DT, dt_buffer);
    }

    // psi_d.
    const std::string psi_d_string = std::string(APP_PSI_D);
    if (application::GetTranslatorData(job, psi_d_string, da, &pdv, application::WRITE_ACCESS)
        && data_config.GetFlag(DataConfig::PSI_D)) {
      translator.WriteScalarArrayBool(
          &array_reg_thin_outer, array_shift, &pdv, &projection.laplace->psi_D);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating psi_d.\n");
    // psi_n.
    const std::string psi_n_string = std::string(APP_PSI_N);
    if (application::GetTranslatorData(job, psi_n_string, da, &pdv, application::WRITE_ACCESS)
        && data_config.GetFlag(DataConfig::PSI_N)) {
      translator.WriteFaceArrayBool(
          &array_reg_thin_outer, array_shift, &pdv, &projection.laplace->psi_N);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating psi_n.\n");
    // pressure.
    const std::string pressure_string = std::string(APP_PRESSURE);
    if (application::GetTranslatorData(job, pressure_string, da, &pdv, application::WRITE_ACCESS)
        && data_config.GetFlag(DataConfig::PRESSURE)) {
      translator.WriteScalarArrayFloat(
          &array_reg_thin_outer, array_shift, &pdv, &projection.p);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating pressure.\n");
    // filled_region_colors.
    const std::string filled_region_colors_string =
        std::string(APP_FILLED_REGION_COLORS);
    if (application::GetTranslatorData(job, filled_region_colors_string, da, &pdv, application::WRITE_ACCESS)
        && data_config.GetFlag(DataConfig::REGION_COLORS)) {
      dbg(APP_LOG, "filled_region_colors is being written to Nimbus.\n");
      translator.WriteScalarArrayInt(
          &array_reg_thin_outer, array_shift, &pdv,
          &projection.laplace->filled_region_colors);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating filled_region_colors.\n");
    // divergence.
    const std::string divergence_string =
        std::string(APP_DIVERGENCE);
    if (application::GetTranslatorData(job, divergence_string, da, &pdv, application::WRITE_ACCESS)
        && data_config.GetFlag(DataConfig::DIVERGENCE)) {
      translator.WriteScalarArrayFloat(
          &array_reg_thin_outer, array_shift, &pdv,
          &projection.laplace->f);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating divergence.\n");

    typedef nimbus::Data Data;
    if (data_config.GetFlag(DataConfig::MATRIX_A)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_MATRIX_A), da, application::WRITE_ACCESS);
      if (data_temp) {
        application::DataSparseMatrix* data_real =
            dynamic_cast<application::DataSparseMatrix*>(data_temp);
        data_real->SaveToNimbus(laplace_solver_wrapper.A_array(1));
        dbg(APP_LOG, "Finish writing MATRIX_A.\n");
      }
    }
    if (data_config.GetFlag(DataConfig::VECTOR_B)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_VECTOR_B), da, application::WRITE_ACCESS);
      if (data_temp) {
        application::DataRawVectorNd* data_real =
            dynamic_cast<application::DataRawVectorNd*>(data_temp);
        data_real->SaveToNimbus(laplace_solver_wrapper.b_array(1));
        dbg(APP_LOG, "Finish writing VECTOR_B.\n");
      }
    }
    if (data_config.GetFlag(DataConfig::INDEX_C2M)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_INDEX_C2M), da, application::WRITE_ACCESS);
      if (data_temp) {
        application::DataRawGridArray* data_real =
            dynamic_cast<application::DataRawGridArray*>(data_temp);
        data_real->SaveToNimbus(
            laplace_solver_wrapper.cell_index_to_matrix_index);
        dbg(APP_LOG, "Finish writing INDEX_C2M.\n");
      }
    }
    if (data_config.GetFlag(DataConfig::INDEX_M2C)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_INDEX_M2C), da, application::WRITE_ACCESS);
      if (data_temp) {
        application::DataRawArrayM2C* data_real =
            dynamic_cast<application::DataRawArrayM2C*>(data_temp);
        data_real->SaveToNimbus(
            laplace_solver_wrapper.matrix_index_to_cell_index_array(1));
        dbg(APP_LOG, "Finish writing INDEX_M2C.\n");
      }
    }
    if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_TOLERANCE)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_PROJECTION_LOCAL_TOLERANCE),
          da, application::WRITE_ACCESS);
      if (data_temp) {
        nimbus::ScalarData<float>* data_real =
            dynamic_cast<nimbus::ScalarData<float>*>(data_temp);
        data_real->set_scalar(projection.elliptic_solver->tolerance);
        dbg(APP_LOG, "Finish writing PROJECTION_LOCAL_TOLERANCE.\n");
      }
    }
    if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_N)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_PROJECTION_LOCAL_N),
          da, application::WRITE_ACCESS);
      if (data_temp) {
        nimbus::ScalarData<int>* data_real =
            dynamic_cast<nimbus::ScalarData<int>*>(data_temp);
        data_real->set_scalar(laplace_solver_wrapper.local_n);
        dbg(APP_LOG, "Finish writing PROJECTION_LOCAL_N.\n");
      }
    }
    if (data_config.GetFlag(DataConfig::PROJECTION_INTERIOR_N)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_PROJECTION_INTERIOR_N),
          da, application::WRITE_ACCESS);
      if (data_temp) {
        nimbus::ScalarData<int>* data_real =
            dynamic_cast<nimbus::ScalarData<int>*>(data_temp);
        data_real->set_scalar(laplace_solver_wrapper.interior_n);
        dbg(APP_LOG, "Finish writing PROJECTION_INTERIOR_N.\n");
      }
    }
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Save_To_Nimbus(const nimbus::Job *job, const nimbus::DataArray &da, const int frame)
{
    if (!(use_cache && application::kUseCache)) {
      Save_To_Nimbus_No_Cache(job, da, frame);
      return;
    }

    // nimbus::int_dimension_t array_shift[3] = {
    //     local_region.x() - 1, local_region.y() - 1, local_region.z() - 1};
    nimbus::PdiVector pdv;

    GeometricRegion array_reg_central(local_region);
    GeometricRegion array_reg_outer(array_reg_central.NewEnlarged(application::kGhostNum));
    GeometricRegion array_reg_thin_outer(array_reg_central.NewEnlarged(1));
    GeometricRegion array_reg_outer_7(array_reg_central.NewEnlarged(7));
    GeometricRegion array_reg_outer_8(array_reg_central.NewEnlarged(8));

    GeometricRegion enlarge(1-application::kGhostNum,
                            1-application::kGhostNum,
                            1-application::kGhostNum,
                            local_region.dx()+2*application::kGhostNum,
                            local_region.dy()+2*application::kGhostNum,
                            local_region.dz()+2*application::kGhostNum);

    nimbus::CacheManager *cm = job->GetCacheManager();
    // mac velocities
    if (cache_fv) {
        T_FACE_ARRAY *fv = cache_fv->data();
        T_FACE_ARRAY::Exchange_Arrays(*fv, face_velocities);
        nimbus::DataArray write;
        cm->ReleaseAccess(cache_fv);
        cache_fv = NULL;
    }

    // mac velocities ghost
    if (cache_fvg) {
        T_FACE_ARRAY *fvg = cache_fvg->data();
        T_FACE_ARRAY::Exchange_Arrays(*fvg, face_velocities_ghost);
        nimbus::DataArray write;
        cm->ReleaseAccess(cache_fvg);
        cache_fvg = NULL;
    }

    {
      // particle leveset quantities
      T_PARTICLE_LEVELSET& particle_levelset = particle_levelset_evolution.particle_levelset;
        nimbus::DataArray write;
        application::GetWriteData(*job, APP_PHI, da, &write);
      // levelset
      if (cache_phi3) {
          T_SCALAR_ARRAY *phi3 = cache_phi3->data();
          T_SCALAR_ARRAY::Exchange_Arrays(*phi3, particle_levelset.levelset.phi);
          cm->ReleaseAccess(cache_phi3);
          cache_phi3 = NULL;
      }
      if (cache_phi7) {
          T_SCALAR_ARRAY *phi7 = cache_phi7->data();
          T_SCALAR_ARRAY::Exchange_Arrays(*phi7, phi_ghost_bandwidth_seven);
          cm->ReleaseAccess(cache_phi7);
          cache_phi7 = NULL;
      }
      if (cache_phi8) {
          T_SCALAR_ARRAY *phi8 = cache_phi8->data();
          T_SCALAR_ARRAY::Exchange_Arrays(*phi8, phi_ghost_bandwidth_eight);
          cm->ReleaseAccess(cache_phi8);
          cache_phi8 = NULL;
      }
      // last unique particle id
      const std::string lupistring = std::string(APP_LAST_UNIQUE_PARTICLE_ID);
      if (Data *d = application::GetFirstData(lupistring, da)) {
          nimbus::ScalarData<int> *sd = static_cast<nimbus::ScalarData<int> * >(d);
          sd->set_scalar(particle_levelset.last_unique_particle_id);
      }

      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_DT), da, application::WRITE_ACCESS);
      if (data_temp) {
        nimbus::ScalarData<float>* data_real =
            dynamic_cast<nimbus::ScalarData<float>*>(data_temp);
        data_real->set_scalar(dt_buffer);
        dbg(APP_LOG, "[Data Saving]%s: %0.9f\n", APP_DT, dt_buffer);
      }

      // ** there should not be any accesses to particle levelset after this **
      if (cache_ple) {
          if (data_config.GetFlag(DataConfig::SHARED_PARTICLES_FLUSH)) {
              nimbus::cache::type_id_t vars[] = { application::POS,
                                                  application::NEG,
                                                  application::POS_REM,
                                                  application::NEG_REM };
              std::vector<nimbus::cache::type_id_t> var_type(
                  vars, vars + sizeof(vars)/sizeof(nimbus::cache::type_id_t));
              std::vector<nimbus::DataArray> write(4), shared(4);
              std::string dtype[] = { APP_POS_PARTICLES,
                                      APP_NEG_PARTICLES,
                                      APP_POS_REM_PARTICLES,
                                      APP_NEG_REM_PARTICLES
                                    };
              for (size_t t = 0; t < application::NUM_PARTICLE_TYPES; ++t)
                  application::GetWriteData(*job, dtype[t], da, &write[t], false);
              nimbus::GeometricRegion inner_reg(
                  array_reg_central.NewEnlarged(-application::kGhostNum));
              for (size_t t = 0; t < application::NUM_PARTICLE_TYPES; ++t) {
                  nimbus::DataArray &write_t = write[t];
                  nimbus::DataArray &shared_t = shared[t];
                  for (size_t i = 0; i < write_t.size(); ++i) {
                      nimbus::Data *d = write_t[i];
                      nimbus::GeometricRegion dr = d->region();
                      if (!inner_reg.Covers(&dr)) {
                          shared_t.push_back(d);
                      }
                  }
              }
              cache_ple->WriteImmediately(var_type, shared);
          }
          cm->ReleaseAccess(cache_ple);
          cache_ple = NULL;
      }
    }

    // psi_d.
    if (cache_psi_d) {
        BOOL_SCALAR_ARRAY *psi_d = cache_psi_d->data();
        BOOL_SCALAR_ARRAY::Exchange_Arrays(*psi_d, projection.laplace->psi_D);
        nimbus::DataArray write;
        cm->ReleaseAccess(cache_psi_d);
        cache_psi_d = NULL;
    }

    // psi_n.
    if (cache_psi_n) {
        BOOL_FACE_ARRAY *psi_n = cache_psi_n->data();
        BOOL_FACE_ARRAY::Exchange_Arrays(*psi_n, projection.laplace->psi_N);
        nimbus::DataArray write;
        cm->ReleaseAccess(cache_psi_n);
        cache_psi_n = NULL;
    }

    // pressure.
    if (cache_pressure) {
      T_SCALAR_ARRAY* pressure = cache_pressure->data();
      T_SCALAR_ARRAY::Exchange_Arrays(*pressure, projection.p);
      cm->ReleaseAccess(cache_pressure);
      cache_pressure = NULL;
    }
    // colors.
    if (cache_colors) {
      typedef typename PhysBAM::ARRAY<int, TV_INT> INT_SCALAR_ARRAY;
      INT_SCALAR_ARRAY* colors = cache_colors->data();
      INT_SCALAR_ARRAY::Exchange_Arrays(
          *colors, projection.laplace->filled_region_colors);
      cm->ReleaseAccess(cache_colors);
      cache_colors = NULL;
    }
    // divergence.
    if (cache_divergence) {
      T_SCALAR_ARRAY* divergence = cache_divergence->data();
      T_SCALAR_ARRAY::Exchange_Arrays(*divergence, projection.laplace->f);
      cm->ReleaseAccess(cache_divergence);
      cache_divergence = NULL;
    }

    // TODO(addcache) the following data translation is implemented by memcpy,
    // caching might not be needed.
    typedef nimbus::Data Data;
    if (data_config.GetFlag(DataConfig::MATRIX_A)) {
      // TODO(quhang) swap rather than add.
      *(cache_matrix_a->data()) = laplace_solver_wrapper.A_array(1);
      cm->ReleaseAccess(cache_matrix_a);
      cache_matrix_a = NULL;
    }
    if (data_config.GetFlag(DataConfig::VECTOR_B)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_VECTOR_B), da, application::WRITE_ACCESS);
      if (data_temp) {
        application::DataRawVectorNd* data_real =
            dynamic_cast<application::DataRawVectorNd*>(data_temp);
        data_real->SaveToNimbus(laplace_solver_wrapper.b_array(1));
        dbg(APP_LOG, "Finish writing VECTOR_B.\n");
      }
    }
    if (data_config.GetFlag(DataConfig::INDEX_C2M)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_INDEX_C2M), da, application::WRITE_ACCESS);
      if (data_temp) {
        application::DataRawGridArray* data_real =
            dynamic_cast<application::DataRawGridArray*>(data_temp);
        data_real->SaveToNimbus(
            laplace_solver_wrapper.cell_index_to_matrix_index);
        dbg(APP_LOG, "Finish writing INDEX_C2M.\n");
      }
    }
    if (data_config.GetFlag(DataConfig::INDEX_M2C)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_INDEX_M2C), da, application::WRITE_ACCESS);
      if (data_temp) {
        application::DataRawArrayM2C* data_real =
            dynamic_cast<application::DataRawArrayM2C*>(data_temp);
        data_real->SaveToNimbus(
            laplace_solver_wrapper.matrix_index_to_cell_index_array(1));
        dbg(APP_LOG, "Finish writing INDEX_M2C.\n");
      }
    }
    if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_TOLERANCE)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_PROJECTION_LOCAL_TOLERANCE),
          da, application::WRITE_ACCESS);
      if (data_temp) {
        nimbus::ScalarData<float>* data_real =
            dynamic_cast<nimbus::ScalarData<float>*>(data_temp);
        data_real->set_scalar(projection.elliptic_solver->tolerance);
        dbg(APP_LOG, "Finish writing PROJECTION_LOCAL_TOLERANCE.\n");
      }
    }
    if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_N)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_PROJECTION_LOCAL_N),
          da, application::WRITE_ACCESS);
      if (data_temp) {
        nimbus::ScalarData<int>* data_real =
            dynamic_cast<nimbus::ScalarData<int>*>(data_temp);
        data_real->set_scalar(laplace_solver_wrapper.local_n);
        dbg(APP_LOG, "Finish writing PROJECTION_LOCAL_N.\n");
      }
    }
    if (data_config.GetFlag(DataConfig::PROJECTION_INTERIOR_N)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_PROJECTION_INTERIOR_N),
          da, application::WRITE_ACCESS);
      if (data_temp) {
        nimbus::ScalarData<int>* data_real =
            dynamic_cast<nimbus::ScalarData<int>*>(data_temp);
        data_real->set_scalar(laplace_solver_wrapper.interior_n);
        dbg(APP_LOG, "Finish writing PROJECTION_INTERIOR_N.\n");
      }
    }
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Load_From_Nimbus_No_Cache(const nimbus::Job *job, const nimbus::DataArray &da, const int frame)
{
    nimbus::int_dimension_t array_shift[3] = {
        local_region.x() - 1, local_region.y() - 1, local_region.z() - 1};
    nimbus::PdiVector pdv;

    GeometricRegion array_reg_central(local_region);
    GeometricRegion array_reg_outer(array_reg_central.NewEnlarged(application::kGhostNum));
    GeometricRegion array_reg_thin_outer(array_reg_central.NewEnlarged(1));
    GeometricRegion array_reg_outer_7(array_reg_central.NewEnlarged(7));
    GeometricRegion array_reg_outer_8(array_reg_central.NewEnlarged(8));

    GeometricRegion enlarge(1-application::kGhostNum,
                            1-application::kGhostNum,
                            1-application::kGhostNum,
                            local_region.dx()+2*application::kGhostNum,
                            local_region.dy()+2*application::kGhostNum,
                            local_region.dz()+2*application::kGhostNum);

    // mac velocities
    const std::string fvstring = std::string(APP_FACE_VEL);
    if (application::GetTranslatorData(job, fvstring, da, &pdv, application::READ_ACCESS)
        && data_config.GetFlag(DataConfig::VELOCITY)) {
      translator.ReadFaceArrayFloat(
          &array_reg_central, array_shift, &pdv, &face_velocities);
    }
    application::DestroyTranslatorObjects(&pdv);

    // mac velocities ghost
    const std::string fvgstring = std::string(APP_FACE_VEL_GHOST);
    if (application::GetTranslatorData(job, fvgstring, da, &pdv, application::READ_ACCESS)
        && data_config.GetFlag(DataConfig::VELOCITY_GHOST)) {
      translator.ReadFaceArrayFloat(
          &array_reg_outer, array_shift, &pdv, &face_velocities_ghost);
    }
    application::DestroyTranslatorObjects(&pdv);

    // particle leveset quantities
    T_PARTICLE_LEVELSET& particle_levelset = particle_levelset_evolution.particle_levelset;

    // levelset
    const std::string lsstring = std::string(APP_PHI);
    if (application::GetTranslatorData(job, lsstring, da, &pdv, application::READ_ACCESS)
        && data_config.GetFlag(DataConfig::LEVELSET)) {
      translator.ReadScalarArrayFloat(
          &array_reg_outer,
          array_shift,
          &pdv,
          &particle_levelset.levelset.phi);
      std::cout << "OMID: Read 3.\n";
    }
    application::DestroyTranslatorObjects(&pdv);
    if (application::GetTranslatorData(job, lsstring, da, &pdv, application::READ_ACCESS)
        && data_config.GetFlag(DataConfig::LEVELSET_READ)) {
      translator.ReadScalarArrayFloat(
          &array_reg_outer,
          array_shift,
          &pdv,
          &particle_levelset.levelset.phi);
      std::cout << "OMID: Read 3.\n";
    }
    application::DestroyTranslatorObjects(&pdv);
    if (application::GetTranslatorData(job, lsstring, da, &pdv, application::READ_ACCESS)
        && data_config.GetFlag(DataConfig::LEVELSET_BW_SEVEN_READ)) {
      translator.ReadScalarArrayFloat(
          &array_reg_outer_7,
          array_shift,
          &pdv,
          &phi_ghost_bandwidth_seven);
      std::cout << "OMID: Read 7.\n";
    }
    application::DestroyTranslatorObjects(&pdv);
    if (application::GetTranslatorData(job, lsstring, da, &pdv, application::READ_ACCESS)
        && data_config.GetFlag(DataConfig::LEVELSET_BW_EIGHT_READ)) {
      translator.ReadScalarArrayFloat(
          &array_reg_outer_8,
          array_shift,
          &pdv,
          &phi_ghost_bandwidth_eight);
      std::cout << "OMID: Read 8.\n";
    }
    application::DestroyTranslatorObjects(&pdv);

    // positive particles
    const std::string ppstring = std::string(APP_POS_PARTICLES);
    if (application::GetTranslatorData(job, ppstring, da, &pdv, application::READ_ACCESS)
        && data_config.GetFlag(DataConfig::POSITIVE_PARTICLE)) {
      translator.ReadParticles(
          &enlarge, array_shift,
          &pdv, particle_levelset, kScale, true);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating positive particles.\n");

    // negative particles
    const std::string npstring = std::string(APP_NEG_PARTICLES);
    if (application::GetTranslatorData(job, npstring, da, &pdv, application::READ_ACCESS)
        && data_config.GetFlag(DataConfig::NEGATIVE_PARTICLE)) {
      translator.ReadParticles(
          &enlarge, array_shift,
          &pdv, particle_levelset, kScale, false);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating negative particles.\n");

    // Removed positive particles.
    const std::string prpstring = std::string(APP_POS_REM_PARTICLES);
    if (application::GetTranslatorData(job, prpstring, da, &pdv, application::READ_ACCESS)
        && data_config.GetFlag(DataConfig::REMOVED_POSITIVE_PARTICLE)) {
      translator.ReadRemovedParticles(
          &enlarge, array_shift,
          &pdv, particle_levelset, kScale, true);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating remove positive particles.\n");

    // Removed negative particles.
    const std::string nrpstring = std::string(APP_NEG_REM_PARTICLES);
    if (application::GetTranslatorData(job, nrpstring, da, &pdv, application::READ_ACCESS)
        && data_config.GetFlag(DataConfig::REMOVED_NEGATIVE_PARTICLE)) {
      translator.ReadRemovedParticles(
          &enlarge, array_shift,
          &pdv, particle_levelset, kScale, false);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating remove negative particles.\n");

    // last unique particle id
    const std::string lupistring = std::string(APP_LAST_UNIQUE_PARTICLE_ID);
    if (Data *d = application::GetFirstData(lupistring, da)) {
        nimbus::ScalarData<int> *sd = static_cast<nimbus::ScalarData<int> * >(d);
        particle_levelset.last_unique_particle_id = sd->scalar();
    }
    dbg(APP_LOG, "Finish translating particle id.\n");

    // Calculate dt.
    if (data_config.GetFlag(DataConfig::DT)) {
        if (application::GetTranslatorData(
                job, std::string(APP_DT), da, &pdv,
                application::READ_ACCESS)) {
        dbg(APP_LOG, "Reducing DT min(");
        // TODO(quhang) maybe not safe. To be put the MAX float value.
        dt_buffer = 1e6;
        nimbus::PdiVector::const_iterator iter = pdv.begin();
        for (; iter != pdv.end(); ++iter) {
          const nimbus::PhysicalDataInstance* instance = *iter;
          nimbus::ScalarData<float>* data_real =
              dynamic_cast<nimbus::ScalarData<float>*>(instance->data());
          float value = data_real->scalar();
          dbg(APP_LOG, "%f ", value);
          if (value < dt_buffer) {
            dt_buffer = value;
          }
        }
        dbg(APP_LOG, ") = %f.\n", dt_buffer);
        }
    }
    dbg(APP_LOG, "Finish reduce dt.\n");

    // psi_d.
    const std::string psi_d_string = std::string(APP_PSI_D);
    if (application::GetTranslatorData(job, psi_d_string, da, &pdv, application::READ_ACCESS)
        && data_config.GetFlag(DataConfig::PSI_D)) {
      translator.ReadScalarArrayBool(
          &array_reg_thin_outer, array_shift, &pdv, &projection.laplace->psi_D);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating psi_d.\n");
    // psi_n.
    const std::string psi_n_string = std::string(APP_PSI_N);
    if (application::GetTranslatorData(job, psi_n_string, da, &pdv, application::READ_ACCESS)
        && data_config.GetFlag(DataConfig::PSI_N)) {
      translator.ReadFaceArrayBool(
          &array_reg_thin_outer, array_shift, &pdv, &projection.laplace->psi_N);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating psi_n.\n");
    // pressure.
    const std::string pressure_string = std::string(APP_PRESSURE);
    if (application::GetTranslatorData(job, pressure_string, da, &pdv, application::READ_ACCESS)
        && data_config.GetFlag(DataConfig::PRESSURE)) {
      translator.ReadScalarArrayFloat(
          &array_reg_thin_outer, array_shift, &pdv, &projection.p);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating pressure.\n");
    // filled_region_colors.
    const std::string filled_region_colors_string =
        std::string(APP_FILLED_REGION_COLORS);
    if (application::GetTranslatorData(job, filled_region_colors_string, da, &pdv, application::READ_ACCESS)
        && data_config.GetFlag(DataConfig::REGION_COLORS)) {
      dbg(APP_LOG, "filled_region_colors is being read from Nimbus.\n");
      translator.ReadScalarArrayInt(
          &array_reg_thin_outer, array_shift, &pdv,
          &projection.laplace->filled_region_colors);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating filled_region_colors.\n");
    // divergence.
    const std::string divergence_string =
        std::string(APP_DIVERGENCE);
    if (application::GetTranslatorData(job, divergence_string, da, &pdv, application::READ_ACCESS)
        && data_config.GetFlag(DataConfig::DIVERGENCE)) {
      translator.ReadScalarArrayFloat(
          &array_reg_thin_outer, array_shift, &pdv,
          &projection.laplace->f);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating divergence.\n");

    typedef nimbus::Data Data;
    if (data_config.GetFlag(DataConfig::MATRIX_A)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_MATRIX_A), da, application::READ_ACCESS);
      if (data_temp) {
        application::DataSparseMatrix* data_real =
            dynamic_cast<application::DataSparseMatrix*>(data_temp);
        data_real->LoadFromNimbus(&laplace_solver_wrapper.A_array(1));
        dbg(APP_LOG, "Finish reading MATRIX_A.\n");
      }
    }
    if (data_config.GetFlag(DataConfig::VECTOR_B)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_VECTOR_B), da, application::READ_ACCESS);
      if (data_temp) {
        application::DataRawVectorNd* data_real =
            dynamic_cast<application::DataRawVectorNd*>(data_temp);
        data_real->LoadFromNimbus(&laplace_solver_wrapper.b_array(1));
        dbg(APP_LOG, "Finish reading VECTOR_B.\n");
      }
    }
    if (data_config.GetFlag(DataConfig::INDEX_C2M)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_INDEX_C2M), da, application::READ_ACCESS);
      if (data_temp) {
        application::DataRawGridArray* data_real =
            dynamic_cast<application::DataRawGridArray*>(data_temp);
        data_real->LoadFromNimbus(
            &laplace_solver_wrapper.cell_index_to_matrix_index);
        dbg(APP_LOG, "Finish reading INDEX_C2M.\n");
      }
    }
    if (data_config.GetFlag(DataConfig::INDEX_M2C)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_INDEX_M2C), da, application::READ_ACCESS);
      if (data_temp) {
        application::DataRawArrayM2C* data_real =
            dynamic_cast<application::DataRawArrayM2C*>(data_temp);
        data_real->LoadFromNimbus(
            &laplace_solver_wrapper.matrix_index_to_cell_index_array(1));
        dbg(APP_LOG, "Finish reading INDEX_M2C.\n");
      }
    }
    if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_TOLERANCE)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_PROJECTION_LOCAL_TOLERANCE),
          da, application::READ_ACCESS);
      if (data_temp) {
        nimbus::ScalarData<float>* data_real =
            dynamic_cast<nimbus::ScalarData<float>*>(data_temp);
        projection.elliptic_solver->tolerance = data_real->scalar();
        dbg(APP_LOG, "Finish reading PROJECTION_LOCAL_TOLERANCE.\n");
      }
    }
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Load_From_Nimbus(const nimbus::Job *job, const nimbus::DataArray &da, const int frame)
{
    if (!(use_cache && application::kUseCache)) {
      Load_From_Nimbus_No_Cache(job, da, frame);
      return;
    }

    // nimbus::int_dimension_t array_shift[3] = {
    //    local_region.x() - 1, local_region.y() - 1, local_region.z() - 1};
    nimbus::PdiVector pdv;

    GeometricRegion array_reg_central(local_region);
    GeometricRegion array_reg_outer(array_reg_central.NewEnlarged(application::kGhostNum));
    GeometricRegion array_reg_thin_outer(array_reg_central.NewEnlarged(1));
    GeometricRegion array_reg_outer_7(array_reg_central.NewEnlarged(7));
    GeometricRegion array_reg_outer_8(array_reg_central.NewEnlarged(8));

    GeometricRegion enlarge(1-application::kGhostNum,
                            1-application::kGhostNum,
                            1-application::kGhostNum,
                            local_region.dx()+2*application::kGhostNum,
                            local_region.dy()+2*application::kGhostNum,
                            local_region.dz()+2*application::kGhostNum);

    // mac velocities
    if (cache_fv)
    {
        T_FACE_ARRAY *fv = cache_fv->data();
        T_FACE_ARRAY::Exchange_Arrays(*fv, face_velocities);
    }

    // mac velocities ghost
    if (cache_fvg)
    {
        T_FACE_ARRAY *fvg = cache_fvg->data();
        T_FACE_ARRAY::Exchange_Arrays(*fvg, face_velocities_ghost);
    }

    // particle leveset quantities
    T_PARTICLE_LEVELSET& particle_levelset = particle_levelset_evolution.particle_levelset;

    // levelset
    if (cache_phi3)
    {
        T_SCALAR_ARRAY *phi3 = cache_phi3->data();
        T_SCALAR_ARRAY::Exchange_Arrays(*phi3, particle_levelset.levelset.phi);
    }
    if (cache_phi7)
    {
        T_SCALAR_ARRAY *phi7 = cache_phi7->data();
        T_SCALAR_ARRAY::Exchange_Arrays(*phi7, phi_ghost_bandwidth_seven);
    }
    if (cache_phi8)
    {
        T_SCALAR_ARRAY *phi8 = cache_phi8->data();
        T_SCALAR_ARRAY::Exchange_Arrays(*phi8, phi_ghost_bandwidth_eight);
    }

    // last unique particle id
    const std::string lupistring = std::string(APP_LAST_UNIQUE_PARTICLE_ID);
    if (Data *d = application::GetFirstData(lupistring, da)) {
        nimbus::ScalarData<int> *sd = static_cast<nimbus::ScalarData<int> * >(d);
        particle_levelset.last_unique_particle_id = sd->scalar();
    }
    dbg(APP_LOG, "Finish translating particle id.\n");

    // Calculate dt.
    if (data_config.GetFlag(DataConfig::DT)) {
        if (application::GetTranslatorData(
                job, std::string(APP_DT), da, &pdv,
                application::READ_ACCESS)) {
        dbg(APP_LOG, "Reducing DT min(");
        // TODO(quhang) maybe not safe. To be put the MAX float value.
        dt_buffer = 1e6;
        nimbus::PdiVector::const_iterator iter = pdv.begin();
        for (; iter != pdv.end(); ++iter) {
          const nimbus::PhysicalDataInstance* instance = *iter;
          nimbus::ScalarData<float>* data_real =
              dynamic_cast<nimbus::ScalarData<float>*>(instance->data());
          float value = data_real->scalar();
          dbg(APP_LOG, "%f ", value);
          if (value < dt_buffer) {
            dt_buffer = value;
          }
        }
        dbg(APP_LOG, ") = %f.\n", dt_buffer);
        }
    }
    dbg(APP_LOG, "Finish reduce dt.\n");

    // psi_d.
    if (cache_psi_d)
    {
        BOOL_SCALAR_ARRAY *psi_d = cache_psi_d->data();
        BOOL_SCALAR_ARRAY::Exchange_Arrays(*psi_d, projection.laplace->psi_D);
    }

    // psi_n.
    if (cache_psi_n)
    {
        BOOL_FACE_ARRAY *psi_n = cache_psi_n->data();
        BOOL_FACE_ARRAY::Exchange_Arrays(*psi_n, projection.laplace->psi_N);
    }

    // pressure.
    if (cache_pressure) {
      T_SCALAR_ARRAY* pressure = cache_pressure->data();
      T_SCALAR_ARRAY::Exchange_Arrays(*pressure, projection.p);
    }
    // colors.
    if (cache_colors) {
      typedef typename PhysBAM::ARRAY<int, TV_INT> INT_SCALAR_ARRAY;
      INT_SCALAR_ARRAY* colors = cache_colors->data();
      INT_SCALAR_ARRAY::Exchange_Arrays(
          *colors, projection.laplace->filled_region_colors);
    }
    // divergence.
    if (cache_divergence) {
      T_SCALAR_ARRAY* divergence = cache_divergence->data();
      T_SCALAR_ARRAY::Exchange_Arrays(*divergence, projection.laplace->f);
    }

    typedef nimbus::Data Data;
    if (cache_matrix_a) {
      // TODO(quhang) matrix_a reading is not handled.
      // laplace_solver_wrapper.A_array(1) = *(cache_matrix_a->data());
    }
    // TODO(addcache), the following data uses memcpy, maybe doesn't need to be
    // cached.
    if (data_config.GetFlag(DataConfig::VECTOR_B)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_VECTOR_B), da, application::READ_ACCESS);
      if (data_temp) {
        application::DataRawVectorNd* data_real =
            dynamic_cast<application::DataRawVectorNd*>(data_temp);
        data_real->LoadFromNimbus(&laplace_solver_wrapper.b_array(1));
        dbg(APP_LOG, "Finish reading VECTOR_B.\n");
      }
    }
    if (data_config.GetFlag(DataConfig::INDEX_C2M)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_INDEX_C2M), da, application::READ_ACCESS);
      if (data_temp) {
        application::DataRawGridArray* data_real =
            dynamic_cast<application::DataRawGridArray*>(data_temp);
        data_real->LoadFromNimbus(
            &laplace_solver_wrapper.cell_index_to_matrix_index);
        dbg(APP_LOG, "Finish reading INDEX_C2M.\n");
      }
    }
    if (data_config.GetFlag(DataConfig::INDEX_M2C)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_INDEX_M2C), da, application::READ_ACCESS);
      if (data_temp) {
        application::DataRawArrayM2C* data_real =
            dynamic_cast<application::DataRawArrayM2C*>(data_temp);
        data_real->LoadFromNimbus(
            &laplace_solver_wrapper.matrix_index_to_cell_index_array(1));
        dbg(APP_LOG, "Finish reading INDEX_M2C.\n");
      }
    }
    if (data_config.GetFlag(DataConfig::PROJECTION_LOCAL_TOLERANCE)) {
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_PROJECTION_LOCAL_TOLERANCE),
          da, application::READ_ACCESS);
      if (data_temp) {
        nimbus::ScalarData<float>* data_real =
            dynamic_cast<nimbus::ScalarData<float>*>(data_temp);
        projection.elliptic_solver->tolerance = data_real->scalar();
        dbg(APP_LOG, "Finish reading PROJECTION_LOCAL_TOLERANCE.\n");
      }
    }
}
//#####################################################################
template class WATER_EXAMPLE<VECTOR<float,3> >;

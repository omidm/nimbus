//#####################################################################
// Copyright 2009, Michael Lentine, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <set>
#include <PhysBAM_Tools/Parallel_Computation/DOMAIN_ITERATOR_THREADED.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <pthread.h>
#include "application/smoke/app_utils.h"
#include "application/smoke/cache_prototypes.h"
#include "application/smoke/cache_options.h"
#include "application/smoke/data_include.h"
#include "application/smoke/options.h"
#include "application/smoke/parameters.h"
#include "application/smoke/reg_def.h"
#include "application/smoke/smoke_example.h"
#include "data/physbam/translator_physbam_old.h"
#include "data/scalar_data.h"
#include "shared/nimbus.h"
#include "worker/physical_data_instance.h"

// TODO(quhang) In three places where nimbus_thread_queue is introduced.

using namespace PhysBAM;
//#####################################################################
// SMOKE_EXAMPLE
//#####################################################################
template<class TV> SMOKE_EXAMPLE<TV>::
SMOKE_EXAMPLE(const STREAM_TYPE stream_type_input,
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
    output_directory(application::kOutputDir),
    restart(0),
    write_debug_data(false),
    number_of_ghost_cells(application::kGhostNum),
    cfl(1),
    mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),//incompressible_fluid_collection(mac_grid),
    projection(mac_grid, false, false, nimbus_thread_queue),
    advection_scalar(nimbus_thread_queue),
    boundary(0)
{
    for (int i = 1; i <= TV::dimension; i++) {
      domain_boundary(i)(1) = true;
      domain_boundary(i)(2) = true;
    }
    pthread_mutex_init(&lock, 0);
    use_cache   = false;
    cache_fv    = NULL;
    cache_fvg   = NULL;
    cache_psi_n = NULL;
    cache_dens = NULL;
    cache_dens_ghost = NULL;
    // cache_phi3  = NULL;
    // cache_phi7  = NULL;
    // cache_phi8  = NULL;
    cache_psi_d = NULL;
    // cache_ple   = NULL;
    cache_pressure = NULL;
    cache_colors = NULL;
    cache_divergence = NULL;
    // create_destroy_ple = true;
}

//#####################################################################
// SMOKE_EXAMPLE
//#####################################################################
template<class TV> SMOKE_EXAMPLE<TV>::
SMOKE_EXAMPLE(const STREAM_TYPE stream_type_input,
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
    output_directory(application::kOutputDir),
    restart(0),
    write_debug_data(false),
    number_of_ghost_cells(application::kGhostNum),
    cfl(1),
    mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),//incompressible_fluid_collection(mac \grid),
    projection(mac_grid, false, false, nimbus_thread_queue),
    advection_scalar(nimbus_thread_queue),
    boundary(0)

{
    for (int i = 1; i <= TV::dimension; i++) {
      domain_boundary(i)(1) = true;
      domain_boundary(i)(2) = true;
    }
    pthread_mutex_init(&lock, 0);
    use_cache   = false;
    cache_fv    = cache->fv;
    cache_fvg   = cache->fvg;
    cache_psi_n = cache->psi_n;
    cache_dens = cache->dens;
    cache_dens_ghost = cache->dens_ghost;
    // cache_phi3  = cache->phi3;
    // cache_phi7  = cache->phi7;
    // cache_phi8  = cache->phi8;
    cache_psi_d = cache->psi_d;
    // cache_ple   = cache->ple;
    cache_pressure = cache->pressure;
    cache_colors = cache->color;
    cache_divergence = cache->divergence;
    // create_destroy_ple = true;
}
//#####################################################################
// ~SMOKE_EXAMPLE
//#####################################################################
template<class TV> SMOKE_EXAMPLE<TV>::
~SMOKE_EXAMPLE()
{
    if (nimbus_thread_queue) {
      delete nimbus_thread_queue;
    }
}

//#####################################################################
// CFL
//#####################################################################
template<class TV> typename TV::SCALAR SMOKE_EXAMPLE<TV>::
CFL(ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities)
{
  T dt=FLT_MAX;
  DOMAIN_ITERATOR_THREADED_ALPHA<SMOKE_EXAMPLE<TV>,TV>(
    mac_grid.Domain_Indices(),
    nimbus_thread_queue).template Run<ARRAY<T,FACE_INDEX<TV::dimension> >&,T&>(*this,&SMOKE_EXAMPLE::CFL_Threaded,face_velocities,dt);
  return dt;
}

//#####################################################################
// CFL_Threaded
//#####################################################################
template<class TV> void SMOKE_EXAMPLE<TV>::
CFL_Threaded(RANGE<TV_INT>& domain, ARRAY<T, FACE_INDEX<TV::dimension> >& face_velocities, T& dt) {
  T dt_convection = 0;
  for (typename GRID<TV>::CELL_ITERATOR iterator(mac_grid, domain);
       iterator.Valid();
       iterator.Next()) {
    TV_INT cell = iterator.Cell_Index(); 
    T local_V_norm = 0;
    for (int axis = 1; axis <= GRID<TV>::dimension; axis++) {
      local_V_norm += mac_grid.one_over_dX[axis]*maxabs(
        face_velocities(axis, mac_grid.First_Face_Index_In_Cell(axis, cell)), 
	face_velocities(axis, mac_grid.Second_Face_Index_In_Cell(axis, cell)));
    }
    dt_convection = max(dt_convection, local_V_norm);
  }
  pthread_mutex_lock(&lock);
  dt = min(dt, (T)1.0/dt_convection);
  pthread_mutex_unlock(&lock);
}


// Sets the boundary conditions before projection. It might read levelset and
// velocity near the boundary. For velocity, it only reads to check if an index
// is valid. It will try to write to psi_D and psi_N of the whole region, and
// write to pressure and velocity in the boundary region.
// From what I read. --quhang
template<class TV> void SMOKE_EXAMPLE<TV>::
Set_Boundary_Conditions(const T time) {
  projection.elliptic_solver->psi_D.Fill(false);
  projection.elliptic_solver->psi_N.Fill(false);
  for (int axis = 1; axis <= TV::dimension; axis++) {
    for (int axis_side = 1; axis_side <= 2; axis_side++) {
      int side = 2 * (axis-1) + axis_side;
      if (domain_boundary(axis)(axis_side)) {
	TV_INT interior_cell_offset =
	  axis_side==1 ? TV_INT() : -TV_INT::Axis_Vector(axis);
        for (typename GRID<TV>::FACE_ITERATOR iterator(
                mac_grid, 1, GRID<TV>::BOUNDARY_REGION, side);
            iterator.Valid();
            iterator.Next()) {
          TV_INT cell = iterator.Face_Index() + interior_cell_offset;
	  TV_INT boundary_face = axis_side == 1 ? 
	    iterator.Face_Index()+TV_INT::Axis_Vector(axis) :
	    iterator.Face_Index() - TV_INT::Axis_Vector(axis);
	  projection.elliptic_solver->psi_D(cell)=true;projection.p(cell)=0;
	}
      }
    }
    for (typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);
            iterator.Valid();
            iterator.Next()) {
      if(source.Lazy_Inside(iterator.Location())){
	projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
	if(iterator.Axis()==2)face_velocities(iterator.Full_Index())=1;
	else face_velocities(iterator.Full_Index())=0;
      }
    }
  }
}

//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void SMOKE_EXAMPLE<TV>::
Write_Output_Files(const int frame, int rank)
{
    //if(!write_output_files) return;
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
      FILE_UTILITIES::Write_To_File(stream_type,rank_name+"/"+f+"/density",density);
      FILE_UTILITIES::Write_To_File(stream_type,rank_name+"/"+f+"/pressure",projection.p);
      FILE_UTILITIES::Write_To_File(stream_type,rank_name+"/"+f+"/psi_N",projection.elliptic_solver->psi_N);
      FILE_UTILITIES::Write_To_File(stream_type,rank_name+"/"+f+"/psi_D",projection.elliptic_solver->psi_D);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
#else
      PHYSBAM_FATAL_ERROR("Cannot read doubles");
#endif
    } else {
      std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
      FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities);
      FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",mac_grid);
      FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/density",density);
      FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/pressure",projection.p);
      FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_N",projection.elliptic_solver->psi_N);
      FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_D",projection.elliptic_solver->psi_D);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
#else
      PHYSBAM_FATAL_ERROR("Cannot read doubles");
#endif
    }
}

//#####################################################################
// Read_Output_Files
//#####################################################################
template<class TV> void SMOKE_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
  std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
  FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/"+f+"/density",density);
  std::string filename;
  filename=output_directory+"/"+f+"/mac_velocities";
  if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading mac_velocities "<<filename<<std::endl;LOG::filecout(ss.str());FILE_UTILITIES::Read_From_File(stream_type,filename,face_velocities);}
  filename=output_directory+"/"+f+"/pressure";
  if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading pressure "<<filename<<std::endl;LOG::filecout(ss.str());FILE_UTILITIES::Read_From_File(stream_type,filename,projection.p);}
}

//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void SMOKE_EXAMPLE<TV>::
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

    //TODO: save density and density ghost

    // density
    const std::string dstring = std::string(APP_DENSITY);
    if (application::GetTranslatorData(job, dstring, da, &pdv, application::WRITE_ACCESS)
	&& data_config.GetFlag(DataConfig::DENSITY)) {
      //TODO: translator stuff ???
      translator.WriteScalarArrayFloat(
          &array_reg_outer, array_shift, &pdv, &density);
    }
    
    //density ghost
    const std::string dgstring = std::string(APP_DENSITY_GHOST);
    if (application::GetTranslatorData(job, dgstring, da, &pdv, application::WRITE_ACCESS)
	&& data_config.GetFlag(DataConfig::DENSITY_GHOST)) {
      //TODO: translator stuff ???
      translator.WriteScalarArrayFloat(
          &array_reg_outer, array_shift, &pdv, &density_ghost);
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
template<class TV> void SMOKE_EXAMPLE<TV>::
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

    if (cache_dens) {
      T_SCALAR_ARRAY *dens = cache_dens->data();
      T_SCALAR_ARRAY::Exchange_Arrays(*dens, density);
      nimbus::DataArray write;
      cm->ReleaseAccess(cache_dens);
      cache_dens = NULL;
    }

    if (cache_dens_ghost) {
      T_SCALAR_ARRAY *dens_ghost = cache_dens_ghost->data();
      T_SCALAR_ARRAY::Exchange_Arrays(*dens_ghost, density_ghost);
      nimbus::DataArray write;
      cm->ReleaseAccess(cache_dens_ghost);
      cache_dens_ghost = NULL;
    }

    {
      // ??? what is GetTheOnlyData ???
      Data* data_temp = application::GetTheOnlyData(
          job, std::string(APP_DT), da, application::WRITE_ACCESS);
      if (data_temp) {
        nimbus::ScalarData<float>* data_real =
            dynamic_cast<nimbus::ScalarData<float>*>(data_temp);
        data_real->set_scalar(dt_buffer);
        dbg(APP_LOG, "[Data Saving]%s: %0.9f\n", APP_DT, dt_buffer);
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

    // TODO(addcache).
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
template<class TV> void SMOKE_EXAMPLE<TV>::
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

    // TODO: load density and density ghost

    // density
    const std::string dstring = std::string(APP_DENSITY);
    if (application::GetTranslatorData(job, dstring, da, &pdv, application::READ_ACCESS) 
	&& data_config.GetFlag(DataConfig::DENSITY)) {
      //TODO: translator stuff ???
      translator.ReadScalarArrayFloat(
          &array_reg_outer, array_shift, &pdv, &density);
				      
    }

    // density ghost
    const std::string dgstring = std::string(APP_DENSITY_GHOST);
    if (application::GetTranslatorData(job, dgstring, da, &pdv, application::READ_ACCESS) 
	&& data_config.GetFlag(DataConfig::DENSITY_GHOST)) {
      //TODO: translator stuff ???
      translator.ReadScalarArrayFloat(
          &array_reg_outer, array_shift, &pdv, &density);
    }

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
template<class TV> void SMOKE_EXAMPLE<TV>::
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

    if (cache_dens)
    {
	T_SCALAR_ARRAY *dens = cache_dens->data();
	T_SCALAR_ARRAY::Exchange_Arrays(*dens, density);
    }

    if (cache_dens_ghost)
    {
	T_SCALAR_ARRAY *dens_ghost = cache_dens_ghost->data();
	T_SCALAR_ARRAY::Exchange_Arrays(*dens_ghost, density_ghost);
    }

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

    // TODO(addcache), the following data uses memcpy, maybe doesn't need to be
    // cached.
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
template class SMOKE_EXAMPLE<VECTOR<float,3> >;

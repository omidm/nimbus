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
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 * Modifier for smoke: Andrew Lim <alim16@stanford.edu> 
 */

#include <string>
#include <vector>

#include "applications/physbam/smoke/app_utils.h"
#include "applications/physbam/smoke/app_data_include.h"
#include "applications/physbam/smoke/app_data_options.h"
#include "applications/physbam/smoke/app_data_prototypes.h"
#include "applications/physbam/smoke/data_names.h"
#include "applications/physbam/smoke/options.h"
#include "applications/physbam/smoke/parameters.h"
#include "applications/physbam/smoke/physbam_utils.h"
#include "applications/physbam/smoke/smoke_driver.h"
#include "applications/physbam/smoke/smoke_example.h"
#include "src/shared/geometric_region.h"
#include "src/shared/nimbus.h"

#define PHYSBAM_INIT_LOG

namespace application {

Range GridToRange(
    const GeometricRegion& global_region,
    const GeometricRegion& local_region) {
  TV start, end;
  start(1) = (float)(local_region.x() - 1) / (float)global_region.dx();
  start(2) = (float)(local_region.y() - 1) / (float)global_region.dy();
  start(3) = (float)(local_region.z() - 1) / (float)global_region.dz();
  end(1) =  (float)(local_region.x() + local_region.dx() - 1) / (float)global_region.dx();
  end(2) =  (float)(local_region.y() + local_region.dy() - 1) / (float)global_region.dy();
  end(3) =  (float)(local_region.z() + local_region.dz() - 1) / (float)global_region.dz();
  return Range(start, end);
}

void GetAppAppObjects(
    const InitConfig &init_config,
    const DataConfig &data_config,
    const nimbus::Job &job,
    const nimbus::DataArray &da,
    AppAppObjects *cache) {
  nimbus::GeometricRegion local_region = init_config.local_region;
  nimbus::GeometricRegion array_reg(local_region);
  nimbus::GeometricRegion array_reg_outer_1(array_reg.NewEnlarged(1));
  nimbus::GeometricRegion array_reg_outer_3(array_reg.NewEnlarged(kGhostNum));

  nimbus::AppDataManager *cm = job.GetAppDataManager();

  // vector_b.
  if (data_config.GetFlag(DataConfig::VECTOR_B)) {
    nimbus::DataArray read, write;
    const std::string vector_b_string = std::string(APP_VECTOR_B);
    GetReadData(job, vector_b_string, da, &read);
    GetWriteData(job, vector_b_string, da, &write);
    nimbus::AppVar* cache_var =
      cm->GetAppVar(
	  read, array_reg,
	  write, array_reg,
	  kAppDataVectorB, array_reg,
	  nimbus::app_data::EXCLUSIVE);
    cache->vector_b = dynamic_cast<AppDataVector*>(cache_var);
    assert(cache->vector_b != NULL);
  }
  // matrix_a.
  if (data_config.GetFlag(DataConfig::MATRIX_A)) {
    nimbus::DataArray read, write;
    const std::string matrix_a_string = std::string(APP_MATRIX_A);
    GetReadData(job, matrix_a_string, da, &read);
    GetWriteData(job, matrix_a_string, da, &write);
    nimbus::AppVar* cache_var =
      cm->GetAppVar(
		    read, array_reg,
		    write, array_reg,
		    kAppDataSparseMatrixA, array_reg,
		    nimbus::app_data::EXCLUSIVE);
    cache->matrix_a = dynamic_cast<AppDataSparseMatrix*>(cache_var);
    assert(cache->matrix_a != NULL);
  }
  // index_m2c.
  if (data_config.GetFlag(DataConfig::INDEX_M2C)) {
    nimbus::DataArray read, write;
    const std::string index_m2c_string = std::string(APP_INDEX_M2C);
    GetReadData(job, index_m2c_string, da, &read);
    GetWriteData(job, index_m2c_string, da, &write);
    nimbus::AppVar* cache_var =
      cm->GetAppVar(
		    read, array_reg,
		    write, array_reg,
		    kAppDataArrayM2C, array_reg,
		    nimbus::app_data::EXCLUSIVE);
    cache->index_m2c = dynamic_cast<AppDataArrayM2C*>(cache_var);
    assert(cache->index_m2c != NULL);
  }
  // pressure.
  if (data_config.GetFlag(DataConfig::PRESSURE)) {
    nimbus::DataArray read, write;
    const std::string pressure_string = std::string(APP_PRESSURE);
    GetReadData(job, pressure_string, da, &read);
    GetWriteData(job, pressure_string, da, &write);
    nimbus::AppVar* cache_var =
        cm->GetAppVar(
            read, array_reg_outer_1,
            write, array_reg_outer_1,
            kAppDataPressure, array_reg_outer_1,
            nimbus::app_data::EXCLUSIVE);
    cache->pressure = dynamic_cast<AppDataScalarArray<T>*>(cache_var);
    assert(cache->pressure != NULL);
  }
  // index_c2m.
  if (data_config.GetFlag(DataConfig::INDEX_C2M)) {
    nimbus::DataArray read, write;
    const std::string index_c2m_string = std::string(APP_INDEX_C2M);
    GetReadData(job, index_c2m_string, da, &read);
    GetWriteData(job, index_c2m_string, da, &write);
    nimbus::AppVar* cache_var =
      cm->GetAppVar(
	  read, array_reg,
	  write, array_reg,
	  kAppDataIndexC2M, array_reg,
	  nimbus::app_data::EXCLUSIVE);
    cache->index_c2m = dynamic_cast<AppDataRawGridArray*>(cache_var);
    assert(cache->index_c2m != NULL);
  }
  // filled_region_colors.
  if (data_config.GetFlag(DataConfig::REGION_COLORS)) {
    nimbus::DataArray read, write;
    const std::string color_string = std::string(APP_FILLED_REGION_COLORS);
    GetReadData(job, color_string, da, &read);
    GetWriteData(job, color_string, da, &write);
    nimbus::AppVar* cache_var =
        cm->GetAppVar(
            read, array_reg_outer_1,
            write, array_reg_outer_1,
            kAppDataColors, array_reg_outer_1,
            nimbus::app_data::EXCLUSIVE);
    cache->color = dynamic_cast<AppDataScalarArray<int>*>(cache_var);
    assert(cache->color != NULL);
  }
  // divergence.
  if (data_config.GetFlag(DataConfig::DIVERGENCE)) {
    nimbus::DataArray read, write;
    const std::string divergence_string = std::string(APP_DIVERGENCE);
    GetReadData(job, divergence_string, da, &read);
    GetWriteData(job, divergence_string, da, &write);
    nimbus::AppVar* cache_var =
        cm->GetAppVar(
            read, array_reg_outer_1,
            write, array_reg_outer_1,
            kAppDataDivergence, array_reg_outer_1,
            nimbus::app_data::EXCLUSIVE);
    cache->divergence = dynamic_cast<AppDataScalarArray<T>*>(cache_var);
    assert(cache->divergence != NULL);
  }
  // mac velocities
  if (data_config.GetFlag(DataConfig::VELOCITY))
  {
    nimbus::DataArray read, write;
    const std::string fvstring = std::string(APP_FACE_VEL);
    GetReadData(job, fvstring, da, &read);
    GetWriteData(job, fvstring, da, &write);
    nimbus::AppVar *cache_var =
      cm->GetAppVar(
          read, array_reg,
          write, array_reg,
          kAppDataFaceVel, array_reg,
          nimbus::app_data::EXCLUSIVE);
    cache->fv = dynamic_cast<AppDataFaceArray<T> *>(cache_var);
    assert(cache->fv != NULL);
  }
  // mac velocities ghost
  if (data_config.GetFlag(DataConfig::VELOCITY_GHOST))
  {
    nimbus::DataArray read, write;
    const std::string fvgstring = std::string(APP_FACE_VEL_GHOST);
    GetReadData(job, fvgstring, da, &read);
    GetWriteData(job, fvgstring, da, &write);
    nimbus::AppObject *cache_var =
      cm->GetAppVar(
          read, array_reg_outer_3,
          write, array_reg_outer_3,
          kAppDataFaceVelGhost, array_reg_outer_3,
          nimbus::app_data::EXCLUSIVE);
    cache->fvg = dynamic_cast<AppDataFaceArray<T> *>(cache_var);
    assert(cache->fvg != NULL);
  }  
  //density
  if (data_config.GetFlag(DataConfig::DENSITY))
  {
    nimbus::DataArray read, write;
    const std::string dstring = std::string(APP_DENSITY);
    GetReadData(job, dstring, da, &read);
    GetWriteData(job, dstring, da, &write);
    nimbus::AppObject *cache_var = 
      cm->GetAppVar(read, array_reg,
		    write, array_reg,
		    kAppDataDensity, array_reg, 
		    nimbus::app_data::EXCLUSIVE);
    cache->dens = dynamic_cast<AppDataScalarArray<T> *>(cache_var);
    assert(cache->dens != NULL);
  }
  //density ghost
  if (data_config.GetFlag(DataConfig::DENSITY_GHOST))
  {
    nimbus::DataArray read, write;
    const std::string dgstring = std::string(APP_DENSITY_GHOST);
    GetReadData(job, dgstring, da, &read);
    GetWriteData(job, dgstring, da, &write);
    nimbus::AppObject *cache_var = 
      cm->GetAppVar(read, array_reg_outer_3, 
		    write, array_reg_outer_3,
		    kAppDataDensityGhost, array_reg_outer_3,
		    nimbus::app_data::EXCLUSIVE);
    cache->dens_ghost = dynamic_cast<AppDataScalarArray<T> *>(cache_var);
    assert(cache->dens_ghost != NULL);
  }
  // psi_d.
  if (data_config.GetFlag(DataConfig::PSI_D))
  {
    nimbus::DataArray read, write;
    const std::string psi_d_string = std::string(APP_PSI_D);
    GetReadData(job, psi_d_string, da, &read);
    GetWriteData(job, psi_d_string, da, &write);
    nimbus::AppVar *cache_var =
      cm->GetAppVar(
          read, array_reg_outer_1,
          write, array_reg_outer_1,
          kAppDataPsiD, array_reg_outer_1,
          nimbus::app_data::EXCLUSIVE);
    cache->psi_d = dynamic_cast<AppDataScalarArray<bool> *>(cache_var);
    assert(cache->psi_d != NULL);
  }
  // psi_n.
  if (data_config.GetFlag(DataConfig::PSI_N))
  {
    nimbus::DataArray read, write;
    const std::string psi_n_string = std::string(APP_PSI_N);
    GetReadData(job, psi_n_string, da, &read);
    GetWriteData(job, psi_n_string, da, &write);
    nimbus::AppVar *cache_var =
      cm->GetAppVar(
          read, array_reg_outer_1,
          write, array_reg_outer_1,
          kAppDataPsiN, array_reg_outer_1,
          nimbus::app_data::EXCLUSIVE);
    cache->psi_n = dynamic_cast<AppDataFaceArray<bool> *>(cache_var);
    assert(cache->psi_n != NULL);
  }
}

bool InitializeExampleAndDriver(
    const InitConfig& init_config,
    const DataConfig& data_config,
    const nimbus::Job* job,
    const nimbus::DataArray& da,
    PhysBAM::SMOKE_EXAMPLE<TV>*& example,
    PhysBAM::SMOKE_DRIVER<TV>*& driver) {
//   static Log init_log(std::string("physbam-init-log"));
//
// #ifdef PHYSBAM_INIT_LOG
//   {
//     std::stringstream msg;
//     msg << "~~~ App InitializeExampleAndDriver start : " << init_log.GetTime();
//     init_log.WriteToFile(msg.str());
//   }
// #endif

  dbg(APP_LOG, "Enter initialize_example_driver.\n");
  dbg(APP_LOG, "Global region: %s\n", init_config.global_region.ToNetworkData().c_str());
  dbg(APP_LOG, "Local region: %s\n", init_config.local_region.ToNetworkData().c_str());
  {
    AppAppObjects cache;
    GetAppAppObjects(init_config, data_config, *job, da, &cache);
    example = new PhysBAM::SMOKE_EXAMPLE<TV>(PhysBAM::STREAM_TYPE((RW())),
				       &cache,
				       init_config.use_threading,
				       init_config.core_quota);
    // parameters for nimbus
    example->local_region = init_config.local_region;
    example->kScale = init_config.global_region.dx();
    example->relative_region.Rebuild(1, 1, 1,
        init_config.local_region.dx(),
        init_config.local_region.dy(),
        init_config.local_region.dz());
    // physbam intiialization
    example->Initialize_Grid(
        TV_INT(init_config.local_region.dx(),
          init_config.local_region.dy(),
          init_config.local_region.dz()),
        GridToRange(init_config.global_region, init_config.local_region));
    // add source
    TV point1=TV::All_Ones_Vector()*.2,point2=TV::All_Ones_Vector()*.3;point1(2)=0;point2(2)=.05;
    example->source.min_corner=point1;example->source.max_corner=point2;

    example->data_config.Set(data_config);
  }
  {
    application::ScopeTimer scope_timer("init_driver");
    driver= new PhysBAM::SMOKE_DRIVER<TV>(*example);
    // parameters
    driver->init_phase = init_config.init_phase;
    driver->current_frame = init_config.frame;
    driver->time = init_config.time;
    dbg(APP_LOG, "Before enter driver->Initialize.\n");
    // physbam initialization
    if (init_config.init_phase)
      driver->InitializeFirstDistributed(job, da);
    else
      driver->InitializeUseCache(job, da);
  }

  dbg(APP_LOG, "Exit initialize_example_driver.\n");

// #ifdef PHYSBAM_INIT_LOG
//   {
//     std::stringstream msg;
//     msg << "~~~ App InitializeExampleAndDriver end : " << init_log.GetTime();
//     init_log.WriteToFile(msg.str());
//   }
// #endif

  return true;
}

void DestroyExampleAndDriver(
    PhysBAM::SMOKE_EXAMPLE<TV>*& example,
    PhysBAM::SMOKE_DRIVER<TV>*& driver) {
  application::ScopeTimer scope_timer("delete_example_and_driver");
  delete example;
  example = NULL;
  delete driver;
  driver = NULL;
}

}  // namespace application

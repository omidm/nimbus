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
 */

#include <string>
#include <vector>

#include "application/smoke/app_utils.h"
#include "application/smoke/cache_data_include.h"
#include "application/smoke/cache_options.h"
#include "application/smoke/cache_prototypes.h"
#include "application/smoke/data_names.h"
#include "application/smoke/options.h"
#include "application/smoke/parameters.h"
#include "application/smoke/physbam_utils.h"
#include "application/smoke/smoke_driver.h"
#include "application/smoke/smoke_example.h"
#include "shared/geometric_region.h"
#include "shared/nimbus.h"

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

void GetAppCacheObjects(
    const InitConfig &init_config,
    const DataConfig &data_config,
    const nimbus::Job &job,
    const nimbus::DataArray &da,
    AppCacheObjects *cache) {
  nimbus::GeometricRegion local_region = init_config.local_region;
  nimbus::GeometricRegion array_reg(local_region);
  nimbus::GeometricRegion array_reg_outer_1(array_reg.NewEnlarged(1));
  nimbus::GeometricRegion array_reg_outer_3(array_reg.NewEnlarged(kGhostNum));
  nimbus::GeometricRegion array_reg_thin_outer(array_reg.NewEnlarged(1));

  nimbus::CacheManager *cm = job.GetCacheManager();

  // pressure.
  if (data_config.GetFlag(DataConfig::PRESSURE)) {
    nimbus::DataArray read, write;
    const std::string pressure_string = std::string(APP_PRESSURE);
    GetReadData(job, pressure_string, da, &read);
    GetWriteData(job, pressure_string, da, &write);
    nimbus::CacheVar* cache_var =
        cm->GetAppVar(
            read, array_reg_thin_outer,
            write, array_reg_thin_outer,
            kCachePressure, array_reg_thin_outer,
            nimbus::cache::EXCLUSIVE);
    cache->pressure = dynamic_cast<CacheScalarArray<T>*>(cache_var);
    assert(cache->pressure != NULL);
  }
  // filled_region_colors.
  if (data_config.GetFlag(DataConfig::REGION_COLORS)) {
    nimbus::DataArray read, write;
    const std::string color_string = std::string(APP_FILLED_REGION_COLORS);
    GetReadData(job, color_string, da, &read);
    GetWriteData(job, color_string, da, &write);
    nimbus::CacheVar* cache_var =
        cm->GetAppVar(
            read, array_reg_thin_outer,
            write, array_reg_thin_outer,
            kCacheColors, array_reg_thin_outer,
            nimbus::cache::EXCLUSIVE);
    cache->color = dynamic_cast<CacheScalarArray<int>*>(cache_var);
    assert(cache->color != NULL);
  }
  // divergence.
  if (data_config.GetFlag(DataConfig::DIVERGENCE)) {
    nimbus::DataArray read, write;
    const std::string divergence_string = std::string(APP_DIVERGENCE);
    GetReadData(job, divergence_string, da, &read);
    GetWriteData(job, divergence_string, da, &write);
    nimbus::CacheVar* cache_var =
        cm->GetAppVar(
            read, array_reg_thin_outer,
            write, array_reg_thin_outer,
            kCacheDivergence, array_reg_thin_outer,
            nimbus::cache::EXCLUSIVE);
    cache->divergence = dynamic_cast<CacheScalarArray<T>*>(cache_var);
    assert(cache->divergence != NULL);
  }
  // mac velocities
  if (data_config.GetFlag(DataConfig::VELOCITY))
  {
    nimbus::DataArray read, write;
    const std::string fvstring = std::string(APP_FACE_VEL);
    GetReadData(job, fvstring, da, &read);
    GetWriteData(job, fvstring, da, &write);
    nimbus::CacheVar *cache_var =
      cm->GetAppVar(
          read, array_reg,
          write, array_reg,
          kCacheFaceVel, array_reg,
          nimbus::cache::EXCLUSIVE);
    cache->fv = dynamic_cast<CacheFaceArray<T> *>(cache_var);
    assert(cache->fv != NULL);
  }
  // mac velocities ghost
  if (data_config.GetFlag(DataConfig::VELOCITY_GHOST))
  {
    nimbus::DataArray read, write;
    const std::string fvgstring = std::string(APP_FACE_VEL_GHOST);
    GetReadData(job, fvgstring, da, &read);
    GetWriteData(job, fvgstring, da, &write);
    nimbus::CacheObject *cache_var =
      cm->GetAppVar(
          read, array_reg_outer_3,
          write, array_reg_outer_3,
          kCacheFaceVelGhost, array_reg_outer_3,
          nimbus::cache::EXCLUSIVE);
    cache->fvg = dynamic_cast<CacheFaceArray<T> *>(cache_var);
    assert(cache->fvg != NULL);
  }  
  //density
  if (data_config.GetFlag(DataConfig::DENSITY))
  {
    nimbus::DataArray read, write;
    const std::string dstring = std::string(APP_DENSITY);
    GetReadData(job, dstring, da, &read);
    GetWriteData(job, dstring, da, &write);
    nimbus::CacheObject *cache_var = 
      cm->GetAppVar(read, array_reg,
		    write, array_reg,
		    kCacheDensity, array_reg, 
		    nimbus::cache::EXCLUSIVE);
    cache->dens = dynamic_cast<CacheScalarArray<T> *>(cache_var);
    assert(cache->dens != NULL);
  }
  //density ghost
  if (data_config.GetFlag(DataConfig::DENSITY_GHOST))
  {
    nimbus::DataArray read, write;
    const std::string dgstring = std::string(APP_DENSITY_GHOST);
    GetReadData(job, dgstring, da, &read);
    GetWriteData(job, dgstring, da, &write);
    nimbus::CacheObject *cache_var = 
      cm->GetAppVar(read, array_reg_outer_3, 
		    write, array_reg_outer_3,
		    kCacheDensityGhost, array_reg_outer_3,
		    nimbus::cache::EXCLUSIVE);
    cache->dens_ghost = dynamic_cast<CacheScalarArray<T> *>(cache_var);
    assert(cache->dens_ghost != NULL);
  }
  // psi_d.
  if (data_config.GetFlag(DataConfig::PSI_D))
  {
    nimbus::DataArray read, write;
    const std::string psi_d_string = std::string(APP_PSI_D);
    GetReadData(job, psi_d_string, da, &read);
    GetWriteData(job, psi_d_string, da, &write);
    nimbus::CacheVar *cache_var =
      cm->GetAppVar(
          read, array_reg_outer_1,
          write, array_reg_outer_1,
          kCachePsiD, array_reg_outer_1,
          nimbus::cache::EXCLUSIVE);
    cache->psi_d = dynamic_cast<CacheScalarArray<bool> *>(cache_var);
    assert(cache->psi_d != NULL);
  }
  // psi_n.
  if (data_config.GetFlag(DataConfig::PSI_N))
  {
    nimbus::DataArray read, write;
    const std::string psi_n_string = std::string(APP_PSI_N);
    GetReadData(job, psi_n_string, da, &read);
    GetWriteData(job, psi_n_string, da, &write);
    nimbus::CacheVar *cache_var =
      cm->GetAppVar(
          read, array_reg_outer_1,
          write, array_reg_outer_1,
          kCachePsiN, array_reg_outer_1,
          nimbus::cache::EXCLUSIVE);
    cache->psi_n = dynamic_cast<CacheFaceArray<bool> *>(cache_var);
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
  static Log init_log(std::string("physbam-init-log"));

#ifdef PHYSBAM_INIT_LOG
  {
    std::stringstream msg;
    msg << "~~~ App InitializeExampleAndDriver start : " << init_log.GetTime();
    init_log.WriteToFile(msg.str());
  }
#endif


  dbg(APP_LOG, "Enter initialize_example_driver.\n");
  dbg(APP_LOG, "Global region: %s\n", init_config.global_region.toString().c_str());
  dbg(APP_LOG, "Local region: %s\n", init_config.local_region.toString().c_str());
  {
    if (init_config.use_cache && kUseCache) {
      AppCacheObjects cache;
      GetAppCacheObjects(init_config, data_config, *job, da, &cache);
      example = new PhysBAM::SMOKE_EXAMPLE<TV>(PhysBAM::STREAM_TYPE((RW())),
					       &cache,
					       init_config.use_threading,
					       init_config.core_quota);
      example->use_cache = true;
    } else {
      example = new PhysBAM::SMOKE_EXAMPLE<TV>(PhysBAM::STREAM_TYPE((RW())),
                                               init_config.use_threading,
                                               init_config.core_quota);
      example->use_cache = false;
    }
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
    driver= new PhysBAM::SMOKE_DRIVER<TV>(*example);
    // parameters
    driver->init_phase = init_config.init_phase;
    driver->current_frame = init_config.frame;
    driver->time = init_config.time;
    dbg(APP_LOG, "Before enter driver->Initialize.\n");
    // physbam initialization
    if (init_config.init_phase)
      driver->InitializeFirstDistributed(job, da);
    else if (init_config.use_cache && kUseCache)
      driver->InitializeUseCache(job, da);
    else
      driver->Initialize(job, da);
  }

  dbg(APP_LOG, "Exit initialize_example_driver.\n");

#ifdef PHYSBAM_INIT_LOG
  {
    std::stringstream msg;
    msg << "~~~ App InitializeExampleAndDriver end : " << init_log.GetTime();
    init_log.WriteToFile(msg.str());
  }
#endif

  return true;
}

void DestroyExampleAndDriver(
    PhysBAM::SMOKE_EXAMPLE<TV>*& example,
    PhysBAM::SMOKE_DRIVER<TV>*& driver) {
  delete example;
  example = NULL;
  delete driver;
  driver = NULL;
}

}  // namespace application

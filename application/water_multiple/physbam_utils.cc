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

#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/cache_data_include.h"
#include "application/water_multiple/cache_options.h"
#include "application/water_multiple/cache_prototypes.h"
#include "application/water_multiple/data_names.h"
#include "application/water_multiple/options.h"
#include "application/water_multiple/parameters.h"
#include "application/water_multiple/physbam_utils.h"
#include "application/water_multiple/water_driver.h"
#include "application/water_multiple/water_example.h"
#include "application/water_multiple/water_sources.h"
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
  nimbus::GeometricRegion array_reg_outer_7(array_reg.NewEnlarged(7));
  nimbus::GeometricRegion array_reg_outer_8(array_reg.NewEnlarged(8));

  nimbus::CacheManager *cm = job.GetCacheManager();

  // matrix_a.
  if (data_config.GetFlag(DataConfig::MATRIX_A)) {
    nimbus::DataArray read, write;
    const std::string matrix_a_string = std::string(APP_MATRIX_A);
    GetReadData(job, matrix_a_string, da, &read);
    GetWriteData(job, matrix_a_string, da, &write);
    nimbus::CacheVar* cache_var =
        cm->GetAppVar(
            read, array_reg,
            write, array_reg,
            kCacheSparseMatrixA, array_reg,
            nimbus::cache::EXCLUSIVE);
    cache->matrix_a = dynamic_cast<CacheSparseMatrix*>(cache_var);
    assert(cache->matrix_a != NULL);
  }
  // pressure.
  if (data_config.GetFlag(DataConfig::PRESSURE)) {
    nimbus::DataArray read, write;
    const std::string pressure_string = std::string(APP_PRESSURE);
    GetReadData(job, pressure_string, da, &read);
    GetWriteData(job, pressure_string, da, &write);
    nimbus::CacheVar* cache_var =
        cm->GetAppVar(
            read, array_reg_outer_1,
            write, array_reg_outer_1,
            kCachePressure, array_reg_outer_1,
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
            read, array_reg_outer_1,
            write, array_reg_outer_1,
            kCacheColors, array_reg_outer_1,
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
            read, array_reg_outer_1,
            write, array_reg_outer_1,
            kCacheDivergence, array_reg_outer_1,
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
  // levelset
  {
    const std::string lsstring = std::string(APP_PHI);
    nimbus::DataArray read3, read7, read8, write3, write7, write8;
    nimbus::DataArray read, write;
    bool l = data_config.GetFlag(DataConfig::LEVELSET);
    bool lr = data_config.GetFlag(DataConfig::LEVELSET_READ);
    bool lw = data_config.GetFlag(DataConfig::LEVELSET_WRITE);
    bool l7r = data_config.GetFlag(DataConfig::LEVELSET_BW_SEVEN_READ);
    bool l7w = data_config.GetFlag(DataConfig::LEVELSET_BW_SEVEN_WRITE);
    bool l8r = data_config.GetFlag(DataConfig::LEVELSET_BW_EIGHT_READ);
    bool l8w = data_config.GetFlag(DataConfig::LEVELSET_BW_EIGHT_WRITE);
    if (l || lr || lw || l7r || l7w || l8r || l8w) {
      GetReadData(job, lsstring, da, &read);
      GetWriteData(job, lsstring, da, &write);
    }
    if (l || lr) read3 = read;
    if (l || lw) write3 = write;
    if (l7r) read7 = read;
    if (l7w) write7 = write;
    if (l8r) read8 = read;
    if (l8w) write8 = write;
    nimbus::DataArray write_empty;
    if (l || lr || lw) {
      nimbus::CacheVar *cache_var =
        cm->GetAppVar(
            read3, array_reg_outer_3,
            write_empty, array_reg_outer_3,
            kCachePhi3, array_reg_outer_3,
            nimbus::cache::EXCLUSIVE);
      cache->phi3 = dynamic_cast<CacheScalarArray<T> *>(cache_var);
      assert(cache->phi3 != NULL);
    }
    if (l7r || l7w) {
      if (cache->phi3)
        cache->phi3->WriteImmediately(write7);
      nimbus::CacheVar *cache_var =
        cm->GetAppVar(
            read7, array_reg_outer_7,
            write7, array_reg_outer_7,
            kCachePhi7, array_reg_outer_7,
            nimbus::cache::EXCLUSIVE);
      cache->phi7 = dynamic_cast<CacheScalarArray<T> *>(cache_var);
      assert(cache->phi7 != NULL);
    }
    if (l8r || l8w) {
      if (cache->phi3)
        cache->phi3->WriteImmediately(write8);
      nimbus::CacheVar *cache_var =
        cm->GetAppVar(
            read8, array_reg_outer_8,
            write8, array_reg_outer_8,
            kCachePhi8, array_reg_outer_8,
            nimbus::cache::EXCLUSIVE);
      cache->phi8 = dynamic_cast<CacheScalarArray<T> *>(cache_var);
      assert(cache->phi8 != NULL);
    }
    if (!write3.empty()) {
      cm->DoSetUpWrite(cache->phi3, write3, array_reg_outer_3);
      // cache->phi3->SetUpWrite(write3, array_reg_outer_3);
    }
    // TODO(chinmayee): comment these later, not needed
    //if (!write7.empty()) {
    //  cm->DoSetUpWrite(cache->phi7, write7, array_reg_outer_7);
    //  // cache->phi7->SetUpWrite(write7, array_reg_outer_7);
    //}
    //if (!write8.empty()) {
    //  cm->DoSetUpWrite(cache->phi8, write8, array_reg_outer_8);
    //  // cache->phi8->SetUpWrite(write8, array_reg_outer_8);
    //}
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
  bool dflag[] = {  data_config.GetFlag(DataConfig::POSITIVE_PARTICLE),
                    data_config.GetFlag(DataConfig::NEGATIVE_PARTICLE),
                    data_config.GetFlag(DataConfig::REMOVED_POSITIVE_PARTICLE),
                    data_config.GetFlag(DataConfig::REMOVED_NEGATIVE_PARTICLE)
                 };
  if (dflag[POS] || dflag[NEG] || dflag[POS_REM] || dflag[NEG_REM]) {
    nimbus::cache::type_id_t vars[] = {POS, NEG, POS_REM, NEG_REM};
    std::vector<nimbus::cache::type_id_t> var_type(
        vars, vars + sizeof(vars)/sizeof(nimbus::cache::type_id_t));
    std::vector<nimbus::DataArray> read(4), write(4);
    std::string dtype[] = { APP_POS_PARTICLES,
                            APP_NEG_PARTICLES,
                            APP_POS_REM_PARTICLES,
                            APP_NEG_REM_PARTICLES
                          };
    for (size_t t = 0; t < NUM_PARTICLE_TYPES; ++t) {
      if (!dflag[t])
        continue;
      GetReadData(job, dtype[t], da, &read[t], false);
      GetWriteData(job, dtype[t], da, &write[t], false);
    }
    nimbus::CacheStruct *cache_struct =
      cm->GetAppStruct(
          var_type,
          read, array_reg_outer_3,
          write, array_reg_outer_3,
          kCachePLE, array_reg_outer_3,
          nimbus::cache::EXCLUSIVE);
    cache->ple = dynamic_cast<CacheParticleLevelsetEvolution<T> *>(cache_struct);
    assert(cache->ple != NULL);
  }
}

bool InitializeExampleAndDriver(
    const InitConfig& init_config,
    const DataConfig& data_config,
    const nimbus::Job* job,
    const nimbus::DataArray& da,
    PhysBAM::WATER_EXAMPLE<TV>*& example,
    PhysBAM::WATER_DRIVER<TV>*& driver) {
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
      if (cache.ple)
        example = new PhysBAM::WATER_EXAMPLE<TV>(PhysBAM::STREAM_TYPE((RW())),
                                                 &cache,
                                                 cache.ple->data(),
                                                 init_config.use_threading,
                                                 init_config.core_quota);
      else
        example = new PhysBAM::WATER_EXAMPLE<TV>(PhysBAM::STREAM_TYPE((RW())),
                                                 &cache,
                                                 init_config.use_threading,
                                                 init_config.core_quota);
      example->use_cache = true;
    } else {
      example = new PhysBAM::WATER_EXAMPLE<TV>(PhysBAM::STREAM_TYPE((RW())),
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
    PhysBAM::WaterSources::Add_Source(example);
    example->data_config.Set(data_config);
  }
  {
    driver= new PhysBAM::WATER_DRIVER<TV>(*example);
    // parameters
    driver->init_phase = init_config.init_phase;
    driver->current_frame = init_config.frame;
    driver->time = init_config.time;
    dbg(APP_LOG, "Before enter driver->Initialize.\n");
    // physbam initialization
    if (init_config.init_phase)
      driver->InitializeFirstDistributed(job, da);
      // driver->InitializeFirst(job, da);
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
    PhysBAM::WATER_EXAMPLE<TV>*& example,
    PhysBAM::WATER_DRIVER<TV>*& driver) {
  if (example->create_destroy_ple)
    delete &example->particle_levelset_evolution;
  delete example;
  example = NULL;
  delete driver;
  driver = NULL;
}

}  // namespace application

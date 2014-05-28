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

#include "application/water_multiple/app_utils.h"
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
  nimbus::GeometricRegion array_reg_outer_3(array_reg.NewEnlarged(application::kGhostNum));
  nimbus::GeometricRegion array_reg_outer_7(array_reg.NewEnlarged(7));
  nimbus::GeometricRegion array_reg_outer_8(array_reg.NewEnlarged(8));

  dbg(DBG_WARN, "\n--- *** --- LOAD \n");
  nimbus::CacheManager *cm = job.GetCacheManager();

  // mac velocities
  if (data_config.GetFlag(DataConfig::VELOCITY))
  {
    nimbus::DataArray read, write;
    const std::string fvstring = std::string(APP_FACE_VEL);
    GetReadData(job, fvstring, da, &read);
    GetWriteData(job, fvstring, da, &write);
    dbg(DBG_WARN, "\n--- Requesting %i elements into face velocity for region %s\n",
        read.size(), array_reg.toString().c_str());
    nimbus::CacheObject *cache_obj =
      cm->GetAppObject(read, write,
          array_reg,
          application::kCacheFaceVel,
          nimbus::EXCLUSIVE, write.empty());
    cache->fv = dynamic_cast<CacheFaceArray<T> *>(cache_obj);
    assert(cache->fv != NULL);
    dbg(APP_LOG, "Finish translating velocity.\n");
  }
  // mac velocities ghost
  if (data_config.GetFlag(DataConfig::VELOCITY_GHOST))
  {
    nimbus::DataArray read, write;
    const std::string fvgstring = std::string(APP_FACE_VEL_GHOST);
    GetReadData(job, fvgstring, da, &read);
    GetWriteData(job, fvgstring, da, &write);
    dbg(DBG_WARN, "\n--- Requesting %i elements into face velocity ghost for region %s\n",
        read.size(), array_reg_outer_3.toString().c_str());
    nimbus::CacheObject *cache_obj =
      cm->GetAppObject(read, write,
          array_reg_outer_3,
          application::kCacheFaceVelGhost,
          nimbus::EXCLUSIVE, write.empty());
    cache->fvg = dynamic_cast<CacheFaceArray<T> *>(cache_obj);
    assert(cache->fvg != NULL);
    dbg(APP_LOG, "Finish translating ghost velocity.\n");
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
    if (l) {
      read3 = read;
      write3 = write;
    }
    if (lr)
      read3 = read;
    if (lw)
      write3 = write;
    if (l7r)
      read7 = read;
    if (l7w)
      write7 = write;
    if (l8r)
      read8 = read;
    if (l8w)
      write8 = write;
    // TODO(Chinmayee): revisit this design later. Multiple cache objects for
    // same physical data.
    int order[3];
    if (!write3.empty()) {
      order[0] = 7; order[1] = 8; order[2] = 3;
    }
    else if (!write7.empty()) {
      order[0] = 3; order[1] = 8; order[2] = 7;
    }
    else {
      order[0] = 3; order[1] = 7; order[2] = 8;
    }
    for (size_t i = 0; i < 3; ++i) {
      if ((!(read3.empty() && write3.empty())) && order[i] == 3)
      {
        dbg(DBG_WARN, "\n--- Requesting %i elements into levelset 3 for region %s\n",
            read3.size(), array_reg_outer_3.toString().c_str());
        nimbus::CacheObject *cache_obj =
          cm->GetAppObject(read3, write3,
              array_reg_outer_3,
              application::kCachePhi3,
              nimbus::EXCLUSIVE, write3.empty() && write7.empty() && write8.empty());
        cache->phi3 = dynamic_cast<CacheScalarArray<T> *>(cache_obj);
        assert(cache->phi3 != NULL);
        dbg(APP_LOG, "Finish translating velocity levelset 3.\n");
      }
      if ((!(read7.empty() && write7.empty())) && order[i] == 7) {
        dbg(DBG_WARN, "\n--- Requesting %i elements into levelset 7 for region %s\n",
            read7.size(), array_reg_outer_7.toString().c_str());
        nimbus::CacheObject *cache_obj =
          cm->GetAppObject(read7, write7,
              array_reg_outer_7,
              application::kCachePhi7,
              nimbus::EXCLUSIVE, false);
        cache->phi7 = dynamic_cast<CacheScalarArray<T> *>(cache_obj);
        assert(cache->phi7 != NULL);
        dbg(APP_LOG, "Finish translating velocity levelset 7.\n");
      }
      if ((!(read8.empty() && write8.empty())) && order[i] == 8) {
        dbg(DBG_WARN, "\n--- Requesting %i elements into levelset 8 for region %s\n",
            read8.size(), array_reg_outer_8.toString().c_str());
        nimbus::CacheObject *cache_obj =
          cm->GetAppObject(read8, write8,
              array_reg_outer_8,
              application::kCachePhi8,
              nimbus::EXCLUSIVE, false);
        cache->phi8 = dynamic_cast<CacheScalarArray<T> *>(cache_obj);
        assert(cache->phi8 != NULL);
        dbg(APP_LOG, "Finish translating levelset 8.\n");
      }
    }
  }
  // psi_d.
  if (data_config.GetFlag(DataConfig::PSI_D))
  {
    nimbus::DataArray read, write;
    const std::string psi_d_string = std::string(APP_PSI_D);
    GetReadData(job, psi_d_string, da, &read);
    GetWriteData(job, psi_d_string, da, &write);
    dbg(DBG_WARN, "\n--- Requesting %i elements into psi_d for region %s\n",
        read.size(), array_reg_outer_1.toString().c_str());
    nimbus::CacheObject *cache_obj =
      cm->GetAppObject(read, write,
          array_reg_outer_1,
          application::kCachePsiD,
          nimbus::EXCLUSIVE, write.empty());
    cache->psi_d = dynamic_cast<CacheScalarArray<bool> *>(cache_obj);
    assert(cache->psi_d != NULL);
    dbg(APP_LOG, "Finish translating psi_d.\n");
  }
  // psi_n.
  if (data_config.GetFlag(DataConfig::PSI_N))
  {
    nimbus::DataArray read, write;
    const std::string psi_n_string = std::string(APP_PSI_N);
    GetReadData(job, psi_n_string, da, &read);
    GetWriteData(job, psi_n_string, da, &write);
    dbg(DBG_WARN, "\n--- Requesting %i elements into psi_n for region %s\n",
        read.size(), array_reg_outer_1.toString().c_str());
    nimbus::CacheObject *cache_obj =
      cm->GetAppObject(read, write,
          array_reg_outer_1,
          application::kCachePsiN,
          nimbus::EXCLUSIVE, write.empty());
    cache->psi_n = dynamic_cast<CacheFaceArray<bool> *>(cache_obj);
    assert(cache->psi_n != NULL);
    dbg(APP_LOG, "Finish translating psi_n.\n");
  }
  if (data_config.GetFlag(DataConfig::POSITIVE_PARTICLE) ||
      data_config.GetFlag(DataConfig::NEGATIVE_PARTICLE) ||
      data_config.GetFlag(DataConfig::REMOVED_POSITIVE_PARTICLE) ||
      data_config.GetFlag(DataConfig::REMOVED_NEGATIVE_PARTICLE))
  {
    nimbus::DataArray read, write;
    const std::string pp_string = std::string(APP_POS_PARTICLES);
    GetReadData(job, pp_string, da, &read, false);
    GetWriteData(job, pp_string, da, &write, false);
    const std::string np_string = std::string(APP_NEG_PARTICLES);
    GetReadData(job, np_string, da, &read, false);
    GetWriteData(job, np_string, da, &write, false);
    const std::string prp_string = std::string(APP_POS_REM_PARTICLES);
    GetReadData(job, prp_string, da, &read, false);
    GetWriteData(job, prp_string, da, &write, false);
    const std::string nrp_string = std::string(APP_NEG_REM_PARTICLES);
    GetReadData(job, nrp_string, da, &read, false);
    GetWriteData(job, nrp_string, da, &write, false);
    dbg(DBG_WARN, "\n--- Requesting %i elements into particles for region %s\n",
        read.size(), array_reg_outer_3.toString().c_str());
    nimbus::CacheObject *cache_obj =
      cm->GetAppObject(read, write,
          array_reg_outer_3,
          application::kCachePLE,
          nimbus::EXCLUSIVE, false, true);
    cache->ple = dynamic_cast<CacheParticleLevelsetEvolution<T> *>(cache_obj);
    assert(cache->ple != NULL);

    if (init_config.clear_read_shared_particles) {
      nimbus::GeometricRegion inner_reg(array_reg.NewEnlarged(-3));
      nimbus::DataArray shared;
      for (size_t k = 0; k < read.size(); ++k) {
        nimbus::Data *d = read[k];
        nimbus::GeometricRegion dr = d->region();
        if (!inner_reg.Covers(&dr))
          shared.push_back(d);
      }
      cache->ple->InvalidateCacheObject(shared);
    }
    dbg(APP_LOG, "Finish translating particles.\n");
  }
}

bool InitializeExampleAndDriver(
    const InitConfig& init_config,
    const DataConfig& data_config,
    const nimbus::Job* job,
    const nimbus::DataArray& da,
    PhysBAM::WATER_EXAMPLE<TV>*& example,
    PhysBAM::WATER_DRIVER<TV>*& driver) {
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
                                                 cache.ple->data());
      else
        example = new PhysBAM::WATER_EXAMPLE<TV>(PhysBAM::STREAM_TYPE((RW())),
                                                 &cache);
      example->use_cache = true;
    } else {
      example = new PhysBAM::WATER_EXAMPLE<TV>(PhysBAM::STREAM_TYPE((RW())));
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
      driver->InitializeFirst(job, da);
      // driver->InitializeFirstDistributed(job, da);
    else if (init_config.use_cache && kUseCache)
      driver->InitializeUseCache(job, da);
    else
      driver->Initialize(job, da);
  }

  dbg(APP_LOG, "Exit initialize_example_driver.\n");
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

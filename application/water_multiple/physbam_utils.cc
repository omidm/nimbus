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
#include "application/water_multiple/app_data_include.h"
#include "application/water_multiple/app_data_options.h"
#include "application/water_multiple/app_data_prototypes.h"
#include "application/water_multiple/data_names.h"
#include "application/water_multiple/options.h"
#include "application/water_multiple/parameters.h"
#include "application/water_multiple/physbam_include.h"
#include "application/water_multiple/physbam_utils.h"
#include "application/water_multiple/water_driver.h"
#include "application/water_multiple/water_example.h"
#include "application/water_multiple/water_sources.h"
#include "shared/fast_log.hh"
#include "shared/geometric_region.h"
#include "shared/nimbus.h"
#include "worker/worker_thread.h"
#include "worker/static_config_manager.h"
#include "worker/task_thread_pool.h"

#define PHYSBAM_INIT_LOG

namespace application {

void PrintGridDbg(const PhysBAM::GRID<PhysBAM::VECTOR<float, 3> >& grid) {
  printf("GRID: %d %d %d %f %f %f %f %f %f\n", grid.counts[1], grid.counts[2], grid.counts[3],
         grid.domain.min_corner[1], grid.domain.min_corner[2], grid.domain.min_corner[3], 
         grid.domain.max_corner[1], grid.domain.max_corner[2], grid.domain.max_corner[3]);
}

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
    AppAppObjects *app_data) {
  nimbus::timer::StartTimer(nimbus::timer::kAssemblingCache);
  nimbus::GeometricRegion local_region = init_config.local_region;
  nimbus::GeometricRegion array_reg(local_region);
  nimbus::GeometricRegion array_reg_outer_1(array_reg.NewEnlarged(1));
  nimbus::GeometricRegion array_reg_outer_3(array_reg.NewEnlarged(kGhostNum));
  nimbus::GeometricRegion array_reg_outer_7(array_reg.NewEnlarged(7));
  nimbus::GeometricRegion array_reg_outer_8(array_reg.NewEnlarged(8));

  nimbus::AppDataManager *cm = job.GetAppDataManager();

  // vector_b.
  if (data_config.GetFlag(DataConfig::VECTOR_B)) {
    nimbus::DataArray read, write;
    const std::string vector_b_string = std::string(APP_VECTOR_B);
    GetReadData(job, vector_b_string, da, &read);
    GetWriteData(job, vector_b_string, da, &write);
    nimbus::AppVar* app_var =
        cm->GetAppVar(
            read, array_reg,
            write, array_reg,
            kAppDataVectorB, array_reg,
            nimbus::app_data::SHARED);
    app_data->vector_b = dynamic_cast<AppDataVector*>(app_var);
    assert(app_data->vector_b != NULL);
  }
  // matrix_a.
  if (data_config.GetFlag(DataConfig::MATRIX_A)) {
    nimbus::DataArray read, write;
    const std::string matrix_a_string = std::string(APP_MATRIX_A);
    GetReadData(job, matrix_a_string, da, &read);
    GetWriteData(job, matrix_a_string, da, &write);
    nimbus::AppVar* app_var =
        cm->GetAppVar(
            read, array_reg,
            write, array_reg,
            kAppDataSparseMatrixA, array_reg,
            nimbus::app_data::SHARED);
    app_data->matrix_a = dynamic_cast<AppDataSparseMatrix*>(app_var);
    assert(app_data->matrix_a != NULL);
  }
  // index_m2c.
  if (data_config.GetFlag(DataConfig::INDEX_M2C)) {
    nimbus::DataArray read, write;
    const std::string index_m2c_string = std::string(APP_INDEX_M2C);
    GetReadData(job, index_m2c_string, da, &read);
    GetWriteData(job, index_m2c_string, da, &write);
    nimbus::AppVar* app_var =
        cm->GetAppVar(
            read, array_reg,
            write, array_reg,
            kAppDataArrayM2C, array_reg,
            nimbus::app_data::SHARED);
    app_data->index_m2c = dynamic_cast<AppDataArrayM2C*>(app_var);
    assert(app_data->index_m2c != NULL);
  }
  // pressure.
  if (data_config.GetFlag(DataConfig::PRESSURE)) {
    nimbus::DataArray read, write;
    const std::string pressure_string = std::string(APP_PRESSURE);
    GetReadData(job, pressure_string, da, &read);
    GetWriteData(job, pressure_string, da, &write);
    nimbus::AppVar* app_var =
        cm->GetAppVar(
            read, array_reg_outer_1,
            write, array_reg_outer_1,
            kAppDataPressure, array_reg_outer_1,
            nimbus::app_data::SHARED);
    app_data->pressure = dynamic_cast<AppDataScalarArray<T>*>(app_var);
    assert(app_data->pressure != NULL);
  }
  // index_c2m.
  if (data_config.GetFlag(DataConfig::INDEX_C2M)) {
    nimbus::DataArray read, write;
    const std::string index_c2m_string = std::string(APP_INDEX_C2M);
    GetReadData(job, index_c2m_string, da, &read);
    GetWriteData(job, index_c2m_string, da, &write);
    nimbus::AppVar* app_var =
        cm->GetAppVar(
            read, array_reg,
            write, array_reg,
            kAppDataIndexC2M, array_reg,
            nimbus::app_data::SHARED);
    app_data->index_c2m = dynamic_cast<AppDataRawGridArray*>(app_var);
    assert(app_data->index_c2m != NULL);
  }
  // filled_region_colors.
  if (data_config.GetFlag(DataConfig::REGION_COLORS)) {
    nimbus::DataArray read, write;
    const std::string color_string = std::string(APP_FILLED_REGION_COLORS);
    GetReadData(job, color_string, da, &read);
    GetWriteData(job, color_string, da, &write);
    nimbus::AppVar* app_var =
        cm->GetAppVar(
            read, array_reg_outer_1,
            write, array_reg_outer_1,
            kAppDataColors, array_reg_outer_1,
            nimbus::app_data::SHARED);
    app_data->color = dynamic_cast<AppDataScalarArray<int>*>(app_var);
    assert(app_data->color != NULL);
  }
  // divergence.
  if (data_config.GetFlag(DataConfig::DIVERGENCE)) {
    nimbus::DataArray read, write;
    const std::string divergence_string = std::string(APP_DIVERGENCE);
    GetReadData(job, divergence_string, da, &read);
    GetWriteData(job, divergence_string, da, &write);
    nimbus::AppVar* app_var =
        cm->GetAppVar(
            read, array_reg_outer_1,
            write, array_reg_outer_1,
            kAppDataDivergence, array_reg_outer_1,
            nimbus::app_data::SHARED);
    app_data->divergence = dynamic_cast<AppDataScalarArray<T>*>(app_var);
    assert(app_data->divergence != NULL);
  }
  // mac velocities
  if (data_config.GetFlag(DataConfig::VELOCITY))
  {
    nimbus::DataArray read, write;
    const std::string fvstring = std::string(APP_FACE_VEL);
    GetReadData(job, fvstring, da, &read);
    GetWriteData(job, fvstring, da, &write);
    nimbus::AppVar *app_var =
      cm->GetAppVar(
          read, array_reg,
          write, array_reg,
          kAppDataFaceVel, array_reg,
          nimbus::app_data::SHARED);
    app_data->fv = dynamic_cast<AppDataFaceArray<T> *>(app_var);
    assert(app_data->fv != NULL);
  }
  // mac velocities ghost
  if (data_config.GetFlag(DataConfig::VELOCITY_GHOST))
  {
    nimbus::DataArray read, write;
    const std::string fvgstring = std::string(APP_FACE_VEL_GHOST);
    GetReadData(job, fvgstring, da, &read);
    GetWriteData(job, fvgstring, da, &write);
    nimbus::AppObject *app_var =
      cm->GetAppVar(
          read, array_reg_outer_3,
          write, array_reg_outer_3,
          kAppDataFaceVelGhost, array_reg_outer_3,
          nimbus::app_data::SHARED);
    app_data->fvg = dynamic_cast<AppDataFaceArray<T> *>(app_var);
    assert(app_data->fvg != NULL);
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

    int order[3] = {3, 7, 8};

    if (l || lw) {
      order[0] = 7;
      order[1] = 8;
      order[2] = 3;
    } else if (l7w) {
      order[0] = 3;
      order[1] = 8;
      order[2] = 7;
    }

    for (int i = 0; i < 3; ++i) {
      if (order[i] == 3) {
        if (l || lr || lw) {
          nimbus::AppVar *app_var =
            cm->GetAppVar(
                read3, array_reg_outer_3,
                write3, array_reg_outer_3,
                kAppDataPhi3, array_reg_outer_3,
                nimbus::app_data::SHARED);
          app_data->phi3 = dynamic_cast<AppDataScalarArray<T> *>(app_var);
          assert(app_data->phi3 != NULL);
        }
      }
      if (order[i] == 7) {
        if (l7r || l7w) {
          nimbus::AppVar *app_var =
            cm->GetAppVar(
                read7, array_reg_outer_7,
                write7, array_reg_outer_7,
                kAppDataPhi7, array_reg_outer_7,
                nimbus::app_data::SHARED);
          app_data->phi7 = dynamic_cast<AppDataScalarArray<T> *>(app_var);
          assert(app_data->phi7 != NULL);
        }
      }
      if (order[i] == 8) {
        if (l8r || l8w) {
          nimbus::AppVar *app_var =
            cm->GetAppVar(
                read8, array_reg_outer_8,
                write8, array_reg_outer_8,
                kAppDataPhi8, array_reg_outer_8,
                nimbus::app_data::SHARED);
          app_data->phi8 = dynamic_cast<AppDataScalarArray<T> *>(app_var);
          assert(app_data->phi8 != NULL);
        }
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
    nimbus::AppVar *app_var =
      cm->GetAppVar(
          read, array_reg_outer_1,
          write, array_reg_outer_1,
          kAppDataPsiD, array_reg_outer_1,
          nimbus::app_data::SHARED);
    app_data->psi_d = dynamic_cast<AppDataScalarArray<bool> *>(app_var);
    assert(app_data->psi_d != NULL);
  }
  // psi_n.
  if (data_config.GetFlag(DataConfig::PSI_N))
  {
    nimbus::DataArray read, write;
    const std::string psi_n_string = std::string(APP_PSI_N);
    GetReadData(job, psi_n_string, da, &read);
    GetWriteData(job, psi_n_string, da, &write);
    nimbus::AppVar *app_var =
      cm->GetAppVar(
          read, array_reg_outer_1,
          write, array_reg_outer_1,
          kAppDataPsiN, array_reg_outer_1,
          nimbus::app_data::SHARED);
    app_data->psi_n = dynamic_cast<AppDataFaceArray<bool> *>(app_var);
    assert(app_data->psi_n != NULL);
  }
  bool dflag[] = {  data_config.GetFlag(DataConfig::POSITIVE_PARTICLE),
                    data_config.GetFlag(DataConfig::NEGATIVE_PARTICLE),
                    data_config.GetFlag(DataConfig::REMOVED_POSITIVE_PARTICLE),
                    data_config.GetFlag(DataConfig::REMOVED_NEGATIVE_PARTICLE)
                 };
  if (dflag[POS] || dflag[NEG] || dflag[POS_REM] || dflag[NEG_REM]) {
    nimbus::app_data::type_id_t vars[] = {POS, NEG, POS_REM, NEG_REM};
    std::vector<nimbus::app_data::type_id_t> var_type(
        vars, vars + sizeof(vars)/sizeof(nimbus::app_data::type_id_t));
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
    nimbus::AppStruct *app_struct =
      cm->GetAppStruct(
          var_type,
          read, array_reg_outer_3,
          write, array_reg_outer_3,
          kAppDataPLE, array_reg_outer_3,
          nimbus::app_data::SHARED);
    app_data->ple = dynamic_cast<AppDataParticleLevelsetEvolution<T> *>(app_struct);
    assert(app_data->ple != NULL);
  }

  nimbus::StaticConfigManager *config_manager = job.GetStaticConfigManager();
  if (data_config.GetFlag(DataConfig::VALID_MASK)) {
    app_data->static_config_valid_mask = dynamic_cast<StaticConfigValidMask*>(
        config_manager->GetStaticConfigVariable(STATIC_CONFIG_VALID_MASK,
                                                local_region));
    assert(app_data->static_config_valid_mask != NULL);
  }
  if (data_config.GetFlag(DataConfig::U_INTERFACE)) {
    app_data->static_config_u_interface = dynamic_cast<StaticConfigUInterface*>(
        config_manager->GetStaticConfigVariable(STATIC_CONFIG_U_INTERFACE,
                                                local_region));
    assert(app_data->static_config_u_interface != NULL);
  }
  app_data->static_config_force = dynamic_cast<StaticConfigForce*>(
      config_manager->GetStaticConfigVariable(STATIC_CONFIG_FORCE,
                                              local_region));
  assert(app_data->static_config_force != NULL);
  nimbus::timer::StopTimer(nimbus::timer::kAssemblingCache);
}

bool InitializeExampleAndDriver(
    const InitConfig& init_config,
    const DataConfig& data_config,
    const nimbus::Job* job,
    const nimbus::DataArray& da,
    PhysBAM::WATER_EXAMPLE<TV>*& example,
    PhysBAM::WATER_DRIVER<TV>*& driver) {
  dbg(APP_LOG, "Enter initialize_example_driver.\n");
  dbg(APP_LOG, "Global region: %s\n", init_config.global_region.ToNetworkData().c_str());
  dbg(APP_LOG, "Local region: %s\n", init_config.local_region.ToNetworkData().c_str());
  {
    application::ScopeTimer* example_scope_timer = NULL;
    // Branch on CACHE/ NO CACHE, depending on the stage
    if (!init_config.init_phase) {
      // CASE CACHE
      application::ScopeTimer* app_data_scope_timer =
          new application::ScopeTimer("app_data_lookup");
      AppAppObjects app_data;
      GetAppAppObjects(init_config, data_config, *job, da, &app_data);
      delete app_data_scope_timer;
      example_scope_timer =
          new application::ScopeTimer("init_example");

      StaticConfigCollisionBody* collision_body
       = dynamic_cast<StaticConfigCollisionBody*>(
          job->GetStaticConfigManager()
          ->GetStaticConfigVariable(STATIC_CONFIG_COLLISION_BODY,
                                    init_config.local_region));
      assert(collision_body != NULL);
      if (app_data.ple)
        example = new PhysBAM::WATER_EXAMPLE<TV>(collision_body,
                                                 PhysBAM::STREAM_TYPE((RW())),
                                                 &app_data,
                                                 app_data.ple->data(),
                                                 &job->worker_thread()->allocated_threads);
      else
        example = new PhysBAM::WATER_EXAMPLE<TV>(collision_body,
                                                 PhysBAM::STREAM_TYPE((RW())),
                                                 &app_data,
                                                 &job->worker_thread()->allocated_threads);
    } else {
      // CASE NO CACHE
      example_scope_timer =
          new application::ScopeTimer("init_example");
      StaticConfigCollisionBody* collision_body
       = dynamic_cast<StaticConfigCollisionBody*>(
          job->GetStaticConfigManager()
          ->GetStaticConfigVariable(STATIC_CONFIG_COLLISION_BODY,
                                    init_config.local_region));
      assert(collision_body != NULL);
      example = new PhysBAM::WATER_EXAMPLE<TV>(collision_body,
                                               PhysBAM::STREAM_TYPE((RW())),
                                               &job->worker_thread()->allocated_threads);
    }
    // parameters for nimbus
    example->local_region = init_config.local_region;
    example->kScale = init_config.global_region.dx();
    example->relative_region.Rebuild(1, 1, 1,
        init_config.local_region.dx(),
        init_config.local_region.dy(),
        init_config.local_region.dz());
    // physbam grid and source intiialization
    example->Initialize_Grid(
        TV_INT(init_config.local_region.dx(),
          init_config.local_region.dy(),
          init_config.local_region.dz()),
        GridToRange(init_config.global_region, init_config.local_region));
    PhysBAM::WaterSources::Add_Source(example);
    example->data_config.Set(data_config);
    delete example_scope_timer;
  }
  {
    application::ScopeTimer scope_timer("init_driver");
    driver= new PhysBAM::WATER_DRIVER<TV>(*example);
    // parameters
    driver->init_phase = init_config.init_phase;
    driver->current_frame = init_config.frame;
    driver->time = init_config.time;
    dbg(APP_LOG, "Before enter driver->Initialize.\n");
    // physbam driver initialization
    if (init_config.init_phase)
      driver->InitializeFirstDistributed(job, da);
    else
      driver->InitializeUseCachedAppData(job, da);
  }
  return true;
}

void DestroyExampleAndDriver(
    PhysBAM::WATER_EXAMPLE<TV>*& example,
    PhysBAM::WATER_DRIVER<TV>*& driver) {
  application::ScopeTimer scope_timer("delete_example_and_driver");
  if (example->create_destroy_ple) {
    example->particle_levelset_evolution.particle_levelset.Set_Thread_Queue(NULL);
    example->particle_levelset_evolution.particle_levelset.levelset.thread_queue=NULL;
    delete &example->particle_levelset_evolution;
  }
  delete example;
  example = NULL;
  delete driver;
  driver = NULL;
}

int NumParticles(PhysBAM::WATER_EXAMPLE<TV> &example) {
  typedef typename PhysBAM::PARTICLE_LEVELSET_PARTICLES<TV> ParticleBucket;
  typedef typename PhysBAM::ARRAY<ParticleBucket*, TV_INT> ParticleArray;
  nimbus::GeometricRegion lr = example.local_region;
  int gw = example.number_of_ghost_cells;
  int count = 0;

  ParticleArray *positive_particles =
    &example.particle_levelset_evolution.particle_levelset.positive_particles;
  ParticleArray *negative_particles =
    &example.particle_levelset_evolution.particle_levelset.negative_particles;
  for (int x = -gw + 1; x <= lr.dx() + gw + 1; ++x) {
    for (int y = -gw + 1; y <= lr.dy() + gw + 1; ++y) {
      for (int z = -gw + 1; z <= lr.dz() + gw + 1; ++z) {
        TV_INT bucket_index(x, y, z);
        ParticleBucket* pos_bucket = (*positive_particles)(bucket_index);
        ParticleBucket* neg_bucket = (*negative_particles)(bucket_index);
	if (pos_bucket != NULL) {
	  count += pos_bucket->Number();
	}
	if (neg_bucket != NULL) {
	  count += neg_bucket->Number();
	}
      }
    }
  }
  return(count);
}

int NumRemovedParticles(PhysBAM::WATER_EXAMPLE<TV> &example) {
  typedef typename PhysBAM::PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> RemovedParticleBucket;
  typedef typename PhysBAM::ARRAY<RemovedParticleBucket*, TV_INT> RemovedParticleArray;
  nimbus::GeometricRegion lr = example.local_region;
  int gw = example.number_of_ghost_cells;
  int count = 0;

  RemovedParticleArray *removed_positive_particles =
    &example.particle_levelset_evolution.particle_levelset.removed_positive_particles;
  RemovedParticleArray *removed_negative_particles =
    &example.particle_levelset_evolution.particle_levelset.removed_negative_particles;
  for (int x = -gw + 1; x <= lr.dx() + gw + 1; ++x) {
    for (int y = -gw + 1; y <= lr.dy() + gw + 1; ++y) {
      for (int z = -gw + 1; z <= lr.dz() + gw + 1; ++z) {
        TV_INT bucket_index(x, y, z);
        RemovedParticleBucket* pos_bucket = (*removed_positive_particles)(bucket_index);
        RemovedParticleBucket* neg_bucket = (*removed_negative_particles)(bucket_index);
	if (pos_bucket != NULL) {
	  count += pos_bucket->Number();
	}
	if (neg_bucket != NULL) {
	  count += neg_bucket->Number();
	}
      }
    }
  }
  return(count);
}

}  // namespace application

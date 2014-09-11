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
 * This file contains a loop iteration job that spawns the sub step jobs to
 * calculate the current frame. It keeps spawning the iteration in a loop as
 * long as frame computation in not complete. When the frame is done it will
 * spawn the loop frame job for the next frame. The granularity of the sub step
 * jobs could be controlled by changing the changing the GRANULARITY_STATE
 * macro in this file as follows:
 * ONE_JOB:              calculate the frame iteration in one job (like water coarse).
 * SUPER_JOBS:           break the frame iteration in to three super jobs.
 * BREAK_SUPER_JOB_1:    further break super job 1 in to its components.
 * BREAK_SUPER_JOB_2:    further break super job 2 in to its components.
 * BREAK_SUPER_JOB_3:    further break super job 3 in to its components.
 * BREAK_ALL_SUPER_JOBS: break all three super jobs in to their components.
 *
 * Author: Omid Mashayekhi <omidm@stanford.edu>
 */

#define GRANULARITY_STATE BREAK_ALL_SUPER_JOBS

#include <sys/time.h>
#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/job_loop_iteration.h"
#include "application/water_multiple/physbam_utils.h"
#include "application/water_multiple/water_driver.h"
#include "application/water_multiple/water_example.h"
#include "application/water_multiple/water_sources.h"
#include "application/water_multiple/job_names.h"
#include "application/water_multiple/data_names.h"
#include "application/water_multiple/reg_def.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"
#include "worker/job_query.h"
#include <sstream>
#include <string>

namespace application {

  JobLoopIteration::JobLoopIteration(nimbus::Application *app) {
    set_application(app);
  };

  nimbus::Job* JobLoopIteration::Clone() {
    return new JobLoopIteration(application());
  }

  void JobLoopIteration::Execute(
      nimbus::Parameter params,
      const nimbus::DataArray& da) {
    dbg(APP_LOG, "Executing loop iteration job\n");

    // Get parameters: frame, time
    InitConfig init_config;
    // Threading settings.
    init_config.use_threading = use_threading();
    init_config.core_quota = core_quota();
    std::string params_str(params.ser_data().data_ptr_raw(),
        params.ser_data().size());
    LoadParameter(params_str, &init_config.frame, &init_config.time,
                  &init_config.global_region);
    init_config.local_region = init_config.global_region;

    const int& frame = init_config.frame;
    const T& time = init_config.time;

    dbg(APP_LOG, "Frame %i and time %f in iteration job\n",
        frame, time);

    // Initialize the state of example and driver.
    PhysBAM::WATER_EXAMPLE<TV>* example =
      new PhysBAM::WATER_EXAMPLE<TV>(PhysBAM::STREAM_TYPE((RW())), false, 1);
    DataConfig data_config;
    data_config.SetFlag(DataConfig::DT);
    example->data_config.Set(data_config);
    example->Load_From_Nimbus(this, da, frame);

    *thread_queue_hook() = example->nimbus_thread_queue;

    // check whether the frame is done or not
    bool done = false;
    T target_time = example->Time_At_Frame(frame + 1);
    T dt = example->dt_buffer;
    if (time + dt >= target_time) {
      dt = target_time - time;
      done = true;
    } else {
      if (time + 2 * dt >= target_time)
        dt = .5 * (target_time - time);
    }

    if (done) {
      dbg(APP_LOG, "[CONTROL FLOW] First part, Loop done.\n");
    } else {
      dbg(APP_LOG, "[CONTROL FLOW] First part, Loop not done.\n");
    }
    dbg(APP_LOG, "[CONTROL FLOW] First part, Frame=%d, Time=%f, dt=%f\n",
        frame, time, dt);

    // done flag no more matters.
    SpawnWithBreakAllGranularity(false, init_config.frame, init_config.time,
                                 dt, da, init_config.global_region);

    // Free resources.
    example->Save_To_Nimbus(this, da, frame+1);
    *thread_queue_hook() = NULL;
    delete example;
  }

  void JobLoopIteration::SpawnWithBreakAllGranularity(
      bool done, int frame, T time, T dt, const nimbus::DataArray& da,
      const GeometricRegion& global_region) {
    nimbus::JobQuery job_query(this);
    dbg(APP_LOG, "Loop frame is spawning super job 1, 2, 3 for frame %i.\n", frame);

    int job_num = 13;
    std::vector<nimbus::job_id_t> job_ids;
    GetNewJobID(&job_ids, job_num);
    nimbus::IDSet<nimbus::logical_data_id_t> read, write;

    // because of Half Region definition this number could be either 1 or 2 for now -omidm
    int update_ghost_velocities_job_num = kAppPartNum;
    std::vector<nimbus::job_id_t> update_ghost_velocities_job_ids;
    GetNewJobID(&update_ghost_velocities_job_ids, update_ghost_velocities_job_num);

    int advect_v_job_num = kAppPartNum;
    std::vector<nimbus::job_id_t> advect_v_job_ids;
    GetNewJobID(&advect_v_job_ids, advect_v_job_num);

    int apply_forces_job_num = kAppPartNum;
    std::vector<nimbus::job_id_t> apply_forces_job_ids;
    GetNewJobID(&apply_forces_job_ids, apply_forces_job_num);

    int adjust_phi_job_num = kAppPartNum;
    std::vector<nimbus::job_id_t> adjust_phi_job_ids;
    GetNewJobID(&adjust_phi_job_ids, adjust_phi_job_num);

    int first_extrapolate_phi_job_num = kAppPartNum;
    std::vector<nimbus::job_id_t> first_extrapolate_phi_job_ids;
    GetNewJobID(&first_extrapolate_phi_job_ids, first_extrapolate_phi_job_num);

    int advect_phi_job_num = kAppPartNum;
    std::vector<nimbus::job_id_t> advect_phi_job_ids;
    GetNewJobID(&advect_phi_job_ids, advect_phi_job_num);

    int projection_job_num = 4;
    std::vector<nimbus::job_id_t> projection_job_ids;
    GetNewJobID(&projection_job_ids, projection_job_num);

    int extrapolation_job_num = kAppPartNum;
    std::vector<nimbus::job_id_t> extrapolation_job_ids;
    GetNewJobID(&extrapolation_job_ids, extrapolation_job_num);

    int extrapolate_phi_job_num = kAppPartNum;
    std::vector<nimbus::job_id_t> extrapolate_phi_job_ids;
    GetNewJobID(&extrapolate_phi_job_ids, extrapolate_phi_job_num);

    // jobs that touch particles
    size_t particle_partitions = kAppPartNum;
    // step particles
    size_t step_particles_job_num = particle_partitions;
    std::vector<nimbus::job_id_t> step_particles_job_ids;
    GetNewJobID(&step_particles_job_ids, step_particles_job_num);
    bool step_particles_single = (step_particles_job_num == 1);
    // step particles sync
    size_t step_particles_sync_job_num = (step_particles_job_num == 1)? 0 : step_particles_job_num;
    std::vector<nimbus::job_id_t> step_particles_sync_job_ids;
    GetNewJobID(&step_particles_sync_job_ids, step_particles_sync_job_num);
    // advect removed particles
    size_t advect_removed_particles_job_num = particle_partitions;
    std::vector<nimbus::job_id_t> advect_removed_particles_job_ids;
    GetNewJobID(&advect_removed_particles_job_ids, advect_removed_particles_job_num);
    bool advect_removed_particles_single = (advect_removed_particles_job_num == 1);
    // modify levelset
    size_t modify_levelset_job_num = particle_partitions;
    std::vector<nimbus::job_id_t> modify_levelset_part_one_job_ids;
    GetNewJobID(&modify_levelset_part_one_job_ids, modify_levelset_job_num);
    std::vector<nimbus::job_id_t> modify_levelset_part_two_job_ids;
    GetNewJobID(&modify_levelset_part_two_job_ids, modify_levelset_job_num);
    bool modify_levelset_single = (modify_levelset_job_num == 1);
    // delete particles
    size_t delete_particles_job_num = particle_partitions;
    std::vector<nimbus::job_id_t> delete_particles_job_ids;
    GetNewJobID(&delete_particles_job_ids, delete_particles_job_num);
    bool delete_particles_single = (delete_particles_job_num == 1);
    // reincorporate particles
    size_t reincorporate_particles_job_num = particle_partitions;
    std::vector<nimbus::job_id_t> reincorporate_particles_job_ids;
    GetNewJobID(&reincorporate_particles_job_ids, reincorporate_particles_job_num);
    bool reincorporate_particles_single = (reincorporate_particles_job_num == 1);

    struct timeval start_time;
    gettimeofday(&start_time, NULL);
    /*
     * Spawning UPDATE_GHOST_VELOCITIES stage over two workrs
     */
    for (int i = 0; i < update_ghost_velocities_job_num; ++i) {
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[i], APP_FACE_VEL_GHOST, NULL);

      nimbus::Parameter s11_params;
      std::string s11_str;
      SerializeParameter(frame, time, dt, global_region, kRegY2W3Central[i], &s11_str);
      s11_params.set_ser_data(SerializedData(s11_str));
      job_query.StageJob(UPDATE_GHOST_VELOCITIES,
          update_ghost_velocities_job_ids[i],
          read, write,
          s11_params, true);
      job_query.Hint(update_ghost_velocities_job_ids[i],
                     kRegY2W3Central[i]);
    }

    job_query.CommitStagedJobs();

    /*
     * Spawning first extrapolate phi stage over multiple workers.
     */
    for (int i = 0; i < first_extrapolate_phi_job_num; ++i) {
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegY2W8Outer[i], APP_PHI, NULL);
      LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegY2W8CentralWGB[i], APP_PHI, NULL);

      nimbus::Parameter s_extra_params;
      std::string s_extra_str;
      SerializeParameter(frame, time, dt, global_region, kRegY2W8Central[i], &s_extra_str);
      s_extra_params.set_ser_data(SerializedData(s_extra_str));
      job_query.StageJob(EXTRAPOLATE_PHI,
          first_extrapolate_phi_job_ids[i],
          read, write,
          s_extra_params, true);
      job_query.Hint(first_extrapolate_phi_job_ids[i],
                     kRegY2W8Central[i]);
    }

    job_query.CommitStagedJobs();

    /*
     * Spawning advect phi stage over multiple workers
     */
    for(int i = 0; i < advect_phi_job_num; ++i) {
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL, APP_PHI, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegY2W3Central[i], APP_FACE_VEL, APP_PHI, NULL);
      // LoadLogicalIdsInSet(this, &write, kRegY2W3Central[i], APP_PHI, NULL);

      nimbus::Parameter s12_params;
      std::string s12_str;
      SerializeParameter(frame, time, dt, global_region, kRegY2W8Central[i], &s12_str);
      s12_params.set_ser_data(SerializedData(s12_str));
      job_query.StageJob(ADVECT_PHI,
          advect_phi_job_ids[i],
          read, write,
          s12_params, true);
      job_query.Hint(advect_phi_job_ids[i], kRegY2W8Central[i]);
    }

    job_query.CommitStagedJobs();

    /*
     * Spawning step particles.
     */
    for (size_t sj = 0; sj < step_particles_job_num; sj++) {
        read.clear();
        write.clear();
        std::string step_particles_str;

        // there is just 1 last unique particle id: need to figure out how to
        // handle the case of splitting last unique particle id
        if (step_particles_single) {
            LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL_GHOST, NULL);
            LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_POS_PARTICLES,
                    APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);
            LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_POS_PARTICLES,
                    APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);
            SerializeParameter(frame,
                               time,
                               dt,
                               global_region,
                               kRegW3Central[0],
                               &step_particles_str);
        } else {
            LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[sj], APP_FACE_VEL_GHOST, NULL);
            LoadLogicalIdsInSet(this, &read, kRegY2W3CentralWGB[sj], APP_POS_PARTICLES,
                    APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);
            LoadLogicalIdsInSet(this, &write, kRegY2W3Inner[sj], APP_POS_PARTICLES,
                    APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);
            kScratchPosParticles.GetJobScratchData(this, kRegY2W3Central[sj], &write);
            kScratchNegParticles.GetJobScratchData(this, kRegY2W3Central[sj], &write);
            kScratchPosRemParticles.GetJobScratchData(this, kRegY2W3Central[sj], &write);
            kScratchNegRemParticles.GetJobScratchData(this, kRegY2W3Central[sj], &write);
            SerializeParameter(frame,
                               time,
                               dt,
                               global_region,
                               kRegY2W3Central[sj],
                               &step_particles_str);
        }

        nimbus::Parameter step_particles_params;
        step_particles_params.set_ser_data(SerializedData(step_particles_str));

        job_query.StageJob(STEP_PARTICLES,
                step_particles_job_ids[sj],
                read, write,
                step_particles_params, true);
        job_query.Hint(step_particles_job_ids[sj], kRegY2W3Central[sj]);
    }

    job_query.CommitStagedJobs();

    /*
     * Conditionally spawn synchronize jobs.
     */
    if (!step_particles_single) {
        for (size_t i = 0; i < step_particles_sync_job_num; ++i) {
            std::string sync_particles_str;
            SerializeParameter(frame,
                               time,
                               dt,
                               global_region,
                               kRegY2W3Central[i],
                               &sync_particles_str);
            nimbus::Parameter sync_particles_params;
            sync_particles_params.set_ser_data(SerializedData(sync_particles_str));

            // positive
            read.clear();
            write.clear();
            kScratchPosParticles.GetAllScratchData(this, kRegY2W3CentralWGB[i], &read);
            kScratchPosParticles.GetMainForScratchData(this,
                                                       kRegY2W3CentralWGB[i],
                                                       kRegY2W3Inner[i],
                                                       &write);
            kScratchNegParticles.GetAllScratchData(this, kRegY2W3CentralWGB[i], &read);
            kScratchNegParticles.GetMainForScratchData(this,
                                                       kRegY2W3CentralWGB[i],
                                                       kRegY2W3Inner[i],
                                                       &write);
            kScratchPosRemParticles.GetAllScratchData(this, kRegY2W3CentralWGB[i], &read);
            kScratchPosRemParticles.GetMainForScratchData(this,
                                                          kRegY2W3CentralWGB[i],
                                                          kRegY2W3Inner[i],
                                                          &write);
            kScratchNegRemParticles.GetAllScratchData(this, kRegY2W3CentralWGB[i], &read);
            kScratchNegRemParticles.GetMainForScratchData(this,
                                                          kRegY2W3CentralWGB[i],
                                                          kRegY2W3Inner[i],
                                                          &write);
            LoadLogicalIdsInSet(this, &read, kRegY2W3Inner[i],
                APP_POS_PARTICLES, APP_NEG_PARTICLES,
                APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);

            job_query.StageJob(SYNCHRONIZE_PARTICLES,
                    step_particles_sync_job_ids[i],
                    read, write,
                    sync_particles_params, true);
            job_query.Hint(step_particles_sync_job_ids[i],
                           kRegY2W3Central[i]);
        }
    }

    job_query.CommitStagedJobs();

    /* 
     * Spawning advect removed particles.
     */

    for (size_t sj = 0; sj < advect_removed_particles_job_num; sj++) {
        read.clear();
        write.clear();
        std::string advect_rem_particles_str;

        // there is just 1 last unique particle id: need to figure out how to
        // handle the case of splitting last unique particle id
        if (advect_removed_particles_single) {
            LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL_GHOST,
                    APP_PHI, NULL);
            LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_POS_REM_PARTICLES,
                APP_NEG_REM_PARTICLES, NULL);
            LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_POS_REM_PARTICLES,
                APP_NEG_REM_PARTICLES,  NULL);
            SerializeParameter(frame,
                               time,
                               dt,
                               global_region,
                               kRegW3Central[0],
                               &advect_rem_particles_str);
        } else {
            LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[sj], APP_FACE_VEL_GHOST,
                    APP_PHI, NULL);
            LoadLogicalIdsInSet(this, &read, kRegY2W3CentralWGB[sj], APP_POS_REM_PARTICLES,
                APP_NEG_REM_PARTICLES, NULL);
            LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[sj], APP_POS_REM_PARTICLES,
                APP_NEG_REM_PARTICLES,  NULL);
            SerializeParameter(frame,
                               time,
                               dt,
                               global_region,
                               kRegY2W3Central[sj],
                               &advect_rem_particles_str);
        }

        nimbus::Parameter advect_rem_particles_params;
        advect_rem_particles_params.set_ser_data(SerializedData(advect_rem_particles_str));

        job_query.StageJob(ADVECT_REMOVED_PARTICLES,
            advect_removed_particles_job_ids[sj],
            read, write,
            advect_rem_particles_params, true);
        job_query.Hint(advect_removed_particles_job_ids[sj],
                       kRegY2W3Central[sj]);
    }

    job_query.CommitStagedJobs();

    /* 
     * Spawning multiple jobs for Advect V stage
     */

    for (int i = 0; i < advect_v_job_num; ++i) {
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL_GHOST, APP_PHI, NULL);
      LoadLogicalIdsInSet(this, &read, kRegY2W3Central[i], APP_FACE_VEL, NULL);
      LoadLogicalIdsInSet(this, &read, kRegY2W1Outer[i], APP_PSI_D, APP_PSI_N, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegY2W3Central[i], APP_FACE_VEL, APP_PHI, NULL);

      nimbus::Parameter s15_params;
      std::string s15_str;
      SerializeParameter(frame, time, dt, global_region, kRegY2W3Central[i], &s15_str);
      s15_params.set_ser_data(SerializedData(s15_str));
      job_query.StageJob(ADVECT_V,
          advect_v_job_ids[i],
          read, write,
          s15_params, true);
      job_query.Hint(advect_v_job_ids[i],
                     kRegY2W3Central[i]);
    }

    job_query.CommitStagedJobs();

    {
      std::vector<nimbus::job_id_t> temp_job_ids;
      GetNewJobID(&temp_job_ids, kAppPartNum);
      /*
       * Spawning UPDATE_GHOST_VELOCITIES stage over two workrs
       */
      for (int i = 0; i < update_ghost_velocities_job_num; ++i) {
        read.clear();
        LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL, NULL);
        write.clear();
        LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[i], APP_FACE_VEL_GHOST, NULL);

        nimbus::Parameter temp_params;
        std::string temp_str;
        SerializeParameter(frame, time, dt, global_region, kRegY2W3Central[i],
                           &temp_str);
        temp_params.set_ser_data(SerializedData(temp_str));
        job_query.StageJob(UPDATE_GHOST_VELOCITIES,
                           temp_job_ids[i],
                           read, write,
                           temp_params, true);
        job_query.Hint(temp_job_ids[i],
                       kRegY2W3Central[i]);
      }
      job_query.CommitStagedJobs();
    }


    /* 
     * Spawning multiple jobs for apply forces stage
     */
    for (int i = 0; i < apply_forces_job_num; ++i) {
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL, APP_FACE_VEL_GHOST, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegY2W3Central[i], APP_FACE_VEL, NULL);

      nimbus::Parameter s16_params;
      std::string s16_str;
      SerializeParameter(frame, time, dt, global_region, kRegY2W3Central[i], &s16_str);
      s16_params.set_ser_data(SerializedData(s16_str));
      job_query.StageJob(APPLY_FORCES,
          apply_forces_job_ids[i],
          read, write,
          s16_params, true);
      job_query.Hint(apply_forces_job_ids[i],
                     kRegY2W3Central[i]);
    }

    job_query.CommitStagedJobs();


    {
      std::vector<nimbus::job_id_t> temp_job_ids;
      GetNewJobID(&temp_job_ids, kAppPartNum);
      /*
       * Spawning UPDATE_GHOST_VELOCITIES stage over two workrs
       */
      for (int i = 0; i < update_ghost_velocities_job_num; ++i) {
        read.clear();
        LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL, NULL);
        write.clear();
        LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[i], APP_FACE_VEL_GHOST, NULL);

        nimbus::Parameter temp_params;
        std::string temp_str;
        SerializeParameter(frame, time, dt, global_region, kRegY2W3Central[i],
                           &temp_str);
        temp_params.set_ser_data(SerializedData(temp_str));
        job_query.StageJob(UPDATE_GHOST_VELOCITIES,
                           temp_job_ids[i],
                           read, write,
                           temp_params, true);
        job_query.Hint(temp_job_ids[i],
                       kRegY2W3Central[i]);
      }
      job_query.CommitStagedJobs();
    }

    /* 
     * Spawning modify levelset -- part one.
     */
    for (size_t mj = 0; mj < modify_levelset_job_num; mj++) {
        read.clear();
        write.clear();
        std::string modify_levelset_part_one_str;

        // there is just 1 last unique particle id: need to figure out how to
        // handle the case of splitting last unique particle id
        if (modify_levelset_single) {
            LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL_GHOST,
                    APP_FACE_VEL, APP_PHI, NULL);
            LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_POS_PARTICLES,
                APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);
            LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_PHI, NULL);
            LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_POS_PARTICLES,
                APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);
            SerializeParameter(frame,
                               time,
                               dt,
                               global_region,
                               kRegW3Central[0],
                               &modify_levelset_part_one_str);
        } else {
            LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[mj], APP_FACE_VEL_GHOST,
                    APP_FACE_VEL, APP_PHI, NULL);
            LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[mj], APP_POS_PARTICLES,
                APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);
            LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[mj], APP_PHI, NULL);
            LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[mj], APP_POS_PARTICLES,
                APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);
            SerializeParameter(frame,
                               time,
                               dt,
                               global_region,
                               kRegY2W3Central[mj],
                               &modify_levelset_part_one_str);
        }

        nimbus::Parameter modify_levelset_params;
        modify_levelset_params.set_ser_data(SerializedData(modify_levelset_part_one_str));

        job_query.StageJob(MODIFY_LEVELSET_PART_ONE,
            modify_levelset_part_one_job_ids[mj],
            read, write,
            modify_levelset_params, true);
        job_query.Hint(modify_levelset_part_one_job_ids[mj],
                       kRegY2W3Central[mj]);
    }

    job_query.CommitStagedJobs();

    /* 
     * Spawning modify levelset -- part two.
     */
    for (size_t mj = 0; mj < modify_levelset_job_num; mj++) {
        read.clear();
        write.clear();
        std::string modify_levelset_part_two_str;

        // there is just 1 last unique particle id: need to figure out how to
        // handle the case of splitting last unique particle id
        if (modify_levelset_single) {
            LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL_GHOST,
                    APP_FACE_VEL, APP_PHI, NULL);
            LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_POS_PARTICLES,
                APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);
            LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_PHI, NULL);
            LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_POS_PARTICLES,
                APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);
            SerializeParameter(frame,
                               time,
                               dt,
                               global_region,
                               kRegW3Central[0],
                               &modify_levelset_part_two_str);
        } else {
            LoadLogicalIdsInSet(this, &read, kRegY2W7Outer[mj], APP_PHI, NULL);
            LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[mj], APP_FACE_VEL_GHOST,
                    APP_FACE_VEL, NULL);
            LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[mj], APP_POS_PARTICLES,
                APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);
            LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[mj], APP_PHI, NULL);
            LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[mj], APP_POS_PARTICLES,
                APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);
            SerializeParameter(frame,
                               time,
                               dt,
                               global_region,
                               kRegY2W3Central[mj],
                               &modify_levelset_part_two_str);
        }

        nimbus::Parameter modify_levelset_params;
        modify_levelset_params.set_ser_data(SerializedData(modify_levelset_part_two_str));

        job_query.StageJob(MODIFY_LEVELSET_PART_TWO,
            modify_levelset_part_two_job_ids[mj],
            read, write,
            modify_levelset_params, true);
        job_query.Hint(modify_levelset_part_two_job_ids[mj],
                       kRegY2W3Central[mj]);
    }

    job_query.CommitStagedJobs();

    /* 
     * Spawning adjust phi stage for multiple workers.
     */

    for (int i = 0; i < adjust_phi_job_num; ++i) {
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_PHI, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegY2W3Central[i], APP_PHI, NULL);

      nimbus::Parameter adjust_phi_params;
      std::string adjust_phi_str;
      SerializeParameter(frame, time, dt, global_region, kRegY2W3Central[i], &adjust_phi_str);
      adjust_phi_params.set_ser_data(SerializedData(adjust_phi_str));
      job_query.StageJob(ADJUST_PHI,
          adjust_phi_job_ids[i],
          read, write,
          adjust_phi_params, true);
      job_query.Hint(adjust_phi_job_ids[i],
                     kRegY2W3Central[i]);
    }

    job_query.CommitStagedJobs();

    /* 
     * Spawning delete particles.
     */

    for (size_t dj = 0; dj < delete_particles_job_num; dj++) {
        read.clear();
        write.clear();
        std::string delete_particles_str;

        // there is just 1 last unique particle id: need to figure out how to
        // handle the case of splitting last unique particle id
        if (delete_particles_single) {
            LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL_GHOST,
                    APP_PHI, NULL);
            LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_POS_PARTICLES,
                APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);
            LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_POS_PARTICLES,
                APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);
            SerializeParameter(frame,
                               time,
                               dt,
                               global_region,
                               kRegW3Central[0],
                               &delete_particles_str);
        } else {
            LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[dj], APP_FACE_VEL_GHOST,
                    APP_PHI, NULL);
            LoadLogicalIdsInSet(this, &read, kRegY2W3CentralWGB[dj], APP_POS_PARTICLES,
                APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);
            LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[dj], APP_POS_PARTICLES,
                APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);
            SerializeParameter(frame,
                               time,
                               dt,
                               global_region,
                               kRegY2W3Central[dj],
                               &delete_particles_str);
        }

        nimbus::Parameter delete_particles_params;
        delete_particles_params.set_ser_data(SerializedData(delete_particles_str));

        job_query.StageJob(DELETE_PARTICLES,
            delete_particles_job_ids[dj],
            read, write,
            delete_particles_params, true);
        job_query.Hint(delete_particles_job_ids[dj],
                       kRegY2W3Central[dj]);
    }

    job_query.CommitStagedJobs();

    /* 
     * Spawning reincorporate particles stage over entire block
     */

    for (size_t rj = 0; rj < reincorporate_particles_job_num; rj++) {
        read.clear();
        write.clear();
        std::string reincorporate_particles_str;

        // there is just 1 last unique particle id: need to figure out how to
        // handle the case of splitting last unique particle id
        if (reincorporate_particles_single) {
            LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL,
                    APP_PHI, NULL);
            LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_POS_PARTICLES,
                APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);
            LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_POS_PARTICLES,
                APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);
            SerializeParameter(frame,
                               time,
                               dt,
                               global_region,
                               kRegW3Central[0],
                               &reincorporate_particles_str);
        } else {
            LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[rj], APP_FACE_VEL,
                    APP_PHI, NULL);
            LoadLogicalIdsInSet(this, &read, kRegY2W3CentralWGB[rj], APP_POS_PARTICLES,
                APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);
            LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[rj], APP_POS_PARTICLES,
                APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);
            SerializeParameter(frame,
                               time,
                               dt,
                               global_region,
                               kRegY2W3Central[rj],
                               &reincorporate_particles_str);
        }

        nimbus::Parameter reincorporate_particles_params;
        reincorporate_particles_params.set_ser_data(SerializedData(reincorporate_particles_str));

        job_query.StageJob(REINCORPORATE_PARTICLES,
            reincorporate_particles_job_ids[rj],
            read, write,
            reincorporate_particles_params, true);
        job_query.Hint(reincorporate_particles_job_ids[rj],
                       kRegY2W3Central[rj]);
    }

    job_query.CommitStagedJobs();

    /*
    std::vector<nimbus::job_id_t> barrier_job_ids;
    GetNewJobID(&barrier_job_ids, )1;
    read.clear();
    write.clear();
    job_query.StageJob(BARRIER_JOB,
                       barrier_job_ids[0],
                       read, write,
                       params,
                       false, true);
    job_query.CommitStagedJobs();
    */

    /*
     * Loop iteration part two.
     */

    read.clear();
    write.clear();

    nimbus::Parameter projection_main_params;
    std::string projection_main_str;
    SerializeParameter(frame, time, dt, global_region, global_region,
                       &projection_main_str);
    projection_main_params.set_ser_data(
        SerializedData(projection_main_str));
    job_query.StageJob(PROJECTION_MAIN,
                    job_ids[10],
                    read, write,
                    projection_main_params,
                    false, true);
    job_query.Hint(job_ids[10], kRegW3Central[0], true);
    job_query.CommitStagedJobs();
    job_query.PrintTimeProfile();
    {
      struct timeval t;
      gettimeofday(&t, NULL);
      double time  = (static_cast<double>(t.tv_sec - start_time.tv_sec)) +
          .000001 * (static_cast<double>(t.tv_usec - start_time.tv_usec));
      dbg(APP_LOG, "\nThe query time spent in job LOOP_ITERATION is %f seconds.\n",
          time);
    }
    if (time == 0) {
      dbg(APP_LOG, "Print job dependency figure.\n");
      job_query.GenerateDotFigure("loop_iteration.dot");
    }
  }
} // namespace application

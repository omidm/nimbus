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
    PhysBAM::WATER_EXAMPLE<TV>* example;
    PhysBAM::WATER_DRIVER<TV>* driver;
    init_config.set_boundary_condition = true;
    DataConfig data_config;
    data_config.SetAll();
    InitializeExampleAndDriver(init_config, data_config,
                               this, da, example, driver);

    // check whether the frame is done or not
    bool done = false;
    T target_time = example->Time_At_Frame(driver->current_frame+1);
    T dt = example->cfl * example->incompressible.CFL(example->face_velocities);
    T temp_dt =
        example->particle_levelset_evolution.cfl_number
        * example->particle_levelset_evolution.particle_levelset.levelset.CFL(
            example->face_velocities);
    if (temp_dt < dt) {
      dt = temp_dt;
    }
    if (time + dt >= target_time) {
      dt = target_time - time;
      done = true;
    } else {
      if (time + 2*dt >= target_time)
        dt = .5 * (target_time-time);
    }

    dbg(APP_LOG, "Frame=%d, Time=%f, dt=%f\n", frame, time, dt);

    // spawn the jobs to compute the frame, depending on the
    // level of granularity we will have different sub jobs.
    switch (GRANULARITY_STATE) {
      /*
      case ONE_JOB:
        SpawnWithOneJobGranularity(done, frame, time, dt, da,
                                   init_config.global_region);
        break;
        */
      case BREAK_ALL_SUPER_JOBS:
        SpawnWithBreakAllGranularity(done, frame, time, dt, da,
                                     init_config.global_region);
        break;
      default:
        dbg(APP_LOG, "ERROR: The granularity state is not defined.");
        exit(-1);
        break;
    }

    // Free resources.
    DestroyExampleAndDriver(example, driver);
  }

  /*
  void JobLoopIteration::SpawnWithOneJobGranularity(
      bool done, int frame, T time, T dt, const nimbus::DataArray& da) {
    if (!done) {
      //Spawn one iteration of frame computation and then one loop
      //iteration job, for next iteration.
      dbg(APP_LOG, "Loop frame is spawning calculate frame job for frame %i.\n", frame);

      int job_num = 2;
      std::vector<nimbus::job_id_t> job_ids;
      GetNewJobID(&job_ids, job_num);
      nimbus::IDSet<nimbus::logical_data_id_t> read, write;
      nimbus::IDSet<nimbus::job_id_t> before, after;
      nimbus::Parameter cal_params;
      nimbus::Parameter iter_params;

      LoadReadWriteSets(this, &read, &write);

      std::string cal_str;
      SerializeParameter(frame, time, dt, &cal_str);
      cal_params.set_ser_data(SerializedData(cal_str));
      after.clear();
      after.insert(job_ids[1]);
      before.clear();
      SpawnComputeJob(CALCULATE_FRAME,
          job_ids[0],
          read, write,
          before, after,
          cal_params);

      std::string iter_str;
      SerializeParameter(frame, time + dt, &iter_str);
      iter_params.set_ser_data(SerializedData(iter_str));
      after.clear();
      before.clear();
      before.insert(job_ids[0]);
      SpawnComputeJob(LOOP_ITERATION,
          job_ids[1],
          read, write,
          before, after,
          iter_params);
    } else {
      // compute one last iteration, then spawn write frame job. Finally,
      // spawn loop frame job for next frame computation.

      int job_num = 3;
      std::vector<nimbus::job_id_t> job_ids;
      GetNewJobID(&job_ids, job_num);
      nimbus::IDSet<nimbus::logical_data_id_t> read, write;
      nimbus::IDSet<nimbus::job_id_t> before, after;
      nimbus::Parameter cal_params;
      nimbus::Parameter write_params;
      nimbus::Parameter frame_params;

      LoadReadWriteSets(this, &read, &write);

      std::string cal_str;
      SerializeParameter(frame, time, dt, &cal_str);
      cal_params.set_ser_data(SerializedData(cal_str));
      after.clear();
      after.insert(job_ids[1]);
      before.clear();
      SpawnComputeJob(CALCULATE_FRAME,
          job_ids[0],
          read, write,
          before, after,
          cal_params);

      std::string write_str;
      SerializeParameter(frame, time + dt, 0, &write_str);
      write_params.set_ser_data(SerializedData(write_str));
      after.clear();
      after.insert(job_ids[2]);
      before.clear();
      before.insert(job_ids[0]);
      SpawnComputeJob(WRITE_FRAME,
          job_ids[1],
          read, write,
          before, after,
          write_params);

      std::string frame_str;
      SerializeParameter(frame + 1, &frame_str);
      frame_params.set_ser_data(SerializedData(frame_str));
      after.clear();
      before.clear();
      before.insert(job_ids[1]);
      SpawnComputeJob(LOOP_FRAME,
          job_ids[2],
          read, write,
          before, after,
          frame_params);
    }
  }
  */

  void JobLoopIteration::SpawnWithBreakAllGranularity(
      bool done, int frame, T time, T dt, const nimbus::DataArray& da,
      const GeometricRegion& global_region) {
    dbg(APP_LOG, "Loop frame is spawning super job 1, 2, 3 for frame %i.\n", frame);

    int job_num = 13;
    std::vector<nimbus::job_id_t> job_ids;
    GetNewJobID(&job_ids, job_num);
    nimbus::IDSet<nimbus::logical_data_id_t> read, write;
    nimbus::IDSet<nimbus::job_id_t> before, after;

    // because of Half Region definition this number could be either 1 or 2 for now -omidm
    int update_ghost_velocities_job_num = 2;
    std::vector<nimbus::job_id_t> update_ghost_velocities_job_ids;
    GetNewJobID(&update_ghost_velocities_job_ids, update_ghost_velocities_job_num);

    int advect_v_job_num = 2;
    std::vector<nimbus::job_id_t> advect_v_job_ids;
    GetNewJobID(&advect_v_job_ids, advect_v_job_num);

    int apply_forces_job_num = 2;
    std::vector<nimbus::job_id_t> apply_forces_job_ids;
    GetNewJobID(&apply_forces_job_ids, apply_forces_job_num);

    int adjust_phi_job_num = 2;
    std::vector<nimbus::job_id_t> adjust_phi_job_ids;
    GetNewJobID(&adjust_phi_job_ids, adjust_phi_job_num);

    int advect_phi_job_num = 2;
    std::vector<nimbus::job_id_t> advect_phi_job_ids;
    GetNewJobID(&advect_phi_job_ids, advect_phi_job_num);

    int projection_job_num = 4;
    std::vector<nimbus::job_id_t> projection_job_ids;
    GetNewJobID(&projection_job_ids, projection_job_num);

    int extrapolation_job_num = 1;
    std::vector<nimbus::job_id_t> extrapolation_job_ids;
    GetNewJobID(&extrapolation_job_ids, extrapolation_job_num);

    int extrapolate_phi_job_num = 1;
    std::vector<nimbus::job_id_t> extrapolate_phi_job_ids;
    GetNewJobID(&extrapolate_phi_job_ids, extrapolate_phi_job_num);

    // jobs that touch particles
    size_t particle_partitions = 1;
    // step particles
    size_t step_particles_job_num = particle_partitions;
    std::vector<nimbus::job_id_t> step_particles_job_ids;
    GetNewJobID(&step_particles_job_ids, step_particles_job_num);
    bool step_particles_single = (step_particles_job_num == 1);
    // step particles sync
    size_t step_particles_sync_job_num = (step_particles_job_num == 1)? 0 : 4 * kRegY2W3Scratch_len;
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

    // Original adjust phi with objects that operates over entire block.
/*
    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL, APP_PHI, NULL);
    write.clear();
    LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_FACE_VEL_GHOST, NULL);

    nimbus::Parameter s11_params;
    std::string s11_str;
    SerializeParameter(frame, time, dt, global_region, global_region, &s11_str);
    s11_params.set_ser_data(SerializedData(s11_str));
    before.clear();
    after.clear();
    // after.insert(extrapolate_phi_job_ids[0]);
    after.insert(job_ids[1]);
    SpawnComputeJob(UPDATE_GHOST_VELOCITIES,
        update_ghost_velocities_job_ids[0],
        read, write,
        before, after,
        s11_params);
*/

    /* 
     * Spawning adjust phi with objects stage over two workrs
     */
    for (int i = 0; i < update_ghost_velocities_job_num; ++i) {
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL, APP_PHI, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[i], APP_FACE_VEL_GHOST, NULL);

      nimbus::Parameter s11_params;
      std::string s11_str;
      SerializeParameter(frame, time, dt, global_region, kRegY2W3Central[i], &s11_str);
      s11_params.set_ser_data(SerializedData(s11_str));
      before.clear();
      after.clear();
      after.insert(job_ids[1]);
      // after.insert(extrapolate_phi_job_ids[0]);
      SpawnComputeJob(UPDATE_GHOST_VELOCITIES,
          update_ghost_velocities_job_ids[i],
          read, write,
          before, after,
          s11_params, true);
    }


    // Original ADVECT_PHI job spawning.
    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL, APP_PHI, NULL);
    write.clear();
    LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_FACE_VEL, APP_PHI, NULL);

    nimbus::Parameter s12_params;
    std::string s12_str;
    SerializeParameter(frame, time, dt, global_region, global_region, &s12_str);
    s12_params.set_ser_data(SerializedData(s12_str));
    before.clear();
    for (int j = 0; j < update_ghost_velocities_job_num; ++j) {
        before.insert(update_ghost_velocities_job_ids[j]);
    }
    // before.insert(job_ids[0]);
    // before.insert(update_ghost_velocities_job_ids[0]);
    // before.insert(update_ghost_velocities_job_ids[1]);
    after.clear();
    for (size_t j = 0; j < step_particles_job_num; ++j) {
        after.insert(step_particles_job_ids[j]);
    }
    SpawnComputeJob(ADVECT_PHI,
        job_ids[1],
        read, write,
        before, after,
        s12_params, true);


    /* 
     * Spawning extrapolate phi stage over entire block
     */
/*
    read.clear();
    // TODO(quhang): read set should be the central region if it is right.
    LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_PHI, APP_FACE_VEL,
                        NULL);
    write.clear();
    LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_PHI, NULL);

    nimbus::Parameter s_extra_params;
    std::string s_extra_str;
    SerializeParameter(frame, time, dt, global_region, global_region,
                       &s_extra_str);
    s_extra_params.set_ser_data(SerializedData(s_extra_str));
    before.clear();
    before.insert(job_ids[0]);
    // before.insert(update_ghost_velocities_job_ids[0]);
    // before.insert(update_ghost_velocities_job_ids[1]);
    after.clear();
    after.insert(advect_phi_job_ids[0]);
    after.insert(advect_phi_job_ids[1]);
    SpawnComputeJob(EXTRAPOLATE_PHI,
                    extrapolate_phi_job_ids[0],
                    read, write,
                    before, after,
                    s_extra_params);

*/
    /* 
     * Start, Running ADVECT_PHI on two workers.
     */
/*
    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[0],
        APP_FACE_VEL, APP_PHI, NULL);
    write.clear();
    LoadLogicalIdsInSet(this, &write, kRegY2W3Central[0], APP_PHI, NULL);

    nimbus::Parameter s12_params;
    std::string s12_str;
    SerializeParameter(frame, time, dt, global_region,
        kRegY2W3Central[0], &s12_str);
    s12_params.set_ser_data(SerializedData(s12_str));
    before.clear();
    before.insert(extrapolate_phi_job_ids[0]);
    after.clear();
    for (int i = 0; i < step_particles_job_num; i++)
        after.insert(step_particles_job_ids[i]);
    SpawnComputeJob(ADVECT_PHI,
                    advect_phi_job_ids[0],
                    read, write,
                    before, after,
                    s12_params);

    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[1],
        APP_FACE_VEL, APP_PHI, NULL);
    write.clear();
    LoadLogicalIdsInSet(this, &write, kRegY2W3Central[1], APP_PHI, NULL);

    nimbus::Parameter s12r_params;
    std::string s12r_str;
    SerializeParameter(frame, time, dt, global_region,
        kRegY2W3Central[1], &s12r_str);
    s12r_params.set_ser_data(SerializedData(s12r_str));
    before.clear();
    before.insert(extrapolate_phi_job_ids[0]);
    after.clear();
    for (int i = 0; i < step_particles_job_num; i++)
        after.insert(step_particles_job_ids[i]);
    SpawnComputeJob(ADVECT_PHI,
                    advect_phi_job_ids[1],
                    read, write,
                    before, after,
                    s12r_params);

    // End Running ADVECT_PHI on two workers.
*/

    /* 
     * Spawning step particles.
     */
    for (size_t sj = 0; sj < step_particles_job_num; sj++) {
        read.clear();
        write.clear();
        std::string step_particles_str;

        // there is just 1 last unique particle id: need to figure out how to
        // handle the case of splitting last unique particle id
        LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_LAST_UNIQUE_PARTICLE_ID , NULL);

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

        before.clear();
        after.clear();
        // before.insert(advect_phi_job_ids[0]);
        // before.insert(advect_phi_job_ids[1]);
        before.insert(job_ids[1]);
        if (step_particles_single) {
            for (size_t i = 0; i < advect_removed_particles_job_num; i++)
                after.insert(advect_removed_particles_job_ids[i]);
        } else {
            for (size_t i = 0; i < step_particles_sync_job_num; i++)
                after.insert(step_particles_sync_job_ids[i]);
        }

        SpawnComputeJob(STEP_PARTICLES,
                step_particles_job_ids[sj],
                read, write,
                before, after,
                step_particles_params, true);
    }

    /*
     * Conditionally spawn synchronize jobs.
     */
    if (!step_particles_single) {
        nimbus::Parameter step_particles_sync_params;
        before.clear();
        after.clear();
        for (size_t i = 0; i < step_particles_job_num; i++)
            before.insert(step_particles_job_ids[i]);
        for (size_t i = 0; i < advect_removed_particles_job_num; i++)
            after.insert(advect_removed_particles_job_ids[i]);
        for (size_t i = 0; i < kRegY2W3Scratch_len; i++) {
            size_t ii = 4*i;
            // positive
            read.clear();
            write.clear();
            kScratchPosParticles.GetAllScratchData(this, kRegY2W3Scratch[i], &read);
            LoadLogicalIdsInSet(this, &write, kRegY2W3Scratch[i], APP_POS_PARTICLES, NULL);
            SpawnComputeJob(SYNCHRONIZE_PARTICLES,
                    step_particles_sync_job_ids[ii],
                    read, write,
                    before, after,
                    step_particles_sync_params, true);
            // negative
            read.clear();
            write.clear();
            kScratchNegParticles.GetAllScratchData(this, kRegY2W3Scratch[i], &read);
            LoadLogicalIdsInSet(this, &write, kRegY2W3Scratch[i], APP_NEG_PARTICLES, NULL);
            SpawnComputeJob(SYNCHRONIZE_PARTICLES,
                    step_particles_sync_job_ids[ii+1],
                    read, write,
                    before, after,
                    step_particles_sync_params, true);
            // positive removed
            read.clear();
            write.clear();
            kScratchPosRemParticles.GetAllScratchData(this, kRegY2W3Scratch[i], &read);
            LoadLogicalIdsInSet(this, &write, kRegY2W3Scratch[i], APP_POS_REM_PARTICLES, NULL);
            SpawnComputeJob(SYNCHRONIZE_PARTICLES,
                    step_particles_sync_job_ids[ii+2],
                    read, write,
                    before, after,
                    step_particles_sync_params, true);
            // negative removed
            read.clear();
            write.clear();
            kScratchNegRemParticles.GetAllScratchData(this, kRegY2W3Scratch[i], &read);
            LoadLogicalIdsInSet(this, &write, kRegY2W3Scratch[i], APP_NEG_REM_PARTICLES, NULL);
            SpawnComputeJob(SYNCHRONIZE_PARTICLES,
                    step_particles_sync_job_ids[ii+3],
                    read, write,
                    before, after,
                    step_particles_sync_params, true);
        }
    }


    /* 
     * Spawning advect removed particles.
     */

    for (size_t sj = 0; sj < advect_removed_particles_job_num; sj++) {
        read.clear();
        write.clear();
        std::string advect_rem_particles_str;

        // there is just 1 last unique particle id: need to figure out how to
        // handle the case of splitting last unique particle id
        LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_LAST_UNIQUE_PARTICLE_ID , NULL);

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
            LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[sj], APP_POS_REM_PARTICLES,
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

        before.clear();
        after.clear();
        if (step_particles_single) {
            for (size_t i = 0; i < step_particles_job_num; i++)
                before.insert(step_particles_job_ids[i]);
        } else {
            for (size_t i = 0; i < step_particles_sync_job_num; i++)
                before.insert(step_particles_sync_job_ids[i]);
        }
        // after.insert(job_ids[4]);
        for (int j = 0; j < advect_v_job_num; ++j) {
          after.insert(advect_v_job_ids[j]);
        }
        SpawnComputeJob(ADVECT_REMOVED_PARTICLES,
            advect_removed_particles_job_ids[sj],
            read, write,
            before, after,
            advect_rem_particles_params, true);
    }


/*
    // Spawning Advect V stage over entire block.

    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL, APP_FACE_VEL_GHOST, APP_PHI, NULL);
    write.clear();
    LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_FACE_VEL, APP_FACE_VEL_GHOST, APP_PHI, NULL);

   nimbus::Parameter s15_params;
    std::string s15_str;
    SerializeParameter(frame, time, dt, global_region, global_region, &s15_str);
    s15_params.set_ser_data(SerializedData(s15_str));
    before.clear();
    for (size_t j = 0; j < advect_removed_particles_job_num; j++)
      before.insert(advect_removed_particles_job_ids[j]);
    after.clear();
    after.insert(job_ids[5]);
    SpawnComputeJob(ADVECT_V,
        advect_v_job_ids[0],
        read, write,
        before, after,
        s15_params);
*/

    /* 
     * Spawning multiple jobs for Advect V stage
     */

    for (int i = 0; i < advect_v_job_num; ++i) {
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL_GHOST, APP_PHI, NULL);
      LoadLogicalIdsInSet(this, &read, kRegY2W3Central[i], APP_FACE_VEL, NULL);
      LoadLogicalIdsInSet(this, &read, kRegY2W1Outer[i], APP_PSI_D, NULL);
      LoadLogicalIdsInSet(this, &read, kRegY2W0Central[i], APP_PSI_N, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegY2W3Central[i], APP_FACE_VEL, APP_PHI, NULL);

      nimbus::Parameter s15_params;
      std::string s15_str;
      SerializeParameter(frame, time, dt, global_region, kRegY2W3Central[i], &s15_str);
      s15_params.set_ser_data(SerializedData(s15_str));
      before.clear();
      for (size_t j = 0; j < advect_removed_particles_job_num; j++)
        before.insert(advect_removed_particles_job_ids[j]);
      after.clear();
      // after.insert(job_ids[5]);
      for (int j = 0; j < apply_forces_job_num; ++j) {
        after.insert(apply_forces_job_ids[j]);
      }
      SpawnComputeJob(ADVECT_V,
          advect_v_job_ids[i],
          read, write,
          before, after,
          s15_params, true);
    }

    /* 
     * Spawning apply forces stage over entire block
     */
/*
    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL, APP_PHI, NULL);
    LoadLogicalIdsInSet(this, &read, kRegW1Outer[0], APP_PSI_D, NULL);
    LoadLogicalIdsInSet(this, &read, kRegW0Central[0], APP_PSI_N, NULL);
    write.clear();
    LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_FACE_VEL, APP_PHI, NULL);

    nimbus::Parameter s16_params;
    std::string s16_str;
    SerializeParameter(frame, time, dt, global_region, global_region, &s16_str);
    s16_params.set_ser_data(SerializedData(s16_str));
    before.clear();
    // before.insert(job_ids[4]);
    for (int j = 0; j < advect_v_job_num; ++j) {
      before.insert(advect_v_job_ids[j]);
    }
    after.clear();
    for (size_t j = 0; j < modify_levelset_job_num; j++)
        after.insert(modify_levelset_part_one_job_ids[j]);
    SpawnComputeJob(APPLY_FORCES,
        apply_forces_job_ids[0],
        read, write,
        before, after,
        s16_params);
*/

    /* 
     * Spawning multiple jobs for apply forces stage
     */
    for (int i = 0; i < apply_forces_job_num; ++i) {
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL, APP_FACE_VEL_GHOST, APP_PHI, NULL);
      LoadLogicalIdsInSet(this, &read, kRegY2W1Outer[i], APP_PSI_D, NULL);
      LoadLogicalIdsInSet(this, &read, kRegY2W0Central[i], APP_PSI_N, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegY2W3Central[i], APP_FACE_VEL, APP_PHI, NULL);

      nimbus::Parameter s16_params;
      std::string s16_str;
      SerializeParameter(frame, time, dt, global_region, kRegY2W3Central[i], &s16_str);
      s16_params.set_ser_data(SerializedData(s16_str));
      before.clear();
      // before.insert(job_ids[4]);
      for (int j = 0; j < advect_v_job_num; ++j) {
        before.insert(advect_v_job_ids[j]);
      }
      after.clear();
      for (size_t j = 0; j < modify_levelset_job_num; j++)
        after.insert(modify_levelset_part_one_job_ids[j]);
      SpawnComputeJob(APPLY_FORCES,
          apply_forces_job_ids[i],
          read, write,
          before, after,
          s16_params, true);
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
        LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_LAST_UNIQUE_PARTICLE_ID , NULL);

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

        before.clear();
        after.clear();
        // before.insert(job_ids[5]);
        for (int j = 0; j < apply_forces_job_num; ++j) {
          before.insert(apply_forces_job_ids[j]);
        }
        // after.insert(job_ids[7]);
        for (size_t j = 0; j < modify_levelset_job_num; ++j) {
          after.insert(modify_levelset_part_two_job_ids[j]);
        }
        SpawnComputeJob(MODIFY_LEVELSET_PART_ONE,
            modify_levelset_part_one_job_ids[mj],
            read, write,
            before, after,
            modify_levelset_params, true);
    }

    /* 
     * Spawning modify levelset -- part two.
     */
    for (size_t mj = 0; mj < modify_levelset_job_num; mj++) {
        read.clear();
        write.clear();
        std::string modify_levelset_part_two_str;

        // there is just 1 last unique particle id: need to figure out how to
        // handle the case of splitting last unique particle id
        LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_LAST_UNIQUE_PARTICLE_ID , NULL);

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

        before.clear();
        after.clear();
        // before.insert(job_ids[5]);
        for (size_t j = 0; j < modify_levelset_job_num; ++j) {
          before.insert(modify_levelset_part_one_job_ids[j]);
        }
        // after.insert(job_ids[7]);
        for (int j = 0; j < adjust_phi_job_num; ++j) {
          after.insert(adjust_phi_job_ids[j]);
        }
        SpawnComputeJob(MODIFY_LEVELSET_PART_TWO,
            modify_levelset_part_two_job_ids[mj],
            read, write,
            before, after,
            modify_levelset_params, true);
    }

    /* 
     * Spawning adjust phi stage over entire block
     */
/*
    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_PHI, NULL);
    write.clear();
    LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_PHI, NULL);

    nimbus::Parameter adjust_phi_params;
    std::string adjust_phi_str;
    SerializeParameter(frame, time, dt, global_region, global_region, &adjust_phi_str);
    adjust_phi_params.set_ser_data(SerializedData(adjust_phi_str));
    after.clear();
    after.insert(job_ids[8]);
    before.clear();
    for (size_t j = 0; j < modify_levelset_job_num; j++)
        before.insert(modify_levelset_part_two_job_ids[j]);
    SpawnComputeJob(ADJUST_PHI,
        adjust_phi_job_ids[0],
        read, write,
        before, after,
        adjust_phi_params);
*/

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
      after.clear();
      for (size_t dj = 0; dj < delete_particles_job_num; dj++)
          after.insert(delete_particles_job_ids[dj]);
      before.clear();
      for (size_t j = 0; j < modify_levelset_job_num; j++)
        before.insert(modify_levelset_part_two_job_ids[j]);
      SpawnComputeJob(ADJUST_PHI,
          adjust_phi_job_ids[i],
          read, write,
          before, after,
          adjust_phi_params, true);
    }


    /* 
     * Spawning delete particles.
     */

    for (size_t dj = 0; dj < delete_particles_job_num; dj++) {
        read.clear();
        write.clear();
        std::string delete_particles_str;

        // there is just 1 last unique particle id: need to figure out how to
        // handle the case of splitting last unique particle id
        LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_LAST_UNIQUE_PARTICLE_ID , NULL);
        LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_LAST_UNIQUE_PARTICLE_ID , NULL);

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
            LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[dj], APP_POS_PARTICLES,
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

        before.clear();
        for (int j = 0; j < adjust_phi_job_num; ++j) 
            before.insert(adjust_phi_job_ids[j]);
        after.clear();
        for (size_t rj = 0; rj < reincorporate_particles_job_num; rj++)
            after.insert(reincorporate_particles_job_ids[rj]);
        // before.insert(job_ids[7]);
        SpawnComputeJob(DELETE_PARTICLES,
            delete_particles_job_ids[dj],
            read, write,
            before, after,
            delete_particles_params, true);
    }


    /* 
     * Spawning reincorporate particles stage over entire block
     */

    for (size_t rj = 0; rj < reincorporate_particles_job_num; rj++) {
        read.clear();
        write.clear();
        std::string reincorporate_particles_str;

        // there is just 1 last unique particle id: need to figure out how to
        // handle the case of splitting last unique particle id
        LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_LAST_UNIQUE_PARTICLE_ID , NULL);
        LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_LAST_UNIQUE_PARTICLE_ID , NULL);

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
            LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[rj], APP_POS_PARTICLES,
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

        before.clear();
        after.clear();
        after.insert(projection_job_ids[0]);
        for (int j = 0; j < adjust_phi_job_num; ++j) {
          before.insert(delete_particles_job_ids[j]);
        }
        // before.insert(job_ids[7]);
        SpawnComputeJob(REINCORPORATE_PARTICLES,
            reincorporate_particles_job_ids[rj],
            read, write,
            before, after,
            reincorporate_particles_params, true);
    }


    /*
     * TAG_PROJECTION
     * std::vector<nimbus::job_id_t> projection_job_ids;
     * Spawning projection stage over entire block
     */

    {
      // Read velocity, pressure, levelset.
      // Write velocity, pressure, psi_D, psi_N, filled_region_colors,
      //     divergence.
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL, APP_PHI, NULL);
      LoadLogicalIdsInSet(this, &read, kRegW1Outer[0],
                          APP_DIVERGENCE, APP_PSI_D, APP_FILLED_REGION_COLORS,
                          APP_PRESSURE, NULL);
      LoadLogicalIdsInSet(this, &read, kRegW1Central[0], APP_PSI_N,
                          APP_U_INTERFACE, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_FACE_VEL, APP_PHI, NULL);
      LoadLogicalIdsInSet(this, &write, kRegW1Outer[0],
                          APP_DIVERGENCE, APP_PSI_D, APP_FILLED_REGION_COLORS,
                          APP_PRESSURE, NULL);
      LoadLogicalIdsInSet(this, &write, kRegW1Central[0], APP_PSI_N,
                          APP_U_INTERFACE, NULL);

      nimbus::Parameter projection_calculate_boundary_condition_params;
      std::string projection_calculate_boundary_condition_str;
      SerializeParameter(frame, time, dt, global_region, global_region,
                         &projection_calculate_boundary_condition_str);
      projection_calculate_boundary_condition_params.set_ser_data(
          SerializedData(projection_calculate_boundary_condition_str));

      before.clear();
      for (size_t rj = 0; rj < reincorporate_particles_job_num; rj++)
        before.insert(reincorporate_particles_job_ids[rj]);
      after.clear();
      after.insert(projection_job_ids[1]);

      SpawnComputeJob(PROJECTION_CALCULATE_BOUNDARY_CONDITION,
                      projection_job_ids[0],
                      read, write,
                      before, after,
                      projection_calculate_boundary_condition_params, true);

      // Read psi_D, psi_N, filled_region_colors, divergence, pressure.
      // Write pressure.
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL, APP_PHI, NULL);
      LoadLogicalIdsInSet(this, &read, kRegW1Outer[0],
                          APP_DIVERGENCE, APP_PSI_D, APP_FILLED_REGION_COLORS,
                          APP_PRESSURE, NULL);
      LoadLogicalIdsInSet(this, &read, kRegW1Central[0], APP_PSI_N,
                          APP_U_INTERFACE, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_FACE_VEL, APP_PHI, NULL);
      LoadLogicalIdsInSet(this, &write, kRegW1Outer[0],
                          APP_DIVERGENCE, APP_PSI_D, APP_FILLED_REGION_COLORS,
                          APP_PRESSURE, NULL);
      LoadLogicalIdsInSet(this, &write, kRegW1Central[0], APP_PSI_N,
                          APP_U_INTERFACE, APP_MATRIX_A, APP_VECTOR_X,
                          APP_VECTOR_B, APP_PROJECTION_LOCAL_TOLERANCE,
                          APP_INDEX_M2C, APP_INDEX_C2M, NULL);

      nimbus::Parameter projection_construct_matrix_params;
      std::string projection_construct_matrix_str;
      SerializeParameter(frame, time, dt, global_region, global_region,
                         &projection_construct_matrix_str);
      projection_construct_matrix_params.set_ser_data(
          SerializedData(projection_construct_matrix_str));

      before.clear();
      before.insert(projection_job_ids[0]);
      after.clear();
      after.insert(projection_job_ids[2]);
      SpawnComputeJob(PROJECTION_CONSTRUCT_MATRIX,
                      projection_job_ids[1],
                      read, write,
                      before, after,
                      projection_construct_matrix_params, true);

      // Read psi_D, psi_N, filled_region_colors, divergence, pressure.
      // Write pressure.
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL, APP_PHI, NULL);
      LoadLogicalIdsInSet(this, &read, kRegW1Outer[0],
                          APP_DIVERGENCE, APP_PSI_D, APP_FILLED_REGION_COLORS,
                          APP_PRESSURE, NULL);
      LoadLogicalIdsInSet(this, &read, kRegW1Central[0], APP_PSI_N,
                          APP_U_INTERFACE, APP_MATRIX_A, APP_VECTOR_X,
                          APP_VECTOR_B, APP_PROJECTION_LOCAL_TOLERANCE,
                          APP_INDEX_M2C, APP_INDEX_C2M, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_FACE_VEL, APP_PHI, NULL);
      LoadLogicalIdsInSet(this, &write, kRegW1Outer[0],
                          APP_DIVERGENCE, APP_PSI_D, APP_FILLED_REGION_COLORS,
                          APP_PRESSURE, NULL);
      LoadLogicalIdsInSet(this, &write, kRegW1Central[0], APP_PSI_N,
                          APP_U_INTERFACE, APP_VECTOR_X, NULL);

      nimbus::Parameter projection_core_params;
      std::string projection_core_str;
      SerializeParameter(frame, time, dt, global_region, global_region,
                         &projection_core_str);
      projection_core_params.set_ser_data(SerializedData(projection_core_str));

      before.clear();
      before.insert(projection_job_ids[1]);
      after.clear();
      after.insert(projection_job_ids[3]);
      SpawnComputeJob(PROJECTION_CORE,
                      projection_job_ids[2],
                      read, write,
                      before, after,
                      projection_core_params);


      // ? u_interface
      // Read pressure, levelset, psi_D, psi_N, velocity.
      // Write velocity, write pressure.
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL, APP_PHI, NULL);
      LoadLogicalIdsInSet(this, &read, kRegW1Outer[0],
                          APP_DIVERGENCE, APP_PSI_D, APP_FILLED_REGION_COLORS,
                          APP_PRESSURE, NULL);
      LoadLogicalIdsInSet(this, &read, kRegW1Central[0], APP_PSI_N,
                          APP_U_INTERFACE, APP_INDEX_M2C, APP_VECTOR_X, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_FACE_VEL, APP_PHI, NULL);
      LoadLogicalIdsInSet(this, &write, kRegW1Outer[0],
                          APP_DIVERGENCE, APP_PSI_D, APP_FILLED_REGION_COLORS,
                          APP_PRESSURE, NULL);
      LoadLogicalIdsInSet(this, &write, kRegW1Central[0], APP_PSI_N,
                          APP_U_INTERFACE, NULL);

      nimbus::Parameter projection_wrapup_params;
      std::string projection_wrapup_str;
      SerializeParameter(frame, time, dt, global_region, global_region,
                         &projection_wrapup_str);
      projection_wrapup_params.set_ser_data(SerializedData(projection_wrapup_str));

      before.clear();
      before.insert(projection_job_ids[2]);
      after.clear();
      SpawnComputeJob(PROJECTION_WRAPUP,
                      projection_job_ids[3],
                      read, write,
                      before, after,
                      projection_wrapup_params, true);
    }

    /* 
     * Spawning extrapolate phi stage over entire block
     */
    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_PHI, APP_FACE_VEL,
                        NULL);
    write.clear();
    LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_PHI, NULL);

    nimbus::Parameter s_extra_params;
    std::string s_extra_str;
    SerializeParameter(frame, time, dt, global_region, global_region,
                       &s_extra_str);
    s_extra_params.set_ser_data(SerializedData(s_extra_str));
    before.clear();
    before.insert(projection_job_ids[3]);
    after.clear();
    SpawnComputeJob(EXTRAPOLATE_PHI,
                    extrapolate_phi_job_ids[0],
                    read, write,
                    before, after,
                    s_extra_params);


    /* 
     * Spawning extrapolate phi stage over multiple workers.
     */
/*    
    for (int i = 0; i < extrapolate_phi_job_num; ++i) {
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegY2W8Outer[i], APP_PHI, APP_FACE_VEL, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegY2W8CentralWGB[i], APP_PHI, NULL);

      nimbus::Parameter s_extra_params;
      std::string s_extra_str;
      SerializeParameter(frame, time, dt, global_region, kRegY2W8Central[i], &s_extra_str);
      s_extra_params.set_ser_data(SerializedData(s_extra_str));
      before.clear();
      before.insert(projection_job_ids[3]);
      after.clear();
      SpawnComputeJob(EXTRAPOLATE_PHI,
          extrapolate_phi_job_ids[i],
          read, write,
          before, after,
          s_extra_params);
    }
*/

    /*
     * Spawning extrapolation stage over entire block
     */
    {
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL, APP_PHI, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_FACE_VEL, APP_PHI, NULL);

      int index = 11;
      nimbus::Parameter extrapolation_params;
      std::string extrapolation_str;
      SerializeParameter(frame, time, dt, global_region, global_region, &extrapolation_str);
      extrapolation_params.set_ser_data(SerializedData(extrapolation_str));
      after.clear();
      after.insert(job_ids[index+1]);
      before.clear();
      for (int j = 0; j < extrapolate_phi_job_num; ++j) {
        before.insert(extrapolate_phi_job_ids[j]);
      }
      SpawnComputeJob(EXTRAPOLATION,
                      extrapolation_job_ids[0],
                      read, write,
                      before, after,
                      reincorporate_particles_params, true);
    }

    /*
     * Spawning extrapolation stage over multiple workers
     */
/*
    for (int i = 0; i < extrapolation_job_num; ++i) {
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL, APP_PHI, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegY2W3Central[i], APP_FACE_VEL, APP_PHI, NULL);

      int index = 11;
      nimbus::Parameter extrapolation_params;
      std::string extrapolation_str;
      SerializeParameter(frame, time, dt, global_region, kRegY2W3Central[i], &extrapolation_str);
      extrapolation_params.set_ser_data(SerializedData(extrapolation_str));
      after.clear();
      after.insert(job_ids[index+1]);
      before.clear();
      before.insert(projection_job_ids[3]);
      SpawnComputeJob(EXTRAPOLATION,
                      extrapolation_job_ids[i],
                      read, write,
                      before, after,
                      extrapolation_params, true);
    }

*/


    if (!done) {

    /* 
     * Spawning loop iteration fopr next iteration.
     */

      {
        read.clear();
        LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL, APP_FACE_VEL_GHOST, APP_PHI, NULL);
        LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_POS_PARTICLES,
            APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES,
            APP_LAST_UNIQUE_PARTICLE_ID , NULL);
        write.clear();
        LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_FACE_VEL, APP_FACE_VEL_GHOST, APP_PHI, NULL);
        LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_POS_PARTICLES,
            APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES,
            APP_LAST_UNIQUE_PARTICLE_ID , NULL);

        int index = 12;
        nimbus::Parameter iter_params;
        std::string iter_str;
        SerializeParameter(frame, time + dt, global_region, &iter_str);
        iter_params.set_ser_data(SerializedData(iter_str));
        after.clear();
        before.clear();
        for (int j = 0; j < extrapolation_job_num; ++j) {
          before.insert(extrapolation_job_ids[j]);
        }
        //before.insert(job_ids[index-1]);
        SpawnComputeJob(LOOP_ITERATION,
                        job_ids[index],
                        read, write,
                        before, after,
                        iter_params);
      }
    } else {

      std::vector<nimbus::job_id_t> loop_job_id;
      GetNewJobID(&loop_job_id, 1);

    /* 
     * Spawning write frame over entire block
     */

      {
        read.clear();
        LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL, APP_FACE_VEL_GHOST, APP_PHI, NULL);
        LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_POS_PARTICLES,
            APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES,
            APP_LAST_UNIQUE_PARTICLE_ID , NULL);
        write.clear();
        LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_FACE_VEL, APP_FACE_VEL_GHOST, APP_PHI, NULL);
        LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_POS_PARTICLES,
            APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES,
            APP_LAST_UNIQUE_PARTICLE_ID , NULL);

        int index = 12;
        nimbus::Parameter write_params;
        std::string write_str;
        SerializeParameter(frame, time + dt, 0,
                           global_region, global_region, &write_str);
        write_params.set_ser_data(SerializedData(write_str));
        after.clear();
        after.insert(loop_job_id[0]);
        before.clear();
        for (int j = 0; j < extrapolation_job_num; ++j) {
          before.insert(extrapolation_job_ids[j]);
        }
        // before.insert(job_ids[index-1]);
        SpawnComputeJob(WRITE_FRAME,
                        job_ids[index],
                        read, write,
                        before, after,
                        write_params, true);
      }

    /* 
     * Spawning loop frame to compute next frame.
     */

      {
        read.clear();
        write.clear();

        int index = 13;
        nimbus::Parameter frame_params;
        std::string frame_str;
        SerializeParameter(frame + 1, global_region, &frame_str);
        frame_params.set_ser_data(SerializedData(frame_str));
        after.clear();
        before.clear();
        before.insert(job_ids[index-1]);
        SpawnComputeJob(LOOP_FRAME,
                        loop_job_id[0],
                        read, write,
                        before, after,
                        frame_params);
      }
    }
  }


} // namespace application

/* Copyright 2013 Stanford University.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 vd* are met:
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
 * This file contains the "main" job that Nimbus launches after loading an
 * application. All subsequent jobs are spawned from here.
 *
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/data_def.h"
#include "application/water_multiple/data_names.h"
#include "application/water_multiple/job_main.h"
#include "application/water_multiple/job_names.h"
#include "application/water_multiple/reg_def.h"
#include "data/scratch_data_helper.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"
#include "worker/job_query.h"
#include <vector>

namespace application {

JobMain::JobMain(nimbus::Application *app) {
  set_application(app);
};

nimbus::Job* JobMain::Clone() {
  return new JobMain(application());
}

void JobMain::Execute(nimbus::Parameter params, const nimbus::DataArray& da) {
  nimbus::JobQuery job_query(this);
  dbg(APP_LOG, "Executing main job\n");
  DefineNimbusData(this);

  // Job setup
  int init_job_num = kAppPartNum;
  std::vector<nimbus::job_id_t> init_job_ids;
  GetNewJobID(&init_job_ids, init_job_num);

  int make_signed_distance_job_num = kAppPartNum;
  std::vector<nimbus::job_id_t> make_signed_distance_job_ids;
  GetNewJobID(&make_signed_distance_job_ids, make_signed_distance_job_num);

  int extrapolate_phi_job_num = kAppPartNum;
  std::vector<nimbus::job_id_t> extrapolate_phi_job_ids;
  GetNewJobID(&extrapolate_phi_job_ids, extrapolate_phi_job_num);

  int extrapolate_phi_2_job_num = kAppPartNum;
  std::vector<nimbus::job_id_t> extrapolate_phi_2_job_ids;
  GetNewJobID(&extrapolate_phi_2_job_ids, extrapolate_phi_2_job_num);

  int extrapolate_v_job_num = kAppPartNum;
  std::vector<nimbus::job_id_t> extrapolate_v_job_ids;
  GetNewJobID(&extrapolate_v_job_ids, extrapolate_v_job_num);

  int reseed_particles_job_num = kAppPartNum;
  std::vector<nimbus::job_id_t> reseed_particles_job_ids;
  GetNewJobID(&reseed_particles_job_ids, reseed_particles_job_num);

  int write_output_job_num = kAppPartNum;
  std::vector<nimbus::job_id_t> write_output_job_ids;
  GetNewJobID(&write_output_job_ids, write_output_job_num);

  int loop_frame_job_num = 1;
  std::vector<nimbus::job_id_t> loop_frame_job_ids;
  GetNewJobID(&loop_frame_job_ids, loop_frame_job_num);


  nimbus::IDSet<nimbus::logical_data_id_t> read, write;

  int frame = 0;
  T time = 0;
  T dt = 0;

  /*
   * Spawning initialize stage over multiple workers
   */
  for (int i = 0; i < init_job_num; ++i) {
    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL, APP_FACE_VEL_GHOST, APP_PHI, NULL);
    LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_POS_PARTICLES,
                        APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES,
                        APP_LAST_UNIQUE_PARTICLE_ID , NULL);
    LoadLogicalIdsInSet(this, &read, kRegY2W1Outer[i], APP_PSI_D, APP_PSI_N, NULL);

    write.clear();
    LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[i], APP_FACE_VEL, APP_FACE_VEL_GHOST, APP_PHI, NULL);
    LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[i], APP_POS_PARTICLES,
                        APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES,
                        APP_LAST_UNIQUE_PARTICLE_ID , NULL);
    LoadLogicalIdsInSet(this, &write, kRegY2W1CentralWGB[i],
                        APP_PRESSURE,APP_PSI_D, APP_PSI_N,  NULL);

    nimbus::Parameter init_params;
    std::string init_str;
    SerializeParameter(frame, time, dt, kDefaultRegion, kRegY2W3Central[i], &init_str);
    init_params.set_ser_data(SerializedData(init_str));
    job_query.StageJob(INITIALIZE,
                       init_job_ids[i],
                       read, write,
                       init_params, true);
  }
  job_query.CommitStagedJobs();

  /*
   * Spawning extrapolate phi stage over multiple workers.
   */
  for (int i = 0; i < extrapolate_phi_job_num; ++i) {
    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegY2W8Outer[i], APP_PHI,
                        APP_FACE_VEL, NULL);

    write.clear();
    LoadLogicalIdsInSet(this, &write,
                        kRegY2W8CentralWGB[i], APP_PHI, NULL);

    nimbus::Parameter phi_params;
    std::string phi_str;
    SerializeParameter(frame, time, dt, kDefaultRegion, kRegY2W8Central[i], &phi_str);
    phi_params.set_ser_data(SerializedData(phi_str));

    job_query.StageJob(EXTRAPOLATE_PHI,
                       extrapolate_phi_job_ids[i],
                       read, write,
                       phi_params, true);
  }
  job_query.CommitStagedJobs();


  /*
   * Spawning make signed distance.
   */
  for (int i = 0; i < make_signed_distance_job_num; ++i) {
    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegY2W7Outer[i], APP_PHI, NULL);
    LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL_GHOST,
                        APP_FACE_VEL, NULL);
    LoadLogicalIdsInSet(this, &read, kRegY2W1Outer[i], APP_PSI_D, APP_PSI_N, NULL);

    write.clear();
    LoadLogicalIdsInSet(this, &write, kRegY2W7CentralWGB[i], APP_PHI, NULL);
    LoadLogicalIdsInSet(this, &write, kRegY2W1CentralWGB[i], APP_PSI_D, APP_PSI_N,  NULL);

    std::string make_signed_distance_str;
    SerializeParameter(frame, time, dt, kDefaultRegion, kRegY2W3Central[i], &make_signed_distance_str);
    nimbus::Parameter make_signed_distance_params;
    make_signed_distance_params.set_ser_data(SerializedData(make_signed_distance_str));

    job_query.StageJob(MAKE_SIGNED_DISTANCE,
                       make_signed_distance_job_ids[i],
                       read, write,
                       make_signed_distance_params, true);
  }
  job_query.CommitStagedJobs();



  for (int i = 0; i < reseed_particles_job_num; ++i) {
    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL,
                        APP_PHI, NULL);
    LoadLogicalIdsInSet(this, &read, kRegY2W1Outer[i], APP_PSI_D,
                        APP_PSI_N, NULL);
    LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_POS_PARTICLES,
                        APP_NEG_PARTICLES, APP_POS_REM_PARTICLES,
                        APP_NEG_REM_PARTICLES, APP_LAST_UNIQUE_PARTICLE_ID,
                        NULL);
    write.clear();
    LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[i], APP_POS_PARTICLES,
                        APP_NEG_PARTICLES, APP_POS_REM_PARTICLES,
                        APP_NEG_REM_PARTICLES, APP_LAST_UNIQUE_PARTICLE_ID,
                        NULL);

    nimbus::Parameter temp_params;
    std::string temp_str;
    SerializeParameter(frame, time, dt, kDefaultRegion, kRegY2W3Central[i],
                       &temp_str);
    temp_params.set_ser_data(SerializedData(temp_str));
    job_query.StageJob(RESEED_PARTICLES,
                       reseed_particles_job_ids[i],
                       read, write,
                       temp_params, true);
  }
  job_query.CommitStagedJobs();

  /*
   * Spawning extrapolate phi stage over multiple workers.
   */
  for (int i = 0; i < extrapolate_phi_2_job_num; ++i) {
    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegY2W8Outer[i], APP_PHI,
                        APP_FACE_VEL, NULL);

    write.clear();
    LoadLogicalIdsInSet(this, &write,
                        kRegY2W8CentralWGB[i], APP_PHI, NULL);

    nimbus::Parameter phi_params;
    std::string phi_str;
    SerializeParameter(frame, time, dt, kDefaultRegion, kRegY2W8Central[i], &phi_str);
    phi_params.set_ser_data(SerializedData(phi_str));

    job_query.StageJob(EXTRAPOLATE_PHI,
                       extrapolate_phi_2_job_ids[i],
                       read, write,
                       phi_params, true);
  }
  job_query.CommitStagedJobs();


  /*
   * Spawning extrapolate v stage over multiple workers
   */
  for (int i = 0; i < extrapolate_v_job_num; ++i) {
    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegY2W8Outer[i],
                        APP_FACE_VEL, APP_PHI, NULL);

    write.clear();
    LoadLogicalIdsInSet(this, &write, kRegY2W8Central[i],
                        APP_FACE_VEL, NULL);

    nimbus::Parameter v_params;
    std::string v_str;
    SerializeParameter(frame, time, dt, kDefaultRegion, kRegY2W3Central[i], &v_str);
    v_params.set_ser_data(SerializedData(v_str));

    job_query.StageJob(EXTRAPOLATION,
                       extrapolate_v_job_ids[i],
                       read, write,
                       v_params, true);
  }
  job_query.CommitStagedJobs();


  if (kUseGlobalWrite) {
    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL,
                        APP_PHI, NULL);
    LoadLogicalIdsInSet(this, &read, kRegW1Outer[0], APP_PSI_D,
                        APP_PSI_N, NULL);
    LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_POS_PARTICLES,
                        APP_NEG_PARTICLES, APP_POS_REM_PARTICLES,
                        APP_NEG_REM_PARTICLES, APP_LAST_UNIQUE_PARTICLE_ID,
                        NULL);
    write.clear();

    nimbus::Parameter temp_params;
    std::string temp_str;
    SerializeParameter(frame - 1, time + dt, 0, -1, kDefaultRegion, kDefaultRegion,
                       &temp_str);
    temp_params.set_ser_data(SerializedData(temp_str));
    job_query.StageJob(WRITE_OUTPUT,
                       write_output_job_ids[0],
                       read, write,
                       temp_params, true);
    job_query.CommitStagedJobs();
  } else {
    for (int i = 0; i < write_output_job_num; ++i) {
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[0], APP_FACE_VEL,
                          APP_PHI, NULL);
      LoadLogicalIdsInSet(this, &read, kRegY2W1Outer[0], APP_PSI_D,
                          APP_PSI_N, NULL);
      LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[0], APP_POS_PARTICLES,
                          APP_NEG_PARTICLES, APP_POS_REM_PARTICLES,
                          APP_NEG_REM_PARTICLES, APP_LAST_UNIQUE_PARTICLE_ID,
                          NULL);
      write.clear();

      nimbus::Parameter temp_params;
      std::string temp_str;
      SerializeParameter(frame - 1, time + dt, 0, i, kDefaultRegion, kRegY2W3Central[i],
                         &temp_str);
      temp_params.set_ser_data(SerializedData(temp_str));
      job_query.StageJob(WRITE_OUTPUT,
                         write_output_job_ids[i],
                         read, write,
                         temp_params, true);
    }
    job_query.CommitStagedJobs();
  }

  /*
   * Spawning loop frame job.
   */
  read.clear();
  write.clear();

  nimbus::Parameter loop_params;
  std::string loop_str;
  SerializeParameter(frame, kDefaultRegion, &loop_str);
  loop_params.set_ser_data(SerializedData(loop_str));

  job_query.StageJob(LOOP_FRAME,
                     loop_frame_job_ids[0],
                     read, write,
                     loop_params, false, true);
  job_query.CommitStagedJobs();

  dbg(APP_LOG, "Completed executing main job\n");

  dbg(APP_LOG, "Print job dependency figure.\n");
  job_query.GenerateDotFigure("job_main.dot");
}

} // namespace application

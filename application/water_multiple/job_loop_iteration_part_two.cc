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
 * Author: Hang Qu <quhang@stanford.edu>
 */


#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/physbam_utils.h"
#include "application/water_multiple/job_names.h"
#include "application/water_multiple/data_names.h"
#include "application/water_multiple/reg_def.h"
#include "application/water_multiple/water_example.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"
#include "worker/job_query.h"
#include "worker/worker_thread.h"
#include <sstream>
#include <string>

#include "application/water_multiple/job_loop_iteration_part_two.h"

namespace application {

JobLoopIterationPartTwo::JobLoopIterationPartTwo(nimbus::Application *app) {
  set_application(app);
};

nimbus::Job* JobLoopIterationPartTwo::Clone() {
  return new JobLoopIterationPartTwo(application());
}

void JobLoopIterationPartTwo::Execute(
    nimbus::Parameter params,
    const nimbus::DataArray& da) {
  dbg(APP_LOG, "Executing LOOP_ITERATION_PART_TWO job\n");

  InitConfig init_config;
  // Threading settings.
  std::string params_str(params.ser_data().data_ptr_raw(),
                         params.ser_data().size());
  LoadParameter(params_str, &init_config);
  T dt = init_config.dt;
  dbg(APP_LOG, "Frame %i in LOOP_ITERATION_PART_TWO job\n", init_config.frame);

  const int& frame = init_config.frame;
  const T& time = init_config.time;

  dbg(APP_LOG, "Frame %i and time %f in iteration job\n",
      frame, time);

  // Initialize the state of example and driver.
  PhysBAM::WATER_EXAMPLE<TV>* example =
      new PhysBAM::WATER_EXAMPLE<TV>(NULL,
                                     PhysBAM::STREAM_TYPE((RW())),
                                     &worker_thread()->allocated_threads);

  // check whether the frame is done or not
  bool done = false;
  if (time + dt >= example->Time_At_Frame(frame + 1)) {
    done = true;
  }

  // done = true;

  // This job does not use application app_data.
  delete example;

  SpawnJobs(
      done, frame, time, dt, da, init_config.global_region);

}

void JobLoopIterationPartTwo::SpawnJobs(
    bool done, int frame, T time, T dt, const nimbus::DataArray& da,
    const nimbus::GeometricRegion& global_region) {
  struct timeval start_time;
  gettimeofday(&start_time, NULL);

  // nimbus::JobQuery job_query(this);

  int job_num = 2;
  std::vector<nimbus::job_id_t> job_ids;
  GetNewJobID(&job_ids, job_num);
  nimbus::IDSet<nimbus::logical_data_id_t> read, write;
  nimbus::IDSet<nimbus::job_id_t> before, after;

  int extrapolation_job_num = kAppPartNum;
  std::vector<nimbus::job_id_t> extrapolation_job_ids;
  GetNewJobID(&extrapolation_job_ids, extrapolation_job_num);

  int extrapolate_phi_job_num = kAppPartNum;
  std::vector<nimbus::job_id_t> extrapolate_phi_job_ids;
  GetNewJobID(&extrapolate_phi_job_ids, extrapolate_phi_job_num);

  int reseed_particles_job_num = kAppPartNum;
  std::vector<nimbus::job_id_t> reseed_particles_job_ids;
  GetNewJobID(&reseed_particles_job_ids, reseed_particles_job_num);

  int write_output_job_num = kAppPartNum;
  std::vector<nimbus::job_id_t> write_output_job_ids;
  GetNewJobID(&write_output_job_ids, write_output_job_num);

  int calculate_dt_job_num = kAppPartNum;
  std::vector<nimbus::job_id_t> calculate_dt_job_ids;
  GetNewJobID(&calculate_dt_job_ids, calculate_dt_job_num);

  if (!done) {
    StartTemplate("loop_iteration_part_two");
  } else {
    StartTemplate("loop_iteration_part_two_end");
  }

  /*
   * Spawning extrapolate phi stage over multiple workers.
   */
  for (int i = 0; i < extrapolate_phi_job_num; ++i) {
    read.clear();
    LoadLdoIdsInSet(&read, kRegY2W3Outer[i], APP_PHI,
                        APP_FACE_VEL, NULL);
    write.clear();
    LoadLdoIdsInSet(&write,
                        kRegY2W3CentralWGB[i], APP_PHI,
                        NULL);

    nimbus::Parameter s_extra_params;
    std::string s_extra_str;
    SerializeParameter(frame, time, dt, kPNAInt,
                       global_region, kRegY2W3Central[i],
                       kPNAInt, &s_extra_str);
    s_extra_params.set_ser_data(SerializedData(s_extra_str));
    before.clear();
    StageJobAndLoadBeforeSet(&before, EXTRAPOLATE_PHI, extrapolate_phi_job_ids[i],
                       read, write);

    SpawnComputeJob(EXTRAPOLATE_PHI, extrapolate_phi_job_ids[i],
                       read, write, before, after,
                       s_extra_params, true,
                       kRegY2W3Central[i]);
  }

  MarkEndOfStage();

  /*
   * Spawning extrapolation stage over multiple workers
   */
  for (int i = 0; i < extrapolation_job_num; ++i) {
    read.clear();
    LoadLdoIdsInSet(&read, kRegY2W3Outer[i],
                        APP_FACE_VEL, APP_PHI, NULL);
    write.clear();
    LoadLdoIdsInSet(&write, kRegY2W3Central[i],
                        APP_FACE_VEL, NULL);

    nimbus::Parameter extrapolation_params;
    std::string extrapolation_str;
    SerializeParameter(frame, time, dt, kPNAInt,
                       global_region, kRegY2W3Central[i],
                       kPNAInt, &extrapolation_str);
    extrapolation_params.set_ser_data(SerializedData(extrapolation_str));
    before.clear();
    StageJobAndLoadBeforeSet(&before, EXTRAPOLATION, extrapolation_job_ids[i],
                       read, write);

    SpawnComputeJob(EXTRAPOLATION, extrapolation_job_ids[i],
                       read, write, before, after,
                       extrapolation_params, true,
                       kRegY2W3Central[i]);
  }

  MarkEndOfStage();

  if (done) {
    dbg(APP_LOG, "[CONTROL FLOW] Second part, Loop done.\n");
  } else {
    dbg(APP_LOG, "[CONTROL FLOW] Second part, Loop not done.\n");
  }
  dbg(APP_LOG, "[CONTROL FLOW] Second part, Frame=%d, Time=%f, dt=%f\n",
      frame, time, dt);

  if (!done) {

    // Spawning loop iteration for next iteration.

    for (int i = 0; i < calculate_dt_job_num; ++i) {
      read.clear();
      LoadLdoIdsInSet(&read, kRegY2W3Outer[i], APP_FACE_VEL,
                          APP_PHI, NULL);
      write.clear();
      LoadLdoIdsInSet(&write, kRegY2W3Central[i], APP_DT, NULL);

      nimbus::Parameter dt_params;
      std::string dt_str;
      SerializeParameter(frame, time, 0, kPNAInt,
                         global_region, kRegY2W3Central[i],
                         kPNAInt, &dt_str);
      dt_params.set_ser_data(SerializedData(dt_str));
      before.clear();
      StageJobAndLoadBeforeSet(&before, CALCULATE_DT, calculate_dt_job_ids[i],
                         read, write);

      SpawnComputeJob(CALCULATE_DT, calculate_dt_job_ids[i],
                         read, write, before, after,
                         dt_params, true,
                         kRegY2W3Central[i]);
    }
    MarkEndOfStage();

    read.clear();
    LoadLdoIdsInSet(&read, kRegW3Central[0], APP_DT, NULL);
    write.clear();
    nimbus::Parameter iter_params;
    std::string iter_str;
    SerializeParameter(frame, time + dt, kPNAFloat, kPNAInt,
                       global_region,
                       kPNAReg,
                       kPNAInt, &iter_str);
    iter_params.set_ser_data(SerializedData(iter_str));
    before.clear();
    StageJobAndLoadBeforeSet(&before, LOOP_ITERATION, job_ids[1],
                       read, write,
                       true);

    SpawnComputeJob(LOOP_ITERATION, job_ids[1],
                       read, write, before, after,
                       iter_params, false,
                       kRegW3Central[0],
                       true);
    MarkEndOfStage();

    EndTemplate("loop_iteration_part_two");
  } else {

    std::vector<nimbus::job_id_t> loop_job_id;
    GetNewJobID(&loop_job_id, 1);

    for (int i = 0; i < reseed_particles_job_num; ++i) {
      read.clear();
      LoadLdoIdsInSet(&read, kRegY2W3Outer[i], APP_FACE_VEL,
                          APP_PHI, NULL);
      LoadLdoIdsInSet(&read, kRegY2W1Outer[i], APP_PSI_D,
                          APP_PSI_N, NULL);
      LoadLdoIdsInSet(&read, kRegY2W3Outer[i], APP_POS_PARTICLES,
                          APP_NEG_PARTICLES, APP_POS_REM_PARTICLES,
                          APP_NEG_REM_PARTICLES, NULL);
      LoadLdoIdsInSet(&read, kRegY2W3CentralWGB[i], APP_LAST_UNIQUE_PARTICLE_ID , NULL);
      write.clear();
      LoadLdoIdsInSet(&write, kRegY2W3CentralWGB[i], APP_FACE_VEL,
                          APP_PHI, NULL);
      LoadLdoIdsInSet(&write, kRegY2W3CentralWGB[i], APP_POS_PARTICLES,
                          APP_NEG_PARTICLES, APP_POS_REM_PARTICLES,
                          APP_NEG_REM_PARTICLES, APP_LAST_UNIQUE_PARTICLE_ID,
                          NULL);

      nimbus::Parameter temp_params;
      std::string temp_str;
      SerializeParameter(frame, time, dt, kPNAInt,
                         global_region, kRegY2W3Central[i],
                         kPNAInt, &temp_str);
      temp_params.set_ser_data(SerializedData(temp_str));
      before.clear();
      StageJobAndLoadBeforeSet(&before, RESEED_PARTICLES,
          reseed_particles_job_ids[i],
          read, write);

      SpawnComputeJob(RESEED_PARTICLES,
          reseed_particles_job_ids[i],
          read, write, before, after,
          temp_params, true,
          kRegY2W3Central[i]);
    }
    MarkEndOfStage();

    if (kUseGlobalWrite) {
      read.clear();
      LoadLdoIdsInSet(&read, kRegW3Outer[0], APP_FACE_VEL,
                          APP_PHI, NULL);
      LoadLdoIdsInSet(&read, kRegW1Outer[0], APP_PSI_D,
                          APP_PSI_N, NULL);
      LoadLdoIdsInSet(&read, kRegW3Outer[0], APP_POS_PARTICLES,
                          APP_NEG_PARTICLES, APP_POS_REM_PARTICLES,
                          APP_NEG_REM_PARTICLES, NULL);
      LoadLdoIdsInSet(&read, kRegW3Central[0], APP_LAST_UNIQUE_PARTICLE_ID , NULL);
      write.clear();

      nimbus::Parameter temp_params;
      std::string temp_str;
      SerializeParameter(frame, time, dt, -1,
                         global_region, global_region,
                         kPNAInt, &temp_str);
      temp_params.set_ser_data(SerializedData(temp_str));
      before.clear();
      StageJobAndLoadBeforeSet(&before, WRITE_OUTPUT,
                         write_output_job_ids[0],
                         read, write);

      SpawnComputeJob(WRITE_OUTPUT,
                         write_output_job_ids[0],
                         read, write, before, after,
                         temp_params, true,
                         kRegW3Central[0]);
      MarkEndOfStage();
    } else {
      for (int i = 0; i < write_output_job_num; ++i) {
        read.clear();
        LoadLdoIdsInSet(&read, kRegY2W3Outer[i], APP_FACE_VEL,
                            APP_PHI, NULL);
        LoadLdoIdsInSet(&read, kRegY2W1Outer[i], APP_PSI_D,
                            APP_PSI_N, NULL);
        LoadLdoIdsInSet(&read, kRegY2W3Outer[i], APP_POS_PARTICLES,
                            APP_NEG_PARTICLES, APP_POS_REM_PARTICLES,
                            APP_NEG_REM_PARTICLES, NULL);
        LoadLdoIdsInSet(&read, kRegY2W3CentralWGB[i], APP_LAST_UNIQUE_PARTICLE_ID , NULL);
        write.clear();

        nimbus::Parameter temp_params;
        std::string temp_str;
        SerializeParameter(frame, time, dt, i+1,
                           global_region, kRegY2W3Central[i],
                           kPNAInt, &temp_str);
        temp_params.set_ser_data(SerializedData(temp_str));
        before.clear();
        StageJobAndLoadBeforeSet(&before, WRITE_OUTPUT,
                           write_output_job_ids[i],
                           read, write);

        SpawnComputeJob(WRITE_OUTPUT,
                           write_output_job_ids[i],
                           read, write, before, after,
                           temp_params, true,
                           kRegY2W3Central[i]);
      }
      MarkEndOfStage();
    }

    // Spawning loop frame to compute next frame.

    read.clear();
    write.clear();

    nimbus::Parameter frame_params;
    std::string frame_str;
    SerializeParameter(frame + 1, kPNAFloat, kPNAFloat, kPNAInt,
                       global_region,
                       kPNAReg,
                       kPNAInt, &frame_str);
    frame_params.set_ser_data(SerializedData(frame_str));
    before.clear();
    StageJobAndLoadBeforeSet(&before, LOOP_FRAME, loop_job_id[0],
                       read, write,
                       true);

    SpawnComputeJob(LOOP_FRAME, loop_job_id[0],
                       read, write, before, after,
                       frame_params, false,
                       kRegW3Central[0],
                       true);
    MarkEndOfStage();


    EndTemplate("loop_iteration_part_two_end");
  }  // end loop "if (done)".
  // job_query.PrintTimeProfile();
  if (time == 0) {
    dbg(APP_LOG, "Print job dependency figure.\n");
    // job_query.GenerateDotFigure("loop_iteration_part_two.dot");
  }
  {
    struct timeval t;
    gettimeofday(&t, NULL);
    double time  = (static_cast<double>(t.tv_sec - start_time.tv_sec)) +
        .000001 * (static_cast<double>(t.tv_usec - start_time.tv_usec));
    dbg(APP_LOG, "\nThe query time spent in job LOOP_ITERATION_PART_TWO is %f seconds.\n",
        time);
  }
}


}  // namespace application

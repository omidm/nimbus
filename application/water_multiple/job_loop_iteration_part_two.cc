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
  T dt;
  std::string params_str(params.ser_data().data_ptr_raw(),
                         params.ser_data().size());
  LoadParameter(params_str, &init_config.frame, &init_config.time, &dt,
                &init_config.global_region, &init_config.local_region);
  dbg(APP_LOG, "Frame %i in LOOP_ITERATION_PART_TWO job\n", init_config.frame);

  const int& frame = init_config.frame;
  const T& time = init_config.time;

  dbg(APP_LOG, "Frame %i and time %f in iteration job\n",
      frame, time);

  // Initialize the state of example and driver.
  PhysBAM::WATER_EXAMPLE<TV>* example =
      new PhysBAM::WATER_EXAMPLE<TV>(PhysBAM::STREAM_TYPE((RW())));

  // check whether the frame is done or not
  bool done = false;
  if (time + dt >= example->Time_At_Frame(frame + 1)) {
    done = true;
  }

  delete example;

  SpawnJobs(
      done, frame, time, dt, da, init_config.global_region);

}

void JobLoopIterationPartTwo::SpawnJobs(
    bool done, int frame, T time, T dt, const nimbus::DataArray& da,
    const nimbus::GeometricRegion& global_region) {

  int job_num = 2;
  std::vector<nimbus::job_id_t> job_ids;
  GetNewJobID(&job_ids, job_num);
  nimbus::IDSet<nimbus::logical_data_id_t> read, write;
  nimbus::IDSet<nimbus::job_id_t> before, after;

  // Spawning extrapolation stage over entire block

  read.clear();
  LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL, APP_PHI, NULL);
  write.clear();
  LoadLogicalIdsInSet(this, &write,
                      kRegW3Outer[0], APP_FACE_VEL, APP_PHI, NULL);

  nimbus::Parameter extrapolation_params;
  std::string extrapolation_str;
  SerializeParameter(frame, time, dt, global_region, global_region,
                     &extrapolation_str);
  extrapolation_params.set_ser_data(SerializedData(extrapolation_str));
  after.clear();
  after.insert(job_ids[1]);
  before.clear();
  SpawnComputeJob(EXTRAPOLATION,
                  job_ids[0],
                  read, write,
                  before, after,
                  extrapolation_params);

  if (!done) {

    // Spawning loop iteration for next iteration.

    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL,
                        APP_FACE_VEL_GHOST, APP_PHI, NULL);
    LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_POS_PARTICLES,
                        APP_NEG_PARTICLES, APP_POS_REM_PARTICLES,
                        APP_NEG_REM_PARTICLES, APP_LAST_UNIQUE_PARTICLE_ID,
                        NULL);
    write.clear();
    LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_FACE_VEL,
                        APP_FACE_VEL_GHOST, APP_PHI, NULL);
    LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_POS_PARTICLES,
                        APP_NEG_PARTICLES, APP_POS_REM_PARTICLES,
                        APP_NEG_REM_PARTICLES, APP_LAST_UNIQUE_PARTICLE_ID,
                        NULL);

    nimbus::Parameter iter_params;
    std::string iter_str;
    SerializeParameter(frame, time + dt, global_region, &iter_str);
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

    std::vector<nimbus::job_id_t> loop_job_id;
    GetNewJobID(&loop_job_id, 1);

    // Spawning write frame over entire block

    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL,
                        APP_FACE_VEL_GHOST, APP_PHI, NULL);
    LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_POS_PARTICLES,
                        APP_NEG_PARTICLES, APP_POS_REM_PARTICLES,
                        APP_NEG_REM_PARTICLES, APP_LAST_UNIQUE_PARTICLE_ID,
                        NULL);
    write.clear();
    LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_FACE_VEL,
                        APP_FACE_VEL_GHOST, APP_PHI, NULL);
    LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_POS_PARTICLES,
                        APP_NEG_PARTICLES, APP_POS_REM_PARTICLES,
                        APP_NEG_REM_PARTICLES, APP_LAST_UNIQUE_PARTICLE_ID,
                        NULL);

    nimbus::Parameter write_params;
    std::string write_str;
    SerializeParameter(frame, time + dt, 0,
                       global_region, global_region, &write_str);
    write_params.set_ser_data(SerializedData(write_str));
    after.clear();
    after.insert(loop_job_id[0]);
    before.clear();
    before.insert(job_ids[0]);
    SpawnComputeJob(WRITE_FRAME,
                    job_ids[1],
                    read, write,
                    before, after,
                    write_params);

    // Spawning loop frame to compute next frame.

    read.clear();
    write.clear();

    nimbus::Parameter frame_params;
    std::string frame_str;
    SerializeParameter(frame + 1, global_region, &frame_str);
    frame_params.set_ser_data(SerializedData(frame_str));
    after.clear();
    before.clear();
    before.insert(job_ids[1]);
    SpawnComputeJob(LOOP_FRAME,
                    loop_job_id[0],
                    read, write,
                    before, after,
                    frame_params);

  }  // end loop "if (done)".
}


}  // namespace application

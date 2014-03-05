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
#include "shared/dbg.h"
#include "shared/nimbus.h"
#include <sstream>
#include <string>

#include "application/water_multiple/projection/job_projection_main.h"

namespace application {

JobProjectionMain::JobProjectionMain(nimbus::Application *app) {
  set_application(app);
};

nimbus::Job* JobProjectionMain::Clone() {
  return new JobProjectionMain(application());
}

void JobProjectionMain::Execute(
    nimbus::Parameter params,
    const nimbus::DataArray& da) {
  dbg(APP_LOG, "Executing PROJECTION_MAIN job\n");

  // Get parameters: frame, time
  InitConfig init_config;
  T dt;
  std::string params_str(params.ser_data().data_ptr_raw(),
                         params.ser_data().size());
  LoadParameter(params_str, &init_config.frame, &init_config.time, &dt,
                &init_config.global_region, &init_config.local_region);

  const int& frame = init_config.frame;
  const T& time = init_config.time;

  dbg(APP_LOG, "Frame %i and time %f in PROJECTION_MAIN job\n",
      frame, time);
  SpawnJobs(frame, time, dt, da, init_config.global_region);
}

void JobProjectionMain::SpawnJobs(
    int frame, T time, T dt, const nimbus::DataArray& da,
    const nimbus::GeometricRegion& global_region) {
  int projection_job_num = 5;
  std::vector<nimbus::job_id_t> projection_job_ids;
  GetNewJobID(&projection_job_ids, projection_job_num);
  nimbus::IDSet<nimbus::logical_data_id_t> read, write;
  nimbus::IDSet<nimbus::job_id_t> before, after;
  nimbus::Parameter default_params;
  std::string default_params_str;
  SerializeParameter(
      frame, time, dt, global_region, global_region,
      &default_params_str);
  default_params.set_ser_data(SerializedData(default_params_str));

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

  before.clear();
  after.clear();
  after.insert(projection_job_ids[1]);

  SpawnComputeJob(PROJECTION_CALCULATE_BOUNDARY_CONDITION,
                  projection_job_ids[0],
                  read, write,
                  before, after,
                  default_params);

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
                      APP_U_INTERFACE, APP_MATRIX_A,
                      APP_VECTOR_B, APP_PROJECTION_LOCAL_TOLERANCE,
                      APP_INDEX_M2C, APP_INDEX_C2M,
                      APP_PROJECTION_LOCAL_N, APP_PROJECTION_INTERIOR_N,
                      NULL);

  before.clear();
  before.insert(projection_job_ids[0]);
  after.clear();
  after.insert(projection_job_ids[2]);
  SpawnComputeJob(PROJECTION_CONSTRUCT_MATRIX,
                  projection_job_ids[1],
                  read, write,
                  before, after,
                  default_params);

  // Loop job LOCAL_INITIALIZE/GLOBAL_INITIALIZE, and then spawns
  // PROJECTION_LOOP_ITERATION.
  read.clear();
  LoadLogicalIdsInSet(this, &read, kRegW0Central[0],
                      APP_PROJECTION_LOCAL_N, APP_PROJECTION_INTERIOR_N,
                      APP_PRESSURE, APP_INDEX_M2C,
                      APP_VECTOR_B,
                      APP_MATRIX_A, NULL);
  write.clear();
  LoadLogicalIdsInSet(this, &write, kRegW0Central[0],
                      APP_VECTOR_B, APP_PROJECTION_LOCAL_RESIDUAL, APP_MATRIX_C,
                      APP_VECTOR_TEMP, APP_VECTOR_P, APP_VECTOR_Z, NULL);
  before.clear();
  before.insert(projection_job_ids[1]);
  after.clear();
  after.insert(projection_job_ids[3]);
  SpawnComputeJob(PROJECTION_LOCAL_INITIALIZE,
                  projection_job_ids[2],
                  read, write,
                  before, after,
                  default_params);

  read.clear();
  LoadLogicalIdsInSet(this, &read, kRegW0Central[0],
                      APP_PROJECTION_INTERIOR_N, APP_PROJECTION_LOCAL_TOLERANCE,
                      NULL);
  write.clear();
  LoadLogicalIdsInSet(this, &write, kRegW0Central[0],
                      APP_PROJECTION_GLOBAL_N,
                      APP_PROJECTION_GLOBAL_TOLERANCE,
                      APP_PROJECTION_DESIRED_ITERATIONS, NULL);
  before.clear();
  before.insert(projection_job_ids[2]);
  after.clear();
  after.insert(projection_job_ids[4]);
  SpawnComputeJob(PROJECTION_GLOBAL_INITIALIZE,
                  projection_job_ids[3],
                  read, write,
                  before, after,
                  default_params);

  read.clear();
  LoadLogicalIdsInSet(this, &read, kRegW0Central[0],
                      APP_PROJECTION_LOCAL_RESIDUAL,
                      APP_PROJECTION_GLOBAL_TOLERANCE,
                      APP_PROJECTION_DESIRED_ITERATIONS,
                      NULL);
  write.clear();
  before.clear();
  before.insert(projection_job_ids[3]);
  after.clear();
  nimbus::Parameter projection_loop_iteration_params;
  std::string projection_loop_iteration_str;
  SerializeParameter(
      frame, time, dt, global_region, global_region, 1,
      &projection_loop_iteration_str);
  projection_loop_iteration_params.set_ser_data(
      SerializedData(projection_loop_iteration_str));
  SpawnComputeJob(PROJECTION_LOOP_ITERATION,
                  projection_job_ids[4],
                  read, write,
                  before, after,
                  projection_loop_iteration_params);
}


}  // namespace application

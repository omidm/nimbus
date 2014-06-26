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

#include <sstream>
#include <string>
#include <vector>

#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/data_names.h"
#include "application/water_multiple/job_names.h"
#include "application/water_multiple/physbam_utils.h"
#include "application/water_multiple/reg_def.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"
#include "worker/job_query.h"

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
  init_config.use_cache = true;
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
  struct timeval start_time;
  gettimeofday(&start_time, NULL);
  nimbus::JobQuery job_query(this);
  int projection_job_num = 5;
  std::vector<nimbus::job_id_t> projection_job_ids;
  GetNewJobID(&projection_job_ids, projection_job_num);
  nimbus::IDSet<nimbus::logical_data_id_t> read, write;

  nimbus::Parameter default_params;
  std::string default_params_str;
  SerializeParameter(
      frame, time, dt, global_region, global_region,
      &default_params_str);
  default_params.set_ser_data(SerializedData(default_params_str));

  std::vector<nimbus::Parameter> default_part_params;
  default_part_params.resize(kAppPartNum);
  for (int i = 0; i < kAppPartNum; ++i) {
    std::string default_params_str;
    SerializeParameter(
        frame, time, dt, global_region, kRegY2W0Central[i],
        &default_params_str);
    default_part_params[i].set_ser_data(SerializedData(default_params_str));
  }

  int construct_matrix_job_num = kAppPartNum;
  std::vector<nimbus::job_id_t> construct_matrix_job_ids;
  GetNewJobID(&construct_matrix_job_ids, construct_matrix_job_num);

  int local_initialize_job_num = kAppPartNum;
  std::vector<nimbus::job_id_t> local_initialize_job_ids;
  GetNewJobID(&local_initialize_job_ids, local_initialize_job_num);

  int calculate_boundary_condition_part_one_job_num = kAppPartNum;
  std::vector<nimbus::job_id_t> calculate_boundary_condition_part_one_job_ids;
  GetNewJobID(&calculate_boundary_condition_part_one_job_ids,
              calculate_boundary_condition_part_one_job_num);

  int calculate_boundary_condition_part_two_job_num = kAppPartNum;
  std::vector<nimbus::job_id_t> calculate_boundary_condition_part_two_job_ids;
  GetNewJobID(&calculate_boundary_condition_part_two_job_ids,
              calculate_boundary_condition_part_two_job_num);

  for (int index = 0;
       index < calculate_boundary_condition_part_one_job_num;
       ++index) {
    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[index], APP_FACE_VEL, APP_PHI, NULL);
    LoadLogicalIdsInSet(this, &read, kRegY2W1Outer[index],
                        APP_DIVERGENCE, APP_PSI_D, APP_PSI_N,
                        APP_FILLED_REGION_COLORS, APP_PRESSURE, NULL);
    LoadLogicalIdsInSet(this, &read, kRegY2W0Central[index],
                        APP_U_INTERFACE, NULL);
    write.clear();
    LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[index], APP_FACE_VEL, APP_PHI, NULL);
    LoadLogicalIdsInSet(this, &write, kRegY2W1CentralWGB[index],
                        APP_DIVERGENCE, APP_PSI_D, APP_PSI_N,
                        APP_FILLED_REGION_COLORS, APP_PRESSURE, NULL);
    LoadLogicalIdsInSet(this, &write, kRegY2W0Central[index], APP_U_INTERFACE, NULL);

    job_query.StageJob(PROJECTION_CALCULATE_BOUNDARY_CONDITION_PART_ONE,
                       calculate_boundary_condition_part_one_job_ids[index],
                       read, write,
                       default_part_params[index], true);
  }
  job_query.CommitStagedJobs();

  for (int index = 0;
       index < calculate_boundary_condition_part_two_job_num;
       ++index) {
    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[index], APP_FACE_VEL, APP_PHI, NULL);
    LoadLogicalIdsInSet(this, &read, kRegY2W1Outer[index],
                        APP_DIVERGENCE, APP_PSI_D, APP_PSI_N,
                        APP_FILLED_REGION_COLORS, APP_PRESSURE, NULL);
    LoadLogicalIdsInSet(this, &read, kRegY2W0Central[index], APP_U_INTERFACE, NULL);
    write.clear();
    LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[index], APP_FACE_VEL, APP_PHI, NULL);
    LoadLogicalIdsInSet(this, &write, kRegY2W1CentralWGB[index],
                        APP_DIVERGENCE, APP_PSI_D, APP_PSI_N,
                        APP_FILLED_REGION_COLORS, APP_PRESSURE, NULL);
    LoadLogicalIdsInSet(this, &write, kRegY2W0Central[index],
                        APP_U_INTERFACE, NULL);

    job_query.StageJob(PROJECTION_CALCULATE_BOUNDARY_CONDITION_PART_TWO,
                       calculate_boundary_condition_part_two_job_ids[index],
                       read, write,
                       default_part_params[index], true);
  }
  job_query.CommitStagedJobs();

  // Construct matrix.
  for (int index = 0; index < construct_matrix_job_num; ++index) {
    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[index], APP_FACE_VEL, APP_PHI, NULL);
    LoadLogicalIdsInSet(this, &read, kRegY2W1Outer[index],
                        APP_DIVERGENCE, APP_PSI_D, APP_PSI_N,
                        APP_FILLED_REGION_COLORS, APP_PRESSURE, NULL);
    LoadLogicalIdsInSet(this, &read, kRegY2W0Central[index], APP_U_INTERFACE, NULL);
    write.clear();
    LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[index], APP_FACE_VEL, APP_PHI, NULL);
    LoadLogicalIdsInSet(this, &write, kRegY2W1CentralWGB[index],
                        APP_DIVERGENCE, APP_PSI_D, APP_PSI_N,
                        APP_FILLED_REGION_COLORS, APP_PRESSURE, NULL);
    LoadLogicalIdsInSet(this, &write, kRegY2W0Central[index],
                        APP_U_INTERFACE, APP_MATRIX_A,
                        APP_VECTOR_B, APP_PROJECTION_LOCAL_TOLERANCE,
                        APP_INDEX_M2C, APP_INDEX_C2M,
                        APP_PROJECTION_LOCAL_N, APP_PROJECTION_INTERIOR_N,
                        NULL);

    job_query.StageJob(PROJECTION_CONSTRUCT_MATRIX,
                       construct_matrix_job_ids[index],
                       read, write,
                       default_part_params[index], true);
  }
  job_query.CommitStagedJobs();

  // Local initialize.
  for (int index = 0; index < local_initialize_job_num; ++index) {
    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegY2W0Central[index],
                        APP_PROJECTION_LOCAL_N, APP_PROJECTION_INTERIOR_N,
                        APP_INDEX_M2C,
                        APP_VECTOR_B,
                        APP_MATRIX_A, NULL);
    LoadLogicalIdsInSet(this, &read, kRegY2W1Outer[index], APP_PRESSURE,
                        NULL);
    write.clear();
    LoadLogicalIdsInSet(this, &write, kRegY2W0Central[index],
                        APP_VECTOR_B, APP_PROJECTION_LOCAL_RESIDUAL, APP_MATRIX_C,
                        APP_VECTOR_TEMP, APP_VECTOR_Z,
                        APP_VECTOR_P_GRID_FORMAT, APP_VECTOR_P_LINEAR_FORMAT,
                        NULL);
    job_query.StageJob(PROJECTION_LOCAL_INITIALIZE,
                       local_initialize_job_ids[index],
                       read, write,
                       default_part_params[index], true);
  }
  job_query.CommitStagedJobs();

  // Global initialize.
  read.clear();
  LoadLogicalIdsInSet(this, &read, kRegW0Central[0],
                      APP_PROJECTION_INTERIOR_N, APP_PROJECTION_LOCAL_TOLERANCE,
                      NULL);
  write.clear();
  LoadLogicalIdsInSet(this, &write, kRegW0Central[0],
                      APP_PROJECTION_GLOBAL_N,
                      APP_PROJECTION_GLOBAL_TOLERANCE,
                      APP_PROJECTION_DESIRED_ITERATIONS, NULL);
  job_query.StageJob(PROJECTION_GLOBAL_INITIALIZE,
                     projection_job_ids[3],
                     read, write,
                     default_params, true);
  job_query.CommitStagedJobs();

  // Projection loop.
  read.clear();
  LoadLogicalIdsInSet(this, &read, kRegW0Central[0],
                      APP_PROJECTION_LOCAL_RESIDUAL,
                      APP_PROJECTION_GLOBAL_TOLERANCE,
                      APP_PROJECTION_DESIRED_ITERATIONS,
                      NULL);
  write.clear();
  nimbus::Parameter projection_loop_iteration_params;
  std::string projection_loop_iteration_str;
  SerializeParameter(
      frame, time, dt, global_region, global_region, 1,
      &projection_loop_iteration_str);
  projection_loop_iteration_params.set_ser_data(
      SerializedData(projection_loop_iteration_str));
  job_query.StageJob(PROJECTION_LOOP_ITERATION,
                     projection_job_ids[4],
                     read, write,
                     projection_loop_iteration_params, false, true);
  job_query.CommitStagedJobs();

  if (time == 0) {
    dbg(APP_LOG, "Print job dependency figure.\n");
    job_query.GenerateDotFigure("projection_main.dot");
  }
  {
    struct timeval t;
    gettimeofday(&t, NULL);
    double time  = (static_cast<double>(t.tv_sec - start_time.tv_sec)) +
        .000001 * (static_cast<double>(t.tv_usec - start_time.tv_usec));
    dbg(APP_LOG, "\nThe query time spent in job PROJECTION_LOOP_ITERATION_MAIN is %f seconds.\n",
        time);
  }
}


}  // namespace application

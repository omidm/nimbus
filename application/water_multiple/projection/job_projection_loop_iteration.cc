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
 * Author: Hang Qu <quhang@stanford.edu>
 */

#include <sstream>
#include <string>

#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>

#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/data_names.h"
#include "application/water_multiple/job_names.h"
#include "application/water_multiple/physbam_utils.h"
#include "application/water_multiple/projection/projection_driver.h"
#include "application/water_multiple/reg_def.h"
#include "application/water_multiple/water_driver.h"
#include "application/water_multiple/water_example.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"

#include "application/water_multiple/projection/job_projection_loop_iteration.h"

namespace application {

JobProjectionLoopIteration::JobProjectionLoopIteration(
    nimbus::Application *app) {
  set_application(app);
}

nimbus::Job* JobProjectionLoopIteration::Clone() {
  return new JobProjectionLoopIteration(application());
}

void JobProjectionLoopIteration::Execute(
    nimbus::Parameter params,
    const nimbus::DataArray& da) {
  dbg(APP_LOG, "Executing PROJECTION_LOOP_ITERATION job.\n");

  InitConfig init_config;
  T dt;
  int iteration;
  std::string params_str(params.ser_data().data_ptr_raw(),
                         params.ser_data().size());
  LoadParameter(params_str, &init_config.frame, &init_config.time, &dt,
                &init_config.global_region, &init_config.local_region,
                &iteration);
  const nimbus::GeometricRegion& global_region = init_config.global_region;
  const int& frame = init_config.frame;
  const T& time = init_config.time;

  DataConfig data_config;
  // TODO(quhang), LOCAL_N and INTERIOR_N might not be needed.
  data_config.SetFlag(DataConfig::PROJECTION_LOCAL_N);
  data_config.SetFlag(DataConfig::PROJECTION_INTERIOR_N);
  data_config.SetFlag(DataConfig::PROJECTION_LOCAL_RESIDUAL);
  data_config.SetFlag(DataConfig::PROJECTION_GLOBAL_TOLERANCE);
  data_config.SetFlag(DataConfig::PROJECTION_DESIRED_ITERATIONS);

  PhysBAM::PCG_SPARSE<float> pcg_temp;
  pcg_temp.Set_Maximum_Iterations(40);
  pcg_temp.evolution_solver_type = PhysBAM::krylov_solver_cg;
  pcg_temp.cg_restart_iterations = 0;
  pcg_temp.Show_Results();
  PhysBAM::ProjectionDriver projection_driver(
      pcg_temp, init_config, data_config);

  dbg(APP_LOG, "Job PROJECTION_LOOP_ITERATION starts (dt=%f).\n", dt);

  projection_driver.LoadFromNimbus(this, da);
  projection_driver.projection_data.residual =
      projection_driver.projection_data.local_residual;

  // Decides whether to spawn a new projection loop or finish it.
  if (projection_driver.projection_data.residual <=
      projection_driver.projection_data.global_tolerance ||
      projection_driver.projection_data.iteration ==
      projection_driver.projection_data.desired_iterations) {
    // Finishes projection iterating.
    int projection_job_num = 2;
    std::vector<nimbus::job_id_t> projection_job_ids;
    GetNewJobID(&projection_job_ids, projection_job_num);

    int wrapup_job_num = 2;
    std::vector<nimbus::job_id_t> wrapup_job_ids;
    GetNewJobID(&wrapup_job_ids, wrapup_job_num);

    nimbus::IDSet<nimbus::logical_data_id_t> read, write;
    nimbus::IDSet<nimbus::job_id_t> before, after;

    // Projection wrapup.
    /*
    for (int index = 0; index < 1; ++index) {
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[0], APP_FACE_VEL, APP_PHI,
                          NULL);
      LoadLogicalIdsInSet(this, &read, kRegY2W1Outer[0], APP_PSI_D,
                          APP_PRESSURE, NULL);
      LoadLogicalIdsInSet(this, &read, kRegY2W1Central[0], APP_PSI_N,
                          APP_U_INTERFACE, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegY2W0Central[0], APP_FACE_VEL, NULL);
      LoadLogicalIdsInSet(this, &write, kRegY2W1Outer[0], APP_PRESSURE, NULL);
      nimbus::Parameter projection_wrapup_params;
      std::string projection_wrapup_str;
      SerializeParameter(frame, time, dt, global_region, global_region,
                         &projection_wrapup_str);
      projection_wrapup_params.set_ser_data(
          SerializedData(projection_wrapup_str));
      before.clear();
      after.clear();
      after.insert(projection_job_ids[1]);
      SpawnComputeJob(PROJECTION_WRAPUP, projection_job_ids[0],
                      read, write, before, after, projection_wrapup_params);
    }
    */

      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL, APP_PHI,
                          NULL);
      LoadLogicalIdsInSet(this, &read, kRegW1Outer[0], APP_PSI_D,
                          APP_PRESSURE, NULL);
      LoadLogicalIdsInSet(this, &read, kRegW1Central[0], APP_PSI_N,
                          APP_U_INTERFACE, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegW0Central[0], APP_FACE_VEL, NULL);
      LoadLogicalIdsInSet(this, &write, kRegW1Outer[0], APP_PRESSURE, NULL);
      nimbus::Parameter projection_wrapup_params;
      std::string projection_wrapup_str;
      SerializeParameter(frame, time, dt, global_region, global_region,
                         &projection_wrapup_str);
      projection_wrapup_params.set_ser_data(
          SerializedData(projection_wrapup_str));
      before.clear();
      after.clear();
      after.insert(projection_job_ids[1]);
      SpawnComputeJob(PROJECTION_WRAPUP, projection_job_ids[0],
                      read, write, before, after, projection_wrapup_params);
    // Loop iteration part two.
    read.clear();
    write.clear();
    nimbus::Parameter loop_iteration_part_two_params;
    std::string loop_iteration_part_two_str;
    SerializeParameter(frame, time, dt, global_region, global_region,
                       &loop_iteration_part_two_str);
    loop_iteration_part_two_params.set_ser_data(
        SerializedData(loop_iteration_part_two_str));
    before.clear();
    before.insert(wrapup_job_ids[0]);
    // before.insert(wrapup_job_ids[1]);
    after.clear();
    SpawnComputeJob(LOOP_ITERATION_PART_TWO, projection_job_ids[1],
                    read, write, before, after, loop_iteration_part_two_params);
  } else {
    // Spawns a new projection iteration.
    nimbus::Parameter default_params;
    std::string default_params_str;
    SerializeParameter(frame, time, dt, global_region, global_region, iteration,
                       &default_params_str);
    default_params.set_ser_data(SerializedData(default_params_str));

    nimbus::Parameter default_part_params[2];
    std::string default_left_params_str;
    SerializeParameter(
        frame, time, dt, global_region, kRegY2W0Central[0],
        &default_left_params_str);
    default_part_params[0].set_ser_data(
        SerializedData(default_left_params_str));
    std::string default_right_params_str;
    SerializeParameter(
        frame, time, dt, global_region, kRegY2W0Central[1],
        &default_right_params_str);
    default_part_params[1].set_ser_data(
        SerializedData(default_right_params_str));

    int projection_job_num = 7;
    std::vector<nimbus::job_id_t> projection_job_ids;
    GetNewJobID(&projection_job_ids, projection_job_num);

    int step_one_job_num = 2;
    std::vector<nimbus::job_id_t> step_one_job_ids;
    GetNewJobID(&step_one_job_ids, step_one_job_num);

    int step_two_job_num = 2;
    std::vector<nimbus::job_id_t> step_two_job_ids;
    GetNewJobID(&step_two_job_ids, step_two_job_num);

    int step_three_job_num = 2;
    std::vector<nimbus::job_id_t> step_three_job_ids;
    GetNewJobID(&step_three_job_ids, step_three_job_num);

    int step_four_job_num = 2;
    std::vector<nimbus::job_id_t> step_four_job_ids;
    GetNewJobID(&step_four_job_ids, step_four_job_num);

    nimbus::IDSet<nimbus::logical_data_id_t> read, write;
    nimbus::IDSet<nimbus::job_id_t> before, after;

    // STEP_ONE.
    for (int index = 0; index < step_one_job_num; ++index) {
      read.clear();
      LoadLogicalIdsInSet(
          this, &read, kRegY2W0Central[index], APP_PROJECTION_LOCAL_N,
          APP_PROJECTION_INTERIOR_N, APP_MATRIX_C, APP_VECTOR_B, APP_VECTOR_Z,
          APP_INDEX_M2C, NULL);
      write.clear();
      LoadLogicalIdsInSet(
          this, &write, kRegY2W0Central[index], APP_VECTOR_Z,
          APP_PROJECTION_LOCAL_RHO, NULL);
      before.clear();
      after.clear();
      after.insert(projection_job_ids[1]);
      SpawnComputeJob(PROJECTION_STEP_ONE, step_one_job_ids[index],
                      read, write, before, after, default_part_params[index]);
    }

    // REDUCE_RHO
    read.clear();
    LoadLogicalIdsInSet(
        this, &read, kRegW0Central[0], APP_PROJECTION_LOCAL_N,
        APP_PROJECTION_INTERIOR_N, APP_PROJECTION_LOCAL_RHO,
        APP_PROJECTION_GLOBAL_RHO, NULL);
    write.clear();
    LoadLogicalIdsInSet(
        this, &write, kRegW0Central[0], APP_PROJECTION_LOCAL_RHO,
        APP_PROJECTION_GLOBAL_RHO, APP_PROJECTION_GLOBAL_RHO_OLD,
        APP_PROJECTION_BETA, NULL);
    before.clear();
    before.insert(step_one_job_ids[0]);
    before.insert(step_one_job_ids[1]);
    after.clear();
    after.insert(step_two_job_ids[0]);
    after.insert(step_two_job_ids[1]);
    SpawnComputeJob(PROJECTION_REDUCE_RHO, projection_job_ids[1],
                    read, write, before, after, default_params);

    // STEP_TWO
    for (int index = 0; index < step_two_job_num; ++index) {
      read.clear();
      LoadLogicalIdsInSet(
          this, &read, kRegW0Central[0], APP_PROJECTION_BETA, NULL);
      LoadLogicalIdsInSet(
          this, &read, kRegY2W0Central[index], APP_VECTOR_Z, APP_VECTOR_P,
          APP_INDEX_M2C, APP_PROJECTION_LOCAL_N, APP_PROJECTION_INTERIOR_N,
          NULL);
      write.clear();
      LoadLogicalIdsInSet(
          this, &write, kRegY2W1CentralWGB[index], APP_VECTOR_P, NULL);
      before.clear();
      before.insert(projection_job_ids[1]);
      after.clear();
      after.insert(step_three_job_ids[0]);
      after.insert(step_three_job_ids[1]);
      SpawnComputeJob(PROJECTION_STEP_TWO, step_two_job_ids[index],
                      read, write, before, after, default_part_params[index]);
    }

    // STEP_THREE
    for (int index = 0; index < step_three_job_num; ++index) {
      read.clear();
      LoadLogicalIdsInSet(
          this, &read, kRegY2W1Outer[index], APP_VECTOR_P, NULL);
      LoadLogicalIdsInSet(
          this, &read, kRegY2W0Central[index], APP_MATRIX_A, APP_INDEX_M2C,
          APP_PROJECTION_LOCAL_N, APP_PROJECTION_INTERIOR_N, NULL);
      write.clear();
      LoadLogicalIdsInSet(
          this, &write, kRegY2W0Central[index], APP_VECTOR_TEMP,
          APP_PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA, NULL);
      before.clear();
      before.insert(step_two_job_ids[0]);
      before.insert(step_two_job_ids[1]);
      after.clear();
      after.insert(projection_job_ids[4]);
      SpawnComputeJob(PROJECTION_STEP_THREE, step_three_job_ids[0],
                      read, write, before, after, default_part_params[index]);
    }

    // REDUCE_ALPHA
    read.clear();
    LoadLogicalIdsInSet(
        this, &read, kRegW0Central[0], APP_PROJECTION_LOCAL_N,
        APP_PROJECTION_INTERIOR_N, APP_PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA,
        NULL);
    LoadLogicalIdsInSet(
        this, &read, kRegW0Central[0], APP_PROJECTION_GLOBAL_RHO, NULL);
    write.clear();
    LoadLogicalIdsInSet(
        this, &write, kRegW0Central[0], APP_PROJECTION_ALPHA, NULL);
    before.clear();
    before.insert(step_three_job_ids[0]);
    before.insert(step_three_job_ids[1]);
    after.clear();
    after.insert(step_four_job_ids[0]);
    after.insert(step_four_job_ids[1]);
    SpawnComputeJob(PROJECTION_REDUCE_ALPHA, projection_job_ids[4],
                    read, write, before, after, default_params);

    // STEP_FOUR
    // Only interior p is needed.
    for (int index = 0; index < step_four_job_num; ++index) {
      read.clear();
      LoadLogicalIdsInSet(
          this, &read, kRegY2W0Central[index], APP_PROJECTION_LOCAL_N,
          APP_PROJECTION_INTERIOR_N, APP_INDEX_M2C, APP_VECTOR_P,
          APP_VECTOR_TEMP, APP_VECTOR_B, APP_PRESSURE, APP_INDEX_M2C, NULL);
      LoadLogicalIdsInSet(
          this, &read, kRegW0Central[0], APP_PROJECTION_ALPHA, NULL);
      write.clear();
      LoadLogicalIdsInSet(
          this, &write, kRegY2W0Central[index], APP_VECTOR_B,
          APP_PROJECTION_LOCAL_RESIDUAL, APP_PRESSURE, NULL);
      before.clear();
      before.insert(projection_job_ids[4]);
      after.clear();
      after.insert(projection_job_ids[6]);
      SpawnComputeJob(PROJECTION_STEP_FOUR, step_four_job_ids[0],
                      read, write, before, after, default_part_params[index]);
    }

    // NEXT_ITERATION
    // TODO(quhang), needs a clean way to deal with LOCAL_N and INTERIOR_N here.
    read.clear();
    LoadLogicalIdsInSet(
        this, &read, kRegW0Central[0], APP_PROJECTION_LOCAL_RESIDUAL,
        APP_PROJECTION_GLOBAL_TOLERANCE, APP_PROJECTION_DESIRED_ITERATIONS,
        NULL);
    write.clear();
    before.clear();
    before.insert(step_four_job_ids[0]);
    before.insert(step_four_job_ids[1]);
    after.clear();
    nimbus::Parameter next_iteration_params;
    std::string next_iteration_params_str;
    SerializeParameter(frame, time, dt, global_region, global_region,
                       iteration + 1,
                       &next_iteration_params_str);
    next_iteration_params.set_ser_data(
        SerializedData(next_iteration_params_str));
    SpawnComputeJob(PROJECTION_LOOP_ITERATION, projection_job_ids[6],
                    read, write, before, after,
                    next_iteration_params);
  }

  // TODO(quhang), removes the saving if possible.
  projection_driver.SaveToNimbus(this, da);

  dbg(APP_LOG, "Completed executing PROJECTION_LOOP_ITERATION job\n");
}

}  // namespace application

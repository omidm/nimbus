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
 * Author: Omid Mashayekhi <omidm@stanford.edu>
 */

#include <sstream>
#include <string>
#include <vector>

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
#include "worker/job_query.h"
#include "worker/worker_thread.h"

#include "application/water_multiple/projection/job_projection_loop_bottleneck.h"

namespace application {

JobProjectionLoopBottleneck::JobProjectionLoopBottleneck(
    nimbus::Application *app) {
  set_application(app);
}

nimbus::Job* JobProjectionLoopBottleneck::Clone() {
  return new JobProjectionLoopBottleneck(application());
}

void JobProjectionLoopBottleneck::Execute(
    nimbus::Parameter params,
    const nimbus::DataArray& da) {
  struct timeval start_time;
  gettimeofday(&start_time, NULL);
  dbg(APP_LOG, "Executing PROJECTION_LOOP_BOTTLENECK job.\n");
  // nimbus::JobQuery job_query(this);

  InitConfig init_config;
  std::string params_str(params.ser_data().data_ptr_raw(),
                         params.ser_data().size());
  LoadParameter(params_str, &init_config);
  T dt = init_config.dt;
  int iteration = init_config.projection_iteration;
  // const nimbus::GeometricRegion& global_region = init_config.global_region;
  // const int& frame = init_config.frame;
  // const T& time = init_config.time;

  DataConfig data_config;
  // TODO(quhang), LOCAL_N and INTERIOR_N might not be needed.
  data_config.SetFlag(DataConfig::PROJECTION_LOCAL_N);
  data_config.SetFlag(DataConfig::PROJECTION_INTERIOR_N);
  data_config.SetFlag(DataConfig::PROJECTION_LOCAL_RESIDUAL);
  data_config.SetFlag(DataConfig::PROJECTION_GLOBAL_TOLERANCE);
  data_config.SetFlag(DataConfig::PROJECTION_DESIRED_ITERATIONS);

  PhysBAM::PCG_SPARSE<float> pcg_temp;
  pcg_temp.Set_Maximum_Iterations(application::kMaxIterations);
  pcg_temp.evolution_solver_type = PhysBAM::krylov_solver_cg;
  pcg_temp.cg_restart_iterations = 0;
  pcg_temp.Show_Results();
  PhysBAM::ProjectionDriver projection_driver(
      pcg_temp, init_config, data_config, &worker_thread()->allocated_threads);

  dbg(APP_LOG, "Job PROJECTION_LOOP_BOTTLENECK starts (dt=%f).\n", dt);

  projection_driver.LoadFromNimbus(this, da);
  projection_driver.projection_data.residual =
      projection_driver.projection_data.local_residual;
  projection_driver.projection_data.iteration = iteration;

  dbg(APP_LOG, "[CONTROL FLOW] size = %d, Iteration = %d, "
      "Desired iteration = %d, "
      "Residual = %f, Global tolerance = %f\n",
      projection_driver.projection_data.interior_n,
      projection_driver.projection_data.iteration,
      projection_driver.projection_data.desired_iterations,
      projection_driver.projection_data.residual,
      projection_driver.projection_data.global_tolerance);
  // Decides whether to spawn a new projection loop or finish it.
  bool end_iterations = 
    (projection_driver.projection_data.residual <=
    projection_driver.projection_data.global_tolerance) ||
    (projection_driver.projection_data.iteration ==
    projection_driver.projection_data.desired_iterations);

  if (end_iterations) {
    dbg(APP_LOG, "Iterations could end!\n");
  }

  // TODO(quhang), removes the saving if possible.
  projection_driver.SaveToNimbus(this, da);

  dbg(APP_LOG, "Completed executing PROJECTION_LOOP_BOTTLENECK job\n");
}

}  // namespace application

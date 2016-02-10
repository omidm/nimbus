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

#include "applications/physbam/water//app_utils.h"
#include "applications/physbam/water//physbam_utils.h"
#include "applications/physbam/water//projection/projection_driver.h"
#include "applications/physbam/water//water_driver.h"
#include "applications/physbam/water//water_example.h"
#include "src/shared/dbg.h"
#include "src/shared/nimbus.h"

#include "src/data/scalar_data.h"
#include "applications/physbam/water//data_include.h"
#include "src/worker/worker_thread.h"
#include "applications/physbam/water//projection/job_projection_reduce_rho.h"

namespace application {

JobProjectionReduceRho::JobProjectionReduceRho(nimbus::Application *app) {
  set_application(app);
};

nimbus::Job* JobProjectionReduceRho::Clone() {
  return new JobProjectionReduceRho(application());
}

void JobProjectionReduceRho::Execute(
    nimbus::Parameter params,
    const nimbus::DataArray& da) {
  dbg(APP_LOG, "Executing PROJECTION_REDUCE_RHO job.\n");

  InitConfig init_config;
  std::string params_str(params.ser_data().data_ptr_raw(),
                         params.ser_data().size());
  LoadParameter(params_str, &init_config);
  T dt = init_config.dt;
  int iteration = init_config.projection_iteration;

  DataConfig data_config;
  data_config.SetFlag(DataConfig::PROJECTION_LOCAL_RHO);
  data_config.SetFlag(DataConfig::PROJECTION_GLOBAL_RHO);
  data_config.SetFlag(DataConfig::PROJECTION_GLOBAL_RHO_OLD);
  data_config.SetFlag(DataConfig::PROJECTION_BETA);

  PhysBAM::PCG_SPARSE<float> pcg_temp;
  pcg_temp.Set_Maximum_Iterations(application::kMaxIterations);
  pcg_temp.evolution_solver_type = PhysBAM::krylov_solver_cg;
  pcg_temp.cg_restart_iterations = 0;
  pcg_temp.Show_Results();

  PhysBAM::ProjectionDriver projection_driver(
      pcg_temp, init_config, data_config, &worker_thread()->allocated_threads);
  projection_driver.projection_data.iteration = iteration;
  dbg(APP_LOG, "Job PROJECTION_REDUCE_RHO starts (dt=%f).\n", dt);

  projection_driver.LoadFromNimbus(this, da);

  // Read PROJECTION_LOCAL_RHO, PROJECTION_GLOBAL_RHO.
  // Write PROJECTION_GLOBAL_RHO, PROJECTION_GLOBAL_RHO_OLD, PROJECTION_BETA.
  // TODO(quhang), seems like rho_old is not needed.
  {
    application::ScopeTimer scope_timer(name());
    projection_driver.ReduceRho();
  }

  projection_driver.SaveToNimbus(this, da);

  dbg(APP_LOG, "Completed executing PROJECTION_REDUCE_RHO job\n");
}

}  // namespace application

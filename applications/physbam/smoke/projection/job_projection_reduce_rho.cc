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
 * Modifier for smoke: Andrew Lim <alim16@stanford.edu> 
 */

#include <sstream>
#include <string>

#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>

#include "applications/physbam/smoke/app_utils.h"
#include "applications/physbam/smoke/physbam_utils.h"
#include "applications/physbam/smoke/projection/projection_driver.h"
#include "applications/physbam/smoke/smoke_driver.h"
#include "applications/physbam/smoke/smoke_example.h"
#include "src/shared/dbg.h"
#include "src/shared/nimbus.h"

#include "src/data/scalar_data.h"
#include "applications/physbam/smoke/data_include.h"
#include "applications/physbam/smoke/projection/job_projection_reduce_rho.h"

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
  init_config.use_cache = true;
  T dt;
  int iteration;
  std::string params_str(params.ser_data().data_ptr_raw(),
                         params.ser_data().size());
  LoadParameter(params_str, &init_config.frame, &init_config.time, &dt,
                &init_config.global_region, &init_config.local_region,
                &iteration);

  DataConfig data_config;
  data_config.SetFlag(DataConfig::PROJECTION_LOCAL_N);
  data_config.SetFlag(DataConfig::PROJECTION_INTERIOR_N);
  data_config.SetFlag(DataConfig::PROJECTION_LOCAL_RHO);
  data_config.SetFlag(DataConfig::PROJECTION_GLOBAL_RHO);
  data_config.SetFlag(DataConfig::PROJECTION_GLOBAL_RHO_OLD);
  data_config.SetFlag(DataConfig::PROJECTION_BETA);

  PhysBAM::PCG_SPARSE<float> pcg_temp;
  pcg_temp.Set_Maximum_Iterations(10);
  pcg_temp.evolution_solver_type = PhysBAM::krylov_solver_cg;
  pcg_temp.cg_restart_iterations = 40;

  PhysBAM::ProjectionDriver projection_driver(
      pcg_temp, init_config, data_config);
  projection_driver.projection_data.iteration = iteration;
  dbg(APP_LOG, "Job PROJECTION_REDUCE_RHO starts (dt=%f).\n", dt);

  Log log_timer;

  log_timer.StartTimer();

  {
    application::ScopeTimer scope_timer(name() + "-load");
    projection_driver.LoadFromNimbus(this, da);
  }

  // Read PROJECTION_LOCAL_RHO, PROJECTION_GLOBAL_RHO.
  // Write PROJECTION_GLOBAL_RHO, PROJECTION_GLOBAL_RHO_OLD, PROJECTION_BETA.
  // TODO(quhang), seems like rho_old is not needed.
  {
    application::ScopeTimer scope_timer(name());
    projection_driver.ReduceRho();
  }

  {
    application::ScopeTimer scope_timer(name() + "-save");
    projection_driver.SaveToNimbus(this, da);
  }

  dbg(APP_LOG, "[PROJECTION] PROJECTION_REDUCE_RHO total time:%f.\n",
      log_timer.GetTime());

  dbg(APP_LOG, "Completed executing PROJECTION_REDUCE_RHO job\n");
}

}  // namespace application

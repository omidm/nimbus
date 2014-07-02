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
#include "application/water_multiple/physbam_utils.h"
#include "application/water_multiple/projection/projection_driver.h"
#include "application/water_multiple/water_driver.h"
#include "application/water_multiple/water_example.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"

#include "data/scalar_data.h"
#include "application/water_multiple/data_include.h"
#include "application/water_multiple/projection/job_projection_reduce_alpha.h"

namespace application {

JobProjectionReduceAlpha::JobProjectionReduceAlpha(nimbus::Application *app) {
  set_application(app);
};

nimbus::Job* JobProjectionReduceAlpha::Clone() {
  return new JobProjectionReduceAlpha(application());
}

void JobProjectionReduceAlpha::Execute(
    nimbus::Parameter params,
    const nimbus::DataArray& da) {
  dbg(APP_LOG, "Executing PROJECTION_REDUCE_ALPHA job.\n");

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
  // TODO(quhang), remove.
  data_config.SetFlag(DataConfig::PROJECTION_LOCAL_N);
  data_config.SetFlag(DataConfig::PROJECTION_INTERIOR_N);
  data_config.SetFlag(DataConfig::PROJECTION_GLOBAL_RHO);
  data_config.SetFlag(DataConfig::PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA);
  data_config.SetFlag(DataConfig::PROJECTION_ALPHA);

  PhysBAM::PCG_SPARSE<float> pcg_temp;
  pcg_temp.Set_Maximum_Iterations(40);
  pcg_temp.evolution_solver_type = PhysBAM::krylov_solver_cg;
  pcg_temp.cg_restart_iterations = 0;
  pcg_temp.Show_Results();

  PhysBAM::ProjectionDriver projection_driver(
      pcg_temp, init_config, data_config);
  projection_driver.projection_data.iteration = iteration;
  dbg(APP_LOG, "Job PROJECTION_REDUCE_ALPHA starts (dt=%f).\n", dt);

  Log log_timer;

  log_timer.StartTimer();

  projection_driver.LoadFromNimbus(this, da);

  // Read PROJECTION_GLOBAL_RHO, PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA.
  // Write PROJECTION_ALPHA.
  projection_driver.ReduceAlpha();

  projection_driver.SaveToNimbus(this, da);
  dbg(APP_LOG, "[PROJECTION] PROJECTION_REDUCE_ALPHA total time:%f.\n",
      log_timer.GetTime());

  dbg(APP_LOG, "Completed executing PROJECTION_REDUCE_ALPHA job\n");
}

}  // namespace application

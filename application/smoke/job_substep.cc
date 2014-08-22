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
 * This file contains the substep job, which computes the time step for
 * an iteration during the computation of a simulation frame.
 *
 * Author: Andrew Lim <alim16@stanford.edu> 
 */

#include "application/smoke/app_utils.h"
#include "application/smoke/job_substep.h"
#include "application/smoke/physbam_utils.h"
#include "application/smoke/smoke_driver.h"
#include "application/smoke/smoke_example.h"
#include "application/smoke/job_names.h"
#include "application/smoke/data_names.h"
#include "application/smoke/reg_def.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"
#include "worker/job_query.h"
#include <sstream>
#include <string>

namespace application {

JobSubstep::JobSubstep(nimbus::Application *app) {
  set_application(app);
};

nimbus::Job* JobSubstep::Clone() {
  return new JobSubstep(application());
}

void JobSubstep::Execute(
    nimbus::Parameter params,
    const nimbus::DataArray& da) {
  dbg(APP_LOG, "Executing substep job.\n");
  InitConfig init_config;
  T dt_dumb;
  std::string params_str(params.ser_data().data_ptr_raw(),
                         params.ser_data().size());
  LoadParameter(params_str, &init_config.frame, &init_config.time, &dt_dumb,
                &init_config.global_region, &init_config.local_region);

  // initialize configuration and state
  PhysBAM::SMOKE_EXAMPLE<TV> *example;
  PhysBAM::SMOKE_DRIVER<TV> *driver;

  init_config.set_boundary_condition = false;
  init_config.use_cache = true;
  // Threading settings.
  init_config.use_threading = use_threading();
  init_config.core_quota = core_quota();

  DataConfig data_config;
  data_config.SetFlag(DataConfig::VELOCITY);
  data_config.SetFlag(DataConfig::DENSITY);
  data_config.SetFlag(DataConfig::DT);

  InitializeExampleAndDriver(init_config, data_config,
                             this, da, example, driver);
  *thread_queue_hook() = example->nimbus_thread_queue;

  {
    application::ScopeTimer scope_timer(name());
    T dt = example->cfl * example->CFL(example->face_velocities);
    dbg(APP_LOG, "[CONTROL FLOW] dt=%f\n", dt);
    example->dt_buffer = dt;
  }

  *thread_queue_hook() = NULL;

  example->Save_To_Nimbus(this, da, driver->current_frame + 1);

  // Free resources.
  DestroyExampleAndDriver(example, driver);
  dbg(APP_LOG, "Completed executing substep.\n");
}

} // namespace application

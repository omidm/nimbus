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

#include "application/smoke/app_utils.h"
#include "application/smoke/physbam_utils.h"
#include "application/smoke/smoke_driver.h"
#include "application/smoke/smoke_example.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"

#include "application/smoke/job_write_output.h"

namespace application {

JobWriteOutput::JobWriteOutput(nimbus::Application *app) {
  set_application(app);
};

nimbus::Job* JobWriteOutput::Clone() {
  return new JobWriteOutput(application());
}

void JobWriteOutput::Execute(nimbus::Parameter params,
                            const nimbus::DataArray& da) {
  dbg(APP_LOG, "Executing WRITE_OUTPUT job.\n");

  InitConfig init_config;
  // Threading settings.
  init_config.use_threading = use_threading();
  init_config.core_quota = core_quota();
  init_config.use_cache = true;
  init_config.set_boundary_condition = false;
  T dt;
  std::string params_str(params.ser_data().data_ptr_raw(),
                         params.ser_data().size());
  int rank;
  LoadParameter(params_str, &init_config.frame, &init_config.time, &dt, &rank,
                &init_config.global_region, &init_config.local_region);

  dbg(APP_LOG,
      "In WRITE_OUTPUT: Initialize SMOKE_DRIVER/SMOKE_EXAMPLE"
      "(Output=%d, Time=%f).\n",
      init_config.frame, init_config.time, dt);

  PhysBAM::SMOKE_EXAMPLE<TV> *example;
  PhysBAM::SMOKE_DRIVER<TV> *driver;
  DataConfig data_config;
  data_config.SetFlag(DataConfig::VELOCITY);
  data_config.SetFlag(DataConfig::DENSITY);
  data_config.SetFlag(DataConfig::PSI_N);
  data_config.SetFlag(DataConfig::PSI_D);
  {
    application::ScopeTimer scope_timer(name() + "-load");
    InitializeExampleAndDriver(init_config, data_config,
			       this, da, example, driver);
  }
  *thread_queue_hook() = example->nimbus_thread_queue;

  dbg(APP_LOG, "Job WRITE_OUTPUT starts.\n");
  // Write frame.
  driver->WriteOutputSplitImpl(this, da, true, dt, rank);

  *thread_queue_hook() = NULL;
  {
    application::ScopeTimer scope_timer(name() + "-save");
    example->Save_To_Nimbus(this, da, driver->current_frame + 1);
  }
   
  // Free resources.
  DestroyExampleAndDriver(example, driver);

  dbg(APP_LOG, "Completed executing WRITE_OUTPUT job\n");
}

}  // namespace application

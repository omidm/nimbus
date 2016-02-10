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

#include "applications/physbam/water//app_utils.h"
#include "applications/physbam/water//job_calculate_dt.h"
#include "applications/physbam/water//physbam_utils.h"
#include "applications/physbam/water//water_driver.h"
#include "applications/physbam/water//water_example.h"
#include "applications/physbam/water//water_sources.h"
#include "applications/physbam/water//job_names.h"
#include "applications/physbam/water//data_names.h"
#include "applications/physbam/water//reg_def.h"
#include "src/shared/dbg.h"
#include "src/shared/nimbus.h"
#include "src/shared/fast_log.h"
#include "src/worker/job_query.h"
#include <sstream>
#include <string>

namespace application {

JobCalculateDt::JobCalculateDt(nimbus::Application *app) {
  set_application(app);
};

nimbus::Job* JobCalculateDt::Clone() {
  return new JobCalculateDt(application());
}

void JobCalculateDt::Execute(
    nimbus::Parameter params,
    const nimbus::DataArray& da) {
  nimbus::timer::AddCounter(nimbus::timer::kCalculateDt);
  dbg(APP_LOG, "Executing CALCULATE_DT job\n");
  InitConfig init_config;
  init_config.set_boundary_condition = false;

  std::string params_str(params.ser_data().data_ptr_raw(),
                         params.ser_data().size());
  LoadParameter(params_str, &init_config);

  // initialize configuration and state
  PhysBAM::WATER_EXAMPLE<TV> *example;
  PhysBAM::WATER_DRIVER<TV> *driver;

  DataConfig data_config;
  data_config.SetFlag(DataConfig::VELOCITY);
  data_config.SetFlag(DataConfig::DT);
  data_config.SetFlag(DataConfig::LEVELSET);
  nimbus::timer::StartTimer(nimbus::timer::kInitializeCalculateDt);
  InitializeExampleAndDriver(init_config, data_config,
                             this, da, example, driver);
  nimbus::timer::StopTimer(nimbus::timer::kInitializeCalculateDt);

  {
    application::ScopeTimer scope_timer(name());
    T dt = example->cfl * example->incompressible.CFL(example->face_velocities);
    T temp_dt =
        example->particle_levelset_evolution.cfl_number
        * example->particle_levelset_evolution.particle_levelset.levelset.CFL(
            example->face_velocities);
    dbg(APP_LOG, "[CONTROL FLOW] dt=%f\n", temp_dt);
    if (temp_dt < dt) {
      dt = temp_dt;
    }
    example->dt_buffer = dt;
  }

  // Write, free resources.
  example->Save_To_Nimbus(this, da, driver->current_frame + 1);
  DestroyExampleAndDriver(example, driver);

  dbg(APP_LOG, "Completed executing CALCULATE_DT job\n");
}

} // namespace application

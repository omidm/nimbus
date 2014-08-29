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

#include "application/smoke/app_utils.h"
#include "application/smoke/physbam_utils.h"
#include "application/smoke/smoke_driver.h"
#include "application/smoke/smoke_example.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"

#include "application/smoke/projection/job_projection_calculate_boundary_condition_part_two.h"

namespace application {

JobProjectionCalculateBoundaryConditionPartTwo::
    JobProjectionCalculateBoundaryConditionPartTwo(nimbus::Application *app) {
  set_application(app);
};

nimbus::Job* JobProjectionCalculateBoundaryConditionPartTwo::Clone() {
  return new JobProjectionCalculateBoundaryConditionPartTwo(application());
}

void JobProjectionCalculateBoundaryConditionPartTwo::Execute(
    nimbus::Parameter params,
    const nimbus::DataArray& da) {
  dbg(APP_LOG,
      "Executing PROJECTION_CALCULATE_BOUNDARY_CONDITION_PART_TWO job.\n");

  InitConfig init_config;
  init_config.use_cache = true;
  init_config.set_boundary_condition = false;

  T dt;
  std::string params_str(params.ser_data().data_ptr_raw(),
                         params.ser_data().size());
  LoadParameter(params_str, &init_config.frame, &init_config.time, &dt,
                &init_config.global_region, &init_config.local_region);

  // Assume time, dt, frame is ready from here.
  dbg(APP_LOG,
      "In PROJECTION_CALCULATE_BOUNDARY_CONDITION_PART_TWO: "
      "Initialize SMOKE_DRIVER/SMOKE_EXAMPLE"
      "(Frame=%d, Time=%f).\n",
      init_config.frame, init_config.time);

  PhysBAM::SMOKE_EXAMPLE<TV> *example;
  PhysBAM::SMOKE_DRIVER<TV> *driver;

  DataConfig data_config;
  data_config.SetFlag(DataConfig::VELOCITY);
  data_config.SetFlag(DataConfig::DIVERGENCE);
  data_config.SetFlag(DataConfig::PSI_N);
  data_config.SetFlag(DataConfig::PSI_D);
  data_config.SetFlag(DataConfig::REGION_COLORS);
  data_config.SetFlag(DataConfig::PRESSURE);
  data_config.SetFlag(DataConfig::U_INTERFACE);

  {
    application::ScopeTimer scope_timer(name() + "-load");
    InitializeExampleAndDriver(init_config, data_config,
			       this, da, example, driver);
  }

  dbg(APP_LOG,
      "Job PROJECTION_CALCULATE_BOUNDARY_CONDITION_PART_TWO starts (dt=%f).\n",
      dt);

  {
    application::ScopeTimer scope_timer(name());
    driver->ProjectionCalculateBoundaryConditionPartTwoImpl(this, da, dt);
  }

  {
    application::ScopeTimer scope_timer(name() + "-save");
    example->Save_To_Nimbus(this, da, driver->current_frame + 1);
    // Free resources.
    DestroyExampleAndDriver(example, driver);
  }

  dbg(APP_LOG,
      "Completed executing"
      " PROJECTION_CALCULATE_BOUNDARY_CONDITION_PART_TWO job\n");
}

}  // namespace application

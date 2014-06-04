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

#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/physbam_utils.h"
#include "application/water_multiple/water_driver.h"
#include "application/water_multiple/water_example.h"
#include "application/water_multiple/water_sources.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"

#include "application/water_multiple/job_write_output.h"

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
  init_config.use_cache = true;
  init_config.set_boundary_condition = false;
  T dt;
  std::string params_str(params.ser_data().data_ptr_raw(),
                         params.ser_data().size());
  int rank;
  LoadParameter(params_str, &init_config.frame, &init_config.time, &dt, &rank,
                &init_config.global_region, &init_config.local_region);

  dbg(APP_LOG,
      "In WRITE_OUTPUT: Initialize WATER_DRIVER/WATER_EXAMPLE"
      "(Output=%d, Time=%f).\n",
      init_config.frame, init_config.time, dt);

  PhysBAM::WATER_EXAMPLE<TV> *example;
  PhysBAM::WATER_DRIVER<TV> *driver;
  DataConfig data_config;
  data_config.SetFlag(DataConfig::VELOCITY);
  data_config.SetFlag(DataConfig::LEVELSET_READ);
  data_config.SetFlag(DataConfig::POSITIVE_PARTICLE);
  data_config.SetFlag(DataConfig::NEGATIVE_PARTICLE);
  data_config.SetFlag(DataConfig::REMOVED_POSITIVE_PARTICLE);
  data_config.SetFlag(DataConfig::REMOVED_NEGATIVE_PARTICLE);
  data_config.SetFlag(DataConfig::PSI_N);
  data_config.SetFlag(DataConfig::PSI_D);
  InitializeExampleAndDriver(init_config, data_config,
                             this, da, example, driver);

  dbg(APP_LOG, "Job WRITE_OUTPUT starts.\n");
  // Reseed particles and write frame.
  driver->WriteOutputSplitImpl(this, da, true, dt, rank);

  // Free resources.
  DestroyExampleAndDriver(example, driver);

  dbg(APP_LOG, "Completed executing WRITE_OUTPUT job\n");
}

}  // namespace application

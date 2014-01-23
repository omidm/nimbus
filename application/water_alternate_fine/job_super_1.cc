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
 * This file contains super job 1 which is the first portion of computations in
 * each iteration of computing the frame, it includes 
 *  Compute Occupied Blocks
 *  Adjust Phi With Objects
 *  Advect Phi
 *  Step Particles
 *  Advect Removed Particles
 *  Advect V
 *  Forces
 *
 * Author: Omid Mashayekhi <omidm@stanford.edu>
 */

#include <sstream>
#include <string>

#include "application/water_alternate_fine/app_utils.h"
#include "application/water_alternate_fine/physbam_utils.h"
#include "application/water_alternate_fine/water_driver.h"
#include "application/water_alternate_fine/water_example.h"
#include "application/water_alternate_fine/water_sources.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"

#include "application/water_alternate_fine/job_super_1.h"

namespace application {

JobSuper1::JobSuper1(nimbus::Application *app) {
  set_application(app);
};

nimbus::Job* JobSuper1::Clone() {
  return new JobSuper1(application());
}

void JobSuper1::Execute(nimbus::Parameter params,
                        const nimbus::DataArray& da) {
  dbg(APP_LOG, "Executing super job 1.\n");

  // get time, dt, frame from the parameters.
  InitConfig init_config;
  T dt;
  std::string params_str(params.ser_data().data_ptr_raw(),
                         params.ser_data().size());
  LoadParameter(params_str, &init_config.frame, &init_config.time, &dt);
  dbg(APP_LOG, " Loaded parameters (Frame=%d, Time=%f, dt=%f).\n",
      init_config.frame, init_config.time, dt);

  // Initializing the example and driver with state and configuration variables.
  PhysBAM::WATER_EXAMPLE<TV> *example;
  PhysBAM::WATER_DRIVER<TV> *driver;

  DataConfig data_config;
  data_config.SetAll();
  InitializeExampleAndDriver(init_config, data_config,
                             this, da, example, driver);

  // Run the steps in the super job.
  dbg(APP_LOG, "Execute the steps is super job 1.");
  driver->SuperJob1Impl(this, da, dt);

  // Free resources.
  DestroyExampleAndDriver(example, driver);

  dbg(APP_LOG, "Completed executing super job 1.\n");
}

}  // namespace application

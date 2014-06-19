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
 * This file contains step particles job which is one of the sub
 * jobs in the iteration of computing a simulation frame.
 *
 * Author: Omid Mashayekhi <omidm@stanford.edu>
 */

#include <sstream>
#include <string>

#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/physbam_utils.h"
#include "application/water_multiple/water_driver.h"
#include "application/water_multiple/water_example.h"
#include "application/water_multiple/water_sources.h"
#include "data/physbam/translator_physbam.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "shared/nimbus.h"

#include "application/water_multiple/job_step_particles.h"

namespace application {

JobStepParticles::JobStepParticles(nimbus::Application *app) {
  set_application(app);
};

nimbus::Job* JobStepParticles::Clone() {
  return new JobStepParticles(application());
}

void JobStepParticles::Execute(nimbus::Parameter params,
                        const nimbus::DataArray& da) {
  dbg(APP_LOG, "--- Executing step particles job.\n");

  // get time, dt, frame from the parameters.
  InitConfig init_config;
  // Threading settings.
  init_config.use_threading = use_threading();
  init_config.core_quota = core_quota();
  init_config.use_cache = true;
  init_config.set_boundary_condition = false;
  T dt;
  std::string params_str(params.ser_data().data_ptr_raw(),
                         params.ser_data().size());
  LoadParameter(params_str, &init_config.frame, &init_config.time, &dt,
                &init_config.global_region, &init_config.local_region);
  dbg(APP_LOG, " Loaded parameters (Frame=%d, Time=%f, dt=%f).\n",
      init_config.frame, init_config.time, dt);

  // Initializing the example and driver with state and configuration variables.
  PhysBAM::WATER_EXAMPLE<TV> *example;
  PhysBAM::WATER_DRIVER<TV> *driver;

  typedef application::DataConfig DataConfig;
  DataConfig data_config;
  data_config.SetFlag(DataConfig::VELOCITY_GHOST);
  data_config.SetFlag(DataConfig::POSITIVE_PARTICLE);
  data_config.SetFlag(DataConfig::NEGATIVE_PARTICLE);
  data_config.SetFlag(DataConfig::REMOVED_POSITIVE_PARTICLE);
  data_config.SetFlag(DataConfig::REMOVED_NEGATIVE_PARTICLE);
  data_config.SetFlag(DataConfig::SHARED_PARTICLES_FLUSH);
  InitializeExampleAndDriver(init_config, data_config,
                             this, da, example, driver);
  *thread_queue_hook() = example->nimbus_thread_queue;

  // Run the computation in the job.
  dbg(APP_LOG, "Execute the step in step particles job.\n");
  {
    //nimbus::Timer timer(std::string("step_particle_") + id().toString());
    typedef typename nimbus::TranslatorPhysBAM<T> Translator;
    GeometricRegion lr = init_config.local_region;
    GeometricRegion gr = init_config.global_region;
    nimbus::Coord shift(lr.x() - gr.x(), lr.y() - gr.y(), lr.z() - gr.z());
    GeometricRegion er = lr.NewEnlarged(kGhostNum);
    Translator::DeleteParticles(shift, lr, er,
        &example->particle_levelset_evolution.particle_levelset, gr.dx(), true);
    Translator::DeleteParticles(shift, lr, er,
        &example->particle_levelset_evolution.particle_levelset, gr.dx(), false);
    Translator::DeleteRemovedParticles(shift, lr, er,
        &example->particle_levelset_evolution.particle_levelset, gr.dx(), true);
    Translator::DeleteRemovedParticles(shift, lr, er,
        &example->particle_levelset_evolution.particle_levelset, gr.dx(), false);
    driver->StepParticlesImpl(this, da, dt);
  }

  *thread_queue_hook() = NULL;
  example->Save_To_Nimbus(this, da, driver->current_frame + 1);
  // Free resources.
  DestroyExampleAndDriver(example, driver);

  dbg(APP_LOG, "Completed executing step particles.\n");
}

}  // namespace application


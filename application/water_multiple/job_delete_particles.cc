/* Copyright 2013 Stanford University.
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
 * This file contains the job that modifies levelset and particles. This job
 * should be spawned after calling advance time step forces and before
 * adjusting phi with sources.
 *
 * This job needs all data (levelset, particles and face velocities).
 *
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/job_delete_particles.h"
#include "application/water_multiple/physbam_utils.h"
#include "application/water_multiple/water_driver.h"
#include "application/water_multiple/water_example.h"
#include "data/physbam/physbam_data.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"
#include <sstream>
#include <string>

namespace application {

JobDeleteParticles::JobDeleteParticles(nimbus::Application *app) {
    set_application(app);
};

nimbus::Job* JobDeleteParticles::Clone() {
    return new JobDeleteParticles(application());
}

void JobDeleteParticles::Execute(nimbus::Parameter params, const nimbus::DataArray& da) {
    dbg(APP_LOG, "--- Executing delete particles job\n");

    InitConfig init_config;
    init_config.use_cache = true;
    // Threading settings.
    init_config.use_threading = use_threading();
    init_config.core_quota = core_quota();
    T dt;
    std::string params_str(params.ser_data().data_ptr_raw(),
                           params.ser_data().size());
    LoadParameter(params_str, &init_config.frame, &init_config.time, &dt,
                  &init_config.global_region, &init_config.local_region);
    dbg(APP_LOG, "Frame %i in delete particles job\n", init_config.frame);

    // initialize configuration and state
    PhysBAM::WATER_EXAMPLE<TV> *example;
    PhysBAM::WATER_DRIVER<TV> *driver;

    init_config.set_boundary_condition = false;
    DataConfig data_config;
    data_config.SetFlag(DataConfig::VELOCITY_GHOST);
    data_config.SetFlag(DataConfig::LEVELSET);
    data_config.SetFlag(DataConfig::POSITIVE_PARTICLE);
    data_config.SetFlag(DataConfig::NEGATIVE_PARTICLE);
    data_config.SetFlag(DataConfig::REMOVED_POSITIVE_PARTICLE);
    data_config.SetFlag(DataConfig::REMOVED_NEGATIVE_PARTICLE);
    // TODO(quhang), why this is needed?
    data_config.SetAll();
    InitializeExampleAndDriver(init_config, data_config,
                               this, da, example, driver);
    *thread_queue_hook() = example->nimbus_thread_queue;

    {
      nimbus::Timer timer(std::string("delete_particles_") + id().toString());
      driver->DeleteParticlesImpl(this, da, dt);
    }

    *thread_queue_hook() = NULL;
    // free resources
    DestroyExampleAndDriver(example, driver);

    dbg(APP_LOG, "Completed executing delete particles job\n");
}

} // namespace application

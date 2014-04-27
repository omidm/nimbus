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
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/data_particle_array.h"
#include "application/water_multiple/job_synchronize_particles.h"
#include "application/water_multiple/job_names.h"
#include "application/water_multiple/options.h"
#include "application/water_multiple/parameters.h"
#include "application/water_multiple/physbam_utils.h"
#include "application/water_multiple/water_driver.h"
#include "application/water_multiple/water_example.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"

namespace application {

JobSynchronizeParticles::JobSynchronizeParticles(nimbus::Application *app) {
    set_application(app);
};

nimbus::Job* JobSynchronizeParticles::Clone() {
    return new JobSynchronizeParticles(application());
}

void JobSynchronizeParticles::Execute(nimbus::Parameter params, const nimbus::DataArray& da) {
    dbg(APP_LOG, "--- Executing synchronize particles job\n");

    // get time, dt, frame from the parameters.
    InitConfig init_config;
    init_config.use_cache = true;
    init_config.clear_shared_particles_read = true;
    init_config.set_boundary_condition = false;
    T dt;
    std::string params_str(params.ser_data().data_ptr_raw(),
                           params.ser_data().size());
    LoadParameter(params_str, &init_config.frame, &init_config.time, &dt,
                  &init_config.global_region, &init_config.local_region);
    dbg(APP_LOG, " Loaded parameters (Frame=%d, Time=%f, dt=%f).\n",
        init_config.frame, init_config.time, dt);

    // initializing the example and driver with state and configuration variables
    PhysBAM::WATER_EXAMPLE<TV> *example;
    PhysBAM::WATER_DRIVER<TV> *driver;
    typedef application::DataConfig DataConfig;
    DataConfig data_config;
    assert(da.size() > 0);
    dbg(APP_LOG, "Syncing over region = %s\n", init_config.local_region.toString().c_str());
    data_config.SetFlag(DataConfig::POSITIVE_PARTICLE);
    dbg(APP_LOG, "Syncing pos particles...\n");
    data_config.SetFlag(DataConfig::NEGATIVE_PARTICLE);
    dbg(APP_LOG, "Syncing neg particles...\n");
    data_config.SetFlag(DataConfig::REMOVED_POSITIVE_PARTICLE);
    dbg(APP_LOG, "Syncing pos removed particles...\n");
    data_config.SetFlag(DataConfig::REMOVED_NEGATIVE_PARTICLE);
    dbg(APP_LOG, "Syncing neg removed particles...\n");
    InitializeExampleAndDriver(init_config, data_config,
                               this, da, example, driver);

    example->Save_To_Nimbus(this, da, init_config.frame + 1);
    DestroyExampleAndDriver(example, driver);
    dbg(APP_LOG, "Completed executing synchronize particles job\n");
}

} // namespace application

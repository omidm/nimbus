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
    dbg(APP_LOG, "Executing synchronize particles job\n");

    bool sync_new = false;
    if (sync_new) {
        // get time, dt, frame from the parameters.
        InitConfig init_config;
        init_config.use_cache = false;
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
        if (da[0]->name() == APP_POS_PARTICLES) {
            data_config.SetFlag(DataConfig::POSITIVE_PARTICLE);
        } else if (da[0]->name() == APP_NEG_PARTICLES) {
            data_config.SetFlag(DataConfig::NEGATIVE_PARTICLE);
        } else if(da[0]->name() == APP_POS_REM_PARTICLES) {
            data_config.SetFlag(DataConfig::REMOVED_POSITIVE_PARTICLE);
        } else if(da[0]->name() == APP_NEG_REM_PARTICLES) {
            data_config.SetFlag(DataConfig::REMOVED_NEGATIVE_PARTICLE);
        } else {
            dbg(DBG_ERROR, "Merge particles fail!!!\n");
            exit(-1);
        }
        InitializeExampleAndDriver(init_config, data_config,
                                   this, da, example, driver);
        example->Save_To_Nimbus(this, da, init_config.frame + 1);
        DestroyExampleAndDriver(example, driver);
    } else {
        DataVec main_copy;
        DataSetVec scratch_copies;
        if (!GroupSyncData(this, da, &main_copy, &scratch_copies)) {
            dbg(DBG_WARN, "Nothing to synchronize\n");
            return;
        }

        for (size_t i = 0; i < main_copy.size(); i++) {
            DataParticleArray *merge_to = dynamic_cast<DataParticleArray *>(main_copy[i]);
            DataVec *scratch = scratch_copies[i];
            if (merge_to != NULL && !scratch->empty())
                merge_to->MergeParticles(*scratch);
            else
                dbg(DBG_WARN, "Passed object is not a particle array\n");
            delete scratch;
        }
    }
    dbg(APP_LOG, "Completed executing synchronize particles job\n");
}

} // namespace application

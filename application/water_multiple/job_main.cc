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
 * This file contains the "main" job that Nimbus launches after loading an
 * application. All subsequent jobs are spawned from here.
 *
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/data_def.h"
#include "application/water_multiple/data_names.h"
#include "application/water_multiple/job_main.h"
#include "application/water_multiple/job_names.h"
#include "application/water_multiple/reg_def.h"
#include "data/scratch_data_helper.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"
#include <vector>

namespace application {

    JobMain::JobMain(nimbus::Application *app) {
        set_application(app);
    };

    nimbus::Job* JobMain::Clone() {
        return new JobMain(application());
    }

    void JobMain::Execute(nimbus::Parameter params, const nimbus::DataArray& da) {
        dbg(APP_LOG, "Executing main job\n");

        DefineNimbusData(this);

        // Job setup
        int job_num = 1;
        std::vector<nimbus::job_id_t> job_ids;
        GetNewJobID(&job_ids, job_num);
        nimbus::IDSet<nimbus::logical_data_id_t> read, write;
        nimbus::IDSet<nimbus::job_id_t> before, after;
        nimbus::IDSet<nimbus::param_id_t> id_set;
        nimbus::Parameter init_params;

        // Init job
        read.clear();
        LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL, APP_FACE_VEL_GHOST, APP_PHI, NULL);
        LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_POS_PARTICLES,
            APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES,
            APP_LAST_UNIQUE_PARTICLE_ID , NULL);
        write.clear();
        LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_FACE_VEL, APP_FACE_VEL_GHOST, APP_PHI, NULL);
        LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_POS_PARTICLES,
            APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES,
            APP_LAST_UNIQUE_PARTICLE_ID , NULL);
        id_set.clear();
        kScratchPosParticles.GetJobScratchData(this, kRegW1Central[0], &id_set);
        kScratchNegParticles.GetJobScratchData(this, kRegW1Central[0], &id_set);
        kScratchPosRemParticles.GetJobScratchData(this, kRegW1Central[0], &id_set);
        kScratchNegRemParticles.GetJobScratchData(this, kRegW1Central[0], &id_set);

        init_params.set_ser_data(SerializedData(""));
        dbg(APP_LOG, "Spawning initialize\n");
        SpawnComputeJob(INITIALIZE,
                        job_ids[0],
                        read, write,
                        before, after,
                        init_params);

        dbg(APP_LOG, "Completed executing main job\n");
    }

} // namespace application

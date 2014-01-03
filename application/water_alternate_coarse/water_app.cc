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

#include "application/water_alternate_coarse/app_utils.h"
#include "application/water_alternate_coarse/data_app.h"
#include "application/water_alternate_coarse/job_initialize.h"
#include "application/water_alternate_coarse/job_iteration.h"
#include "application/water_alternate_coarse/job_loop.h"
#include "application/water_alternate_coarse/job_main.h"
#include "application/water_alternate_coarse/water_app.h"
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include "shared/dbg.h"
#include "shared/nimbus.h"
#include "stdio.h"

namespace application {

    WaterApp::WaterApp() {
    };

    /* Register data and job types and initialize constant quantities used by
     * application jobs. */
    void WaterApp::Load() {

        dbg_add_mode(APP_LOG_STR);
        dbg_add_mode(TRANSLATE_STR);

        dbg(APP_LOG, "Loading water application\n");

        // PhysBAM logging and R/W
        PhysBAM::LOG::Initialize_Logging(false, false, 1<<30, true, kThreadsNum);
        PhysBAM::FILE_UTILITIES::Create_Directory(kOutputDir+"/common");
        PhysBAM::LOG::Instance()->Copy_Log_To_File(kOutputDir+"/common/log.txt", false);

        dbg(APP_LOG, "Registering %s\n", APP_FACE_VEL);
        RegisterData(APP_FACE_VEL, new DataApp(APP_FACE_VEL, kFaceVelBufSize));
        dbg(APP_LOG, "Registering %s\n", APP_PHI);
        RegisterData(APP_PHI, new DataApp(APP_PHI, kPhiBufSize));
        dbg(APP_LOG, "Registering %s\n", APP_PRESSURE);
        RegisterData(APP_PRESSURE, new DataApp(APP_PRESSURE, kPressureBufSize));
        dbg(APP_LOG, "Registering %s\n", APP_POS_PARTICLES);
        RegisterData(APP_POS_PARTICLES, new DataApp(APP_POS_PARTICLES, kParticlesBufSize));
        dbg(APP_LOG, "Registering %s\n", APP_NEG_PARTICLES);
        RegisterData(APP_NEG_PARTICLES, new DataApp(APP_NEG_PARTICLES, kParticlesBufSize));

        RegisterJob(MAIN, new JobMain(this));
        RegisterJob(LOOP, new JobLoop(this));
        RegisterJob(ITERATION, new JobIteration(this));
        RegisterJob(INITIALIZE, new JobInitialize(this));

        dbg(APP_LOG, "Completed loading water application\n");
    }

} // namespace application

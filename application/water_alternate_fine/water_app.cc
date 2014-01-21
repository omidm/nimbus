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

#include "application/water_alternate_fine/app_utils.h"
#include "application/water_alternate_fine/data_app.h"
#include "application/water_alternate_fine/job_adjust_phi.h"
#include "application/water_alternate_fine/job_calculate_frame.h"
#include "application/water_alternate_fine/job_delete_particles.h"
#include "application/water_alternate_fine/job_initialize.h"
#include "application/water_alternate_fine/job_loop_frame.h"
#include "application/water_alternate_fine/job_loop_iteration.h"
#include "application/water_alternate_fine/job_main.h"
#include "application/water_alternate_fine/job_names.h"
#include "application/water_alternate_fine/job_modify_levelset.h"
#include "application/water_alternate_fine/job_reincorporate_removed_particles.h"
#include "application/water_alternate_fine/job_super_1.h"
#include "application/water_alternate_fine/job_super_2.h"
#include "application/water_alternate_fine/job_super_3.h"
#include "application/water_alternate_fine/job_adjust_phi_with_objects.h"
#include "application/water_alternate_fine/job_advect_phi.h"
#include "application/water_alternate_fine/job_step_particles.h"
#include "application/water_alternate_fine/job_advect_removed_particles.h"
#include "application/water_alternate_fine/job_advect_v.h"
#include "application/water_alternate_fine/job_apply_forces.h"
#include "application/water_alternate_fine/job_write_frame.h"
#include "application/water_alternate_fine/water_app.h"
#include "data/scalar_data.h"
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
        dbg(APP_LOG, "Registering %s\n", APP_FACE_VEL_GHOST);
        RegisterData(APP_FACE_VEL_GHOST, new DataApp(APP_FACE_VEL_GHOST, kFaceVelGhostBufSize));
        dbg(APP_LOG, "Registering %s\n", APP_PHI);
        RegisterData(APP_PHI, new DataApp(APP_PHI, kPhiBufSize));
        dbg(APP_LOG, "Registering %s\n", APP_PRESSURE);
        RegisterData(APP_PRESSURE, new DataApp(APP_PRESSURE, kPressureBufSize));
        dbg(APP_LOG, "Registering %s\n", APP_POS_PARTICLES);
        RegisterData(APP_POS_PARTICLES, new DataApp(APP_POS_PARTICLES, kParticlesBufSize));
        dbg(APP_LOG, "Registering %s\n", APP_NEG_PARTICLES);
        RegisterData(APP_NEG_PARTICLES, new DataApp(APP_NEG_PARTICLES, kParticlesBufSize));
        dbg(APP_LOG, "Registering %s\n", APP_POS_REM_PARTICLES);
        RegisterData(APP_POS_REM_PARTICLES, new DataApp(APP_POS_REM_PARTICLES, kParticlesBufSize));
        dbg(APP_LOG, "Registering %s\n", APP_NEG_REM_PARTICLES);
        RegisterData(APP_NEG_REM_PARTICLES, new DataApp(APP_NEG_REM_PARTICLES, kParticlesBufSize));
        dbg(APP_LOG, "Registering %s\n", APP_LAST_UNIQUE_PARTICLE_ID);
        RegisterData(APP_LAST_UNIQUE_PARTICLE_ID, new nimbus::ScalarData<int>(APP_LAST_UNIQUE_PARTICLE_ID));

        RegisterJob(MAIN, new JobMain(this));
        RegisterJob(INITIALIZE, new JobInitialize(this));
        RegisterJob(SUPER_1, new JobSuper1(this));
        RegisterJob(SUPER_2, new JobSuper2(this));
        RegisterJob(SUPER_3, new JobSuper3(this));
        RegisterJob(ADJUST_PHI_WITH_OBJECTS, new JobAdjustPhiWithObjects(this));
        RegisterJob(ADVECT_PHI, new JobAdvectPhi(this));
        RegisterJob(STEP_PARTICLES, new JobStepParticles(this));
        RegisterJob(ADVECT_REMOVED_PARTICLES, new JobAdvectRemovedParticles(this));
        RegisterJob(ADVECT_V, new JobAdvectV(this));
        RegisterJob(APPLY_FORCES, new JobApplyForces(this));
        RegisterJob(LOOP_ITERATION, new JobLoopIteration(this));
        RegisterJob(LOOP_FRAME, new JobLoopFrame(this));
        RegisterJob(CALCULATE_FRAME, new JobCalculateFrame(this));
        RegisterJob(WRITE_FRAME, new JobWriteFrame(this));
        RegisterJob(MODIFY_LEVELSET, new JobModifyLevelset(this));
        RegisterJob(ADJUST_PHI, new JobAdjustPhi(this));
        RegisterJob(DELETE_PARTICLES, new JobDeleteParticles(this));
        RegisterJob(REINCORPORATE_PARTICLES, new JobReincorporateRemovedParticles(this));

        dbg(APP_LOG, "Completed loading water application\n");
    }

} // namespace application

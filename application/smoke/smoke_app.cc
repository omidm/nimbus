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

#include "application/smoke/app_utils.h"
#include "application/smoke/cache_data_include.h"
#include "application/smoke/cache_prototypes.h"
#include "application/smoke/data_include.h"
#include "application/smoke/job_include.h"
#include "application/smoke/reg_def.h"
#include "application/smoke/smoke_app.h"
#include "data/physbam/translator_physbam.h"
#include "data/physbam/translator_physbam_old.h"
#include "data/scalar_data.h"
#include "data/scratch_data_helper.h"
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "shared/log.h"
#include "shared/nimbus.h"
#include "shared/timer.h"
#include "stdio.h"

namespace application {

    SmokeApp::SmokeApp() {};

    /* Register data and job types and initialize constant quantities used by
     * application jobs. */
    void SmokeApp::Load() {
        //nimbus::Timer::Initialize();

        dbg_add_mode(APP_LOG_STR);
        dbg_add_mode(TRANSLATE_STR);

        dbg(APP_LOG, "Loading smoke application\n");

        // Region constants
        InitializeRegions();

        // PhysBAM logging and R/W
        PhysBAM::LOG::Initialize_Logging(false, false, 1<<30, true, kThreadsNum);
        PhysBAM::FILE_UTILITIES::Create_Directory(kOutputDir+"/common");
        PhysBAM::LOG::Instance()->Copy_Log_To_File(kOutputDir+"/common/log.txt", false);

        dbg(APP_LOG, "Registering %s\n", APP_FACE_VEL);
        RegisterData(APP_FACE_VEL, new DataFaceArray<float>(APP_FACE_VEL));
        dbg(APP_LOG, "Registering %s\n", APP_FACE_VEL_GHOST);
        RegisterData(APP_FACE_VEL_GHOST, new DataFaceArray<float>(APP_FACE_VEL_GHOST));
	dbg(APP_LOG, "Registering %s\n", APP_DENSITY);
        RegisterData(APP_DENSITY, new DataScalarArray<float>(APP_DENSITY));
	dbg(APP_LOG, "Registering %s\n", APP_DENSITY_GHOST);
	RegisterData(APP_DENSITY_GHOST, new DataScalarArray<float>(APP_DENSITY_GHOST));
	dbg(APP_LOG, "Registering %s\n", APP_DT);
        RegisterData(APP_DT, new nimbus::ScalarData<float>(APP_DT));

        // These Nimbus data types are used in projection but not used in
        // internal projection loop. They are generally used to generate the
        // matrixes or vectors used in internal projection loop.
        // PSI_D.
        dbg(APP_LOG, "Registering %s\n", APP_PSI_D);
        RegisterData(APP_PSI_D, new DataScalarArray<bool>(APP_PSI_D));
        // PSI_N.
        dbg(APP_LOG, "Registering %s\n", APP_PSI_N);
        RegisterData(APP_PSI_N, new DataFaceArray<bool>(APP_PSI_N));
        // PRESSURE.
        dbg(APP_LOG, "Registering %s\n", APP_PRESSURE);
        RegisterData(APP_PRESSURE, new DataScalarArray<float>(APP_PRESSURE));
        // FILLED_REGION_COLORS.
        dbg(APP_LOG, "Registering %s\n", APP_FILLED_REGION_COLORS);
        RegisterData(APP_FILLED_REGION_COLORS,
            new DataScalarArray<int>(APP_FILLED_REGION_COLORS));
        // DIVERGENCE.
        dbg(APP_LOG, "Registering %s\n", APP_DIVERGENCE);
        RegisterData(APP_DIVERGENCE, new DataScalarArray<float>(APP_DIVERGENCE));
        // U_INTERFACE.
        dbg(APP_LOG, "Registering %s\n", APP_U_INTERFACE);
        RegisterData(APP_U_INTERFACE, new DataFaceArray<float>(APP_U_INTERFACE));

        // These Nimbus data types are used in internal projection loop. They
        // are derived from boundary conditions and act as the linkage between
        // inside projection and outside projection.
        // MATRIX_A.
        dbg(APP_LOG, "Registering %s\n", APP_MATRIX_A);
        RegisterData(APP_MATRIX_A, new DataSparseMatrix(APP_MATRIX_A));
        // VECTOR_B.
        dbg(APP_LOG, "Registering %s\n", APP_VECTOR_B);
        RegisterData(APP_VECTOR_B, new DataRawVectorNd(APP_VECTOR_B));
        // VECTOR_X.
        // dbg(APP_LOG, "Registering %s\n", APP_VECTOR_X);
        // RegisterData(APP_VECTOR_X, new DataRawVectorNd(APP_VECTOR_X));
        // INDEX_C2M.
        dbg(APP_LOG, "Registering %s\n", APP_INDEX_C2M);
        RegisterData(APP_INDEX_C2M, new DataRawGridArray(APP_INDEX_C2M));
        // INDEX_M2C.
        dbg(APP_LOG, "Registering %s\n", APP_INDEX_M2C);
        RegisterData(APP_INDEX_M2C, new DataRawArrayM2C(APP_INDEX_M2C));
        // PROJECTION_LOCAL_N.
        dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_LOCAL_N);
        RegisterData(APP_PROJECTION_LOCAL_N,
            new nimbus::ScalarData<int>(APP_PROJECTION_LOCAL_N));
        // PROJECTION_INTERIOR_N.
        dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_INTERIOR_N);
        RegisterData(APP_PROJECTION_INTERIOR_N,
            new nimbus::ScalarData<int>(APP_PROJECTION_INTERIOR_N));

        // These Nimbus data types are used only in internal projection loop.
        // They are mostly static configurations set for projeciton.
        // PROJECTION_LOCAL_TOLERANCE.
        dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_LOCAL_TOLERANCE);
        RegisterData(APP_PROJECTION_LOCAL_TOLERANCE,
            new nimbus::ScalarData<float>(APP_PROJECTION_LOCAL_TOLERANCE));
        // PROJECTION_GLOBAL_TOLERANCE.
        dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_GLOBAL_TOLERANCE);
        RegisterData(APP_PROJECTION_GLOBAL_TOLERANCE,
            new nimbus::ScalarData<float>(APP_PROJECTION_GLOBAL_TOLERANCE));
        // PROJECTION_GLOBAL_N.
        dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_GLOBAL_N);
        RegisterData(APP_PROJECTION_GLOBAL_N,
            new nimbus::ScalarData<int>(APP_PROJECTION_GLOBAL_N));
        // PROJECTION_DESIRED_N.
        dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_DESIRED_ITERATIONS);
        RegisterData(APP_PROJECTION_DESIRED_ITERATIONS,
            new nimbus::ScalarData<int>(APP_PROJECTION_DESIRED_ITERATIONS));

        // These Nimbus data types are used in internal projection loop.
        // PROJECTION_LOCAL_RESIDUAL.
        dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_LOCAL_RESIDUAL);
        RegisterData(APP_PROJECTION_LOCAL_RESIDUAL,
            new nimbus::ScalarData<double>(APP_PROJECTION_LOCAL_RESIDUAL));
        // PROJECTION_LOCAL_RHO.
        dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_LOCAL_RHO);
        RegisterData(APP_PROJECTION_LOCAL_RHO,
            new nimbus::ScalarData<double>(APP_PROJECTION_LOCAL_RHO));
        // PROJECTION_GLOBAL_RHO.
        dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_GLOBAL_RHO);
        RegisterData(APP_PROJECTION_GLOBAL_RHO,
            new nimbus::ScalarData<double>(APP_PROJECTION_GLOBAL_RHO));
        // PROJECTION_GLOBAL_RHO_OLD.
        dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_GLOBAL_RHO_OLD);
        RegisterData(APP_PROJECTION_GLOBAL_RHO_OLD,
            new nimbus::ScalarData<double>(APP_PROJECTION_GLOBAL_RHO_OLD));
        // PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA.
        dbg(APP_LOG, "Registering %s\n",
            APP_PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA);
        RegisterData(APP_PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA,
            new nimbus::ScalarData<double>(
                APP_PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA));
        // PROJECTION_ALPHA.
        dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_ALPHA);
        RegisterData(APP_PROJECTION_ALPHA,
            new nimbus::ScalarData<float>(APP_PROJECTION_ALPHA));
        // PROJECTION_BETA.
        dbg(APP_LOG, "Registering %s\n", APP_PROJECTION_BETA);
        RegisterData(APP_PROJECTION_BETA,
            new nimbus::ScalarData<float>(APP_PROJECTION_BETA));
        // MATRIX_C.
        dbg(APP_LOG, "Registering %s\n", APP_MATRIX_C);
        RegisterData(APP_MATRIX_C, new DataSparseMatrix(APP_MATRIX_C));
        // VECTOR_Z.
        dbg(APP_LOG, "Registering %s\n", APP_VECTOR_Z);
        RegisterData(APP_VECTOR_Z, new DataRawVectorNd(APP_VECTOR_Z));
        // VECTOR_P.
        dbg(APP_LOG, "Registering %s\n", APP_VECTOR_P);
        RegisterData(APP_VECTOR_P, new DataScalarArray<float>(APP_VECTOR_P));
        // VECTOR_TEMP.
        dbg(APP_LOG, "Registering %s\n", APP_VECTOR_TEMP);
        RegisterData(APP_VECTOR_TEMP, new DataRawVectorNd(APP_VECTOR_TEMP));


	/*
        dbg(APP_LOG, "Registering scratch %s\n", APP_POS_PARTICLES);
        kScratchPosParticles.RegisterScratchNames(this, new DataParticleArray(APP_POS_PARTICLES));
        dbg(APP_LOG, "Registering scratch %s\n", APP_NEG_PARTICLES);
        kScratchNegParticles.RegisterScratchNames(this, new DataParticleArray(APP_NEG_PARTICLES));
        dbg(APP_LOG, "Registering scratch %s\n", APP_POS_REM_PARTICLES);
        kScratchPosRemParticles.RegisterScratchNames(this, new DataParticleArray(APP_POS_REM_PARTICLES));
        dbg(APP_LOG, "Registering scratch %s\n", APP_NEG_REM_PARTICLES);
        kScratchNegRemParticles.RegisterScratchNames(this, new DataParticleArray(APP_NEG_REM_PARTICLES));
	*/

        RegisterJob(MAIN, new JobMain(this));
        RegisterJob(INITIALIZE, new JobInitialize(this));
	
	RegisterJob(SUBSTEP, new JobSubstep(this));
	RegisterJob(SCALAR_ADVANCE, new JobScalarAdvance(this));
	RegisterJob(CONVECT, new JobConvect(this));

        RegisterJob(LOOP_ITERATION, new JobLoopIteration(this));
        RegisterJob(LOOP_ITERATION_PART_TWO, new JobLoopIterationPartTwo(this));
        RegisterJob(LOOP_FRAME, new JobLoopFrame(this));
        RegisterJob(WRITE_OUTPUT, new JobWriteOutput(this));

        RegisterJob(PROJECTION_MAIN, new JobProjectionMain(this));
        RegisterJob(PROJECTION_CALCULATE_BOUNDARY_CONDITION_PART_ONE,
                    new JobProjectionCalculateBoundaryConditionPartOne(this));
        RegisterJob(PROJECTION_CALCULATE_BOUNDARY_CONDITION_PART_TWO,
                    new JobProjectionCalculateBoundaryConditionPartTwo(this));
        RegisterJob(PROJECTION_CONSTRUCT_MATRIX,
                    new JobProjectionConstructMatrix(this));
        RegisterJob(PROJECTION_WRAPUP, new JobProjectionWrapup(this));
        
        RegisterJob(PROJECTION_GLOBAL_INITIALIZE,
                    new JobProjectionGlobalInitialize(this));
        RegisterJob(PROJECTION_LOCAL_INITIALIZE,
                    new JobProjectionLocalInitialize(this));
        RegisterJob(PROJECTION_LOOP_ITERATION, new JobProjectionLoopIteration(this));
        RegisterJob(PROJECTION_STEP_ONE, new JobProjectionStepOne(this));
        RegisterJob(PROJECTION_REDUCE_RHO, new JobProjectionReduceRho(this));
        RegisterJob(PROJECTION_STEP_TWO, new JobProjectionStepTwo(this));
        RegisterJob(PROJECTION_STEP_THREE, new JobProjectionStepThree(this));
        RegisterJob(PROJECTION_REDUCE_ALPHA,
                    new JobProjectionReduceAlpha(this));
        RegisterJob(PROJECTION_STEP_FOUR, new JobProjectionStepFour(this));

        nimbus::TranslatorPhysBAM<float>::log = translator_log;
        nimbus::TranslatorPhysBAMOld<PhysBAM::VECTOR<float, 3> >::log = translator_log;
        dbg(APP_LOG, "Completed loading smoke application\n");
    }

} // namespace application

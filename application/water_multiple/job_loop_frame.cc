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
 * This file contains a loop frame job that spawns loop iteration jobs to
 * calculate a frame of the simulation. It checks whether the last required
 * frame has compututed or not, if so it terminates the application.
 *
 * Author: Omid Mashayekhi <omidm@stanford.edu>
 */


#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/job_loop_frame.h"
#include "application/water_multiple/water_driver.h"
#include "application/water_multiple/water_example.h"
#include "application/water_multiple/job_names.h"
#include "application/water_multiple/data_names.h"
#include "application/water_multiple/reg_def.h"
#include "data/physbam/physbam_data.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"
#include <sstream>
#include <string>

namespace application {

    JobLoopFrame::JobLoopFrame(nimbus::Application *app) {
        set_application(app);
    };

    nimbus::Job* JobLoopFrame::Clone() {
        return new JobLoopFrame(application());
    }

    void JobLoopFrame::Execute(nimbus::Parameter params, const nimbus::DataArray& da) {
        dbg(APP_LOG, "Executing loop frame job.\n");

        // get frame from parameters
        int frame;
        GeometricRegion global_region;
        std::string params_str(params.ser_data().data_ptr_raw(),
                               params.ser_data().size());
        LoadParameter(params_str, &frame, &global_region);

        // get time from frame
        PhysBAM::WATER_EXAMPLE<TV>* example =
          new PhysBAM::WATER_EXAMPLE<TV>(PhysBAM::STREAM_TYPE((RW())));
        T time = example->Time_At_Frame(frame);
        delete example;

        if (frame < kLastFrame) {
          //Spawn the loop iteration job to start computing the frame.
          dbg(APP_LOG, "Loop frame is spawning loop iteration job for frame %i.\n", frame);

          int job_num = 1;
          std::vector<nimbus::job_id_t> job_ids;
          GetNewJobID(&job_ids, job_num);
          nimbus::IDSet<nimbus::logical_data_id_t> read, write;
          nimbus::IDSet<nimbus::job_id_t> before, after;
          nimbus::Parameter iter_params;

          std::string str;
          SerializeParameter(frame, time, global_region, &str);
          iter_params.set_ser_data(SerializedData(str));

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


          SpawnComputeJob(LOOP_ITERATION,
              job_ids[0],
              read, write,
              before, after,
              iter_params);
        } else {
          // Last job has been computed, just terminate the application.
          dbg(APP_LOG, "Simulation is complete, last frame computed.\n");

          TerminateApplication();
        }
    }

} // namespace application

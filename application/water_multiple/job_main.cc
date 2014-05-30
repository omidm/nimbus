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



//        // Job setup
//        int job_num = 1;
//        std::vector<nimbus::job_id_t> job_ids;
//        GetNewJobID(&job_ids, job_num);
//        nimbus::IDSet<nimbus::logical_data_id_t> read, write;
//        nimbus::IDSet<nimbus::job_id_t> before, after;
//        nimbus::IDSet<nimbus::param_id_t> id_set;
//        nimbus::Parameter init_params;
//
//        // Init job
//        read.clear();
//        LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL, APP_FACE_VEL_GHOST, APP_PHI, NULL);
//        LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_POS_PARTICLES,
//            APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES,
//            APP_LAST_UNIQUE_PARTICLE_ID , NULL);
//        LoadLogicalIdsInSet(this, &read, kRegW1Outer[0], APP_PSI_D, APP_PSI_N, NULL);
//        write.clear();
//        LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_FACE_VEL, APP_FACE_VEL_GHOST, APP_PHI, NULL);
//        LoadLogicalIdsInSet(this, &write, kRegW3Outer[0], APP_POS_PARTICLES,
//            APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES,
//            APP_LAST_UNIQUE_PARTICLE_ID , NULL);
//        LoadLogicalIdsInSet(this, &write, kRegW1Outer[0], APP_PRESSURE, NULL);
//        id_set.clear();
//        kScratchPosParticles.GetJobScratchData(this, kRegW3Central[0], &id_set);
//        kScratchNegParticles.GetJobScratchData(this, kRegW3Central[0], &id_set);
//        kScratchPosRemParticles.GetJobScratchData(this, kRegW3Central[0], &id_set);
//        kScratchNegRemParticles.GetJobScratchData(this, kRegW3Central[0], &id_set);
//
//        init_params.set_ser_data(SerializedData(""));
//        dbg(APP_LOG, "Spawning initialize\n");
//        SpawnComputeJob(INITIALIZE,
//                        job_ids[0],
//                        read, write,
//                        before, after,
//                        init_params);




        // Job setup
        int init_job_num = kAppPartNum;
        std::vector<nimbus::job_id_t> init_job_ids;
        GetNewJobID(&init_job_ids, init_job_num);

        int make_signed_distance_job_num = kAppPartNum;
        std::vector<nimbus::job_id_t> make_signed_distance_job_ids;
        GetNewJobID(&make_signed_distance_job_ids, make_signed_distance_job_num);

        int extrapolate_phi_job_num = kAppPartNum;
        std::vector<nimbus::job_id_t> extrapolate_phi_job_ids;
        GetNewJobID(&extrapolate_phi_job_ids, extrapolate_phi_job_num);

        int extrapolate_v_job_num = kAppPartNum;
        std::vector<nimbus::job_id_t> extrapolate_v_job_ids;
        GetNewJobID(&extrapolate_v_job_ids, extrapolate_v_job_num);

        int write_frame_job_num = 1;
        std::vector<nimbus::job_id_t> write_frame_job_ids;
        GetNewJobID(&write_frame_job_ids, write_frame_job_num);

        int loop_frame_job_num = 1;
        std::vector<nimbus::job_id_t> loop_frame_job_ids;
        GetNewJobID(&loop_frame_job_ids, loop_frame_job_num);

        nimbus::IDSet<nimbus::logical_data_id_t> read, write;
        nimbus::IDSet<nimbus::job_id_t> before, after;

        int frame = 0;
        T time = 0;
        T dt = 0;

        /*
         * Spawning initialize stage over multiple workers
         */
        for (int i = 0; i < init_job_num; ++i) {
          read.clear();
          LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL, APP_FACE_VEL_GHOST, APP_PHI, NULL);
          LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_POS_PARTICLES,
              APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES,
              APP_LAST_UNIQUE_PARTICLE_ID , NULL);
          LoadLogicalIdsInSet(this, &read, kRegY2W1Outer[i], APP_PSI_D, APP_PSI_N, NULL);

          write.clear();
          LoadLogicalIdsInSet(this, &write, kRegY2W3Central[i], APP_FACE_VEL, APP_FACE_VEL_GHOST, APP_PHI, NULL);
          LoadLogicalIdsInSet(this, &write, kRegY2W3Central[i], APP_POS_PARTICLES,
              APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES,
              APP_LAST_UNIQUE_PARTICLE_ID , NULL);
          LoadLogicalIdsInSet(this, &write, kRegY2W1Central[i], APP_PRESSURE, NULL);

          before.clear();

          nimbus::Parameter init_params;
          std::string init_str;
          SerializeParameter(frame, time, dt, kDefaultRegion, kRegY2W3Central[i], &init_str);
          init_params.set_ser_data(SerializedData(init_str));

          dbg(APP_LOG, "Spawning initialize\n");
          SpawnComputeJob(INITIALIZE,
              init_job_ids[i],
              read, write,
              before, after,
              init_params, true);
        }


        /* 
         * Spawning make signed distance.
         */
        for (int i = 0; i < make_signed_distance_job_num; ++i) {
          read.clear();
          LoadLogicalIdsInSet(this, &read, kRegY2W7Outer[i], APP_PHI, NULL);
          LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL_GHOST,
              APP_FACE_VEL, NULL);
          LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_POS_PARTICLES,
              APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);

          write.clear();
          LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[i], APP_PHI, NULL);
          LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[i], APP_POS_PARTICLES,
              APP_NEG_PARTICLES, APP_POS_REM_PARTICLES, APP_NEG_REM_PARTICLES, NULL);

          before.clear();
          for (int j = 0; j < init_job_num; ++j) {
            before.insert(init_job_ids[j]);
          }

          std::string make_signed_distance_str;
          SerializeParameter(frame, time, dt, kDefaultRegion, kRegY2W3Central[i], &make_signed_distance_str);
          nimbus::Parameter make_signed_distance_params;
          make_signed_distance_params.set_ser_data(SerializedData(make_signed_distance_str));

          dbg(APP_LOG, "Spawning initialize\n");
          SpawnComputeJob(MAKE_SIGNED_DISTANCE,
              make_signed_distance_job_ids[i],
              read, write,
              before, after,
              make_signed_distance_params, true);
        }


        /*
         * Spawning extrapolate phi stage over multiple workers.
         */
        for (int i = 0; i < extrapolate_phi_job_num; ++i) {
          read.clear();
          LoadLogicalIdsInSet(this, &read, kRegY2W8Outer[i], APP_PHI,
              APP_FACE_VEL, NULL);

          write.clear();
          LoadLogicalIdsInSet(this, &write,
              kRegY2W8CentralWGB[i], APP_PHI, NULL);

          before.clear();
          for (int j = 0; j < make_signed_distance_job_num; ++j) {
            before.insert(make_signed_distance_job_ids[j]);
          }

          nimbus::Parameter phi_params;
          std::string phi_str;
          SerializeParameter(frame, time, dt, kDefaultRegion, kRegY2W8Central[i], &phi_str);
          phi_params.set_ser_data(SerializedData(phi_str));

          dbg(APP_LOG, "Spawning Extrapolate Phi\n");
          SpawnComputeJob(EXTRAPOLATE_PHI,
              extrapolate_phi_job_ids[i],
              read, write,
              before, after,
              phi_params, true);
        }


        /*
         * Spawning extrapolate v stage over multiple workers
         */
        for (int i = 0; i < extrapolate_v_job_num; ++i) {
          read.clear();
          LoadLogicalIdsInSet(this, &read, kRegY2W8Outer[i],
              APP_FACE_VEL, APP_PHI, NULL);

          write.clear();
          LoadLogicalIdsInSet(this, &write, kRegY2W8Central[i],
              APP_FACE_VEL, NULL);

          before.clear();
          for (int j = 0; j < extrapolate_phi_job_num; ++j) {
            before.insert(extrapolate_phi_job_ids[j]);
          }

          nimbus::Parameter v_params;
          std::string v_str;
          SerializeParameter(frame, time, dt, kDefaultRegion, kRegY2W3Central[i], &v_str);
          v_params.set_ser_data(SerializedData(v_str));

          dbg(APP_LOG, "Spawning Extrapolate V\n");
          SpawnComputeJob(EXTRAPOLATION,
              extrapolate_v_job_ids[i],
              read, write,
              before, after,
              v_params, true);
        }



        /*
         * Spawning write frame over entire block.
         */
        read.clear();
        LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL,
            APP_FACE_VEL_GHOST, APP_PHI, NULL);
        LoadLogicalIdsInSet(this, &read, kRegW1Outer[0], APP_PSI_D,
            APP_PSI_N, NULL);
        LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_POS_PARTICLES,
            APP_NEG_PARTICLES, APP_POS_REM_PARTICLES,
            APP_NEG_REM_PARTICLES, APP_LAST_UNIQUE_PARTICLE_ID,
            NULL);

        write.clear();

          before.clear();
          for (int j = 0; j < extrapolate_v_job_num; ++j) {
            before.insert(extrapolate_v_job_ids[j]);
          }

        nimbus::Parameter write_params;
        std::string write_str;
        int write_frame = frame - 1;
        SerializeParameter(write_frame, time + dt, 0, kDefaultRegion, kDefaultRegion, &write_str);
        write_params.set_ser_data(SerializedData(write_str));

          dbg(APP_LOG, "Spawning Write Frame\n");
          SpawnComputeJob(WRITE_FRAME,
              write_frame_job_ids[0],
              read, write,
              before, after,
              write_params, true);



        /*
         * Spawning loop frame job.
         */
        read.clear();
        write.clear();

        before.clear();
        for (int j = 0; j < write_frame_job_num; ++j) {
          before.insert(write_frame_job_ids[j]);
        }

        nimbus::Parameter loop_params;
        std::string loop_str;
        SerializeParameter(frame, kDefaultRegion, &loop_str);
        loop_params.set_ser_data(SerializedData(loop_str));

        dbg(APP_LOG, "Spawning loop frame job for frame %i in main\n", frame);
        SpawnComputeJob(LOOP_FRAME,
            loop_frame_job_ids[0],
            read, write,
            before, after,
            loop_params);






        dbg(APP_LOG, "Completed executing main job\n");
    }

} // namespace application

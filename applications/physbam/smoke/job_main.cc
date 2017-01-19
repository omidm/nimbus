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
 * This file contains the "main" job that Nimbus launches after loading an
 * application. All subsequent jobs are spawned from here.
 *
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 * Modifier for smoke: Andrew Lim <alim16@stanford.edu> 
 */

#include "applications/physbam/smoke/app_utils.h"
#include "applications/physbam/smoke/data_def.h"
#include "applications/physbam/smoke/data_names.h"
#include "applications/physbam/smoke/job_main.h"
#include "applications/physbam/smoke/job_names.h"
#include "applications/physbam/smoke/reg_def.h"
#include "src/data/scratch_data_helper.h"
#include "src/shared/dbg.h"
#include "src/shared/nimbus.h"
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
  int init_job_num = kAppPartNum;
  std::vector<nimbus::job_id_t> init_job_ids;
  GetNewJobID(&init_job_ids, init_job_num);

  int write_output_job_num = kAppPartNum;
  std::vector<nimbus::job_id_t> write_output_job_ids;
  GetNewJobID(&write_output_job_ids, write_output_job_num);

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
    LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL, APP_FACE_VEL_GHOST, 
			APP_DENSITY, APP_DENSITY_GHOST, NULL);
    LoadLogicalIdsInSet(this, &read, kRegY2W1Outer[i], APP_PSI_D, APP_PSI_N, NULL);

    write.clear();
    LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[i], APP_FACE_VEL, APP_FACE_VEL_GHOST, 
			APP_DENSITY, APP_DENSITY_GHOST, NULL);
    LoadLogicalIdsInSet(this, &write, kRegY2W1CentralWGB[i],
                        APP_PRESSURE,APP_PSI_D, APP_PSI_N,  NULL);

    nimbus::Parameter init_params;
    std::string init_str;
    SerializeParameter(frame, time, dt, kDefaultRegion, kRegY2W3Central[i], &init_str);
    init_params.set_ser_data(SerializedData(init_str));

    after.clear();
    before.clear();
    StageJobAndLoadBeforeSet(&before, INITIALIZE,
                             init_job_ids[i],
                             read, write);
    SpawnComputeJob(INITIALIZE,
                    init_job_ids[i],
                    read, write, before, after,
                    init_params, true,
                    kRegY2W3Central[i]);
  }

  MarkEndOfStage();



  if (kUseGlobalWrite) {
    read.clear();
    LoadLogicalIdsInSet(this, &read, kRegW3Outer[0], APP_FACE_VEL, APP_DENSITY, NULL);
    LoadLogicalIdsInSet(this, &read, kRegW1Outer[0], APP_PSI_D,
                        APP_PSI_N, NULL);
    write.clear();

    nimbus::Parameter temp_params;
    std::string temp_str;
    SerializeParameter(frame - 1, time + dt, 0, -1, kDefaultRegion, kDefaultRegion,
                       &temp_str);
    temp_params.set_ser_data(SerializedData(temp_str));

    after.clear();
    before.clear();
    StageJobAndLoadBeforeSet(&before, WRITE_OUTPUT,
                             write_output_job_ids[0],
                             read, write);
    SpawnComputeJob(WRITE_OUTPUT,
                    write_output_job_ids[0],
                    read, write, before, after,
                    temp_params, true,
                    kRegW3Central[0]);
  } else {
    for (int i = 0; i < write_output_job_num; ++i) {
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL, APP_DENSITY, NULL);
      LoadLogicalIdsInSet(this, &read, kRegY2W1Outer[i], APP_PSI_D,
                          APP_PSI_N, NULL);
      write.clear();

      nimbus::Parameter temp_params;
      std::string temp_str;
      SerializeParameter(frame - 1, time + dt, 0, i+1, kDefaultRegion, kRegY2W3Central[i],
                         &temp_str);
      temp_params.set_ser_data(SerializedData(temp_str));

      after.clear();
      before.clear();
      StageJobAndLoadBeforeSet(&before, WRITE_OUTPUT,
                               write_output_job_ids[i],
                               read, write);
      SpawnComputeJob(WRITE_OUTPUT,
                      write_output_job_ids[i],
                      read, write, before, after,
                      temp_params, true,
                      kRegY2W3Central[i]);
    }
  }

  MarkEndOfStage();

  /*
   * Spawning loop frame job.
   */
  read.clear();
  write.clear();

  nimbus::Parameter loop_params;
  std::string loop_str;
  SerializeParameter(frame, kDefaultRegion, &loop_str);
  loop_params.set_ser_data(SerializedData(loop_str));



    after.clear();
    before.clear();
    StageJobAndLoadBeforeSet(&before, LOOP_FRAME,
                             loop_frame_job_ids[0],
                             read, write,
                             true);
    SpawnComputeJob(LOOP_FRAME,
                    loop_frame_job_ids[0],
                    read, write, before, after,
                    loop_params, false,
                    kRegW3Central[0]);

  MarkEndOfStage();

  dbg(APP_LOG, "Completed executing main job\n");

  // dbg(APP_LOG, "Print job dependency figure.\n");
  // job_query.GenerateDotFigure("job_main.dot");
}

} // namespace application

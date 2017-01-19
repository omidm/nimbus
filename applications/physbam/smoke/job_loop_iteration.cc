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
 * This file contains a loop iteration job that spawns the sub step jobs to
 * calculate the current frame. It keeps spawning the iteration in a loop as
 * long as frame computation in not complete. When the frame is done it will
 * spawn the loop frame job for the next frame. 
 *
 * Author: Omid Mashayekhi <omidm@stanford.edu>
 * Modifier for smoke: Andrew Lim <alim16@stanford.edu> 
 */

#include <sys/time.h>
#include "applications/physbam/smoke/app_utils.h"
#include "applications/physbam/smoke/job_loop_iteration.h"
#include "applications/physbam/smoke/physbam_utils.h"
#include "applications/physbam/smoke/smoke_driver.h"
#include "applications/physbam/smoke/smoke_example.h"
#include "applications/physbam/smoke/job_names.h"
#include "applications/physbam/smoke/data_names.h"
#include "applications/physbam/smoke/reg_def.h"
#include "src/shared/dbg.h"
#include "src/shared/nimbus.h"
#include "src/worker/dependency_query.h"
#include <sstream>
#include <string>

namespace application {

  JobLoopIteration::JobLoopIteration(nimbus::Application *app) {
    set_application(app);
  };

  nimbus::Job* JobLoopIteration::Clone() {
    return new JobLoopIteration(application());
  }

  void JobLoopIteration::Execute(
      nimbus::Parameter params,
      const nimbus::DataArray& da) {
    dbg(APP_LOG, "Executing loop iteration job\n");

    // Get parameters: frame, time
    InitConfig init_config;
    // Threading settings.
    init_config.use_threading = use_threading();
    init_config.core_quota = core_quota();
    std::string params_str(params.ser_data().data_ptr_raw(),
        params.ser_data().size());
    LoadParameter(params_str, &init_config.frame, &init_config.time,
                  &init_config.global_region);
    init_config.local_region = init_config.global_region;

    const int& frame = init_config.frame;
    const T& time = init_config.time;

    dbg(APP_LOG, "Frame %i and time %f in iteration job\n",
        frame, time);

    // Initialize the state of example and driver.
    PhysBAM::SMOKE_EXAMPLE<TV>* example = 
      new PhysBAM::SMOKE_EXAMPLE<TV>(PhysBAM::STREAM_TYPE((RW())), false, 1);
    DataConfig data_config;
    data_config.SetFlag(DataConfig::VELOCITY);
    data_config.SetFlag(DataConfig::DENSITY);
    data_config.SetFlag(DataConfig::DT);
    {
      application::ScopeTimer scope_timer(name() + "-load");
      example->data_config.Set(data_config);
      example->Load_From_Nimbus(this, da, frame);
    }
    // *thread_queue_hook() = example->nimbus_thread_queue;

    // check whether the frame is done or not
    bool done = false;
    T target_time = example->Time_At_Frame(frame + 1);
    T dt = example->dt_buffer;
    if (time + dt >= target_time) {
      dt = target_time - time;
      done = true;
    } else if (time + 2 * dt >= target_time) {
        dt = .5 * (target_time - time);
    }

    if (done) {
      dbg(APP_LOG, "[CONTROL FLOW] First part, Loop done.\n");
    } else {
      dbg(APP_LOG, "[CONTROL FLOW] First part, Loop not done.\n");
    }
    dbg(APP_LOG, "[CONTROL FLOW] First part, Frame=%d, Time=%f, dt=%f\n", frame, time, dt);

    // done flag no more matters.
    SpawnWithBreakAllGranularity(false, init_config.frame, init_config.time,
                                 dt, da, init_config.global_region);

    // Free resources.
    {
      application::ScopeTimer scope_timer(name() + "-save");
      example->Save_To_Nimbus(this, da, frame+1);
    }
    // *thread_queue_hook() = NULL;
  }

  void JobLoopIteration::SpawnWithBreakAllGranularity(
      bool done, int frame, T time, T dt, const nimbus::DataArray& da,
      const GeometricRegion& global_region) {
    dbg(APP_LOG, "Loop frame is spawning super job 1, 2, 3 for frame %i.\n", frame);

    nimbus::IDSet<nimbus::logical_data_id_t> read, write;
    nimbus::IDSet<nimbus::job_id_t> before, after;

    int update_ghost_densities_job_num = kAppPartNum;
    std::vector<nimbus::job_id_t> update_ghost_densities_job_ids;
    GetNewJobID(&update_ghost_densities_job_ids, update_ghost_densities_job_num);

    int scalar_advance_job_num = kAppPartNum;
    std::vector<nimbus::job_id_t> scalar_advance_job_ids;
    GetNewJobID(&scalar_advance_job_ids, scalar_advance_job_num);

    int update_ghost_velocities_job_num = kAppPartNum;
    std::vector<nimbus::job_id_t> update_ghost_velocities_job_ids;
    GetNewJobID(&update_ghost_velocities_job_ids, update_ghost_velocities_job_num);

    int convect_job_num = kAppPartNum;
    std::vector<nimbus::job_id_t> convect_job_ids;
    GetNewJobID(&convect_job_ids, convect_job_num);

    std::vector<nimbus::job_id_t> projection_main_job_ids;
    GetNewJobID(&projection_main_job_ids, 1);

    struct timeval start_time;
    gettimeofday(&start_time, NULL);

    for (int i = 0; i < update_ghost_densities_job_num; ++i) {
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_DENSITY, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[i], APP_DENSITY_GHOST, NULL);
      LoadLogicalIdsInSet(this, &write, kRegY2W3Central[i], APP_DENSITY, NULL);

      nimbus::Parameter update_ghost_densities_params;
      std::string update_ghost_densities_str;
      SerializeParameter(frame, time, dt, global_region, kRegY2W3Central[i], &update_ghost_densities_str);
      update_ghost_densities_params.set_ser_data(SerializedData(update_ghost_densities_str));

      after.clear();
      before.clear();
      StageJobAndLoadBeforeSet(&before, UPDATE_GHOST_DENSITIES,
                               update_ghost_densities_job_ids[i],
                               read, write);
      SpawnComputeJob(UPDATE_GHOST_DENSITIES,
                      update_ghost_densities_job_ids[i],
                      read, write, before, after,
                      update_ghost_densities_params, true,
                      kRegY2W3Central[i]);
    }      

    MarkEndOfStage();

    for (int i = 0; i < scalar_advance_job_num; ++i) {
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL, APP_DENSITY, APP_DENSITY_GHOST, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegY2W3Central[i], APP_DENSITY, NULL);

      nimbus::Parameter scalar_params;
      std::string scalar_str;
      SerializeParameter(frame, time, dt, global_region, kRegY2W3Central[i], &scalar_str);
      scalar_params.set_ser_data(SerializedData(scalar_str));


      after.clear();
      before.clear();
      StageJobAndLoadBeforeSet(&before, SCALAR_ADVANCE,
                               scalar_advance_job_ids[i],
                               read, write);
      SpawnComputeJob(SCALAR_ADVANCE,
                      scalar_advance_job_ids[i],
                      read, write, before, after,
                      scalar_params, true,
                      kRegY2W3Central[i]);
    }

    MarkEndOfStage();

    for (int i = 0; i < update_ghost_velocities_job_num; ++i) {
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[i], APP_FACE_VEL_GHOST, NULL);

      nimbus::Parameter update_ghost_velocities_params;
      std::string update_ghost_velocities_str;
      SerializeParameter(frame, time, dt, global_region, kRegY2W3Central[i], &update_ghost_velocities_str);
      update_ghost_velocities_params.set_ser_data(SerializedData(update_ghost_velocities_str));

      after.clear();
      before.clear();
      StageJobAndLoadBeforeSet(&before, UPDATE_GHOST_VELOCITIES,
                               update_ghost_velocities_job_ids[i],
                               read, write);
      SpawnComputeJob(UPDATE_GHOST_VELOCITIES,
                      update_ghost_velocities_job_ids[i],
                      read, write, before, after,
                      update_ghost_velocities_params, true,
                      kRegY2W3Central[i]);
    }

    MarkEndOfStage();

    for (int i = 0; i < convect_job_num; ++i) {
      read.clear();
      LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL_GHOST, NULL);
      LoadLogicalIdsInSet(this, &read, kRegY2W3Central[i], APP_FACE_VEL, NULL);
      LoadLogicalIdsInSet(this, &read, kRegY2W1Outer[i], APP_PSI_D, APP_PSI_N, NULL);
      write.clear();
      LoadLogicalIdsInSet(this, &write, kRegY2W3Central[i], APP_FACE_VEL, NULL);

      nimbus::Parameter convect_params;
      std::string convect_str;
      SerializeParameter(frame, time, dt, global_region, kRegY2W3Central[i],
          &convect_str);
      convect_params.set_ser_data(SerializedData(convect_str));

      after.clear();
      before.clear();
      StageJobAndLoadBeforeSet(&before, CONVECT,
                               convect_job_ids[i],
                               read, write);
      SpawnComputeJob(CONVECT,
                      convect_job_ids[i],
                      read, write, before, after,
                      convect_params, true,
                      kRegY2W3Central[i]);
    }

    MarkEndOfStage();

   /*
    * Projection main.
    */
    read.clear();
    write.clear();

    nimbus::Parameter projection_main_params;
    std::string projection_main_str;
    SerializeParameter(frame, time, dt, global_region, global_region,
        &projection_main_str);
    projection_main_params.set_ser_data(SerializedData(projection_main_str));

    after.clear();
    before.clear();
    StageJobAndLoadBeforeSet(&before, PROJECTION_MAIN,
                             projection_main_job_ids[0],
                             read, write,
                             true);
    SpawnComputeJob(PROJECTION_MAIN,
                    projection_main_job_ids[0],
                    read, write, before, after,
                    projection_main_params, false,
                    kRegW3Central[0]);
 

    MarkEndOfStage();

    // job_query.PrintTimeProfile();
    // {
    //   struct timeval t;
    //   gettimeofday(&t, NULL);
    //   double time  = (static_cast<double>(t.tv_sec - start_time.tv_sec)) +
    //     .000001 * (static_cast<double>(t.tv_usec - start_time.tv_usec));
    //   dbg(APP_LOG, "\nThe query time spent in job LOOP_ITERATION is %f seconds.\n",
    //       time);
    // }
    // if (time == 0) {
    //   dbg(APP_LOG, "Print job dependency figure.\n");
    //   job_query.GenerateDotFigure("loop_iteration.dot");
    // }
  }
} // namespace application

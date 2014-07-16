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
 * spawn the loop frame job for the next frame. The granularity of the sub step
 * jobs could be controlled by changing the changing the GRANULARITY_STATE
 * macro in this file as follows:
 * ONE_JOB:              calculate the frame iteration in one job (like water coarse).
 * SUPER_JOBS:           break the frame iteration in to three super jobs.
 * BREAK_SUPER_JOB_1:    further break super job 1 in to its components.
 * BREAK_SUPER_JOB_2:    further break super job 2 in to its components.
 * BREAK_SUPER_JOB_3:    further break super job 3 in to its components.
 * BREAK_ALL_SUPER_JOBS: break all three super jobs in to their components.
 *
 * Author: Omid Mashayekhi <omidm@stanford.edu>
 */

#define GRANULARITY_STATE BREAK_ALL_SUPER_JOBS

#include <sys/time.h>
#include "application/smoke/app_utils.h"
#include "application/smoke/job_loop_iteration.h"
#include "application/smoke/physbam_utils.h"
#include "application/smoke/smoke_driver.h"
#include "application/smoke/smoke_example.h"
#include "application/smoke/job_names.h"
#include "application/smoke/data_names.h"
#include "application/smoke/reg_def.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"
#include "worker/job_query.h"
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
    PhysBAM::SMOKE_EXAMPLE<TV>* example;
    PhysBAM::SMOKE_DRIVER<TV>* driver;
    init_config.set_boundary_condition = true;
    init_config.use_cache = true;
    DataConfig data_config;
    data_config.SetFlag(DataConfig::VELOCITY);
    data_config.SetFlag(DataConfig::DENSITY);
    data_config.SetFlag(DataConfig::DT);
    InitializeExampleAndDriver(init_config, data_config,
                               this, da, example, driver);
    *thread_queue_hook() = example->nimbus_thread_queue;

    // check whether the frame is done or not
    bool done = false;
    T target_time = example->Time_At_Frame(driver->current_frame+1);
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
    dbg(APP_LOG, "[CONTROL FLOW] First part, Frame=%d, Time=%f, dt=%f\n",
	frame, time, dt);

    // done flag no more matters.
    SpawnWithBreakAllGranularity(false, init_config.frame, init_config.time,
                                 dt, da, init_config.global_region);

    // Free resources.
    example->Save_To_Nimbus(this, da, frame+1);
    *thread_queue_hook() = NULL;
    DestroyExampleAndDriver(example, driver);
  }

  void JobLoopIteration::SpawnWithBreakAllGranularity(
      bool done, int frame, T time, T dt, const nimbus::DataArray& da,
      const GeometricRegion& global_region) {
    nimbus::JobQuery job_query(this);
    dbg(APP_LOG, "Loop frame is spawning super job 1, 2, 3 for frame %i.\n", frame);

    // these job ids are no longer used
    int job_num = 13;
    std::vector<nimbus::job_id_t> job_ids;
    GetNewJobID(&job_ids, job_num);
    nimbus::IDSet<nimbus::logical_data_id_t> read, write;

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

    int projection_job_num = 4;
    std::vector<nimbus::job_id_t> projection_job_ids;
    GetNewJobID(&projection_job_ids, projection_job_num);

    struct timeval start_time;
    gettimeofday(&start_time, NULL);

    {
      for (int i = 0; i < update_ghost_densities_job_num; ++i) {
	read.clear();
	LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_DENSITY, NULL);
	write.clear();
	LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[i], APP_DENSITY, APP_DENSITY_GHOST, NULL);
	// LoadLogicalIdsInSet(this, &write, kRegY2W3Central[i], APP_DENSITY, NULL);
	
	nimbus::Parameter update_ghost_densities_params;
	std::string update_ghost_densities_str;
	SerializeParameter(frame, time, dt, global_region, kRegY2W3Central[i], &update_ghost_densities_str);
	update_ghost_densities_params.set_ser_data(SerializedData(update_ghost_densities_str));
	job_query.StageJob(UPDATE_GHOST_DENSITIES,
			   update_ghost_densities_job_ids[i],
			   read, write,
			   update_ghost_densities_params, true);
	job_query.Hint(update_ghost_densities_job_ids[i],
		       kRegY2W3Central[i]);
      }      
      job_query.CommitStagedJobs();
    }

    {
      for (int i = 0; i < scalar_advance_job_num; ++i) {
	read.clear();
	LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL, APP_DENSITY, APP_DENSITY_GHOST, NULL);
	write.clear();
	LoadLogicalIdsInSet(this, &write, kRegY2W3Central[i], APP_DENSITY, NULL);

	nimbus::Parameter scalar_params;
	std::string scalar_str;
	SerializeParameter(frame, time, dt, global_region, kRegY2W3Central[i], &scalar_str);
	scalar_params.set_ser_data(SerializedData(scalar_str));
	job_query.StageJob(SCALAR_ADVANCE,
			   scalar_advance_job_ids[i],
			   read, write,
			   scalar_params, true);
	job_query.Hint(scalar_advance_job_ids[i],
		       kRegY2W3Central[i]);
      }
      job_query.CommitStagedJobs();
    }

    {
      for (int i = 0; i < update_ghost_velocities_job_num; ++i) {
	read.clear();
	LoadLogicalIdsInSet(this, &read, kRegY2W3Outer[i], APP_FACE_VEL, NULL);
	write.clear();
	LoadLogicalIdsInSet(this, &write, kRegY2W3CentralWGB[i], APP_FACE_VEL_GHOST, NULL);
	
	nimbus::Parameter update_ghost_velocities_params;
	std::string update_ghost_velocities_str;
	SerializeParameter(frame, time, dt, global_region, kRegY2W3Central[i], &update_ghost_velocities_str);
	update_ghost_velocities_params.set_ser_data(SerializedData(update_ghost_velocities_str));
	job_query.StageJob(UPDATE_GHOST_VELOCITIES,
			   update_ghost_velocities_job_ids[i],
			   read, write,
			   update_ghost_velocities_params, true);
	job_query.Hint(update_ghost_velocities_job_ids[i],
		       kRegY2W3Central[i]);
      }
      job_query.CommitStagedJobs();
    }

    {
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
	job_query.StageJob(CONVECT,
			   convect_job_ids[i],
			   read, write,
			   convect_params, true);
	job_query.Hint(convect_job_ids[i],
		       kRegY2W3Central[i]);
      }
      job_query.CommitStagedJobs();
    }


    /*
     * Loop iteration part two.
     */

    read.clear();
    write.clear();

    nimbus::Parameter projection_main_params;
    std::string projection_main_str;
    SerializeParameter(frame, time, dt, global_region, global_region,
                       &projection_main_str);
    projection_main_params.set_ser_data(
        SerializedData(projection_main_str));
    job_query.StageJob(PROJECTION_MAIN,
                    job_ids[10],
                    read, write,
                    projection_main_params,
                    false, true);
    job_query.CommitStagedJobs();
    job_query.PrintTimeProfile();
    {
      struct timeval t;
      gettimeofday(&t, NULL);
      double time  = (static_cast<double>(t.tv_sec - start_time.tv_sec)) +
          .000001 * (static_cast<double>(t.tv_usec - start_time.tv_usec));
      dbg(APP_LOG, "\nThe query time spent in job LOOP_ITERATION is %f seconds.\n",
          time);
    }
    if (time == 0) {
      dbg(APP_LOG, "Print job dependency figure.\n");
      job_query.GenerateDotFigure("loop_iteration.dot");
    }
  }
} // namespace application

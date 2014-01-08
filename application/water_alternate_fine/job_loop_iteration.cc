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
 */

#include "application/water_alternate_fine/app_utils.h"
#include "application/water_alternate_fine/job_loop_iteration.h"
#include "application/water_alternate_fine/job_loop_frame.h"
#include "application/water_alternate_fine/water_driver.h"
#include "application/water_alternate_fine/water_example.h"
#include "application/water_alternate_fine/water_sources.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"
#include <sstream>
#include <string>

namespace application {

    JobLoopIteration::JobLoopIteration(nimbus::Application *app) {
        set_application(app);
    };

    nimbus::Job* JobLoopIteration::Clone() {
        return new JobLoopIteration(application());
    }

    bool JobLoopIteration::InitializeExampleAndDriver(
        const nimbus::DataArray& da,
        const int current_frame,
        const T time,
        const int last_unique_particle_id,
        PhysBAM::WATER_EXAMPLE<TV>*& example,
        PhysBAM::WATER_DRIVER<TV>*& driver) {
      example = new PhysBAM::WATER_EXAMPLE<TV>(PhysBAM::STREAM_TYPE((RW())));
      example->Initialize_Grid(
          TV_INT::All_Ones_Vector()*kScale,
          PhysBAM::RANGE<TV>(TV(), TV::All_Ones_Vector()));
      PhysBAM::WaterSources::Add_Source(example);
      driver= new PhysBAM::WATER_DRIVER<TV>(*example);
      driver->init_phase = false;
      driver->current_frame = current_frame;
      // The returning result is deleted.
      driver->Initialize(this, da, last_unique_particle_id);
      driver->time = time;
      return true;
    }

    void JobLoopIteration::Execute(nimbus::Parameter params, const nimbus::DataArray& da) {
        dbg(APP_LOG, "Executing loop iteration job\n");

        // Get parameters: frame, time, last_unique_particle.
        int frame, last_unique_particle;
        std::stringstream in_ss;
        std::string params_str(params.ser_data().data_ptr_raw(),
                               params.ser_data().size());
        in_ss.str(params_str);
        in_ss >> frame;
        in_ss >> last_unique_particle;
        dbg(APP_LOG, "Frame %i, last unique particle %i in iteration job\n",
                     frame, last_unique_particle);

        // T time = example->Time_At_Frame(driver.current_frame);
        // TODO(quhang): time should be passed.
        T time = 0;
        bool done = true;
        PhysBAM::WATER_EXAMPLE<TV>* example;
        PhysBAM::WATER_DRIVER<TV>* driver;
        InitializeExampleAndDriver(da, frame, time, last_unique_particle,
                                   example, driver);

        T target_time = example->Time_At_Frame(driver->current_frame+1);
        T dt = example->cfl *
            example->incompressible.CFL(example->face_velocities);
        T temp_dt = example->particle_levelset_evolution.CFL(false,false);
        if (temp_dt < dt) {
          dt = temp_dt;
        }
        if (time + dt >= target_time) {
          dt = target_time - time;
          done = true;
        } else {
          if (time + 2*dt >= target_time)
            dt = .5 * (target_time-time);
        }

        // check whether the frame is done or not
        // TODO(omidm): get the right logic for done.

        if (!done) {
          // TODO(omidm): spawn the jobs to compute the frame, depending on the
          // level of granularity we will have different sub jobs.
          // TODO(quhang): Needs to pass "time","dt", "frame_number", and
          // "last_uniqute_particle_id" correctly.
        } else {
          // TODO(omidm): spawn the write frame job, or maybe compute frame job
          // for last time before write frame, and then spawn loop frame job
          // for next fram computation."
          // TODO(quhang): Needs to pass "time", "frame_number", and
          // "last_uniqute_particle_id" correctly.
        }
}

} // namespace application

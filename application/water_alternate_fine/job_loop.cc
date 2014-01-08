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
 * This file contains a loop job that spawns iteration jobs at a coarse
 * granularity.
 *
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#include "application/water_alternate_fine/app_utils.h"
#include "application/water_alternate_fine/job_iteration.h"
#include "application/water_alternate_fine/job_loop.h"
#include "data/physbam/physbam_data.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"
#include <sstream>
#include <string>

namespace application {

    JobLoop::JobLoop(nimbus::Application *app) {
        set_application(app);
    };

    nimbus::Job* JobLoop::Clone() {
        return new JobLoop(application());
    }

    void JobLoop::Execute(nimbus::Parameter params, const nimbus::DataArray& da) {
        dbg(APP_LOG, "Executing loop job\n");

        int frame;
        std::stringstream frame_ss;
        std::string params_str(params.ser_data().data_ptr_raw(),
                               params.ser_data().size());
        frame_ss.str(params_str);
        frame_ss >> frame;

        dbg(APP_LOG, "Loop is spawning iteration for %i frame\n", frame);

        if (frame < kLastFrame) {
            // Job setup
            int job_num = 1;
            std::vector<nimbus::job_id_t> job_ids;
            GetNewJobID(&job_ids, job_num);
            nimbus::IDSet<nimbus::logical_data_id_t> read, write;
            nimbus::IDSet<nimbus::job_id_t> before, after;

            // Iteration job
            for (size_t i = 0; i < da.size(); ++i) {
                nimbus::Data *d = da[i];
                logical_data_id_t id = d->logical_id();
                if (!application::Contains(read, id))
                    read.insert(id);
                if (!application::Contains(write, id))
                    write.insert(id);
            }
            nimbus::Parameter iter_params;
            iter_params.set_ser_data(SerializedData(params_str));
            dbg(APP_LOG, "Spawning iteration after frame %i\n", frame);
            SpawnComputeJob(ITERATION,
                    job_ids[0],
                    read, write,
                    before, after,
                    iter_params);
        }

        dbg(APP_LOG, "Completed executing loop job\n");
    }

} // namespace application

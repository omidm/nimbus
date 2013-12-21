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

#include "application/water_alternate_coarse/app_utils.h"
#include "application/water_alternate_coarse/job_loop.h"
#include "application/water_alternate_coarse/job_main.h"
#include "shared/nimbus.h"
#include <vector>

namespace application {

    JobMain::JobMain(Application *app) {
        set_application(app);
    };

    nimbus::Job* JobMain::Clone() {
        return new JobMain(application());
    }

    void JobMain::Execute(Parameter params, const DataArray& da) {
        dbg(APP_LOG, "Executing main job\n");

        int job_num = 1;
        std::vector<nimbus::job_id_t> job_ids;
        GetNewJobID(&job_ids, job_num);

        nimbus::IDSet<nimbus::logical_data_id_t> read, write;
        nimbus::IDSet<nimbus::job_id_t> before, after;

        nimbus::Parameter loop_params;
        std::stringstream out_frame_ss;
        int frame = 0;
        out_frame_ss << frame;
        loop_params.set_ser_data(SerializedData(out_frame_ss.str()));

        dbg(APP_LOG, "Spawning loop for frame %i in main\n", frame);
        SpawnComputeJob(LOOP,
                        job_ids[0],
                        read, write,
                        before, after,
                        loop_params);

        dbg(APP_LOG, "Completed executing main job\n");
    }

} // namespace application
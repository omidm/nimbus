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
 * This file contains a job corresponding to one iteration consisting of all
 * different simulation stages (advection, projection, extrapolation etc).
 *
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#include "application/water_alternate_coarse/app_utils.h"
#include "application/water_alternate_coarse/job_iteration.h"
#include "application/water_alternate_coarse/job_loop.h"
#include "application/water_alternate_coarse/water_driver.h"
#include "application/water_alternate_coarse/water_example.h"
#include "application/water_alternate_coarse/water_sources.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"
#include <sstream>
#include <string>

namespace application {

    JobIteration::JobIteration(Application *app) {
        set_application(app);
    };

    nimbus::Job* JobIteration::Clone() {
        return new JobIteration(application());
    }

    void JobIteration::Execute(Parameter params, const DataArray& da) {
        dbg(APP_LOG, "Executing iteration job\n");

        // get parameters
        int frame;
        std::stringstream in_frame_ss;
        std::string params_str(params.ser_data().data_ptr_raw(),
                               params.ser_data().size());
        in_frame_ss.str(params_str);
        in_frame_ss >> frame;

        // initialize configuration and state
        PhysBAM::WATER_EXAMPLE<TV> *example =
            new PhysBAM::WATER_EXAMPLE<TV>(PhysBAM::STREAM_TYPE((RW())));

        example->restart = frame;
        example->Initialize_Grid(TV_INT::All_Ones_Vector()*kScale,
                                 PhysBAM::RANGE<TV>(TV(),
                                                    TV::All_Ones_Vector())
                                 );
        PhysBAM::WaterSources::Add_Source(example);
        PhysBAM::WATER_DRIVER<TV> driver(*example);
        driver.Initialize(this, da);

        // simulate - advance a time step
        PhysBAM::LOG::SCOPE scope("FRAME", "Frame %d",
                                  driver.current_frame+1, 1);
        driver.Advance_To_Target_Time(example->Time_At_Frame(driver.current_frame+1));
        PhysBAM::LOG::Time("Reseed");
        example->particle_levelset_evolution.Reseed_Particles(driver.time);
        example->particle_levelset_evolution.Delete_Particles_Outside_Grid();
        driver.Write_Output_Files(++driver.output_number);
        frame++;

        // free resources
        delete example;

        int job_num = 1;
        std::vector<nimbus::job_id_t> job_ids;
        GetNewJobID(&job_ids, job_num);

        nimbus::IDSet<nimbus::logical_data_id_t> read, write;
        nimbus::IDSet<nimbus::job_id_t> before, after;

        nimbus::Parameter loop_params;
        std::stringstream out_frame_ss;
        out_frame_ss << frame;
        loop_params.set_ser_data(SerializedData(out_frame_ss.str()));

        dbg(APP_LOG, "Spawning loop after simulating frame %i\n", frame);
        SpawnComputeJob(LOOP,
                job_ids[0],
                read, write,
                before, after,
                loop_params);

        dbg(APP_LOG, "Completed executing iteration job\n");
}

} // namespace application

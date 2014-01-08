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

#include "application/water_alternate_fine/app_utils.h"
#include "application/water_alternate_fine/job_initialize.h"
#include "application/water_alternate_fine/job_loop.h"
#include "application/water_alternate_fine/water_driver.h"
#include "application/water_alternate_fine/water_example.h"
#include "application/water_alternate_fine/water_sources.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"

namespace application {

    JobInitialize::JobInitialize(nimbus::Application *app) {
        set_application(app);
    };

    nimbus::Job* JobInitialize::Clone() {
        return new JobInitialize(application());
    }

    void JobInitialize::Execute(nimbus::Parameter params, const nimbus::DataArray& da) {
        dbg(APP_LOG, "Executing initialize job\n");

        std::string par_str(params.ser_data().data_ptr_raw(),
                            params.ser_data().size());

        // initialize configuration and state
        PhysBAM::WATER_EXAMPLE<TV> *example =
            new PhysBAM::WATER_EXAMPLE<TV>(PhysBAM::STREAM_TYPE((RW())));

        example->Initialize_Grid(TV_INT::All_Ones_Vector()*kScale,
                                 PhysBAM::RANGE<TV>(TV(),
                                                    TV::All_Ones_Vector())
                                 );
        PhysBAM::WaterSources::Add_Source(example);
        PhysBAM::WATER_DRIVER<TV> driver(*example);
        driver.init_phase = true;
        driver.current_frame = 0;
        int last_unique_particle = driver.Initialize(this, da, 0);

        delete example;

        // next loop
        int job_num = 1;
        std::vector<nimbus::job_id_t> job_ids;
        GetNewJobID(&job_ids, job_num);

        nimbus::IDSet<nimbus::logical_data_id_t> read, write;
        nimbus::IDSet<nimbus::job_id_t> before, after;

        for (size_t i = 0; i < da.size(); ++i) {
            nimbus::Data *d = da[i];
            logical_data_id_t id = d->logical_id();
            if (!application::Contains(read, id))
                read.insert(id);
            if (!application::Contains(write, id))
                write.insert(id);
        }

        nimbus::Parameter loop_params;
        std::stringstream out_ss;
        int frame = 0;
        out_ss << frame;
        out_ss << "\n";
        assert(last_unique_particle >= 0);
        out_ss << last_unique_particle;
        loop_params.set_ser_data(SerializedData(out_ss.str()));

        dbg(APP_LOG, "Spawning loop for frame %i in main\n", frame);
        SpawnComputeJob(LOOP,
                        job_ids[0],
                        read, write,
                        before, after,
                        loop_params);

        dbg(APP_LOG, "Completed executing initialize job\n");
    }

} // namespace application

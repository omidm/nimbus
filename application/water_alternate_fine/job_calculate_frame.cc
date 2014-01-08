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
 * Author: Hang Qu <quhang@stanford.edu>
 */

#include <sstream>
#include <string>

#include "application/water_alternate_fine/app_utils.h"
#include "application/water_alternate_fine/job_loop.h"
#include "application/water_alternate_fine/water_driver.h"
#include "application/water_alternate_fine/water_example.h"
#include "application/water_alternate_fine/water_sources.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"

#include "application/water_alternate_fine/job_calculate_frame.h"

namespace application {

    JobCalculateFrame::JobCalculateFrame(nimbus::Application *app) {
        set_application(app);
    };

    nimbus::Job* JobCalculateFrame::Clone() {
        return new JobCalculateFrame(application());
    }

    void JobCalculateFrame::Execute(nimbus::Parameter params,
                                    const nimbus::DataArray& da) {
        dbg(APP_LOG, "Executing CALCULATE_FRAME job.\n");

        // TODO(quhang): Get last_uniqueid, time, dt and frame right.
        // get parameters
        T dt = 0;
        int frame, last_unique_particle;
        std::stringstream in_ss;
        std::string params_str(params.ser_data().data_ptr_raw(),
                               params.ser_data().size());
        in_ss.str(params_str);
        in_ss >> frame;
        in_ss >> last_unique_particle;
        dbg(APP_LOG, "Frame %i, last unique particle %i in iteration job\n",
                     frame, last_unique_particle);

        // Assume time, dt, frame is ready from here.
        T time = 0;
        dbg(APP_LOG,
            "Initialize WATER_DRIVER/WATER_EXAMPLE"
            "(Frame=%d, Time=%f, dt=%f).\n",
            frame, time, dt);
        PhysBAM::WATER_EXAMPLE<TV> *example =
            new PhysBAM::WATER_EXAMPLE<TV>(PhysBAM::STREAM_TYPE((RW())));

        example->Initialize_Grid(
            TV_INT::All_Ones_Vector()*kScale,
            PhysBAM::RANGE<TV>(TV(), TV::All_Ones_Vector()));
        PhysBAM::WaterSources::Add_Source(example);
        PhysBAM::WATER_DRIVER<TV> driver(*example);
        driver.init_phase = false;
        driver.current_frame = frame;
        driver.Initialize(this, da, last_unique_particle);
        driver.time = time;

        dbg(APP_LOG,
            "Simulation starts"
            "(Frame=%d, Time=%f, dt=%f).\n",
            frame, time, dt);
        // Move forward time "dt" without reseeding and writing frames.
        last_unique_particle = driver.CalculateFrameImpl(this, da, true, dt);

        // TODO(quhang/chinmayee): Fix the passing mechanism for
        // last_unique_particle.

        // Free resources.
        delete example;

        dbg(APP_LOG, "Completed executing CALCULATE_FRAME job\n");
}

} // namespace application

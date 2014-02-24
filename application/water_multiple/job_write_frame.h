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
 * This file contains job WRITE_FRAME that:
 *     reseeds particle, and writes simulation variables to file for rendering.
 * The parameters of WRITE_FRAME:
 *     frame number, simulation time, dt.
 * The read set(not sure) of SUPER_3:
 *     velocity, levelset, particle, removed particle, last_unique_particle_id.
 * The write set(not sure) of SUPER_3:
 *     particles.
 *
 * dt is not needed for job WRITE_FRAME now, and might be removed from the
 * parameter list in the future.
 * It is still unclear whether other simulation variables or states are also
 * needed.
 * For now, all the data is transmitted to guarantee correctness.
 *
 * Reseeding operation is included in this job, because this job includes all
 * the operations that are executed once in each frame, as is reseeding
 * operation. So the job name might be a little bit misleading. Reseeding
 * operation is expected to be moved to another job in the future.
 *
 * Author: Hang Qu <quhang@stanford.edu>
 */

#ifndef NIMBUS_APPLICATION_WATER_MULTIPLE_JOB_WRITE_FRAME_H_
#define NIMBUS_APPLICATION_WATER_MULTIPLE_JOB_WRITE_FRAME_H_

#include "shared/nimbus.h"

namespace application {

    class JobWriteFrame : public nimbus::Job {
        public:
            explicit JobWriteFrame(nimbus::Application *app);
            virtual void Execute(nimbus::Parameter params, const nimbus::DataArray& da);
            virtual nimbus::Job* Clone();
    };

} // namespace application

#endif  // NIMBUS_APPLICATION_WATER_MULTIPLE_JOB_WRITE_FRAME_H_

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

#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/job_synchronize_particles.h"
#include "application/water_multiple/job_names.h"
#include "application/water_multiple/data_particle_array.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"

namespace application {

JobSynchronizeParticles::JobSynchronizeParticles(nimbus::Application *app) {
    set_application(app);
};

nimbus::Job* JobSynchronizeParticles::Clone() {
    return new JobSynchronizeParticles(application());
}

void JobSynchronizeParticles::Execute(nimbus::Parameter params, const nimbus::DataArray& da) {
    dbg(APP_LOG, "Executing synchronize particles job\n");

    // need at least 2 elements - a data object to merge to, and a data
    // object to merge from
    size_t dnum = da.size();
    if (dnum < 2) {
        dbg(DBG_WARN, "Nothing to synchronize\n");
        return;
    }

    DataParticleArray *merge_to = dynamic_cast<DataParticleArray * >(da[dnum-1]);
    if (merge_to == NULL) {
        dbg(DBG_WARN, "Passed object is not a particle array\n");
    }
    GeometricRegion merge_region = merge_to->region();

    std::vector<DataParticleArray * > scratch;
    DataParticleArray *scratch_data;
    GeometricRegion scratch_region;
    for (size_t i = 0; i < dnum - 1; i++) {
        scratch_data = dynamic_cast<DataParticleArray * >(da[i]);
        if (scratch_data == NULL) {
            dbg(DBG_WARN, "Ignoring scratch object since it is not a particle array\n");
            continue;
        }
        scratch_region = scratch_data->region();
        if (scratch_region != merge_region) {
            dbg(DBG_WARN,
                "Ignoring scratch object since region does not match with merge region\n");
            continue;
        }
        scratch.push_back(scratch_data);
    }

    merge_to->MergeParticles(scratch);

    dbg(APP_LOG, "Completed executing synchronize particles job\n");
}

} // namespace application

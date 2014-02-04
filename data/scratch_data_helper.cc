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
 * Scratch data helper class contains a mapping between regions and
 * corresponding scratch region "names" as registered with Nimbus. This map
 * should be populated by the application writer, when defining the
 * corresponding scratch regions.
 *
 * The class provides functions that the application writer can use to get
 * scratch data for a job, and scratch data for a region - for synchronization
 * job.
 *
 * Supports only 3d.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#include <string>
#include <vector>

#include "data/scratch_data_helper.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"
#include "worker/job.h"

namespace nimbus {

ScratchDataHelper::ScratchDataHelper() {
    int temp[] = {7, 3, 3, 1, 1, 1};
    for (int i = 0; i < NUM_TYPES; i++) {
        scratch_type_names_[i] = "";
        num_scratch_[i] = temp[i];
    }
}

ScratchDataHelper::~ScratchDataHelper() {}

void ScratchDataHelper::set_domain(const GeometricRegion &d) {
    domain_ = d;
}

void ScratchDataHelper::set_share_boundary(bool flag) {
    share_boundary_ = flag;
}

void ScratchDataHelper::set_ghost_width(int gw) {
    ghost_width_ = gw;
}

void ScratchDataHelper::SetScratchType(const std::vector<ScratchTypes> &st_index,
                                       const std::vector<std::string> &st_names) {
    int ii = st_index.size();
    int ni = st_names.size();
    int mi = ii < ni? ii : ni;
    if (ii < ni)
        dbg(DBG_WARN, "WARNING: index and name aray sizes do not match in SetScratchType.\n");
    for (int i = 0; i < mi; i ++) {
        if (st_index[i] > NUM_TYPES) {
            dbg(DBG_ERROR, "ERROR: scratch region got an index > %i\n", NUM_TYPES);
            exit(1);
        }
        scratch_type_names_[st_index[i]] = st_names[i];
    }
}

void ScratchDataHelper::GetJobScratchData(const Job *j,
                                          const GeometricRegion &cr,
                                          lIDSet *ids) const {}

void ScratchDataHelper::GetAllScratchData(const Job *j,
                                          std::vector<GeometricRegion> *regions,
                                          std::vector<lIDSet> ids_list) const {}

void ScratchDataHelper::GetScratchRegions(const GeometricRegion &cr,
                                         std::vector<GeometricRegion> *regions) const {
}

}  // namespace nimbus

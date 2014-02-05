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

#ifndef NIMBUS_DATA_SCRATCH_DATA_HELPER_H_
#define NIMBUS_DATA_SCRATCH_DATA_HELPER_H_

#include <string>
#include <vector>

#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"
#include "worker/job.h"

namespace nimbus {

class ScratchDataHelper {
    public:
        enum ScratchTypes { VERTEX = 0,
                            EDGE1,
                            EDGE2,
                            FACE1,
                            FACE2,
                            FACE3,
                            NUM_TYPES
                          };

    private:
        GeometricRegion domain_;
        bool share_boundary_;
        int ghost_width_;
        std::string scratch_type_names_[NUM_TYPES];
        int num_scratch_[NUM_TYPES];

        typedef IDSet<logical_data_id_t> lIDSet;

        void GetScratchRegions(const GeometricRegion &cr,
                               std::vector<GeometricRegion> *regions) const;

    public:
        ScratchDataHelper();
        virtual ~ScratchDataHelper();

        void set_domain(const GeometricRegion &d);
        void set_share_boundary(bool flag);
        void set_ghost_width(int gw);

        void SetScratchType(const std::vector<ScratchTypes> &st_index,
                            const std::vector<std::string> &st_names);
        void GetJobScratchData(const Job *j,
                               const GeometricRegion &cr,
                               lIDSet *ids) const;
        void GetAllScratchData(const Job *j,
                               std::vector<GeometricRegion> *regions,
                               std::vector<lIDSet> ids_list) const;
};
}  // namespace nimbus

#endif  // NIMBUS_DATA_SCRATCH_DATA_HELPER_H_

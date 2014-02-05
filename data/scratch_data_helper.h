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
    enum ScratchType { VERTEX,
                       EDGE,
                       FACE,
                       ELEMS
                     };

    private:
        enum { DIMENSION = 3 };
        enum { XCOORD = 0,
               YCOORD,
               ZCOORD
             };
        enum { VERTEX_TYPES = 8 };
        enum { EDGE_TYPES = 12 };
        enum { FACE_TYPES = 6 };

        GeometricRegion domain_;
        int ghost_width_[DIMENSION];
        bool share_boundary_;

        std::string vertex_types_[VERTEX_TYPES];
        std::string edge_types_[EDGE_TYPES];
        std::string face_types_[FACE_TYPES];

        typedef IDSet<logical_data_id_t> lIDSet;

    public:
        ScratchDataHelper();
        virtual ~ScratchDataHelper();

        void set_domain(const GeometricRegion &d);
        void set_ghost_width(int gw[DIMENSION]);
        void set_share_boundary(bool flag);

        void SetScratchType(const std::vector<ScratchType> &st_index,
                            const std::vector<std::string> &st_names);
        void GetJobScratchData(Job *job,
                               const GeometricRegion &cr,
                               lIDSet *ids) const;
        void GetAllScratchData(Job *job,
                               std::vector<GeometricRegion> *regions,
                               std::vector<lIDSet> ids_list) const;
};
}  // namespace nimbus

#endif  // NIMBUS_DATA_SCRATCH_DATA_HELPER_H_

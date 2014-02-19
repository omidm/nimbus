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
 * The class provides functions that the application writer can use to get
 * scratch data for a job, scratch data from a region for synchronization
 * job and scratch name registration utilities.
 *
 * The user should provide base nimbus type name or a list of scratch names.
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
#include "worker/application.h"
#include "worker/data.h"
#include "worker/job.h"

namespace nimbus {

class ScratchDataHelper {
    public:
        /* scratch types in 3d - scratch regions at the corner, ones at the
         * edge, and onces on 2 sides of a face.
         */
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
        enum { EDGE_TYPES = 4 };
        enum { FACE_TYPES = 2 };

        /* ghost width ( = scratch region width)
         */
        int ghost_width_[DIMENSION];

        /* scratch type names.
         * example: there will be 8 vertex type names for the 8 different
         * corners, since 1 corner region is shared by 8 different application
         * partitions.
         */
        std::string vertex_types_[VERTEX_TYPES];
        std::string edge_types_[EDGE_TYPES];
        std::string face_types_[FACE_TYPES];

        typedef IDSet<logical_data_id_t> lIDSet;

    public:
        ScratchDataHelper();
        explicit ScratchDataHelper(const int gw[DIMENSION]);
        ScratchDataHelper(const int gw[DIMENSION], const std::string b_name);
        virtual ~ScratchDataHelper();

        void set_ghost_width(const int gw[DIMENSION]);

        /* pass a scratch type and corresponding list of names.
         */
        void SetScratchNames(const ScratchType st,
                             const std::vector<std::string> &st_names);
        /* use default scratch type names.
         */
        void SetScratchBaseName(const std::string b_name);

        /* resgister scratch types.
         */
        void RegisterScratchNames(Application *app,
                                  Data *data);

        /* given the inner bounding box for an application partition (that
         * includes the inner scratch region, this gives all the scratch
         * regions associated with the job region.
         * inner region: for a 1d domain 1-30, partitions 2, ghost width 3,
         * inner regions are 1-15 and 16-30.
         */
        void GetJobScratchData(Job *job,
                               const GeometricRegion &cr,
                               lIDSet *ids) const;
        /* given a region, this function obtains all the scratch data ids in
         * that region.
         */
        void GetAllScratchData(Job *job,
                               const GeometricRegion &region,
                               ScratchType st,
                               lIDSet *ids) const;
};
}  // namespace nimbus

#endif  // NIMBUS_DATA_SCRATCH_DATA_HELPER_H_

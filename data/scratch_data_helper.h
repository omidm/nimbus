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
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#ifndef NIMBUS_DATA_SCRATCH_DATA_HELPER_H_
#define NIMBUS_DATA_SCRATCH_DATA_HELPER_H_

#include <map>
#include <string>
#include <vector>

#include "shared/geometric_region.h"

namespace nimbus {

class ScratchDataHelper {
    private:
        typedef std::vector<std::string> NameList;
        struct CompareRegions {
            bool operator() (const GeometricRegion &r1,
                             const GeometricRegion &r2) {
                if (r1.x() != r2.x())
                    return(r1.x() < r2.x());
                else if (r1.y() != r2.y())
                    return(r1.y() < r2.y());
                else if (r1.z() != r2.z())
                    return(r1.z() < r2.z());
                else if (r1.dx() != r2.dx())
                    return(r1.dx() < r2.dx());
                else if (r1.dy() != r2.dy())
                    return(r1.dy() < r2.dy());
                else
                    return(r1.dz() < r2.dz());
            }
        };
        typedef std::map<nimbus::GeometricRegion, NameList, CompareRegions>
            RegionNameListMap;
        RegionNameListMap scratch_map_;

    public:
        ScratchDataHelper();
        virtual ~ScratchDataHelper();
};
}  // namespace nimbus

#endif  // NIMBUS_DATA_SCRATCH_DATA_HELPER_H_

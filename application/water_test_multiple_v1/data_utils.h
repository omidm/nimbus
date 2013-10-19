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
 * This file enumerates the type of data, based on location --
 * (corner/ side/ interior) and utility functions.
 * The code works for only 2d regions.
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#ifndef NIMBUS_APPLICATION_WATER_TEST_MULTIPLE_V1_DATA_UTILS_H_
#define NIMBUS_APPLICATION_WATER_TEST_MULTIPLE_V1_DATA_UTILS_H_

#include "PhysBAM_Tools/Vectors/VECTOR.h"
#include "shared/nimbus.h"

namespace water_app_data {

    namespace {
        const int dimension = 2;
        typedef ::PhysBAM::VECTOR<int, dimension> TV_INT;
    }

    // data locations for a partition
    // update GetDataRegionNames and GetDataRegionSizes
    // when you update DataRegion
    enum DataRegion {
        kDataInterior,
        kDataInUpperLeft,
        kDataInUpper,
        kDataInUpperRight,
        kDataInRight,
        kDataInBottomRight,
        kDataInBottom,
        kDataInBottomLeft,
        kDataInLeft,
        kDataNum
    };
    // this function depends on enum DataRegion, and should be changed if
    // DataRegion is updated
    void GetDataRegionNames(std::string names[]);
    // this function depends on enum DataRegion, and should be changed if
    // DataRegion is updated
    void GetDataRegionSizes(
            TV_INT sizes[],
            TV_INT part_center_size,
            TV_INT sim_center_size,
            int ghost_band);

    class SimData : public ::nimbus::Data {
        private:
            DataRegion region_;
        public:
            SimData() : region_(kDataInterior) {}
            DataRegion region() {
                return region_;
            }
            void set_region(DataRegion region) {
                region_ = region;
            }
    };

} // namespace water_app_data

#endif // NIMBUS_APPLICATION_WATER_TEST_MULTIPLE_V1_DATA_UTILS_H_

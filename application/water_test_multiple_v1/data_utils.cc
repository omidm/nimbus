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

#include "data_utils.h"

namespace water_app_data {

    void GetDataRegionNames(std::string names[]) {
        names[kDataInterior] = "DataInterior";
        names[kDataInUpperLeft] = "DataInUpperLeft";
        names[kDataInUpper] = "DataInUpper";
        names[kDataInUpperRight] = "DataInUpperRight";
        names[kDataInRight] = "DataInRight";
        names[kDataInBottomRight] = "DataInBottomRight";
        names[kDataInBottom] = "DataInBottom";
        names[kDataInBottomLeft] = "DataInBottomLeft";
        names[kDataInLeft] = "DataInLeft";
    }

    void GetDataRegionSizes(
            TV_INT sizes[],
            TV_INT part_center_size,
            TV_INT sim_center_size,
            int ghost_band) {
        TV_INT part_ghost_corner_size(ghost_band, ghost_band);
        TV_INT part_ghost_vert_size = part_center_size;
        TV_INT part_ghost_horiz_size = part_center_size;
        part_ghost_vert_size(0) = ghost_band;
        part_ghost_horiz_size(1) = ghost_band;
        sizes[kDataInterior] = sim_center_size;
        sizes[kDataInUpperLeft] = part_ghost_corner_size;
        sizes[kDataInUpper] = part_ghost_horiz_size;
        sizes[kDataInUpperRight] = part_ghost_corner_size;
        sizes[kDataInRight] = part_ghost_vert_size;
        sizes[kDataInBottomRight] = part_ghost_corner_size;
        sizes[kDataInBottom] = part_ghost_horiz_size;
        sizes[kDataInBottomLeft] = part_ghost_corner_size;
        sizes[kDataInLeft] = part_ghost_vert_size;
    }

} // namespace water_app_data

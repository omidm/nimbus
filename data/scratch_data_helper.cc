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

ScratchDataHelper::ScratchDataHelper() {}

ScratchDataHelper::~ScratchDataHelper() {}

void ScratchDataHelper::set_domain(const GeometricRegion &d) {
    domain_ = d;
}

void ScratchDataHelper::set_ghost_width(int gw[DIMENSION]) {
    for (size_t i = 0; i < DIMENSION; i++)
        ghost_width_[i] = gw[i];
}

void ScratchDataHelper::set_share_boundary(bool flag) {
    share_boundary_ = flag;
}

void ScratchDataHelper::SetScratchType(const std::vector<ScratchType> &st_index,
                                       const std::vector<std::string> &st_names) {
}

void ScratchDataHelper::GetJobScratchData(Job *job,
                                          const GeometricRegion &cr,
                                          lIDSet *ids) const {
    ids->clear();

    const int cl[DIMENSION]  = {cr.x(),  cr.y(),  cr.z()};
    const int cld[DIMENSION] = {cr.dx(), cr.dy(), cr.dx()};
    int l[DIMENSION]  = {0, 0, 0};
    int ld[DIMENSION] = {0, 0, 0};
    size_t n;

    // vertex scratch regions
    for (int d = 0; d < DIMENSION; d++) {
        l[d]  = cl[d] - ghost_width_[d];
        ld[d] = 2*ghost_width_[d];
    }
    n  = 0;
    for (size_t i = 0; i < 2; i++) {
        l[XCOORD] = cl[XCOORD] - ghost_width_[XCOORD] +
                    i * cld[XCOORD];
        for (size_t j = 0; j < 2; j++) {
            l[YCOORD] = cl[YCOORD] - ghost_width_[YCOORD] +
                        j * cld[YCOORD];
            for (size_t k = 0; k < 2; k++) {
                l[ZCOORD] = cl[ZCOORD] - ghost_width_[ZCOORD] +
                            k * cld[ZCOORD];
                CLdoVector ldos;
                GeometricRegion region(l[XCOORD], l[YCOORD], l[ZCOORD],
                                       ld[XCOORD], ld[YCOORD], ld[ZCOORD]);
                job->GetCoveredLogicalObjects(&ldos, vertex_types_[n], &region);
                for (size_t s = 0; s < ldos.size(); s++)
                    ids->insert(ldos[s]->id());
                n++;
            }
        }
    }

    // edge scratch regions
    n = 0;
    for (size_t d = 0; d < DIMENSION; d++) {
        l[d]  = cl[d] + ghost_width_[d];
        ld[d] = cld[d] - 2*ghost_width_[d];
        size_t d1 = (d+1)%DIMENSION;
        size_t d2 = (d+2)%DIMENSION;
        ld[d1] = 2*ghost_width_[d1];
        ld[d2] = 2*ghost_width_[d2];
        for (int i = 0; i < 2; i++) {
            l[d1] = cl[d1] - ghost_width_[d1] +
                    i * cld[d1];
            for (int j = 0; j < 2; j++) {
                l[d2] = cl[d2] - ghost_width_[d2] +
                        j * cld[d2];
                CLdoVector ldos;
                GeometricRegion region(l[XCOORD], l[YCOORD], l[ZCOORD],
                                       ld[XCOORD], ld[YCOORD], ld[ZCOORD]);
                job->GetCoveredLogicalObjects(&ldos, edge_types_[n], &region);
                for (size_t s = 0; s < ldos.size(); s++)
                    ids->insert(ldos[s]->id());
                n++;
            }
        }
    }

    // face scratch regions
    n = 0;
    for (size_t d = 0; d < DIMENSION; d++) {
        size_t d1 = (d+1)%DIMENSION;
        size_t d2 = (d+2)%DIMENSION;
        ld[d] = 2*ghost_width_[d];
        l[d1]  = cl[d1] + ghost_width_[d1];
        l[d2]  = cl[d2] + ghost_width_[d2];
        ld[d1] = cld[d1] - 2*ghost_width_[d1];
        ld[d2] = cld[d2] - 2*ghost_width_[d2];
        for (size_t i = 0; i < 2; i++) {
            l[d]  = cl[d] - ghost_width_[d] + i * cld[d];
            CLdoVector ldos;
            GeometricRegion region(l[XCOORD], l[YCOORD], l[ZCOORD],
                                   ld[XCOORD], ld[YCOORD], ld[ZCOORD]);
            job->GetCoveredLogicalObjects(&ldos, face_types_[n], &region);
            for (size_t s = 0; s < ldos.size(); s++)
                ids->insert(ldos[s]->id());
            n++;
        }
    }
}

void ScratchDataHelper::GetJobScratchData(Job *job,
                                          const GeometricRegion &cr,
                                          const GeometricRegion &region,
                                          lIDSet *ids) const {}

void ScratchDataHelper::GetAllScratchData(Job *j,
                                          std::vector<GeometricRegion> *regions,
                                          std::vector<lIDSet> ids_list) const {}
}  // namespace nimbus

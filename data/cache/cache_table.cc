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
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#include <vector>

#include "data/cache/cache_object.h"
#include "data/cache/cache_table.h"
#include "data/cache/utils.h"
#include "shared/geometric_region.h"
#include "worker/data.h"
#include "worker/job.h"

namespace nimbus {

CacheTable::CacheTable() : table_(Table(GeometricRegionLess)) {}

void CacheTable::AddEntry(const GeometricRegion &region,
                          CacheObject *co) {
    CacheObjects *cos = NULL;
    if (table_.find(region) == table_.end()) {
        cos = new CacheObjects();
        table_[region] = cos;
    } else {
        cos = table_[region];
    }
    cos->push_back(co);
}

CacheObject *CacheTable::GetClosestAvailable(const GeometricRegion &region,
                                             const DataSet &read,
                                             CacheAccess access) {
    if (table_.find(region) == table_.end())
        return NULL;
    int min = GetMinDistanceIndex(table_[region], read, access);
    if (min == -1)
        return NULL;
    else
        return table_[region]->at(min);
}

int CacheTable::GetMinDistanceIndex(const CacheObjects *objects,
                                    const DataSet &read,
                                    CacheAccess access) const {
    size_t num_objects = objects->size();
    std::vector<distance_t> distance_vector(num_objects);
    distance_t min_distance = 2*read.size();
    int min_index = -1;
    for (size_t i = 0; i < objects->size(); ++i) {
        distance_vector[i] = objects->at(i)->GetDistance(read, access);
        if (distance_vector[i] == 0) {
            min_distance = 0;
            min_index = i;
            break;
        } else if (distance_vector[i] < min_distance) {
            min_distance = distance_vector[i];
            min_index = i;
        }
    }
    return min_index;
}

}  // namespace nimbus

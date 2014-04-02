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
#include "worker/data.h"
#include "worker/job.h"

namespace nimbus {

CacheTable::CacheTable() : table_(Table(GeometricRegionLess)) {}

void* CacheTable::GetCachedObject(const GeometricRegion &region,
                                  const Job &job,
                                  const DataArray &da) {
    CacheObjects *objects;
    if (table_.find(region) == table_.end()) {
        objects = new CacheObjects();
        table_[region] = objects;
    } else {
        objects = table_[region];
    }

    int num_objects = objects->size();
    DataSet read, write;
    StringSet read_var, write_var;
    GetReadWrite(job, da, &read, &write, &read_var, &write_var);

    if (num_objects == 0) {
        // faster cache object allocation for common cases
        // TODO(chinmayee): create new object
        return NULL;
    } else if (num_objects == 1 && objects->at(0)->IsAvailable()) {
        // faster cache object allocation for common cases
        objects->at(0)->LockData(read, write, read_var, write_var);
        return objects->at(0);
    } else {
        // case with multiple objects in cache
        int min_index =
            GetMinDistanceIndex(objects, read, write, read_var, write_var);
        if (min_index < 0) {
            // no cached object available for reuse
            // TODO(chinmayee): create new object
            return NULL;
        } else {
            objects->at(min_index)->
                LockData(read, write, read_var, write_var);
            return objects->at(min_index);
        }
    }
}

void CacheTable::GetReadWrite(const Job &job,
                              const DataArray &da,
                              DataSet *read,
                              DataSet *write,
                              StringSet *read_var,
                              StringSet *write_var) {
    size_t num_data = da.size();
    PIDSet read_ids = job.read_set();
    PIDSet write_ids = job.write_set();
    for (size_t i = 0; i < num_data; i++) {
        if (read_ids.contains(da[i]->physical_id())) {
            read->insert(da[i]);
            read_var->insert(da[i]->name());
        }
        if (write_ids.contains(da[i]->physical_id())) {
            write->insert(da[i]);
            write_var->insert(da[i]->name());
        }
    }
}

int CacheTable::GetMinDistanceIndex(const CacheObjects *objects,
                                    const DataSet &read,
                                    const DataSet &write,
                                    const StringSet &read_var,
                                    const StringSet &write_var) {
    size_t num_objects = objects->size();
    std::vector<distance_t> distance_vector(num_objects);
    distance_t min_distance = 2*read.size();
    int min_index = -1;

    for (size_t i = 0; i < objects->size(); i++) {
        distance_vector[i] = objects->at(i)->
            GetDistance(read, write, read_var, write_var);
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

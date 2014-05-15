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
 * A CacheStruct is an application object corresponding to a set of nimbus
 * variables cached by the cache manager.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#include <map>
#include <set>
#include <string>
#include <vector>

#include "data/cache/cache_struct.h"
#include "data/cache/utils.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"
#include "worker/data.h"

namespace nimbus {

/**
 * \details
 */
CacheStruct::CacheStruct(size_t num_variables) : num_variables_(num_variables),
                                                 data_maps_(num_variables),
                                                 write_backs_(num_variables) {
}

/**
 * \details UpdateCache(...) finds out what data is in read_sets and not in
 * CacheStruct instance, and calls a ReadToCache(...) on the diff. Before
 * reading the diff, UpdateCache(...) flushes out data that will be replaced
 * (which is found based on geometric region of the data in cache and in
 * read_sets.
 */
void CacheStruct::UpdateCache(const std::vector<type_id_t> &var_type,
                              const std::vector<DataArray> &read_sets,
                              const GeometricRegion &read_region) {
    size_t num_vars = var_type.size();
    if (read_sets.size() != num_vars) {
        dbg(DBG_ERROR, "Mismatch in number of variable types passed to UpdateCache\n");
        exit(-1);
    }
    std::vector<DataArray> diff(num_vars), flush(num_vars);
    for (size_t t = 0; t < num_vars; ++t) {
        type_id_t type = var_type[t];
        if (type > num_variables_) {
            dbg(DBG_WARN, "Invalid type %u passed to UpdateCache, ignoring it\n", type);
            continue;
        }
        DataArray *diff_t = &diff[t];
        DataArray *flush_t = &flush[t];
        DMap *data_map_t = &data_maps_[type];
        DataSet *write_back_t = &write_backs_[type];
        for (size_t i = 0; i < read_sets[t].size(); ++i) {
            Data *d = read_sets[t].at(i);
            GeometricRegion dreg = d->region();
            DMap::iterator it = data_map_t->find(dreg);
            if (it == data_map_t->end()) {
                diff_t->push_back(d);
            } else {
                Data *d_old = it->second;
                if (d_old != d) {
                    diff_t->push_back(d);
                    // d_old->UnsetCacheObjectMapping(this);
                    if (write_back_t->find(d) != write_back_t->end()) {
                        flush_t->push_back(d);
                    }
                }
            }
        }
    }
    FlushCache(var_type, flush);
    ReadToCache(var_type, diff, read_region);
    for (size_t t = 0; t < num_vars; ++t) {
        type_id_t type = var_type[t];
        DataArray *diff_t = &diff[t];
        DMap *data_map_t = &data_maps_[type];
        for (size_t i = 0; i < diff_t->size(); ++i) {
            Data *d = diff_t->at(i);
            GeometricRegion dreg = d->region();
            (*data_map_t)[dreg] = d;
            // d->SetUpCacheObjectMapping(this);
        }
    }
}

/**
 * \details SetUpWrite(...) sets up mapping for write data and the write
 * region for a struct instance. It also sets up the write back data for the
 * instance. Later, when the instance is released, these mappings and write
 * region are used to pull data into nimbus data instances when needed.
 */
void CacheStruct::SetUpWrite(const std::vector<type_id_t> &var_type,
                             const std::vector<DataArray> &write_sets,
                             const GeometricRegion &write_region) {
    size_t num_vars = var_type.size();
    if (write_sets.size() != num_vars) {
        dbg(DBG_ERROR, "Mismatch in number of variable types passed to SetUpWrite\n");
        exit(-1);
    }
    // TODO(Chinmayee): Twin data - imeplementation incomplete
    std::vector<DataArray> diff(num_vars), flush(num_vars);
    for (size_t t = 0; t < num_vars; ++t) {
        type_id_t type = var_type[t];
        if (type > num_variables_) {
            dbg(DBG_WARN, "Invalid type %u passed to SetUpWrite, ignoring it\n", type);
            continue;
        }
        DataArray *diff_t = &diff[t];
        DataArray *flush_t = &flush[t];
        DMap *data_map_t = &data_maps_[type];
        DataSet *write_back_t = &write_backs_[type];
        for (size_t i = 0; i < write_sets[t].size(); ++i) {
            Data *d = write_sets[t].at(i);
            GeometricRegion dreg = d->region();
            DMap::iterator it = data_map_t->find(dreg);
            if (it == data_map_t->end()) {
                diff_t->push_back(d);
            } else {
                Data *d_old = it->second;
                if (d_old != d) {
                    diff_t->push_back(d);
                    // d_old->UnsetCacheObjectMapping(this);
                    if (write_back_t->find(d) != write_back_t->end()) {
                        flush_t->push_back(d);
                    }
                }
            }
        }
    }
    FlushCache(var_type, flush);
    for (size_t t = 0; t < num_vars; ++t) {
        type_id_t type = var_type[t];
        DataArray *diff_t = &diff[t];
        DMap *data_map_t = &data_maps_[type];
        DataSet *write_back_t = &write_backs_[type];
        for (size_t i = 0; i < diff_t->size(); ++i) {
            Data *d = diff_t->at(i);
            GeometricRegion dreg = d->region();
            (*data_map_t)[dreg] = d;
            write_back_t->insert(d);
            // d->SetUpCacheObjectMapping(this);
            // d->SetUpDirtyCacheObjectMapping(this);
        }
    }
}

}  // namespace nimbus

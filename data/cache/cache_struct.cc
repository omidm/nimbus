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

#include "data/cache/cache_defs.h"
#include "data/cache/cache_struct.h"
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
 * \details
 */
CacheStruct::CacheStruct(const GeometricRegion &ob_reg) : CacheObject(ob_reg) {
}

/**
 * \details UpdateCache(...) finds out what data is in read_sets and not in
 * CacheStruct instance, and calls a ReadToCache(...) on the diff. Before
 * reading the diff, UpdateCache(...) flushes out data that will be replaced
 * (which is found based on geometric region of the data in cache and in
 * read_sets.
 */
void CacheStruct::UpdateCache(const std::vector<cache::type_id_t> &var_type,
                              const std::vector<DataArray> &read_sets,
                              const GeometricRegion &read_region,
                              const GeometricRegion &write_region,
                              bool invalidate_read_minus_write) {
    size_t num_vars = var_type.size();
    if (read_sets.size() != num_vars) {
        dbg(DBG_ERROR, "Mismatch in number of variable types passed to UpdateCache\n");
        exit(-1);
    }
    std::vector<DataArray> diff(num_vars), flush(num_vars), to_map(num_vars);
    for (size_t t = 0; t < num_vars; ++t) {
        cache::type_id_t type = var_type[t];
        if (type > num_variables_) {
            dbg(DBG_WARN, "Invalid type %u passed to UpdateCache, ignoring it\n", type);
            continue;
        }
        const DataArray &read_set_t = read_sets[t];
        DataArray &diff_t = diff[t];
        DataArray &flush_t = flush[t];
        DataArray &to_map_t = to_map[t];
        DMap &data_map_t = data_maps_[type];
        const DataSet &write_back_t = write_backs_[type];
        for (size_t i = 0; i < read_set_t.size(); ++i) {
            Data *d = read_set_t.at(i);
            GeometricRegion dreg = d->region();
            DMap::iterator it = data_map_t.find(dreg);
            if (it == data_map_t.end()) {
                d->SyncData();
                diff_t.push_back(d);
                if (!invalidate_read_minus_write ||
                    write_region.Covers(&dreg)) {
                    to_map_t.push_back(d);
                }
            } else {
                Data *d_old = it->second;
                if (d_old != d) {
                    d->SyncData();
                    diff_t.push_back(d);
                    if (write_back_t.find(d_old) != write_back_t.end()) {
                        flush_t.push_back(d_old);
                    }
                    if (!invalidate_read_minus_write ||
                        write_region.Covers(&dreg)) {
                        to_map_t.push_back(d);
                    }
                    data_map_t.erase(it);
                    d_old->UnsetCacheObject(this);
                }
            }
        }
    }
    FlushCache(var_type, flush);
    ReadToCache(var_type, diff, read_region);
    for (size_t t = 0; t < num_vars; ++t) {
        cache::type_id_t type = var_type[t];
        const DataArray &diff_t = diff[t];
        DMap &data_map_t = data_maps_[type];
        for (size_t i = 0; i < diff_t.size(); ++i) {
            Data *d = diff_t.at(i);
            GeometricRegion dreg = d->region();
            if (!invalidate_read_minus_write || write_region.Covers(&dreg)) {
                data_map_t[dreg] = d;
                d->SetUpCacheObject(this);
            }
        }
    }
}

/**
 * \details SetUpWrite(...) sets up mapping for write data and the write
 * region for a struct instance. It also sets up the write back data for the
 * instance. Later, when the instance is released, these mappings and write
 * region are used to pull data into nimbus data instances when needed.
 */
void CacheStruct::SetUpWrite(const std::vector<cache::type_id_t> &var_type,
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
        cache::type_id_t type = var_type[t];
        if (type > num_variables_) {
            dbg(DBG_WARN, "Invalid type %u passed to SetUpWrite, ignoring it\n", type);
            continue;
        }
        const DataArray &write_set_t = write_sets[t];
        DataArray &diff_t = diff[t];
        DataArray &flush_t = flush[t];
        DMap &data_map_t = data_maps_[type];
        DataSet &write_back_t = write_backs_[type];
        for (size_t i = 0; i < write_set_t.size(); ++i) {
            Data *d = write_set_t.at(i);
            GeometricRegion dreg = d->region();
            DMap::iterator it = data_map_t.find(dreg);
            if (it == data_map_t.end()) {
                diff_t.push_back(d);
            } else {
                Data *d_old = it->second;
                if (d_old != d) {
                    diff_t.push_back(d);
                    if (write_back_t.find(d_old) != write_back_t.end()) {
                        flush_t.push_back(d_old);
                    }
                    data_map_t.erase(it);
                    d_old->UnsetCacheObject(this);
                }
            }
        }
    }
    FlushCache(var_type, flush);
    for (size_t t = 0; t < num_vars; ++t) {
        cache::type_id_t type = var_type[t];
        const DataArray &diff_t = diff[t];
        DMap &data_map_t = data_maps_[type];
        DataSet &write_back_t = write_backs_[type];
        for (size_t i = 0; i < diff_t.size(); ++i) {
            Data *d = diff_t.at(i);
            GeometricRegion dreg = d->region();
            data_map_t[dreg] = d;
            write_back_t.insert(d);
            d->SetUpCacheObject(this);
            d->SetUpDirtyCacheObject(this);
        }
    }
}

/**
 * \details UnsetData(...) removes data d from data_maps_. It checks every data
 * map in data_maps_ for the data d, since no prior type information is
 * available (just like PullData).
 */
void CacheStruct::UnsetData(Data *d) {
    GeometricRegion dreg = d->region();
    for (size_t t = 0; t < num_variables_; ++t) {
        DMap &data_map_t = data_maps_[t];
        if (data_map_t.find(dreg) != data_map_t.end()) {
            data_map_t.erase(dreg);
            break;
        }
    }
}

/**
 * \details UnsetDirtyData(...) removes d from write_backs_.
 */
void CacheStruct::UnsetDirtyData(Data *d) {
    for (size_t t = 0; t < num_variables_; ++t) {
        DataSet &write_back_t = write_backs_[t];
        DataSet::iterator it = write_back_t.find(d);
        if (it != write_back_t.end()) {
            write_back_t.erase(it);
            break;
        }
    }
}

/**
 * \details PullData(...) pulls data from cache, after locking the struct.
 * When data needs to be updated from outside CacheStruct, use PullData.
 * PullData(...)
 * checks each write_back set in the list write_backs_, and finds out which
 * type this data corresponds to. This may be a bad idea, but with twin copy
 * implementation, I expect the overhead to be small. -- Chinmayee
 */
void CacheStruct::PullData(Data *d) {
    AcquireAccess(cache::EXCLUSIVE);
    for (size_t t = 0; t < num_variables_; ++t) {
        DataSet &write_back_t = write_backs_[t];
        if (write_back_t.find(d) == write_back_t.end())
            continue;
        std::vector<cache::type_id_t> var_type(1, t);
        std::vector<DataArray> write_sets(1, DataArray(1, d));
        GeometricRegion dreg = d->region();
        GeometricRegion wreg = GeometricRegion::
            GetIntersection(write_region_, dreg);
        WriteFromCache(var_type, write_sets, wreg);
        d->UnsetDirtyCacheObject(this);
        write_back_t.erase(d);
        break;
    }
    ReleaseAccess();
}

/**
 * \details GetDistance(...) gives the cost of using the CacheStruct instance,
 * given the read set. Current cost function is the sum of geometric sizes of
 * data in the read set.
 * The function does not have a separate argument for write set
 * - if you want write set included in the cost function, either add  another
 * argument or append it to read set that is passed to GetDistance.
 */
cache::distance_t CacheStruct::GetDistance(const std::vector<cache::type_id_t> &var_type,
                                           const std::vector<DataArray> &read_sets) const {
    cache::distance_t cur_distance = 0;
    for (size_t t = 0; t < num_variables_; ++t) {
        cache::type_id_t type = var_type[t];
        if (type > num_variables_) {
            dbg(DBG_WARN, "Invalid type %u passed to SetUpWrite, ignoring it\n", type);
            continue;
        }
        const DMap &data_map_t = data_maps_[type];
        for (size_t i = 0; i < read_sets[t].size(); ++i) {
            Data *d = read_sets[t].at(i);
            GeometricRegion dreg = d->region();
            DMap::const_iterator it = data_map_t.find(dreg);
            if (it->second == d)
                continue;
            cur_distance += dreg.dx() * dreg.dy() * dreg.dz();
        }
    }
    return cur_distance;
}

/**
 * \details FlushCache(...) flushes all data passed to it, and unsets
 * corresponding dirty object mappings. This function should be used by
 * methods of CacheVar only. It does not check if data in flush_sets is in
 * write_backs_, making it unsafe. It also provides no locking.
 */
void CacheStruct::FlushCache(const std::vector<cache::type_id_t> &var_type,
                             const std::vector<DataArray> &flush_sets) {
    size_t num_vars = var_type.size();
    if (flush_sets.size() != num_vars) {
        dbg(DBG_ERROR, "Mismatch in number of variable types passed to FlushCache\n");
        exit(-1);
    }
    WriteFromCache(var_type, flush_sets, write_region_);
    for (size_t t = 0; t < num_vars; ++t) {
        const DataArray &flush_t = flush_sets[t];
        cache::type_id_t type = var_type[t];
        DataSet &write_back_t = write_backs_[type];
        for (size_t i = 0; i < flush_t.size(); ++i) {
            Data *d = flush_t[i];
            d->UnsetDirtyCacheObject(this);
            write_back_t.erase(d);
        }
    }
}

}  // namespace nimbus

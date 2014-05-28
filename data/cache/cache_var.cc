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
 * A CacheVar is an application object corresponding to one nimbus
 * variable, cached by the cache manager.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#include <map>
#include <set>
#include <string>
#include <vector>

#include "data/cache/cache_defs.h"
#include "data/cache/cache_var.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"
#include "worker/data.h"

namespace nimbus {

/**
 * \details
 */
CacheVar::CacheVar() {}

/**
 * \details
 */
CacheVar::CacheVar(const GeometricRegion &ob_reg) : CacheObject(ob_reg) {
}

/**
 * \details UpdateCache(...) finds out what data is in read_sets and not in
 * CacheVar instance, and calls a ReadToCache(...) on the diff. Before
 * reading the diff, UpdateCache(...) flushes out data that will be replaced
 * (which is found based on geometric region of the data in cache and in
 * read_set.
 */
void CacheVar::UpdateCache(const DataArray &read_set,
                           const GeometricRegion &read_region,
                           const GeometricRegion &valid_region,
                           bool invalidate_read_minus_valid) {
    DataArray diff, flush, to_map;
    for (size_t i = 0; i < read_set.size(); ++i) {
        Data *d = read_set.at(i);
        GeometricRegion dreg = d->region();
        DMap::iterator it = data_map_.find(dreg);
        if (it == data_map_.end()) {
            d->SyncData();
            diff.push_back(d);
            if (!invalidate_read_minus_valid ||
                valid_region.Covers(&dreg)) {
                to_map.push_back(d);
            }
        } else {
            Data *d_old = it->second;
            if (d_old != d) {
                d->SyncData();
                diff.push_back(d);
                if (write_back_.find(d_old) != write_back_.end()) {
                    flush.push_back(d_old);
                }
                if (!invalidate_read_minus_valid ||
                    valid_region.Covers(&dreg)) {
                    to_map.push_back(d);
                }
                data_map_.erase(it);
                d_old->UnsetCacheObject(this);
            }
        }
    }
    if (!flush.empty())
        FlushCache(flush);
    ReadToCache(diff, read_region);
    for (size_t i = 0; i < to_map.size(); ++i) {
        Data *d = to_map.at(i);
        GeometricRegion dreg = d->region();
        if (!invalidate_read_minus_valid || valid_region.Covers(&dreg)) {
            data_map_[dreg] = d;
            d->SetUpCacheObject(this);
        }
    }
}

/**
 * \details SetUpWrite(...) sets up mapping for write data and the write
 * region for a struct instance. It also sets up the write back data for the
 * instance. Later, when the instance is released, these mappings and write
 * region are used to pull data into nimbus data instances when needed.
 */
void CacheVar::SetUpWrite(const DataArray &write_set,
                          const GeometricRegion &write_region) {
    // TODO(Chinmayee): Twin data - imeplementation incomplete
    DataArray diff, flush;
    for (size_t i = 0; i < write_set.size(); ++i) {
        Data *d = write_set.at(i);
        GeometricRegion dreg = d->region();
        DMap::iterator it = data_map_.find(dreg);
        if (it == data_map_.end()) {
            diff.push_back(d);
        } else {
            Data *d_old = it->second;
            if (d_old != d) {
                diff.push_back(d);
                if (write_back_.find(d_old) != write_back_.end()) {
                    flush.push_back(d_old);
                }
                data_map_.erase(it);
                d_old->UnsetCacheObject(this);
            }
        }
    }
    if (!flush.empty())
        FlushCache(flush);
    for (size_t i = 0; i < write_set.size(); ++i) {
        Data *d = write_set.at(i);
        d->InvalidateCacheData();
        GeometricRegion dreg = d->region();
        data_map_[dreg] = d;
        write_back_.insert(d);
        d->SetUpCacheObject(this);
        d->SetUpDirtyCacheObject(this);
    }
    write_region_ = write_region;
}

/**
 * \details UnsetData(...) removes data d from data_map_.
 */
void CacheVar::UnsetData(Data *d) {
    GeometricRegion dreg = d->region();
    if (data_map_.find(dreg) != data_map_.end()) {
        assert(data_map_[dreg] == d);
        data_map_.erase(dreg);
    }
}

/**
 * \details UnsetDirtyData(...) removes d from write_back_.
 */
void CacheVar::UnsetDirtyData(Data *d) {
    write_back_.erase(d);
}

/**
 * \details PullData(...) pulls data from cache, after locking the struct.
 * When data needs to be updated from outside CacheVar, use PullData.
 */
void CacheVar::PullData(Data *d) {
    AcquireAccess(cache::EXCLUSIVE);
    if (write_back_.find(d) != write_back_.end()) {
        DataArray write_set(1, d);
        GeometricRegion dreg = d->region();
        GeometricRegion wreg = GeometricRegion::
            GetIntersection(write_region_, dreg);
        WriteFromCache(write_set, wreg);
        d->UnsetDirtyCacheObject(this);
        write_back_.erase(d);
    }
    ReleaseAccess();
}

/**
 * \details GetDistance(...) gives the cost of using the CacheVar instance,
 * given the read set. Current cost function is the sum of geometric sizes of
 * data in the read set.
 * The function does not have a separate argument for write set
 * - if you want write set included in the cost function, either add  another
 * argument or append it to read set that is passed to GetDistance.
 */
cache::distance_t CacheVar::GetDistance(const DataArray &read_set) const {
    cache::distance_t cur_distance = 0;
    for (size_t i = 0; i < read_set.size(); ++i) {
        Data *d = read_set.at(i);
        GeometricRegion dreg = d->region();
        DMap::const_iterator it = data_map_.find(dreg);
        if (it->second == d)
            continue;
        cur_distance += dreg.dx() * dreg.dy() * dreg.dz();
    }
    return cur_distance;
}

void CacheVar::WriteImmediately(const DataArray &write_set) {
    DataArray flush_set;
    for (size_t i = 0; i < write_set.size(); ++i) {
        Data *d = write_set[i];
        if (write_back_.find(d) != write_back_.end())
            flush_set.push_back(d);
    }
    WriteFromCache(flush_set, write_region_);
    for (size_t i = 0; i < flush_set.size(); ++i) {
        Data *d = flush_set[i];
        d->UnsetDirtyCacheObject(this);
        write_back_.erase(d);
    }
}

/**
 * \details FlushCache(...) flushes all data passed to it, and unsets
 * corresponding dirty object mappings. This function should be used by
 * methods of CacheVar only. It does not check if data in flush_set is in
 * the write_back_ set, making it unsafe. It also provides no locking.
 */
void CacheVar::FlushCache(const DataArray &flush_set) {
    WriteFromCache(flush_set, write_region_);
    for (size_t i = 0; i < flush_set.size(); ++i) {
        Data *d = flush_set[i];
        d->UnsetDirtyCacheObject(this);
        write_back_.erase(d);
    }
}

}  // namespace nimbus

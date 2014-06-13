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
    assert(IsAvailable(cache::EXCLUSIVE));
    DataArray write_set(1, d);
    GeometricRegion dreg = d->region();
    GeometricRegion wreg = GeometricRegion::
        GetIntersection(write_region_, dreg);
    WriteFromCache(write_set, wreg);
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

/**
 * \details
 */
void CacheVar::WriteImmediately(const DataArray &write_set) {
    DataArray flush_set;
    for (size_t i = 0; i < write_set.size(); ++i) {
        Data *d = write_set[i];
        if (write_back_.find(d) != write_back_.end())
            flush_set.push_back(d);
    }
    for (size_t i = 0; i < flush_set.size(); ++i) {
        Data *d = flush_set[i];
        d->UnsetDirtyCacheObject(this);
        write_back_.erase(d);
    }
    WriteFromCache(flush_set, write_region_);
}

bool CacheVar::CheckWritePendingFlag(const DataArray &write_set,
                                     GeometricRegion &write_region) {
    for (size_t i = 0; i < write_set.size(); ++i) {
        Data *d = write_set.at(i);
        if (d->pending_flag()) {
            return false;
        }
        GeometricRegion dreg = d->region();
        DMap::iterator it = data_map_.find(dreg);
        if (it != data_map_.end()) {
            Data *d_old = it->second;
            if (d_old->pending_flag()) {
                return false;
            }
        }
    }
    return true;
}
/**
 * \details If data is not already in existing data, create the mappings.
 * If it replaces existing data, flush existing data if dirty. Create
 * dirty object mapping with all data in write set and set write region.
 */
// TODO(chinmayee/quhang) add synchronization.
void CacheVar::SetUpWrite(const DataArray &write_set,
                          GeometricRegion &write_region,
                          DataArray* flush) {
    for (size_t i = 0; i < write_set.size(); ++i) {
        Data *d = write_set.at(i);
        d->set_pending_flag();
        GeometricRegion dreg = d->region();
        DMap::iterator it = data_map_.find(dreg);
        if (it != data_map_.end()) {
            Data *d_old = it->second;
            d_old->set_pending_flag();
            if (d_old != d) {
                if (write_back_.find(d_old) != write_back_.end()) {
                    flush->push_back(d_old);
                    write_back_.erase(d_old);
                    d_old->UnsetDirtyCacheObject(this);
                    d_old->UnsetCacheObject(this);
                    assert(d_old->co_size() == 0);
                }
                d_old->UnsetCacheObject(this);
            }
        }
        d->InvalidateMappings();
        data_map_[dreg] = d;
        d->SetUpCacheObject(this);
        write_back_.insert(d);
        d->SetUpDirtyCacheObject(this);
    }
}

void CacheVar::PerformSetUpWrite(const DataArray &write_set,
                                 GeometricRegion &write_region,
                                 const DataArray& flush) {
    GeometricRegion write_region_old = write_region_;
    write_region_ = write_region;
    WriteFromCache(flush, write_region_old);
}


void CacheVar::ReleaseWritePendingFlag(const DataArray &write_set,
                                       const DataArray& flush) {
    for (size_t i = 0; i < write_set.size(); ++i) {
        Data *d = write_set.at(i);
        d->unset_pending_flag();
    }
    for (size_t i = 0; i < flush.size(); ++i) {
        Data *d = flush.at(i);
        d->unset_pending_flag();
    }
}


/**
 * \details For read set, if data is not already in cache object, insert it in
 * diff set to read, and sync it if necessary. If
 * it replaces existing data, flush existing data if dirty. For write set, if
 * data is not already in existing data, just create the mappings.
 * If it replaces existing data, flush existing data if dirty. Finally create
 * dirty object mapping with all data in write set.
 */
void CacheVar::SetUpReadWrite(const DataArray &read_set,
                              const DataArray &write_set,
                              DataArray *flush,
                              DataArray *diff,
                              DataArray *sync,
                              CacheObjects *sync_co) {
    assert(flush != NULL);
    assert(diff != NULL);
    assert(sync != NULL);
    for (size_t i = 0; i < read_set.size(); ++i) {
        Data *d = read_set.at(i);
        GeometricRegion dreg = d->region();
        DMap::iterator it = data_map_.find(dreg);
        if (it == data_map_.end()) {
            if (d->dirty_cache_object()) {
                sync->push_back(d);
                sync_co->push_back(d->dirty_cache_object());
                d->ClearDirtyMappings();
                assert(d->co_size() == 1);
            }
            diff->push_back(d);
            data_map_[dreg] = d;
            d->SetUpCacheObject(this);
        } else {
            Data *d_old = it->second;
            if (d_old != d) {
                if (d->dirty_cache_object()) {
                    sync->push_back(d);
                    sync_co->push_back(d->dirty_cache_object());
                    d->ClearDirtyMappings();
                    assert(d->co_size() == 1);
                }
                if (write_back_.find(d_old) != write_back_.end()) {
                    flush->push_back(d_old);
                    write_back_.erase(d_old);
                    d_old->UnsetDirtyCacheObject(this);
                    d_old->UnsetCacheObject(this);
                    assert(d_old->co_size() == 0);
                }
                d_old->UnsetCacheObject(this);
                diff->push_back(d);
                data_map_[dreg] = d;
                d->SetUpCacheObject(this);
            }
        }
    }
    for (size_t i = 0; i < write_set.size(); ++i) {
        Data *d = write_set.at(i);
        GeometricRegion dreg = d->region();
        DMap::iterator it = data_map_.find(dreg);
        if (it != data_map_.end()) {
            Data *d_old = it->second;
            if (d_old != d) {
                if (write_back_.find(d_old) != write_back_.end()) {
                    flush->push_back(d_old);
                    write_back_.erase(d_old);
                    d_old->UnsetDirtyCacheObject(this);
                    d_old->UnsetCacheObject(this);
                    assert(d_old->co_size() == 0);
                }
                d_old->UnsetCacheObject(this);
            }
        }
        d->InvalidateMappings();
        data_map_[dreg] = d;
        d->SetUpCacheObject(this);
        write_back_.insert(d);
        d->SetUpDirtyCacheObject(this);
    }
    for (size_t i = 0; i < flush->size(); ++i) {
        flush->at(i)->set_pending_flag();
    }
    for (size_t i = 0; i < diff->size(); ++i) {
        diff->at(i)->set_pending_flag();
    }
    for (size_t i = 0; i < sync->size(); ++i) {
        sync->at(i)->set_pending_flag();
    }
    for (size_t i = 0; i < sync_co->size(); ++i) {
        sync_co->at(i)->set_pending_flag();
    }
}

bool CacheVar::CheckPendingFlag(const DataArray &read_set,
                                const DataArray &write_set) {
    if (pending_flag()) {
        return false;
    }
    // TODO(chinmayee/quhang) some checkings are not required.
    for (size_t i = 0; i < read_set.size(); ++i) {
        Data *d = read_set.at(i);
        if (d->pending_flag()) {
            return false;
        }
        GeometricRegion dreg = d->region();
        DMap::iterator it = data_map_.find(dreg);
        if (it != data_map_.end()) {
            Data *d_old = it->second;
            if (d_old->pending_flag()) {
                return false;
            }
        }
        if (d->dirty_cache_object()) {
            if (d->dirty_cache_object()->pending_flag()) {
                return false;
            }
        }
    }
    for (size_t i = 0; i < write_set.size(); ++i) {
        Data *d = write_set.at(i);
        if (d->pending_flag()) {
            return false;
        }
        GeometricRegion dreg = d->region();
        DMap::iterator it = data_map_.find(dreg);
        if (it != data_map_.end()) {
            Data *d_old = it->second;
            if (d_old->pending_flag()) {
                return false;
            }
        }
    }
    return true;
}

void CacheVar::ReleasePendingFlag(DataArray *flush,
                                  DataArray *diff,
                                  DataArray *sync,
                                  CacheObjects *sync_co) {
    assert(flush != NULL);
    assert(diff != NULL);
    assert(sync != NULL);
    for (size_t i = 0; i < flush->size(); ++i) {
        flush->at(i)->unset_pending_flag();
    }
    for (size_t i = 0; i < diff->size(); ++i) {
        diff->at(i)->unset_pending_flag();
    }
    for (size_t i = 0; i < sync->size(); ++i) {
        sync->at(i)->unset_pending_flag();
    }
    for (size_t i = 0; i < sync_co->size(); ++i) {
        sync_co->at(i)->unset_pending_flag();
    }
}

}  // namespace nimbus

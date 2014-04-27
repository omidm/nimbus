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

#include <set>
#include <string>

#include "data/cache/cache_object.h"
#include "data/cache/utils.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "worker/data.h"

namespace nimbus {

CacheObject::CacheObject(std::string type,
                         const GeometricRegion &app_object_region)
     : type_(type),
       app_object_region_(app_object_region),
       users_(0) {
}

void CacheObject::ReadToCache(const DataArray &read_set,
                              const GeometricRegion &reg) {
    dbg(DBG_ERROR, "CacheObject Read method not imlemented\n");
}

void CacheObject::ReadDiffToCache(const DataArray &read_set,
                                  const DataArray &diff,
                                  const GeometricRegion &reg,
                                  bool all_lids_diff) {
    dbg(DBG_ERROR, "CacheObject Read method not imlemented\n");
}

void CacheObject::Read(const DataArray &read_set,
                       const GeometricRegion &reg,
                       bool read_all_or_none) {
    if (users_ > 1) {
        dbg(DBG_ERROR, "Cache object being shared!");
        exit(-1);
    }
    if (!write_back_.empty()) {
        DataArray flush;
        // TODO(Chinmayee): this is terrible. we cannot change to region
        // because of particles (4 different types of particles share a
        // region). To change to region, we need a group type instead of a
        // simple cache object.
        LIDSet read_lids;
        for (size_t k = 0; k < read_set.size(); ++k) {
            Data *dd = read_set[k];
            read_lids.insert(dd->logical_id());
        }
        std::set<Data *>::iterator iter = write_back_.begin();
        for (; iter != write_back_.end(); ++iter) {
            Data *d = *iter;
            if (read_lids.contains(d->logical_id())) {
                flush.push_back(d);
            }
        }
        FlushCacheData(flush);
    }
    DataArray diff;
    bool all_lids_diff = true;
    for (size_t i = 0; i < read_set.size(); ++i) {
        Data *d = read_set[i];
        if (!pids_.contains(d->physical_id())) {
            diff.push_back(d);
            d->UpdateData(false);
        }
        if (element_map_.find(d->logical_id()) != element_map_.end())
            all_lids_diff = false;
    }
    if (read_all_or_none) {
        if (!diff.empty())
            ReadDiffToCache(read_set, diff, reg, all_lids_diff);
    } else {
        dbg(DBG_WARN, "\n--- Reading %i out of %i\n", diff.size(), read_set.size());
        if (!diff.empty())
            ReadToCache(diff, reg);
    }
}

void CacheObject::WriteFromCache(const DataArray &write_set,
                                 const GeometricRegion &reg) const {
    dbg(DBG_ERROR, "CacheObject Write method not imlemented\n");
}

void CacheObject::WriteImmediately(const DataArray &write_set,
                                   const GeometricRegion &reg,
                                   bool release) {
    DataArray final_write;
    for (size_t i = 0; i < write_set.size(); ++i) {
        if (write_back_.find(write_set[i]) != write_back_.end()) {
            final_write.push_back(write_set[i]);
        }
    }
    write_region_ = reg;
    FlushCacheData(final_write);
    if (release)
        ReleaseAccess();
}

void CacheObject::Write(const GeometricRegion &reg, bool release) {
    write_region_ = reg;
    // FlushCache();
    if (release)
        ReleaseAccess();
}

void CacheObject::FlushCacheData(const DataArray &diff) {
    WriteFromCache(diff, write_region_);
    for (size_t i = 0; i < diff.size(); ++i) {
        Data *d = diff[i];
        d->clear_dirty_cache_object();
        write_back_.erase(d);
    }
}

void CacheObject::FlushCache() {
    if (write_back_.empty()) {
        return;
    }
    DataArray write_set(write_back_.begin(), write_back_.end());
    WriteFromCache(write_set, write_region_);
    std::set<Data *>::iterator iter = write_back_.begin();
    for (; iter != write_back_.end(); ++iter) {
        Data *d = *iter;
        d->clear_dirty_cache_object();
    }
    write_back_.clear();
}

void CacheObject::PullIntoData(Data *d, bool lock_co) {
    if (lock_co)
        AcquireAccess(EXCLUSIVE);
    if (write_back_.find(d) == write_back_.end()) {
        dbg(DBG_ERROR, "Write back set does not contain data that needs to be pulled\n");
        exit(-1);
    }
    DataArray write;
    write.push_back(d);
    GeometricRegion dreg = d->region();
    WriteFromCache(write, dreg);
    d->clear_dirty_cache_object();
    write_back_.erase(d);
    if (lock_co)
        ReleaseAccess();
}

CacheObject *CacheObject::CreateNew(const GeometricRegion &app_object_region) const {
    dbg(DBG_ERROR, "CacheObject CreateNew method not imlemented\n");
    return NULL;
}

std::string CacheObject::type() const {
    return type_;
}

GeometricRegion CacheObject::app_object_region() const {
    return app_object_region_;
}

void CacheObject::AcquireAccess(CacheAccess access) {
    if (users_ != 0 && (access == EXCLUSIVE || access_ == EXCLUSIVE)) {
        dbg(DBG_ERROR, "Error acquiring cache object!!\n");
        assert(false);
    }
    access_ = access;
    users_++;
}

void CacheObject::ReleaseAccess() {
    users_--;
}

void CacheObject::SetUpRead(const DataArray &read_set,
                            bool read_keep_valid) {
    if (read_keep_valid) {
        for (size_t i = 0; i < read_set.size(); ++i) {
            Data *d = read_set[i];
            d->SetUpCacheObjectDataMapping(this);
        }
    } else {
        for (size_t i = 0; i < read_set.size(); ++i) {
            // TODO(Chinmayee): this is broken. Use physical data map at the
            // worker if this really needs to be taken care of.
            Data *d = read_set[i];
            d->UnsetCacheObjectDataMapping(this);
        }
    }
}

void CacheObject::SetUpWrite(const DataArray &write_set) {
    for (size_t i = 0; i < write_set.size(); ++i) {
        Data *d = write_set[i];
        d->InvalidateCacheObjectsDataMapping();
        d->SetUpCacheObjectDataMapping(this);
        d->set_dirty_cache_object(this);
        write_back_.insert(d);
    }
    std::string ple = "particle_levelset_evolution";
    if (type_ == ple)
        dbg(DBG_WARN, "---PLE contains %i ids in the end\n", pids_.size());
}

void CacheObject::SetUpData(Data *d) {
    logical_data_id_t lid = d->logical_id();
    physical_data_id_t pid = d->physical_id();
    if (element_map_.find(lid) != element_map_.end())
        pids_.remove(element_map_[lid]);
    element_map_[lid] = pid;
    pids_.insert(pid);
    data_.insert(d);
}

void CacheObject::UnsetData(Data *d) {
    logical_data_id_t lid = d->logical_id();
    physical_data_id_t pid = d->physical_id();
    pids_.remove(pid);
    element_map_.erase(lid);
    data_.erase(d);
}

void CacheObject::InvalidateCacheObject(const DataArray &da) {
    for (size_t i = 0; i < da.size(); ++i) {
        Data *d = da[i];
        d->UnsetCacheObjectDataMapping(this);
    }
    std::string ple = "particle_levelset_evolution";
    if (type_ == ple)
        dbg(DBG_WARN, "---PLE contains %i ids in the end\n", pids_.size());
}

void CacheObject::InvalidateCacheObjectComplete() {
    std::set<Data *> temp = data_;
    std::set<Data *>::iterator iter = temp.begin();
    for (; iter != temp.end(); ++iter) {
        Data *d = *iter;
        d->UnsetCacheObjectDataMapping(this);
    }
    std::string ple = "particle_levelset_evolution";
    if (type_ == ple)
        dbg(DBG_WARN, "---PLE contains %i ids in the end\n", pids_.size());
}

distance_t CacheObject::GetDistance(const DataArray &data_set) const {
    distance_t cur_distance = 0;
    for (size_t i = 0; i < data_set.size(); ++i) {
        Data *d = data_set[i];
        if (!pids_.contains(d->physical_id()))
            cur_distance++;
    }
    return cur_distance;
}

bool CacheObject::IsAvailable(CacheAccess access) const {
    return ((access == EXCLUSIVE && users_ == 0) ||
            (access == SHARED && access_ == SHARED));
}

}  // namespace nimbus

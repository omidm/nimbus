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
#include "worker/job.h"

namespace nimbus {

CacheObject::CacheObject(std::string type,
                         const GeometricRegion &global_region,
                         const GeometricRegion &local_region)
     : type_(type),
       global_region_(global_region),
       local_region_(local_region),
       users_(0) {
}

void CacheObject::ReadToCache(const DataArray &read_set) {
    dbg(DBG_ERROR, "CacheObject Read method not imlemented\n");
}

void CacheObject::Read(const DataArray &read_set) {
    DataArray read;
    for (size_t i = 0; i < read_set.size(); ++i) {
        Data *d = read_set[i];
        if (!pids_.contains(d->physical_id()))
            read.push_back(d);;
    }
    if (!read.empty())
        ReadToCache(read);
}

void CacheObject::WriteFromCache(const DataArray &write_set) const {
    dbg(DBG_ERROR, "CacheObject Write method not imlemented\n");
}

void CacheObject::Write() const {
    // TODO(chinmayee): remove pointer from data to cache object
    WriteFromCache(write_back_);
}

CacheObject *CacheObject::CreateNew(const GeometricRegion &local_region) const {
    dbg(DBG_ERROR, "CacheObject CreateNew method not imlemented\n");
    return NULL;
}

std::string CacheObject::type() const {
    return type_;
}

GeometricRegion CacheObject::local_region() const {
    return local_region_;
}

GeometricRegion CacheObject::global_region() const {
    return global_region_;
}

void CacheObject::AcquireAccess(CacheAccess access) {
    access_ = access;
    if (access_ == SHARED)
        users_++;
    else
        users_ = 1;
}

void CacheObject::ReleaseAccess() {
    users_--;
}

void CacheObject::SetUpRead(const DataArray &read_set,
                            bool read_only_keep_valid) {
    pids_.clear();
    if (read_only_keep_valid) {
        for (size_t i = 0; i < read_set.size(); ++i) {
            Data *d = read_set[i];
            if (read_only_keep_valid) {
                pids_.insert(d->physical_id());
            }
        }
    }
}

void CacheObject::SetUpWrite(const DataArray &write_set) {
    write_back_ = write_set;
    for (size_t i = 0; i < write_set.size(); ++i) {
        Data *d = write_set[i];
        pids_.insert(d->physical_id());
        // TODO(chinmayee): insert pointer from data to cache object
    }
}

distance_t CacheObject::GetDistance(const DataArray &data_set,
                                    CacheAccess access) const {
    distance_t max_distance = 2*data_set.size();
    distance_t cur_distance = 0;
    if (!IsAvailable(access))
        return max_distance;
    for (size_t i = 0; i < data_set.size(); ++i) {
        Data *d = data_set[i];
        if (!pids_.contains(d->physical_id()))
            cur_distance++;
    }
    return cur_distance;
}

void CacheObject::GetReadSet(const Job &job,
                             const DataArray &da,
                             DataArray *read) const {
    dbg(DBG_WARN, "Using base class implementation for GetReadSet, for %s\n",
            type_.c_str());
    PIDSet read_ids = job.read_set();
    for (size_t i = 0; i < da.size(); ++i) {
        if (read_ids.contains(da[i]->physical_id()) &&
                da[i]->name() == type_) {
            read->push_back(da[i]);
        }
    }
}

void CacheObject::GetWriteSet(const Job &job,
                              const DataArray &da,
                              DataArray *write) const {
    dbg(DBG_WARN, "Using base class implementation for GetWriteSet, for %s\n",
            type_.c_str());
    PIDSet write_ids = job.write_set();
    for (size_t i = 0; i < da.size(); ++i) {
        if (write_ids.contains(da[i]->physical_id()) &&
                da[i]->name() == type_) {
            write->push_back(da[i]);
        }
    }
}

bool CacheObject::IsAvailable(CacheAccess access) const {
    return ((access == EXCLUSIVE && users_ == 0) ||
            (access == SHARED && access_ == SHARED));
}

}  // namespace nimbus

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

void CacheObject::Read(const DataArray &read_set,
                       const GeometricRegion &reg) {
    DataArray read;
    for (size_t i = 0; i < read_set.size(); ++i) {
        Data *d = read_set[i];
        if (!pids_.contains(d->physical_id()))
            read.push_back(d);;
    }
    if (!read.empty())
        ReadToCache(read, reg);
}

void CacheObject::WriteFromCache(const DataArray &write_set,
                                 const GeometricRegion &reg) const {
    dbg(DBG_ERROR, "CacheObject Write method not imlemented\n");
}

void CacheObject::Write(const GeometricRegion &reg, bool release) {
    // TODO(chinmayee): remove pointer from data to cache object
    WriteFromCache(write_back_, reg);
    if (release)
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
                            bool read_keep_valid) {
    pids_.clear();
    if (read_keep_valid) {
        for (size_t i = 0; i < read_set.size(); ++i) {
            Data *d = read_set[i];
            pids_.insert(d->physical_id());
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

bool CacheObject::IsAvailable(CacheAccess access) const {
    return ((access == EXCLUSIVE && users_ == 0) ||
            (access == SHARED && access_ == SHARED));
}

}  // namespace nimbus

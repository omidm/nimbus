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
                         const GeometricRegion &region)
     : type_(type), region_(region), users_(0) {
}

void CacheObject::Read(const Data &read) {
    dbg(DBG_ERROR, "CacheObject Read method not imlemented\n");
}

void CacheObject::Read(const DataSet &read_set,
                       bool read_only_valid) {
    for (DataSet::iterator it = read_set.begin();
            it != read_set.end();
            ++it) {
        Data *d = *it;
        Read(*d);
    }
}

void CacheObject::Write(Data *write) const {
    dbg(DBG_ERROR, "CacheObject Write method not imlemented\n");
}

void CacheObject::Write() const {
    for (DataSet::iterator it = write_back_.begin();
            it != write_back_.end();
            ++it) {
        Data *d = *it;
        Write(d);
    }
}

CacheObject *CacheObject::CreateNew() const {
    dbg(DBG_ERROR, "CacheObject CreateNew method not imlemented\n");
    return NULL;
}

std::string CacheObject::type() const {
    return type_;
}

GeometricRegion CacheObject::region() const {
    return region_;
}

void CacheObject::MarkAccess(CacheAccess access) {
    access_ = access;
    if (access_ == SHARED)
        users_++;
    else
        users_ = 1;
}

void CacheObject::MarkWriteBack(const DataSet &write) {
    write_back_ = write;
    for (DataSet::iterator it = write_back_.begin();
            it != write_back_.end();
            ++it) {
        Data *d = *it;
        dbg(DBG_ERROR, "%p\n", d);
    }
}

distance_t CacheObject::GetDistance(const DataSet &data_set,
                                    CacheAccess access) const {
    distance_t max_distance = 2*data_set.size();
    distance_t cur_distance = 0;
    if (!IsAvailable())
        return max_distance;
    for (DataSet::iterator it = data_set.begin();
            it != data_set.end();
            ++it) {
        Data *d = *it;
        if (!pids_.contains(d->physical_id()))
            cur_distance++;
    }
    return cur_distance;
}

bool CacheObject::IsAvailable() const {
    return (access_ == SHARED || users_ == 1);
}

}  // namespace nimbus

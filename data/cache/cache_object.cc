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
#include "worker/data.h"
#include "worker/job.h"

namespace nimbus {

CacheObject::CacheObject(std::string type)
     : type_(type), region_(region), users_(0) {
}

void CacheObject::Read(const DataSet &read_set) {
    for (DataSet::iterator it = read_set.begin();
            it != read.end();
            ++it) {
        Data *d = *it;
        Read(*d);
    }
}

void CacheObject::Write(DataSet *write_set) {
    for (DataSet::iterator it = write_set->begin();
            it != write_set.end();
            ++it) {
        Data *d = *it;
        Write(d);
    }
}

std::string CacheObject::type() {
    return type_;
}

GeometricRegion CacheObject::region() {
    return region;
}

void CacheObject::MarkAccess(CacheAccess access) {
    access_ = access;
    users_++;
}

void CacheObject::MarkWriteBack(DataSet *write) {
    write_back_ = *write;
    for (DataSet::iterator it = write_back_.begin();
            it != write_back_.end();
            ++it) {
        Data *d = *it;
        d->MarkCacheObject(this);
    }
}

distance_t CacheObject::GetDistance(const DataSet &data_set) const {
    distance_t distance = 0;
    // TODO(chinmayee): implement actual min distance algorithm
    return distance;
}

bool CacheObject::IsAvailable() {
    return (users_ == 0);
}

}  // namespace nimbus

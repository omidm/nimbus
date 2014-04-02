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

#include "data/cache/cache_object.h"
#include "data/cache/utils.h"
#include "worker/data.h"
#include "worker/job.h"

namespace nimbus {

CacheObject::CacheObject() : users_(0) {}

bool CacheObject::IsAvailable() {
    return (users_ == 0);
}

distance_t CacheObject::GetDistance(const DataSet &read,
                                    const DataSet &write,
                                    const StringSet &read_var,
                                    const StringSet &write_var) {
    distance_t max_distance = 2 * read.size();
    for (StringSet::iterator it = read_var.begin();
            it != read_var.end(); it++) {
        if (locked_write_.find(*it) != locked_write_.end())
            return max_distance;
    }
    for (StringSet::iterator it = write_var.begin();
            it != write_var.end(); it++) {
        if (locked_read_.find(*it) != locked_read_.end() &&
                locked_write_.find(*it) != locked_write_.end())
            return max_distance;
    }
    distance_t distance = 0;
    // TODO(chinmayee): implement actual min distance algorithm later
    return distance;
}

void CacheObject::LockData(const DataSet &read,
                           const DataSet &write,
                           const StringSet &read_var,
                           const StringSet &write_var) {
    assert(users_ >= 0);
    users_++;
    // TODO(chinmayee): implement locking later
}

}  // namespace nimbus

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

#include <string>

#include "data/cache/cache_object.h"
#include "data/cache/cache_pool.h"
#include "data/cache/cache_table.h"
#include "data/cache/utils.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "worker/data.h"
#include "worker/job.h"

namespace nimbus {

CachePool::CachePool() {}

CacheObject *CachePool::GetCachedObject(const Job &job,
                                        const DataArray &da,
                                        const GeometricRegion &region,
                                        const CacheObject &prototype,
                                        CacheAccess access,
                                        bool read_only_valid) {
    DataSet read,  write;
    GetReadSet(job, da, &read);
    GetWriteSet(job, da, &write);
    CacheObject *co = NULL;
    if (pool_.find(prototype.type()) == pool_.end()) {
        CacheTable *ct = new CacheTable();
        pool_[prototype.type()] = ct;
        co = prototype.CreateNew();
        if (co == NULL) {
            dbg(DBG_ERROR, "Tried to create a cache object for an unimplemented prototype. Exiting ...\n"); // NOLINT
            exit(-1);
        }
        ct->AddEntry(region, co);
    } else {
        CacheTable *ct = pool_[prototype.type()];
        co = ct->GetClosestAvailable(region, read, access);
        if (co == NULL)
            ct->AddEntry(region, co);
    }
    co->MarkAccess(access);
    co->Read(read, read_only_valid);
    co->MarkWriteBack(write);
    return co;
}

void CachePool::GetReadSet(const Job &job,
                           const DataArray &da,
                           DataSet *read) {
    PIDSet read_ids = job.read_set();
    for (size_t i = 0; i < da.size(); ++i) {
        if (read_ids.contains(da[i]->physical_id())) {
            read->insert(da[i]);
        }
    }
}

void CachePool::GetWriteSet(const Job &job,
                            const DataArray &da,
                            DataSet *write) {
    PIDSet write_ids = job.write_set();
    for (size_t i = 0; i < da.size(); ++i) {
        if (write_ids.contains(da[i]->physical_id())) {
            write->insert(da[i]);
        }
    }
}

}  // namespace nimbus

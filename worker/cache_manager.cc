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
#include "data/cache/cache_table.h"
#include "data/cache/utils.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "worker/cache_manager.h"
#include "worker/data.h"

namespace nimbus {

CacheManager::CacheManager() {
    pool_ = new Pool();
}

CacheObject *CacheManager::GetAppObject(const DataArray &read,
                                        const DataArray &write,
                                        const GeometricRegion &region,
                                        const CacheObject &prototype,
                                        CacheAccess access,
                                        bool read_only_keep_valid) {
    CacheObject *co = NULL;
    if (pool_->find(prototype.type()) == pool_->end()) {
        CacheTable *ct = new CacheTable();
        (*pool_)[prototype.type()] = ct;
        co = prototype.CreateNew(region);
        if (co == NULL) {
            dbg(DBG_ERROR, "Tried to create a cache object for an unimplemented prototype. Exiting ...\n"); // NOLINT
            exit(-1);
        }
        ct->AddEntry(region, co);
    } else {
        CacheTable *ct = (*pool_)[prototype.type()];
        co = ct->GetClosestAvailable(region, read, access);
        if (co == NULL)
            ct->AddEntry(region, co);
    }
    co->AcquireAccess(access);
    co->Read(read, region);
    co->SetUpRead(read, read_only_keep_valid | (access != EXCLUSIVE));
    co->SetUpWrite(write);
    return co;
}

CacheObject *CacheManager::GetAppObject(const DataArray &read,
                                        const DataArray &write,
                                        const GeometricRegion &data_region,
                                        const GeometricRegion &read_region,
                                        const CacheObject &prototype,
                                        CacheAccess access,
                                        bool read_only_keep_valid) {
    CacheObject *co = NULL;
    if (pool_->find(prototype.type()) == pool_->end()) {
        CacheTable *ct = new CacheTable();
        (*pool_)[prototype.type()] = ct;
        co = prototype.CreateNew(data_region);
        if (co == NULL) {
            dbg(DBG_ERROR, "Tried to create a cache object for an unimplemented prototype. Exiting ...\n"); // NOLINT
            exit(-1);
        }
        ct->AddEntry(data_region, co);
    } else {
        CacheTable *ct = (*pool_)[prototype.type()];
        co = ct->GetClosestAvailable(data_region, read, access);
        if (co == NULL)
            ct->AddEntry(data_region, co);
    }
    co->AcquireAccess(access);
    co->Read(read, read_region);
    // if access is shared, all data needs to be kept valid!
    co->SetUpRead(read, read_only_keep_valid | (access != EXCLUSIVE));
    co->SetUpWrite(write);
    return co;
}

}  // namespace nimbus

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

#include <map>
#include <vector>

#include "data/cache/cache_defs.h"
#include "data/cache/cache_object.h"
#include "data/cache/cache_table.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "worker/cache_manager.h"
#include "worker/data.h"

namespace nimbus {

CacheManager::CacheManager() {
    pool_ = new Pool();
}

CacheVar *CacheManager::GetAppVar(const DataArray &read_set,
                                  const GeometricRegion &read_region,
                                  const DataArray &write_set,
                                  const GeometricRegion &write_region,
                                  const CacheVar &prototype,
                                  const GeometricRegion &region,
                                  cache::CacheAccess access,
                                  bool invalidate_read_minus_write) {
    CacheVar *cv = NULL;
    if (pool_->find(prototype.id()) == pool_->end()) {
        CacheTable *ct = new CacheTable(cache::VAR);
        (*pool_)[prototype.id()] = ct;
        cv = prototype.CreateNew(region);
        assert(cv != NULL);
        ct->AddEntry(region, cv);
    } else {
        CacheTable *ct = (*pool_)[prototype.id()];
        cv = ct->GetClosestAvailable(region, read_set, access);
        if (cv == NULL) {
            cv = prototype.CreateNew(region);
            assert(cv != NULL);
            ct->AddEntry(region, cv);
        }
    }
    cv->AcquireAccess(access);
    cv->UpdateCache(read_set, read_region,
                    write_region,
                    invalidate_read_minus_write);
    cv->SetUpWrite(write_set, write_region);
    return cv;
}

CacheStruct *CacheManager::GetAppStruct(const std::vector<cache::type_id_t> &var_type,
                                        const std::vector<DataArray> &read_sets,
                                        const GeometricRegion &read_region,
                                        const std::vector<DataArray> &write_sets,
                                        const GeometricRegion &write_region,
                                        const CacheStruct &prototype,
                                        const GeometricRegion &region,
                                        cache::CacheAccess access,
                                        bool invalidate_read_minus_write) {
    CacheStruct *cs = NULL;
    if (pool_->find(prototype.id()) == pool_->end()) {
        CacheTable *ct = new CacheTable(cache::STRUCT);
        (*pool_)[prototype.id()] = ct;
        cs = prototype.CreateNew(region);
        assert(cs != NULL);
        ct->AddEntry(region, cs);
    } else {
        CacheTable *ct = (*pool_)[prototype.id()];
        cs = ct->GetClosestAvailable(region, var_type, read_sets, access);
        if (cs == NULL) {
            cs = prototype.CreateNew(region);
            assert(cs != NULL);
            ct->AddEntry(region, cs);
        }
    }
    cs->AcquireAccess(access);
    cs->UpdateCache(var_type, read_sets, read_region,
                    write_region,
                    invalidate_read_minus_write);
    cs->SetUpWrite(var_type, write_sets, write_region);
    return cs;
}

}  // namespace nimbus

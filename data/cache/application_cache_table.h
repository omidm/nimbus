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
 * Cache table for caching application state, to avoid reconstructing the
 * state and translating data for every job.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#ifndef NIMBUS_DATA_CACHE_APPLICATION_CACHE_TABLE_H_
#define NIMBUS_DATA_CACHE_APPLICATION_CACHE_TABLE_H_

#include <map>

#include "data/cache/application_cache_entry.h"

namespace nimbus {

class ApplicationCacheTable {
    public:
        ApplicationCacheTable();
    private:
        // tyedef CacheMap
        typedef std::map<GeometricRegion,
                         ApplicationCacheEntry *,
                         bool(*)(const GeometricRegion&, // NOLINT
                                 const GeometricRegion&)>
                        CacheMap;
        // key comparator function for CacheMap
        static bool GeometricRegionLess(const GeometricRegion &r1,
                                        const GeometricRegion &r2);
        // table of cache entries
        CacheMap table_;
};  // class ApplicationCacheTable

}  // namespace nimbus

#endif  // NIMBUS_DATA_CACHE_APPLICATION_CACHE_TABLE_H_

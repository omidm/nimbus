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

#ifndef NIMBUS_WORKER_CACHE_MANAGER_H_
#define NIMBUS_WORKER_CACHE_MANAGER_H_

#include <map>
#include <string>

#include "data/cache/cache_object.h"
#include "data/cache/cache_table.h"
#include "data/cache/utils.h"
#include "shared/geometric_region.h"
#include "worker/data.h"

namespace nimbus {

class CacheManager {
    public:
        CacheManager();

        CacheObject *GetAppObject(const DataArray &read,
                                  const DataArray &write,
                                  const GeometricRegion &region,
                                  const CacheObject &prototype,
                                  CacheAccess access = EXCLUSIVE,
                                  bool read_keep_valid = false);

        CacheObject *GetAppObject(const DataArray &read,
                                  const DataArray &write,
                                  const GeometricRegion &region,
                                  const GeometricRegion &read_region,
                                  const CacheObject &prototype,
                                  CacheAccess access = EXCLUSIVE,
                                  bool read_keep_valid = false);

    private:
        typedef std::map<std::string,
                         CacheTable *> Pool;

        Pool *pool_;
};  // class CacheManager

}  // namespace nimbus

#endif  // NIMBUS_WORKER_CACHE_MANAGER_H_

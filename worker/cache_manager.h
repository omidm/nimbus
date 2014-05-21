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
#include <vector>

#include "data/cache/cache_defs.h"
#include "data/cache/cache_struct.h"
#include "data/cache/cache_table.h"
#include "data/cache/cache_var.h"

namespace nimbus {

class Data;
typedef std::vector<Data *> DataArray;
class GeometricRegion;

class CacheManager {
    public:
        CacheManager();

        CacheVar *GetAppVar(const DataArray &read_set,
                            const GeometricRegion &read_region,
                            const DataArray &write_set,
                            const GeometricRegion &write_region,
                            const CacheVar &prototype,
                            const GeometricRegion &region,
                            cache::CacheAccess access,
                            bool invalidate_read_minus_write = false);

        CacheStruct *GetAppStruct(const std::vector<cache::type_id_t> &var_type,
                                  const std::vector<DataArray> &read_sets,
                                  const GeometricRegion &read_region,
                                  const std::vector<DataArray> &write_sets,
                                  const GeometricRegion &write_region,
                                  const CacheStruct &prototype,
                                  const GeometricRegion &region,
                                  cache::CacheAccess access,
                                  bool invalidate_read_minus_write = false);

    private:
        typedef std::map<cache::co_id_t,
                         CacheTable *> Pool;
        Pool *pool_;
};  // class CacheManager

}  // namespace nimbus

#endif  // NIMBUS_WORKER_CACHE_MANAGER_H_

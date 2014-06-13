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
 * CacheManager is the interface for application jobs to cache. Application
 * jobs send their request for application objects to the cache manager.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#ifndef NIMBUS_WORKER_CACHE_MANAGER_H_
#define NIMBUS_WORKER_CACHE_MANAGER_H_

#include <pthread.h>
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

/**
 * \class CacheManager
 * \details CacheManager is the interface through which application jobs request
 * application objects. Internally, CacheManager contains a two level map - the
 * first level key is a prototype id, and the second level key is an
 * application object geometric region. Together, these keys map to a list of
 * CacheVars or CacheStructs.
 */
class CacheManager {
    public:
        pthread_mutex_t cache_lock;
        pthread_cond_t cache_cond;
        /**
         * \brief Creates a CacheManager instance
         * \return Constructed CacheManager instance
         */
        CacheManager();

        /**
         * \brief Requests a CacheVar instance of type prototype, from the
         * CacheManager
         * \param read_set specifies the data that should be read in the
         * CacheVar instance
         * \param read_region indicates read region
         * \param write_set specifies the data that should be marked as being
         * written
         * \param write_region indicates write region
         * \param prototype represents the application object type
         * \param region is the application object region
         * \access indicates whether application object access should be
         * EXCLUSIVE or SHARED
         * \return A pointer to a CacheVar instance that application can use
         */
        CacheVar *GetAppVar(const DataArray &read_set,
                            const GeometricRegion &read_region,
                            const DataArray &write_set,
                            const GeometricRegion &write_region,
                            const CacheVar &prototype,
                            const GeometricRegion &region,
                            cache::CacheAccess access);

        /**
         * \brief Requests a CacheStruct instance of type prototype, from the
         * CacheManager
         * \param var_type specifies the application variable types for the
         * lists in read_sets/ write_sets
         * \param read_sets specifies the data that should be read in the
         * CacheVar instance
         * \param read_region indicates read region
         * \param write_sets specifies the data that should be marked as being
         * written
         * \param write_region indicates write region
         * \param prototype represents the application object type
         * \param region is the application object region
         * \access indicates whether application object access should be
         * EXCLUSIVE or SHARED
         * \return A pointer to a CacheVar instance that application can use
         */
        CacheStruct *GetAppStruct(const std::vector<cache::type_id_t> &var_type,
                                  const std::vector<DataArray> &read_sets,
                                  const GeometricRegion &read_region,
                                  const std::vector<DataArray> &write_sets,
                                  const GeometricRegion &write_region,
                                  const CacheStruct &prototype,
                                  const GeometricRegion &region,
                                  cache::CacheAccess access);

        /**
         * \brief If data is dirty, syncs with corresponding dirty cache
         * object, and clears the dirty mapping
         * \param Data d to sync
         */
        void SyncData(Data *d);

        /**
         * \brief Invalidates mapping between d and all cache objects (dirty or
         * not dirty)
         * \param Data d
         */
        void InvalidateMappings(Data *d);

    private:
        typedef std::map<cache::co_id_t,
                         CacheTable *> Pool;
        Pool *pool_;
};  // class CacheManager

}  // namespace nimbus

#endif  // NIMBUS_WORKER_CACHE_MANAGER_H_

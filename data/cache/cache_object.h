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
 * Application field contains a pointer to the field in the cached object, and
 * a list of field partitions corresponding to logical data partitions for the
 * field.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#ifndef NIMBUS_DATA_CACHE_CACHE_OBJECT_H_
#define NIMBUS_DATA_CACHE_CACHE_OBJECT_H_

#include <map>
#include <string>
#include <vector>

#include "data/cache/utils.h"
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"
#include "worker/data.h"

namespace nimbus {

enum CacheAccess { SHARED, EXCLUSIVE };

class CacheObject {
    public:
        explicit CacheObject(std::string type,
                             const GeometricRegion &app_object_region);

        virtual void ReadToCache(const DataArray &read_set, const GeometricRegion &reg);
        void Read(const DataArray &read_set, const GeometricRegion &reg, bool read_all = false);
        virtual void WriteFromCache(const DataArray &write_set, const GeometricRegion &reg) const;
        void Write(const GeometricRegion &reg, bool release = true);

        virtual CacheObject *CreateNew(const GeometricRegion &app_object_region) const;

        std::string type() const;
        GeometricRegion app_object_region() const;

        void AcquireAccess(CacheAccess access);
        void ReleaseAccess();

        void SetUpRead(const DataArray &read_set,
                       bool read_keep_valid);
        void SetUpWrite(const DataArray &write_set);
        void InvalidateCacheObject(Data *d);

        bool IsAvailable(CacheAccess access) const;
        distance_t GetDistance(const DataArray &data_set,
                               CacheAccess access = EXCLUSIVE) const;

    private:
        std::string type_;
        GeometricRegion app_object_region_;

        CacheAccess access_;
        int users_;
        bool read_valid_;
        bool write_valid_;

        DataArray write_back_;

        /* Currently, cache object contains only physical id information.
         * Distance (cost) information and validity checks are based on
         * physical id only.
         * TODO(chinmayee): change this later to use logical id & version
         * information.*/
        PIDSet pids_;
        std::map<logical_data_id_t, physical_data_id_t> element_map_;
};  // class CacheObject

typedef std::vector<CacheObject *> CacheObjects;

}  // namespace nimbus

#endif  // NIMBUS_DATA_CACHE_CACHE_OBJECT_H_

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

#ifndef NIMBUS_DATA_CACHE_CACHE_STRUCT_H_
#define NIMBUS_DATA_CACHE_CACHE_STRUCT_H_

#include <map>
#include <set>
#include <string>
#include <vector>

#include "data/cache/utils.h"
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"
#include "worker/data.h"

namespace nimbus {

class CacheStruct {
    public:
        /**
         * \fn CacheStruct::CacheStruct(std::string type,
         *                              size_t num_variables,
         *                              const GeometricRegion *struct_region)
         * \param type specifies the application object/ variable name
         * \param num_variables states the number of different application
         * variables that the CacheStruct instance contains
         * \param struct_region specifies the spatial domain of the CacheStruct
         * instance
         * \return Constructed CacheStruct instance
         */
        explicit CacheStruct(std::string type,
                             size_t num_variables,
                             const GeometricRegion &struct_region);

        /**
         * \fn CacheStruct::CreateNew(const GeometricRegion &struct_region) const
         * \param struct_region specifies the spatial domain of the CacheStruct
         * instance
         * \return Returns a pointer to the newly allocated CacheStruct instance
         */
        virtual CacheStruct *CreateNew(const GeometricRegion &struct_region) const;

        /**
         * \fn CacheStruct::ReadToCache(const Data
         */
        virtual void ReadToCache(const DataArray &read_set,
                                 const GeometricRegion &read_region);
        virtual void WriteFromCache(const DataArray &write_set,
                                    const GeometricRegion &reg) const;

        void Read(const DataArray &read_set,
                  const GeometricRegion &read_region,
                  bool release = false);
        void WriteImmediately(const DataArray &write_set,
                              const GeometricRegion &write_region,
                              bool release = true);
        void Write(const GeometricRegion &write_region,
                   bool release = true);

        // update data
        void FlushCache();
        void PullIntoData(Data *d);
        void RemoveFromWriteBack(Data *d);

        std::string type() const;
        GeometricRegion app_object_region() const;

        void AcquireAccess(CacheAccess access);
        void ReleaseAccess();

        void SetUpRead(const DataArray &read_set,
                       bool read_keep_valid);
        void SetUpWrite(const DataArray &write_set);
        void SetUpData(Data *d);
        void UnsetData(Data *d);

        // Use with care
        void InvalidateCacheStruct(const DataArray &da);
        void InvalidateCacheStructComplete();

        bool IsAvailable(CacheAccess access) const;
        distance_t GetDistance(const DataArray &data_set) const;

    private:
        void FlushCacheData(const DataArray &diff);

        std::string type_;
        GeometricRegion app_object_region_;

        CacheAccess access_;
        int users_;

        GeometricRegion write_region_;
        std::set<Data *> write_back_;

        /* Currently, cache object contains only physical id information.
         * Distance (cost) information and validity checks are based on
         * physical id only.
         * information.*/
        // TODO(Chinmayee): change this to region-pid map (because of
        // scratches)
        std::map<logical_data_id_t, physical_data_id_t> element_map_;
        std::map<logical_data_id_t, Data*> data_map_;
        std::set<Data *> data_;
        PIDSet pids_;
};  // class CacheStruct

typedef std::vector<CacheStruct *> CacheStructs;

}  // namespace nimbus

#endif  // NIMBUS_DATA_CACHE_CACHE_STRUCT_H_

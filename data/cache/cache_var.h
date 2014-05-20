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
 * A CacheVar is an application object corresponding to one nimbus
 * variable, cached by the cache manager.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#ifndef NIMBUS_DATA_CACHE_CACHE_VAR_H_
#define NIMBUS_DATA_CACHE_CACHE_VAR_H_

#include <map>
#include <set>
#include <string>
#include <vector>

#include "data/cache/cache_defs.h"
#include "data/cache/cache_object.h"
#include "shared/nimbus_types.h"

namespace nimbus {

class Data;
typedef std::vector<Data *> DataArray;
typedef std::set<Data *> DataSet;
class GeometricRegion;

/**
 * \class CacheVar
 * \details Application object corresponding to one numbus variable cached by
 * the cache manager. A CacheVar can must be defined over a geometric region.
 */
class CacheVar : public CacheObject {
    public:
        /**
         * \brief Creates a CacheVar
         * \return Constructed CacheVar instance
         */
        explicit CacheVar();

        /**
         * \brief Creates a new CacheVar instance using current instance
         * parameters
         * \param struct_region specifies the spatial domain of the CacheVar
         * instance
         * \return Returns a pointer to the newly allocated CacheVar instance
         * \details This is a virtual function that must be over-written by application
         * writer. When CacheManager cannot satisfy an application object request,
         * using the instances it has already cached, it calls CreateNew(..) on the
         * prototype passed in the request.
         */
        virtual CacheVar *CreateNew(const GeometricRegion &struct_region) const = 0;

        /**
         * \brief Updates CacheVar with data from read_set - performs update
         * only for those physical data that have changed since the CacheVar
         * instance was last updated
         * \param read_set is a data array corresponding to nimbus variables
         * \param read_region is the geometric region to read
         * \param write_region is the geometric region to write (used if
         * invalidate_read_minus_write = true)
         * \param invalidate_read_mibus_write indicates whether to keep the
         * data read into cache object, outside requested write_region valid.
         * This is needed because some application functions will write to
         * cache object regions that are not in the application job's write
         * set.
         */
        void UpdateCache(const DataArray &read_set,
                         const GeometricRegion &read_region,
                         const GeometricRegion &write_region,
                         bool invalidate_read_minus_write = false);

        /**
         * \brief Sets up mappings for data in write_sets, and also sets up the
         * write region
         * \param write_set is a data array corresponding to nimbus variables
         * \param write_region is the geometric region to write
         */
        void SetUpWrite(const DataArray &write_sets,
                        const GeometricRegion &write_region);

        /**
         * \brief Pulls data from cache, removes corresponding dirty data
         * mapping. Locks the cache object when pulling the data.
         * \param d is data to flush to
         */
        virtual void PullData(Data *d);

        /**
         * \brief Calculates distance of a CacheVar, from the given read_set.
         * This distance indicates the cost of reconstruction if
         * this CacheVar instance is used.
         * \param read_set is a data array corresponding to nimbus variables
         * \return Returns distance (cost)
         */
        cache::distance_t GetDistance(const DataArray &read_set) const;

        /**
         * \brief Unsets mapping between data and CacheVar instance
         * \param d denotes the data to unmap
         */
        virtual void UnsetData(Data *d);

    private:
        /**
         * \brief Flushes data from cache, removes corresponding dirty data
         * mapping
         * \param d is data to flush to
         */
        void FlushToData(Data *d);

        /**
         * \brief Flushes data from cache to data in flush_set (immediately)
         * \param flush_set is a data array corresponding to nimbus variables
         */
        void FlushCache(const DataArray &flush_set);

        // cache-data mappings
        typedef std::map<GeometricRegion,
                         Data *> DMap;
        DMap data_map_;
        DataSet write_back_;

    protected:
        /**
         * \brief Reads data from read_sets into CacheVar instance
         * \param read_set is a data array corresponding to nimbus variables
         * \param read_region is the geometric region to read
         * \details This function must be overwritten by the application
         * writer. It provides the transformation from a set of nimbus data to
         * (application) cached instance.
         */
        virtual void ReadToCache(const DataArray &read_set,
                                 const GeometricRegion &read_region) = 0;

        /**
         * \brief Writes data from CacheVar instance to write_sets
         * \param write_set is a data array corresponding to nimbus variables
         * \param write_region is the geometric region to write
         * \details This function must be overwritten by the application
         * writer. It provides the transformation from (application) cached
         * instance to nimbus data.
         */
        virtual void WriteFromCache(const DataArray &write_set,
                                    const GeometricRegion &write_region) = 0;
};  // class CacheVar

typedef std::vector<CacheVar *> CacheVars;

}  // namespace nimbus

#endif  // NIMBUS_DATA_CACHE_CACHE_VAR_H_

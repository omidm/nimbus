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
    friend class CacheManager;
    public:
        /**
         * \brief Creates a CacheVar
         * \return Constructed CacheVar instance
         */
        explicit CacheVar();

        /**
         * \brief Creates a CacheVar
         * \param  ob_reg specifies application object (CacheVar) region
         * \return Constructed CacheVar instance
         */
        explicit CacheVar(const GeometricRegion &ob_reg);

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
         * \brief Unsets mapping between data and CacheVar instance
         * \param d denotes the data to unmap
         */
        virtual void UnsetData(Data *d);

        /**
         * \brief Unsets dirty data mapping between data and CacheVar
         * instance
         * \param d denotes the data to unmap
         */
        virtual void UnsetDirtyData(Data *d);

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
         * \brief Writes data from cache to data in write_set immediately
         * \param write_set is a data array corresponding to nimbus variables
         */
        void WriteImmediately(const DataArray &write_set);

        /**
         * \brief Edits the write set mappings (dirty mappins), flushing data
         * if necessary
         * \param write_set is a data array corresponding to nimbus variables
         * \param write_region is region to write
         */
        void SetUpWrite(const DataArray &write_set,
                        GeometricRegion &write_region,
                        DataArray* flush);

        bool CheckWritePendingFlag(const DataArray &write_set,
                                   GeometricRegion &write_region);

        void PerformSetUpWrite(const DataArray &write_set,
                               GeometricRegion &write_region,
                               const DataArray& flush);
        void ReleaseWritePendingFlag(const DataArray &write_set,
                                     const DataArray& flush);

    private:
        /**
         * \brief Edits all mappings between cache objects, given read and
         * write set
         * \param read_set is a data array corresponding to nimbus variables
         * \param write_set is a data array corresponding to nimbus variables
         * \param flush is where SetUpReadWrite puts all data that needs to be
         * flushed
         * \param diff is where SetUpReadWrite puts all data that needs to be
         * read
         * \param sync is where SetUpReadWrite puts all data that needs to be
         * synced from other cache objects
         * \param sync_co is where SetUpReadWrite puts all cache objects that
         * need to be synced with data in sync, in the same order
         */
        void SetUpReadWrite(const DataArray &read_set,
                            const DataArray &write_set,
                            DataArray *flush,
                            DataArray *diff,
                            DataArray *sync,
                            CacheObjects *sync_co);
        bool CheckPendingFlag(const DataArray &read_set,
                              const DataArray &write_set);
        void ReleasePendingFlag(DataArray *flush,
                                DataArray *diff,
                                DataArray *sync,
                                CacheObjects *sync_co);

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
                                    const GeometricRegion &write_region) const = 0;
};  // class CacheVar

typedef std::vector<CacheVar *> CacheVars;

}  // namespace nimbus

#endif  // NIMBUS_DATA_CACHE_CACHE_VAR_H_

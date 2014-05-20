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
 * A CacheStruct is an application object corresponding to a set of nimbus
 * variables cached by the cache manager.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#ifndef NIMBUS_DATA_CACHE_CACHE_STRUCT_H_
#define NIMBUS_DATA_CACHE_CACHE_STRUCT_H_

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
 * \class CacheStruct
 * \details Application object corresponding to a set of nimbus variables cached
 * by the cache manager. A CacheStruct can have num_variables number of
 * different nimbus variables, and must be defined over a geometric region.
 * Read and write calls contain a type, corresponding to each read/ write set/
 * data. This type is an integer, between zero (inclusive) and num_variables_
 * (exclusive).
 */
class CacheStruct : public CacheObject {
    public:
        /**
         * \brief Creates a CacheStruct
         * \param num_variables indicates number of different application
         * variables that the CacheStruct instance contains
         * \return Constructed CacheStruct instance
         */
        explicit CacheStruct(size_t num_variables);

        /**
         * \brief Creates a new CacheStruct instance using current instance
         * parameters
         * \param struct_region specifies the spatial domain of the CacheStruct
         * instance
         * \return Returns a pointer to the newly allocated CacheStruct instance
         * \details This is a virtual function that must be over-written by application
         * writer. When CacheManager cannot satisfy an application object request,
         * using the instances it has already cached, it calls CreateNew(..) on the
         * prototype passed in the request.
         */
        virtual CacheStruct *CreateNew(const GeometricRegion &struct_region) const = 0;

        /**
         * \brief Updates CacheStruct with data from read_sets - performs update
         * only for those physical data that have changed since the CacheStruct
         * instance was created
         * \param var_type is a list of type_ids corresponding to nimbus variables
         * \param read_sets is a list of data arrays corresponding to nimbus variables
         * \param read_region is the geometric region to read
         * \param write_region is the geometric region to write (used if
         * invalidate_read_minus_write = true)
         * \param invalidate_read_mibus_write indicates whether to keep the
         * data read into cache object, outside requested write_region valid.
         * This is needed because some application functions will write to
         * cache object regions that are not in the application job's write
         * set.
         */
        void UpdateCache(const std::vector<cache::type_id_t> &var_type,
                         const std::vector<DataArray> &read_sets,
                         const GeometricRegion &read_region,
                         const GeometricRegion &write_region,
                         bool invalidate_read_minus_write = false);

        /**
         * \brief Sets up mappings for data in write_sets, and also sets up the
         * write region
         * \param var_type is a list of type_ids corresponding to nimbus variables
         * \param write_sets is a list of data arrays corresponding to nimbus variables
         * \param write_region is the geometric region to write
         */
        void SetUpWrite(const std::vector<cache::type_id_t> &var_type,
                        const std::vector<DataArray> &write_sets,
                        const GeometricRegion &write_region);

        /**
         * \brief Pulls data from cache, removes corresponding dirty data
         * mapping. Locks the cache object when pulling the data.
         * \param d is data to flush to
         */
        virtual void PullData(Data *d);

        /**
         * \brief Calculates distance of a CacheStruct, from the given list of
         * read_sets. This distance indicates the cost of reconstruction if
         * this CacheStruct instance is used.
         * \param var_type is a list of type_ids corresponding to nimbus variables
         * \param read_sets is a list of data arrays corresponding to nimbus variables
         * \return Returns distance (cost)
         */
        cache::distance_t GetDistance(const std::vector<cache::type_id_t> &var_type,
                                      const std::vector<DataArray> &read_sets) const;

        /**
         * \brief Unsets mapping between data and CacheStruct instance
         * \param d denotes the data to unmap
         */
        virtual void UnsetData(Data *d);

    private:
        /**
         * \brief Flushes data from cache, removes corresponding dirty data
         * mapping
         * \param d is data to flush to
         * \param t denotes the type of nimbus variable,
         * as explained in the class description
         */
        void FlushToData(Data *d, cache::type_id_t t);

        /**
         * \brief Flushes data from cache to data in flush_sets (immediately)
         * \param var_type is a list of type_ids corresponding to nimbus variables,
         * as explained in the class description
         * \param flush_sets is a list of data arrays corresponding to nimbus variables
         */
        void FlushCache(const std::vector<cache::type_id_t> &var_type,
                        const std::vector<DataArray> &flush_sets);

        // number of nimbus variables
        size_t num_variables_;

        // cache-data mappings
        typedef std::map<GeometricRegion,
                         Data *> DMap;
        std::vector<DMap> data_maps_;
        std::vector<DataSet> write_backs_;

    protected:
        /**
         * \brief Reads data from read_sets into CacheStruct instance
         * \param var_type is a list of type_ids corresponding to nimbus variables,
         * as explained in the class description
         * \param read_sets is a list of data arrays corresponding to nimbus variables
         * \param read_region is the geometric region to read
         * \details This function must be overwritten by the application
         * writer. It provides the transformation from a set of nimbus data to
         * (application) cached instance.
         */
        virtual void ReadToCache(const std::vector<cache::type_id_t> &var_type,
                                 const std::vector<DataArray> &read_sets,
                                 const GeometricRegion &read_region) = 0;

        /**
         * \brief Writes data from CacheStruct instance to write_sets
         * \param var_type is a list of type_ids corresponding to nimbus variables,
         * as explained in the class description
         * \param write_sets is a list of data arrays corresponding to nimbus variables
         * \param write_region is the geometric region to be write
         * \details This function must be overwritten by the application
         * writer. It provides the transformation from (application) cached
         * instance to nimbus data.
         */
        virtual void WriteFromCache(const std::vector<cache::type_id_t> &var_type,
                                    const std::vector<DataArray> &write_sets,
                                    const GeometricRegion &write_region) = 0;
};  // class CacheStruct

typedef std::vector<CacheStruct *> CacheStructs;

}  // namespace nimbus

#endif  // NIMBUS_DATA_CACHE_CACHE_STRUCT_H_

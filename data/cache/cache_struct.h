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
    friend class CacheManager;
    friend class CacheTable;
    public:
        /**
         * \brief Creates a CacheStruct
         * \param num_variables indicates number of different application
         * variables that the CacheStruct instance contains
         * \return Constructed CacheStruct instance
         */
        explicit CacheStruct(size_t num_variables);

        /**
         * \brief Creates a CacheStruct
         * \param num_variables indicates number of different application
         * variables that the CacheStruct instance contains
         * \param  ob_reg specifies application object (CacheStruct) region
         * \return Constructed CacheStruct instance
         */
        explicit CacheStruct(size_t num_variables, const GeometricRegion &ob_reg);

        /**
         * TODO: Make this protected
         * \brief Unsets mapping between data and CacheStruct instance
         * \param d denotes the data to unmap
         */
        virtual void UnsetData(Data *d);

        /**
         * TODO: Make this protected
         * \brief Unsets dirty data mapping between data and CacheStruct
         * instance
         * \param d denotes the data to unmap
         */
        virtual void UnsetDirtyData(Data *d);

    private:
        // number of nimbus variables
        size_t num_variables_;

        // cache-data mappings
        typedef std::map<GeometricRegion,
                         Data *> DMap;
        std::vector<DMap> data_maps_;
        std::vector<DataSet> write_backs_;

        // disable constructor without arguments
        CacheStruct() {}

        // methods for cache manager to manage mappings and control access
        void SetUpReadWrite(const std::vector<cache::type_id_t> &var_type,
                            const std::vector<DataArray> &read_sets,
                            const std::vector<DataArray> &write_sets,
                            std::vector<DataArray> *flush_sets,
                            std::vector<DataArray> *diff_sets,
                            std::vector<DataArray> *sync_sets,
                            std::vector<CacheObjects> *sync_co_sets);
        void SetUpWrite(const std::vector<cache::type_id_t> &var_type,
                        const std::vector<DataArray> &write_sets,
                        GeometricRegion write_region);

        bool CheckPendingFlag(const std::vector<cache::type_id_t> &var_type,
                              const std::vector<DataArray> &read_sets,
                              const std::vector<DataArray> &write_sets);

        void ReleasePendingFlag(const std::vector<cache::type_id_t> &var_type,
                                std::vector<DataArray> *flush_sets,
                                std::vector<DataArray> *diff_sets,
                                std::vector<DataArray> *sync_sets,
                                std::vector<CacheObjects> *sync_co_sets);

        cache::distance_t GetDistance(const std::vector<cache::type_id_t> &var_type,
                                      const std::vector<DataArray> &read_sets) const;

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
                                    const GeometricRegion &write_region) const = 0;

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

    protected:
        /**
         * \brief Pulls data from cache, removes corresponding dirty data
         * mapping. Locks the cache object when pulling the data.
         * \param d is data to flush to
         */
        virtual void PullData(Data *d);

        // method for cache manager for profiling
        virtual size_t memory_size() {
          return sizeof(*this);
        }
};  // class CacheStruct

typedef std::vector<CacheStruct *> CacheStructs;

}  // namespace nimbus

#endif  // NIMBUS_DATA_CACHE_CACHE_STRUCT_H_

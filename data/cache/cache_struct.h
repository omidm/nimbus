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
 * A CacheStruct is an pplication object corresponding to a set of nimbus
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

#include "data/cache/utils.h"
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"
#include "worker/data.h"

namespace nimbus {

typedef int type_id_t;

/**
 * \class CacheStruct
 * \details Application object corresponding to a set of nimbus variables cached
 * by the cache manager. A CacheStruct can have num_variables number of
 * different nimbus variables, and must be defined over a geometric region.
 */
class CacheStruct {
    public:
        /**
         * \fn CacheStruct::CacheStruct(size_t num_variables)
         * \brief Creates a CacheStruct
         * \param num_variables indicates number of different application
         * variables that the CacheStruct instance contains
         * \return Constructed CacheStruct instance
         */
        explicit CacheStruct(size_t num_variables);

        /**
         * \fn CacheStruct::CacheStruct(size_t num_variables,
         *                              const GeometricRegion *struct_region)
         * \brief Constructs a CacheStruct
         * \param struct_region specifies the spatial domain of the CacheStruct
         * instance
         * \return Constructed CacheStruct instance
         */
        explicit CacheStruct(const GeometricRegion &struct_region);

        /**
         * \fn CacheStruct::CreateNew(const GeometricRegion &struct_region) const
         * \brief Creates a new CacheStruct instance using current instance
         * parameters
         * \param struct_region specifies the spatial domain of the CacheStruct
         * instance
         * \return Returns a pointer to the newly allocated CacheStruct instance
         */
        virtual CacheStruct *CreateNew(const GeometricRegion &struct_region) const;

        /**
         * \fn CacheStruct::UpdateCache(const std::vector<type_id_t> var_type,
         *                              const std::vector<DataArray *> read_sets,
         *                              const GeometricRegion &read_region)
         * \brief Updates CacheStruct with data from read_sets - performs update
         * only for those physical data that have changed since the CacheStruct
         * instance was created
         * \param var_type is a list of type_ids corresponding to nimbus variables
         * \param read_sets is a list of read_sets corresponding to nimbus variables
         * \param read_region is the geometric region to read
         */
        void UpdateCache(const std::vector<type_id_t> var_type,
                         const std::vector<DataArray *> read_sets,
                         const GeometricRegion &read_region);
        /**
         * \fn CacheStruct::CompleteWrite(const GeometricRegion &write_region)
         * \param write_region is the geometric region to write
         * \brief Indicates that write to CacheStruct instance is complete
         */
        void Write(const GeometricRegion &write_region);

        void FlushToData(Data *d);

        void AcquireAccess(CacheAccess access);
        void ReleaseAccess();
        bool IsAvailable(CacheAccess access) const;

        distance_t GetDistance(const DataArray &data_set) const;

        void SetUpRead(const std::vector<type_id_t> var_type,
                       const std::vector<DataArray *> read_sets);
        void SetUpWrite(const DataArray &write_set);

        void SetUpData(Data *d);
        void UnsetData(Data *d);

        GeometricRegion struct_region() const;

    private:
        void FlushCache(const DataArray &diff);

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

    protected:
        /**
         * \fn CacheStruct::ReadToCache(const std::vector<type_id_t> var_type,
         *                              const std::vector<DataArray *> read_sets,
         *                              const GeometricRegion &read_region)
         * \brief Reads data from read_sets into CacheStruct instance
         * \param var_type is a list of type_ids corresponding to nimbus variables
         * \param read_sets is a list of read_sets corresponding to nimbus variables
         * \param read_region is the geometric region to read
         */
        virtual void ReadToCache(const std::vector<type_id_t> var_type,
                                 const std::vector<DataArray *> read_sets,
                                 const GeometricRegion &read_region);

        /**
         * \fn CacheStruct::WriteFromCache(const std::vector<type_id_t> var_type,
         *                                 const std::vector<DataArray *> write_sets,
         *                                 const GeometricRegion &write_region)
         * \brief Writes data from CacheStruct instance to write_sets
         * \param var_type is a list of type_ids corresponding to nimbus variables
         * \param write_sets is a list of write_sets corresponding to nimbus variables
         * \param write_region is the geometric region to be write
         */
        virtual void WriteFromCache(const std::vector<type_id_t> var_type,
                                    const std::vector<DataArray *> read_sets,
                                    const GeometricRegion &write_region);
};  // class CacheStruct

typedef std::vector<CacheStruct *> CacheStructs;

}  // namespace nimbus

#endif  // NIMBUS_DATA_CACHE_CACHE_STRUCT_H_

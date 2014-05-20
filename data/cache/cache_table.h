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
 * A CacheTable contains a map from geometric region to application cache
 * objects.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#ifndef NIMBUS_DATA_CACHE_CACHE_TABLE_H_
#define NIMBUS_DATA_CACHE_CACHE_TABLE_H_

#include <map>
#include <vector>

#include "data/cache/cache_object.h"
#include "data/cache/cache_struct.h"
#include "data/cache/cache_var.h"
#include "data/cache/utils.h"
#include "shared/geometric_region.h"
#include "worker/data.h"

namespace nimbus {

/**
 * \class CacheTable
 * \details A CacheTable contains 2 maps - one from geometric regions to
 * CacheStructs, and another from geometric regions to CacheVars. However,
 * since a CacheTable really only corresponds to one CacheObject prototype,
 * only one of these will be populated and used. 
 */
class CacheTable {
    public:
        /**
         * \brief Creates a CacheTable
         * \param CacheTType specifies whether to create a table of
         * CacheStructs or CacheVars
         * \return Constructed CacheTable isntance
         */
        explicit CacheTable(CacheTType ttype);

        /**
         * \brief Adds cache variable to the map of region\-to\-variables
         * \param region specifies a key in the map, which is the region
         * corresponding to the CacheVar object
         * \param cv is the CacheVar value for given region key
         */
        void AddEntry(const GeometricRegion &region,
                      CacheVar *cv);

        /**
         * \brief Adds cache struct to the map of region\-to\-structs
         * \param region specifies a key in the map, which is the region
         * corresponding to the CacheStruct object
         * \param cs is the CacheStruct value for given region key
         */
        void AddEntry(const GeometricRegion &region,
                      CacheStruct *cs);

        /**
         * \brief Returns closest cache var instance from existing cache var
         * instances, if there is one available for the given region, otherwise
         * returns NULL
         * \param region specifies the geometric region of requested cache var
         * \param read_set specifies data that should be read into cache var
         * instance
         * \param access indicates the access mode for the request (SHARED or
         * EXCLUSIVE)
         * \return Returns a cache var instance that can be used, or NULL in
         * case no usable instance is available
         */
        CacheVar *GetClosestAvailable(const GeometricRegion &region,
                                      const DataArray &read_set,
                                      CacheAccess access = EXCLUSIVE);

        /**
         * \brief Returns closest cache struct instance from existing
         * instances, if there is one available for the given region, otherwise
         * returns NULL
         * \param region specifies the geometric region of requested cache
         * struct
         * \param indicates the data types corresponding to read sets in
         * read_sets
         * \param read_sets specifies data that should be read into cache
         * struct instance
         * \param access indicates the access mode for the request (SHARED or
         * EXCLUSIVE)
         * \return Returns a cache struct instance that can be used, or NULL in
         * case no usable instance is available
         */
        CacheStruct *GetClosestAvailable(const GeometricRegion &region,
                                         const std::vector<cache::type_id_t> &var_type,
                                         const std::vector<DataArray> &read_sets,
                                         CacheAccess access = EXCLUSIVE);

    private:
        /**
         * \brief Returns index corresponding to closest cache var instance,
         * from the given list of cache vars \- cvs
         * \param cvs is the list of cache vars to inspect
         * \param read_set specifies data that should be read into cache var
         * instance
         * \param access indicates the access mode for the request (SHARED or
         * EXCLUSIVE)
         * \return Returns a non-negative number (index) into cvs corresponding
         * to the cache var instance, that is closest to the given read set. If
         * no instances in cvs are usable, the method returns -1.
         */
        int GetMinDistanceIndex(const CacheVars &cvs,
                                const DataArray &read_set,
                                CacheAccess access = EXCLUSIVE) const;

        /**
         * \brief Returns index corresponding to closest cache struct instance,
         * from the given list of cache structs \- css
         * \param css is the list of cache structs to inspect
         * \param indicates the data types corresponding to read sets in
         * read_sets
         * \param read_sets specifies data that should be read into cache
         * struct nstance
         * \param access indicates the access mode for the request (SHARED or
         * EXCLUSIVE)
         * \return Returns a non-negative number (index) into css corresponding
         * to the cache struct instance, that is closest to list of given read
         * sets. If no instances in css are usable, the method returns -1.
         */
        int GetMinDistanceIndex(const CacheStructs &css,
                                const std::vector<cache::type_id_t> &var_type,
                                const std::vector<DataArray> &read_sets,
                                CacheAccess access = EXCLUSIVE) const;

        typedef std::map<GeometricRegion, CacheVars *> TVar;
        typedef std::map<GeometricRegion, CacheStructs *> TStruct;
        TVar *tvar_;
        TStruct *tstruct_;
};  // class CacheTable


}  // namespace nimbus

#endif  // NIMBUS_DATA_CACHE_CACHE_TABLE_H_

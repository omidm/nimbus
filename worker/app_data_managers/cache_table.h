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
 * A CacheTable contains a map from geometric region to cached application
 * objects.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#ifndef NIMBUS_WORKER_APP_DATA_MANAGERS_CACHE_TABLE_H_
#define NIMBUS_WORKER_APP_DATA_MANAGERS_CACHE_TABLE_H_

#include <map>
#include <sstream>
#include <vector>

#include "data/app_data/app_data_defs.h"
#include "data/app_data/app_object.h"
#include "data/app_data/app_struct.h"
#include "data/app_data/app_var.h"

namespace nimbus {

class Data;
typedef std::vector<Data *> DataArray;
class GeometricRegion;

/**
 * \class CacheTable
 * \details CacheTable contains the second level map managed by Nimbus -
 * from geometric regions to AppStructs or AppVars. The methods of this
 * class are visible to only CacheManager.
 */
class CacheTable {
    friend class AppDataManager;
    friend class CacheManager;

    private:
        /**
         * \brief Creates a CacheTable
         * \param app_data::AppDataType specifies whether to create a table of
         * AppStructs or AppVars
         * \return Constructed CacheTable isntance
         */
        explicit CacheTable(app_data::AppDataType ttype);

        /**
         * \brief Adds cache variable to the map of region-to-variables
         * \param region specifies a key in the map, which is the region
         * corresponding to the AppVar object
         * \param cv is the AppVar value for given region key
         */
        void AddEntry(const GeometricRegion &region,
                      AppVar *cv);

        /**
         * \brief Adds cache struct to the map of region-to-structs
         * \param region specifies a key in the map, which is the region
         * corresponding to the AppStruct object
         * \param cs is the AppStruct value for given region key
         */
        void AddEntry(const GeometricRegion &region,
                      AppStruct *cs);

        /**
         * \brief Returns closest cache var instance from existing cache var
         * instances, if there is one available for the given region, otherwise
         * returns NULL
         * \param region specifies the geometric region of requested cache var
         * \param read_set specifies data that should be read into cache var
         * instance
         * \param access indicates the access mode for the request (app_data::SHARED or
         * app_data::EXCLUSIVE)
         * \return Returns a cache var instance that can be used, or NULL in
         * case no usable instance is available
         */
        AppVar *GetClosestAvailable(const GeometricRegion &region,
                                      const DataArray &read_set,
                                      app_data::Access access = app_data::EXCLUSIVE);

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
         * \param access indicates the access mode for the request (app_data::SHARED or
         * app_data::EXCLUSIVE)
         * \return Returns a cache struct instance that can be used, or NULL in
         * case no usable instance is available
         */
        AppStruct *GetClosestAvailable(const GeometricRegion &region,
                                         const std::vector<app_data::type_id_t> &var_type,
                                         const std::vector<DataArray> &read_sets,
                                         app_data::Access access = app_data::EXCLUSIVE);

        /**
         * \brief Returns index corresponding to closest cache var instance,
         * from the given list of cache vars - cvs
         * \param cvs is the list of cache vars to inspect
         * \param read_set specifies data that should be read into cache var
         * instance
         * \param access indicates the access mode for the request (app_data::SHARED or
         * app_data::EXCLUSIVE)
         * \return Returns a non-negative number (index) into cvs corresponding
         * to the cache var instance, that is closest to the given read set. If
         * no instances in cvs are usable, the method returns -1.
         */
        int GetMinDistanceIndex(const AppVars &cvs,
                                const DataArray &read_set,
                                app_data::Access access = app_data::EXCLUSIVE) const;

        /**
         * \brief Returns index corresponding to closest cache struct instance,
         * from the given list of cache structs - css
         * \param css is the list of cache structs to inspect
         * \param indicates the data types corresponding to read sets in
         * read_sets
         * \param read_sets specifies data that should be read into cache
         * struct nstance
         * \param access indicates the access mode for the request (app_data::SHARED or
         * app_data::EXCLUSIVE)
         * \return Returns a non-negative number (index) into css corresponding
         * to the cache struct instance, that is closest to list of given read
         * sets. If no instances in css are usable, the method returns -1.
         */
        int GetMinDistanceIndex(const AppStructs &css,
                                const std::vector<app_data::type_id_t> &var_type,
                                const std::vector<DataArray> &read_sets,
                                app_data::Access access = app_data::EXCLUSIVE) const;

        void PrintProfile(std::stringstream* output);
        typedef std::map<GeometricRegion, AppVars *> TVar;
        typedef std::map<GeometricRegion, AppStructs *> TStruct;
        TVar *tvar_;
        TStruct *tstruct_;
};  // class CacheTable


}  // namespace nimbus

#endif  // NIMBUS_WORKER_APP_DATA_MANAGERS_CACHE_TABLE_H_

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
 * Application data manager is the interface for application and copy jobs to
 * application objects and nimbus data objects. The methods in this class are
 * mostly virtual. There are two implementations for this interface at present,
 * one that caches application objects, and one that does not (writes back data
 * from application object to nimbus objects immediately, and deletes the
 * application object).
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#ifndef NIMBUS_WORKER_APP_DATA_MANAGER_H_
#define NIMBUS_WORKER_APP_DATA_MANAGER_H_

#include <vector>
#include <string>

#include "data/app_data/app_data_defs.h"
#include "data/app_data/app_struct.h"
#include "data/app_data/app_var.h"

namespace nimbus {

class Data;
typedef std::vector<Data *> DataArray;
class GeometricRegion;

/**
 * \class AppDataManager
 * \details AppDataManager is the interface through which application and copy
 * jobs access data in application data objects and nimbus data objects.
 */
class AppDataManager {
    public:
        /**
         * \brief Creates an AppDataManager instance
         * \return Constructed AppDataManager instance
         */
        AppDataManager();

        /**
         * \brief Requests an AppVar instance of type prototype, from the
         * AppDataManager
         * \param read_set specifies the data that should be read in the
         * AppVar instance
         * \param read_region indicates read region
         * \param write_set specifies the data that should be marked as being
         * written
         * \param write_region indicates write region
         * \param prototype represents the application object type
         * \param region is the application object region
         * \access indicates whether application object access should be
         * EXCLUSIVE or SHARED, and is ignored when there is no caching
         * \return A pointer to AppVar instance that application can use
         */
        AppVar *GetAppVar(const DataArray &read_set,
                          const GeometricRegion &read_region,
                          const DataArray &write_set,
                          const GeometricRegion &write_region,
                          const AppVar &prototype,
                          const GeometricRegion &region,
                          appdata::Access access = appdata::EXCLUSIVE,
                          void (*aux)(AppVar*, void*) = NULL,
                          void* aux_data = NULL) = 0;

        /**
         * \brief Requests an AppStruct instance of type prototype, from the
         * AppDataManager
         * \param var_type specifies the application variable types for the
         * lists in read_sets/ write_sets
         * \param read_sets specifies the data that should be read in the
         * AppVar instance
         * \param read_region indicates read region
         * \param write_sets specifies the data that should be marked as being
         * written
         * \param write_region indicates write region
         * \param prototype represents the application object type
         * \param region is the application object region
         * \access indicates whether application object access should be
         * EXCLUSIVE or SHARED
         * \return A pointer to AppStruct instance that application can use
         */
        AppStruct *GetAppStruct(const std::vector<appdata::type_id_t> &var_type,
                                const std::vector<DataArray> &read_sets,
                                const GeometricRegion &read_region,
                                const std::vector<DataArray> &write_sets,
                                const GeometricRegion &write_region,
                                const AppStruct &prototype,
                                const GeometricRegion &region,
                                appdata::Access access) = 0;

        /**
         * \brief Informs application data manager that access to app_object is
         * no longer needed
         * \param app_object specified the application object to release
         */
        void ReleaseAccess(AppObject* app_object) = 0;

        /**
         * \brief Writes data to write_set back nimbus objects synchronously
         * \param app_var is the application variable to write from
         * \param write_set is the set of nimbus objects to write to
         */
        void WriteImmediately(AppVar *app_var, const DataArray &write_set) = 0;

        /**
         * \brief Writes data to write_set back nimbus objects synchronously
         * \param app_struct is the application struct to write from
         * \param var_type specifies the type of variables in the struct that
         * are to be written, corresponding to write_sets
         * \param write_sets is the set of nimbus objects to write to
         */
        void WriteImmediately(AppStruct *app_struct,
                              const std::vector<appdata::type_id_t> &var_type,
                              const std::vector<DataArray> &write_sets);

        /**
         * \brief If data is dirty, syncs with corresponding dirty application
         * object, and clears the dirty mapping, no-op in case there is no
         * caching
         * \param Data d to 
         */
        void SyncData(Data *d);

        /**
         * \brief Invalidates mapping between d and all application objects ,
         * dirty or not dirty, no-op in case there is no caching
         * \param Data d
         */
        void InvalidateMappings(Data *d);

        /**
         * \brief Sets log file names for application data manager
         * \param Worker id, to be used in file names
         */
        void SetLogNames(std::string wid_str);

    private:
        FILE* time_log;
        FILE* block_log;
        FILE* alloc_log;
};  // class AppDataManager
}  // namespace nimbus

#endif  // NIMBUS_WORKER_APP_DATA_MANAGER_H_

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
 * An AppVar is an application object corresponding to one nimbus
 * variable, provided by application managers.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#ifndef NIMBUS_DATA_APP_DATA_APP_VAR_H_
#define NIMBUS_DATA_APP_DATA_APP_VAR_H_

#include <map>
#include <set>
#include <string>
#include <vector>

#include "data/app_data/app_data_defs.h"
#include "data/app_data/app_object.h"
#include "shared/nimbus_types.h"

namespace nimbus {

class Data;
typedef std::vector<Data *> DataArray;
typedef std::set<Data *> DataSet;
class GeometricRegion;

/**
 * \class AppVar
 * \details Application object corresponding to one numbus variable provided
 * by app managers. An AppVar can must be defined over a geometric region.
 */
class AppVar : public AppObject {
    friend class CacheManager;
    friend class CacheTable;
    public:
        /**
         * \brief Creates an AppVar
         * \return Constructed AppVar instance
         */
        explicit AppVar();

        /**
         * \brief Creates an AppVar
         * \param  ob_reg specifies application object (AppVar) region
         * \return Constructed AppVar instance
         */
        explicit AppVar(const GeometricRegion &ob_reg);

    private:
        // app-data mappings
        typedef std::map<GeometricRegion,
                         Data *> DMap;
        DMap data_map_;
        DataSet write_back_;

        // methods for app manager to manage mappings and control access
        void SetUpReadWrite(const DataArray &read_set,
                            const DataArray &write_set,
                            DataArray *flush,
                            DataArray *diff,
                            DataArray *sync,
                            AppObjects *sync_co);

        bool CheckPendingFlag(const DataArray &read_set,
                              const DataArray &write_set);
        void ReleasePendingFlag(DataArray *flush,
                                DataArray *diff,
                                DataArray *sync);

        app_data::distance_t GetDistance(const DataArray &read_set) const;

    protected:
        /**
         * \brief Reads data from read_sets into AppVar instance
         * \param read_set is a data array corresponding to nimbus variables
         * \param read_region is the geometric region to read
         * \details This function must be overwritten by the application
         * writer. It provides the transformation from a set of nimbus data to
         * application instance.
         */
        virtual void ReadAppData(const DataArray &read_set,
                                 const GeometricRegion &read_region) = 0;

        /**
         * \brief Writes data from AppVar instance to write_sets
         * \param write_set is a data array corresponding to nimbus variables
         * \param write_region is the geometric region to write
         * \details This function must be overwritten by the application
         * writer. It provides the transformation from application instance to
         * nimbus data.
         */
        virtual void WriteAppData(const DataArray &write_set,
                                  const GeometricRegion &write_region) const = 0;

        /**
         * \brief Creates a new AppVar instance using current instance
         * parameters
         * \param struct_region specifies the spatial domain of the AppVar
         * instance
         * \return Returns a pointer to the newly allocated AppVar instance
         * \details This is a virtual function that must be over-written by application
         * writer. When AppManager cannot satisfy an application object request,
         * using the instances it already has, it calls CreateNew(..) on the
         * prototype passed in the request.
         */
        virtual AppVar *CreateNew(const GeometricRegion &struct_region) const = 0;

    protected:
        /**
         * \brief Pulls data from app, removes corresponding dirty data
         * mapping. Locks the app object when pulling the data.
         * \param d is data to flush to
         */
        virtual void PullData(Data *d);

        /**
         * \brief Unsets mapping between data and AppVar instance
         * \param d denotes the data to unmap
         */
        virtual void UnsetData(Data *d);

        /**
         * \brief Unsets dirty data mapping between data and AppVar instance
         * \param d denotes the data to unmap
         */
        virtual void UnsetDirtyData(Data *d);
};  // class AppVar

typedef std::vector<AppVar *> AppVars;

}  // namespace nimbus

#endif  // NIMBUS_DATA_APP_DATA_APP_VAR_H_

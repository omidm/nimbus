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
 * SimpleAppDataManager implements application data manager with
 * simple_app_data caching of application data across jobs.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#include <sys/syscall.h>
#include <cstdio>
#include <ctime>
#include <vector>
#include <string>

#include "data/app_data/app_data_defs.h"
#include "data/app_data/app_object.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "worker/app_data_manager.h"
#include "worker/app_data_managers/simple_app_data_manager.h"
#include "worker/data.h"

namespace nimbus {

/**
 * \details
 */
SimpleAppDataManager::SimpleAppDataManager() {}

/**
 * \details
 */
SimpleAppDataManager::~SimpleAppDataManager() {}

/**
 * \details This function does not do anything.
 */
void SimpleAppDataManager::WriteImmediately(AppVar *app_var,
                                    const DataArray &write_set) {
  // No Op
}

/**
 * \details This function does not do anything.
 */
void SimpleAppDataManager::WriteImmediately(AppStruct *app_struct,
                                    const std::vector<app_data::type_id_t> &var_type,
                                    const std::vector<DataArray> &write_sets) {
  // No Op
}

/**
 * \details SimpleAppDataManager read set all data in read within read_region,
 * and sets write back and write region.
 */
// TODO(hang/chinmayee): remove *aux, aux_data.
AppVar *SimpleAppDataManager::GetAppVarV(const DataArray &read_set,
                                 const GeometricRegion &read_region,
                                 const DataArray &write_set,
                                 const GeometricRegion &write_region,
                                 const AppVar &prototype,
                                 const GeometricRegion &region,
                                 app_data::Access access,
                                 void (*aux)(AppVar*, void*),
                                 void* aux_data) {
    AppVar *av = prototype.CreateNew(region);
    assert(av != NULL);
    av->write_back_ = DataSet(write_set.begin(), write_set.end());
    if (aux != NULL) {
      aux(av, aux_data);
    }
    av->write_region_ = write_region;
    av->ReadAppData(read_set, read_region);
    return av;
}

/**
 * \details SimpleAppDataManager reads all data in read sets within read_region,
 * and sets write back and write region.
 */
AppStruct *SimpleAppDataManager::GetAppStructV(const std::vector<app_data::type_id_t> &var_type,
                                       const std::vector<DataArray> &read_sets,
                                       const GeometricRegion &read_region,
                                       const std::vector<DataArray> &write_sets,
                                       const GeometricRegion &write_region,
                                       const AppStruct &prototype,
                                       const GeometricRegion &region,
                                       app_data::Access access) {
    // TODO(chinmayee): Remove this when application objects can be shared by compute jobs
    AppStruct *as = prototype.CreateNew(region);
    assert(as != NULL);
    size_t num_var = var_type.size();
    assert(write_sets.size() == num_var);
    for (size_t t = 0; t < num_var; ++t) {
      const DataArray &write_set = write_sets[t];
      as->write_backs_[var_type[t]] = DataSet(write_set.begin(), write_set.end());
    }
    as->write_region_ = write_region;
    as->ReadAppData(var_type, read_sets, read_region);
    return as;
}

/**
 * \details This function does not do anything.
 */
void SimpleAppDataManager::SyncData(Data *d) {
}

/**
 * \details This function does not do anything.
 */
void SimpleAppDataManager::InvalidateMappings(Data *d) {
}

/**
 * \details Writes back all data synchronously to write set, and deletes the
 * application object, freeing corresponding memory.
 */
void SimpleAppDataManager::ReleaseAccess(AppObject* app_object) {
    if (dynamic_cast<AppVar *>(app_object)) { // NOLINT
      AppVar *app_var = static_cast<AppVar *>(app_object);
      DataSet &write_back = app_var->write_back_;
      DataArray flush_set(write_back.begin(), write_back.end());
      write_back.clear();
      app_var->WriteAppData(flush_set, app_var->write_region_);
    } else {
      AppStruct *app_struct = static_cast<AppStruct *>(app_object);
      assert(app_struct != NULL);
      std::vector<app_data::type_id_t> var_type(app_struct->num_variables_);
      std::vector<DataSet> &write_backs = app_struct->write_backs_;
      std::vector<DataArray> flush_sets;
      for (size_t t = 0; t < app_struct->num_variables_; ++t) {
        var_type[t] = t;
        flush_sets.push_back(DataArray(write_backs[t].begin(),
                                       write_backs[t].end()));
        write_backs[t].clear();
      }
      app_struct->WriteAppData(var_type, flush_sets,
                              app_struct->write_region_);
    }
    app_object->Destroy();
    delete app_object;
}

}  // namespace nimbus

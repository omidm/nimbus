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
    DataArray flush_set;
    DataSet &write_back = app_var->write_back_;
    for (size_t i = 0; i < write_set.size(); ++i) {
        Data *d = write_set[i];
        if (write_back.find(d) != write_back.end()) {
            flush_set.push_back(d);
        }
    }
    for (size_t i = 0; i < flush_set.size(); ++i) {
        Data *d = flush_set[i];
        write_back.erase(d);
    }
    app_var->WriteAppData(flush_set, app_var->write_region_);
}

/**
 * \details This function does not do anything.
 */
void SimpleAppDataManager::WriteImmediately(AppStruct *app_struct,
                                    const std::vector<app_data::type_id_t> &var_type,
                                    const std::vector<DataArray> &write_sets) {
    size_t num_vars = var_type.size();
    if (write_sets.size() != num_vars) {
        dbg(DBG_ERROR, "Mismatch in number of variable types passed to FlushCache\n");
        exit(-1);
    }
    std::vector<DataArray> flush_sets(num_vars);
    std::vector<DataSet> &write_backs = app_struct->write_backs_;
    for (size_t t = 0; t < num_vars; ++t) {
        DataArray &flush_t = flush_sets[t];
        const DataArray &write_set_t = write_sets[t];
        app_data::type_id_t type = var_type[t];
        DataSet &write_back_t = write_backs[type];
        for (size_t i = 0; i < write_set_t.size(); ++i) {
            Data *d = write_set_t[i];
            if (write_back_t.find(d) != write_back_t.end()) {
                flush_t.push_back(d);
            }
        }
    }
    for (size_t t = 0; t < num_vars; ++t) {
        const DataArray &flush_t = flush_sets[t];
        app_data::type_id_t type = var_type[t];
        DataSet &write_back_t = write_backs[type];
        for (size_t i = 0; i < flush_t.size(); ++i) {
            Data *d = flush_t[i];
            write_back_t.erase(d);
        }
    }
    app_struct->WriteAppData(var_type, flush_sets,
                             app_struct->write_region_);
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
    GeometricRegion p1(29, 29, 29, 2, 2, 2);
    av->rc = NULL;
    for (size_t i = 0; i < read_set.size(); ++i) {
        Data *d = read_set[i];
        if (d->region().Intersects(&p1))
            av->rc = d;
    }

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
      if (app_var->rc) {
        Data *rc = app_var->rc;
        Data *wc = rc->Clone();
        wc->Copy(rc);
        float ws0 = wc->FloatingHash();
        app_var->WriteAppData(DataArray(1, wc), wc->region());
        float rs = rc->FloatingHash();
        float ws1 = wc->FloatingHash();
        if (rs - ws1 < -0.01 || rs - ws1 > 0.01) {
            printf("***== %s : Read variable sum %f, write variable sum before %f, after %f, check spurious write\n", wc->name().c_str(), rs, ws0, ws1);  // NOLINT
        } else {
            printf("###== %s : Read variable sum %f, write variable sum before %f, after %f, no spurious write\n", wc->name().c_str(), rs, ws0, ws1);  // NOLINT
        }
        app_var->rc = NULL;
        wc->Destroy();
        free(wc);
      } else {
        GeometricRegion p1(29, 29, 29, 2, 2, 2);
        for (size_t i = 0; i < flush_set.size(); ++i) {
            Data *d = flush_set[i];
            if (d->region().Intersects(&p1)) {
                printf("### non-spurious write variable sum %f\n", d->FloatingHash());
            }
        }
      }
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

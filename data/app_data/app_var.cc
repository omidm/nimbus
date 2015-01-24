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
 * variable, provided by app managers.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#include <map>
#include <set>
#include <string>
#include <vector>

#include "data/app_data/app_data_defs.h"
#include "data/app_data/app_var.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"
#include "worker/data.h"

namespace nimbus {

/**
 * \details
 */
AppVar::AppVar() {}

/**
 * \details
 */
AppVar::AppVar(const GeometricRegion &ob_reg) : AppObject(ob_reg) {
}

/**
 * \details UnsetData(...) removes data d from data_map_.
 */
void AppVar::UnsetData(Data *d) {
    GeometricRegion dreg = d->region();
    if (data_map_.find(dreg) != data_map_.end()) {
        assert(data_map_[dreg] == d);
        data_map_.erase(dreg);
    }
}

/**
 * \details UnsetDirtyData(...) removes d from write_back_.
 */
void AppVar::UnsetDirtyData(Data *d) {
    write_back_.erase(d);
}

/**
 * \details PullData(...) pulls data from app, after locking the struct.
 * When data needs to be updated from outside AppVar, use PullData.
 */
void AppVar::PullData(Data *d) {
    DataArray write_set(1, d);
    GeometricRegion dreg = d->region();
    GeometricRegion wreg = GeometricRegion::
        GetIntersection(write_region_, dreg);
    WriteAppData(write_set, wreg);
}

/**
 * \details GetDistance(...) gives the cost of using the AppVar instance,
 * given the read set. Current cost function is the sum of geometric sizes of
 * data in the read set.
 * The function does not have a separate argument for write set
 * - if you want write set included in the cost function, either add  another
 * argument or append it to read set that is passed to GetDistance.
 */
app_data::distance_t AppVar::GetDistance(const DataArray &read_set) const {
    app_data::distance_t cur_distance = 0;
    for (size_t i = 0; i < read_set.size(); ++i) {
        Data *d = read_set.at(i);
        GeometricRegion dreg = d->region();
        DMap::const_iterator it = data_map_.find(dreg);
        if (it != data_map_.end()) {
          if (it->second == d)
            continue;
        }
        cur_distance += dreg.dx() * dreg.dy() * dreg.dz();
    }
    return cur_distance;
}

/**
 * \details For read set, if data is not already in app object, insert it in
 * diff set to read, and sync it if necessary. If
 * it replaces existing data, flush existing data if dirty. For write set, if
 * data is not already in existing data, just create the mappings.
 * If it replaces existing data, flush existing data if dirty. Finally create
 * dirty object mapping with all data in write set.
 */
void AppVar::SetUpReadWrite(const DataArray &read_set,
                              const DataArray &write_set,
                              DataArray *flush,
                              DataArray *diff,
                              DataArray *sync,
                              AppObjects *sync_ob) {
    assert(flush != NULL);
    assert(diff != NULL);
    assert(sync != NULL);
    assert(sync_ob != NULL);
    for (size_t i = 0; i < read_set.size(); ++i) {
        Data *d = read_set.at(i);
        GeometricRegion dreg = d->region();
        DMap::iterator it = data_map_.find(dreg);
        if (it == data_map_.end()) {
            if (d->dirty_app_object()) {
                sync->push_back(d);
                sync_ob->push_back(d->dirty_app_object());
                d->ClearDirtyMappings();
                assert(d->app_ob_size() == 1);
            }
            diff->push_back(d);
            data_map_[dreg] = d;
            d->SetUpAppObject(this);
        } else {
            Data *d_old = it->second;
            if (d_old != d) {
                if (d->dirty_app_object()) {
                    sync->push_back(d);
                    sync_ob->push_back(d->dirty_app_object());
                    d->ClearDirtyMappings();
                    assert(d->app_ob_size() == 1);
                }
                if (write_back_.find(d_old) != write_back_.end()) {
                    flush->push_back(d_old);
                    write_back_.erase(d_old);
                    d_old->UnsetDirtyAppObject(this);
                    d_old->UnsetAppObject(this);
                    assert(d_old->app_ob_size() == 0);
                }
                d_old->UnsetAppObject(this);
                diff->push_back(d);
                data_map_[dreg] = d;
                d->SetUpAppObject(this);
            }
        }
    }
    for (size_t i = 0; i < write_set.size(); ++i) {
        Data *d = write_set.at(i);
        GeometricRegion dreg = d->region();
        DMap::iterator it = data_map_.find(dreg);
        if (it != data_map_.end()) {
            Data *d_old = it->second;
            if (d_old != d) {
                if (write_back_.find(d_old) != write_back_.end()) {
                    flush->push_back(d_old);
                    write_back_.erase(d_old);
                    d_old->UnsetDirtyAppObject(this);
                    d_old->UnsetAppObject(this);
                    assert(d_old->app_ob_size() == 0);
                }
                d_old->UnsetAppObject(this);
            }
        }
        d->InvalidateMappings();
        data_map_[dreg] = d;
        d->SetUpAppObject(this);
        write_back_.insert(d);
        d->SetUpDirtyAppObject(this);
    }
    for (size_t i = 0; i < flush->size(); ++i) {
      flush->at(i)->set_pending_flag(Data::WRITE);
    }
    for (size_t i = 0; i < sync->size(); ++i) {
      /* Data in the sync array may also be present in the flush array. */
      if (sync->at(i)->pending_flag() != -1) {
        assert(sync->at(i)->pending_flag() == 0);
        sync->at(i)->set_pending_flag(Data::WRITE);
      }
    }
    for (size_t i = 0; i < diff->size(); ++i) {
      /* Data in the diff array may have also be present in the sync and flush
         arrays. Only set read flag if the data was not marked in write mode. */
        // TODO(chinmayee): We probably don't need this. We should remove it.
        if (diff->at(i)->pending_flag() != -1) {
          diff->at(i)->set_pending_flag(Data::READ);
        }
    }
}

// method for app manager to manage mappings and control access - checks if
// any data in read/ write set has pending flag set
bool AppVar::CheckPendingFlag(const DataArray &read_set,
                                const DataArray &write_set) {
    for (size_t i = 0; i < read_set.size(); ++i) {
      Data *d = read_set.at(i);
      if (d->pending_flag() == -1) {
        // chinmayee: assert this. this branch should never happen.
        return false;
      }
      GeometricRegion dreg = d->region();
      DMap::iterator it = data_map_.find(dreg);
      if (it == data_map_.end()) {
        if (d->dirty_app_object()) {
          if (d->pending_flag() != 0) {
            return false;
          }
        }
      } else {
        Data *d_old = it->second;
        if (d_old != d) {
          if (d->dirty_app_object()) {
            if (d->pending_flag() != 0) {
              return false;
            }
          }
          if (write_back_.find(d_old) != write_back_.end()) {
            if (d_old->pending_flag() != 0) {
              return false;
            }
          }
        }
      }
    }
    for (size_t i = 0; i < write_set.size(); ++i) {
      Data *d = write_set.at(i);
      if (d->pending_flag() != 0) {
        // chinmayee: assert this. this branch should never happen.
        return false;
      }
      GeometricRegion dreg = d->region();
      DMap::iterator it = data_map_.find(dreg);
      if (it != data_map_.end()) {
        Data *d_old = it->second;
        if (d_old != d) {
          if (write_back_.find(d_old) != write_back_.end()) {
            // chinmayee: pending flag should never be 1 here.
            // assert that it is not 1, and wait if it is -1.
            if (d_old->pending_flag() != 0) {
              return false;
            }
          }
        }
      }
    }
    return true;
}

// method for app manager to manage mappings and control access - unset pending
// flag for all data in flush, diff and sync sets
void AppVar::ReleasePendingFlag(DataArray *flush,
                                  DataArray *diff,
                                  DataArray *sync) {
    assert(flush != NULL);
    assert(diff != NULL);
    assert(sync != NULL);
    for (size_t i = 0; i < flush->size(); ++i) {
      flush->at(i)->unset_pending_flag(Data::WRITE);
      assert(flush->at(i)->pending_flag() == 0);
    }
    for (size_t i = 0; i < sync->size(); ++i) {
      /* Check for case where data in sync array is also present in flush array. */
      assert(sync->at(i)->pending_flag() <= 0);
      if (sync->at(i)->pending_flag() == -1) {
        sync->at(i)->unset_pending_flag(Data::WRITE);
      }
    }
    for (size_t i = 0; i < diff->size(); ++i) {
      /* Check for case where data in diff array is also present in flush or sync array */
      if (diff->at(i)->pending_flag() > 0) {
        diff->at(i)->unset_pending_flag(Data::READ);
      }
    }
}

}  // namespace nimbus

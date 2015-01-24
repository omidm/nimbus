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

#include <map>
#include <set>
#include <string>
#include <vector>

#include "data/cache/cache_defs.h"
#include "data/cache/cache_struct.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"
#include "worker/data.h"

namespace nimbus {

/**
 * \details
 */
CacheStruct::CacheStruct(size_t num_variables) : num_variables_(num_variables),
                                                 data_maps_(num_variables),
                                                 write_backs_(num_variables) {
}

/**
 * \details
 */
CacheStruct::CacheStruct(size_t num_variables, const GeometricRegion &ob_reg)
    : CacheObject(ob_reg),
      num_variables_(num_variables),
      data_maps_(num_variables),
      write_backs_(num_variables) {
}

/**
 * \details UnsetData(...) removes data d from data_maps_. It checks every data
 * map in data_maps_ for the data d, since no prior type information is
 * available (just like PullData).
 */
void CacheStruct::UnsetData(Data *d) {
    GeometricRegion dreg = d->region();
    cache::type_id_t type = d->cache_type();
    assert(type != MAGIC_CACHE_TYPE);
    assert(type < num_variables_);
    DMap &data_map_t = data_maps_[type];
    bool success = false;
    if (data_map_t.find(dreg) != data_map_t.end()) {
        assert(data_map_t[dreg] == d);
        data_map_t.erase(dreg);
        success = true;
    }
    assert(success);
}

/**
 * \details UnsetDirtyData(...) removes d from write_backs_.
 */
void CacheStruct::UnsetDirtyData(Data *d) {
    cache::type_id_t type = d->cache_type();
    assert(type != MAGIC_CACHE_TYPE);
    assert(type < num_variables_);
    DataSet &write_back_t = write_backs_[type];
    assert(write_back_t.find(d) != write_back_t.end());
    write_back_t.erase(d);
}

/**
 * \details PullData(...) pulls data from cache, after locking the struct.
 * When data needs to be updated from outside CacheStruct, use PullData.
 * PullData(...)
 * checks each write_back set in the list write_backs_, and finds out which
 * type this data corresponds to. This may be a bad idea, but with twin copy
 * implementation, I expect the overhead to be small. -- Chinmayee
 */
void CacheStruct::PullData(Data *d) {
    cache::type_id_t type = d->cache_type();
    std::vector<cache::type_id_t> var_type(1, type);
    std::vector<DataArray> write_sets(1, DataArray(1, d));
    GeometricRegion dreg = d->region();
    GeometricRegion wreg = GeometricRegion::
        GetIntersection(write_region_, dreg);
    WriteFromCache(var_type, write_sets, wreg);
}

/**
 * \details GetDistance(...) gives the cost of using the CacheStruct instance,
 * given the read set. Current cost function is the sum of geometric sizes of
 * data in the read set.
 * The function does not have a separate argument for write set
 * - if you want write set included in the cost function, either add  another
 * argument or append it to read set that is passed to GetDistance.
 */
cache::distance_t CacheStruct::GetDistance(const std::vector<cache::type_id_t> &var_type,
                                           const std::vector<DataArray> &read_sets) const {
    cache::distance_t cur_distance = 0;
    for (size_t t = 0; t < num_variables_; ++t) {
        cache::type_id_t type = var_type[t];
        if (type > num_variables_) {
            dbg(DBG_WARN, "Invalid type %u passed to GetDistance, ignoring it\n", type);
            continue;
        }
        const DMap &data_map_t = data_maps_[type];
        for (size_t i = 0; i < read_sets[t].size(); ++i) {
            Data *d = read_sets[t].at(i);
            GeometricRegion dreg = d->region();
            DMap::const_iterator it = data_map_t.find(dreg);
            if (it != data_map_t.end()) {
              if (it->second == d)
                continue;
            }
            cur_distance += dreg.dx() * dreg.dy() * dreg.dz();
        }
    }
    return cur_distance;
}

/**
 * \details If data is not already in existing data, create the mappings.
 * If it replaces existing data, flush existing data if dirty. Create
 * dirty object mapping with all data in write set and set write region.
 */
void CacheStruct::SetUpWrite(const std::vector<cache::type_id_t> &var_type,
                             const std::vector<DataArray> &write_sets,
                             GeometricRegion write_region) {
    // TODO(quhang): synchronization not implemented.
    assert(false);
    size_t num_vars = var_type.size();
    assert(write_sets.size() == num_vars);
    assert(num_vars <= num_variables_);
    std::vector<DataArray> flush_sets;
    for (size_t t = 0; t < num_vars; ++t) {
        const DataArray &write_set_t = write_sets[t];
        DataArray &flush_t = flush_sets[t];
        cache::type_id_t type = var_type[t];
        DMap &data_map_t = data_maps_[type];
        DataSet &write_back_t = write_backs_[type];
        for (size_t i = 0; i < write_set_t.size(); ++i) {
            Data *d = write_set_t.at(i);
            GeometricRegion dreg = d->region();
            DMap::iterator it = data_map_t.find(dreg);
            if (it != data_map_t.end()) {
                Data *d_old = it->second;
                if (d_old != d) {
                    if (write_back_t.find(d_old) != write_back_t.end()) {
                        flush_t.push_back(d_old);
                        write_back_t.erase(d_old);
                        d_old->UnsetDirtyCacheObject(this);
                        d_old->UnsetCacheObject(this);
                        assert(d_old->co_size() == 0);
                    }
                    d_old->UnsetCacheObject(this);
                }
            }
            d->InvalidateMappings();
            data_map_t[dreg] = d;
            d->SetUpCacheObject(this);
            write_back_t.insert(d);
            d->SetUpDirtyCacheObject(this);
            d->set_cache_type(t);
        }
    }
    GeometricRegion write_region_old = write_region_;
    write_region_ = write_region;
    WriteFromCache(var_type, flush_sets, write_region_old);
}

/**
 * \details For read set, if data is not already in cache object, insert it in
 * diff set to read, and sync it if necessary. If
 * it replaces existing data, flush existing data if dirty. For write set, if
 * data is not already in existing data, just create the mappings.
 * If it replaces existing data, flush existing data if dirty. Finally create
 * dirty object mapping with all data in write set.
 */
void CacheStruct::SetUpReadWrite(const std::vector<cache::type_id_t> &var_type,
                                 const std::vector<DataArray> &read_sets,
                                 const std::vector<DataArray> &write_sets,
                                 std::vector<DataArray> *flush_sets,
                                 std::vector<DataArray> *diff_sets,
                                 std::vector<DataArray> *sync_sets,
                                 std::vector<CacheObjects> *sync_co_sets) {
    assert(flush_sets != NULL);
    assert(diff_sets != NULL);
    assert(sync_sets != NULL);
    assert(sync_co_sets != NULL);
    size_t num_vars = var_type.size();
    assert(read_sets.size() == num_vars &&
           write_sets.size() == num_vars &&
           flush_sets->size() == num_vars &&
           diff_sets->size() == num_vars &&
           sync_sets->size() == num_vars &&
           sync_co_sets->size() == num_vars);
    assert(num_vars <= num_variables_);
    for (size_t t = 0; t < num_vars; ++t) {
        const DataArray &read_set_t = read_sets[t];
        const DataArray &write_set_t = write_sets[t];
        DataArray &flush_t = flush_sets->at(t);
        DataArray &diff_t = diff_sets->at(t);
        DataArray &sync_t = sync_sets->at(t);
        CacheObjects &sync_co_t = sync_co_sets->at(t);
        cache::type_id_t type = var_type[t];
        DMap &data_map_t = data_maps_[type];
        DataSet &write_back_t = write_backs_[type];
        for (size_t i = 0; i < read_set_t.size(); ++i) {
            Data *d = read_set_t.at(i);
            GeometricRegion dreg = d->region();
            DMap::iterator it = data_map_t.find(dreg);
            if (it == data_map_t.end()) {
                if (d->dirty_cache_object()) {
                    sync_t.push_back(d);
                    sync_co_t.push_back(d->dirty_cache_object());
                    d->ClearDirtyMappings();
                    assert(d->co_size() == 1);
                }
                diff_t.push_back(d);
                data_map_t[dreg] = d;
                d->SetUpCacheObject(this);
                d->set_cache_type(type);
            } else {
                Data *d_old = it->second;
                if (d_old != d) {
                    if (d->dirty_cache_object()) {
                        sync_t.push_back(d);
                        sync_co_t.push_back(d->dirty_cache_object());
                        d->ClearDirtyMappings();
                        assert(d->co_size() == 1);
                    }
                    if (write_back_t.find(d_old) != write_back_t.end()) {
                        flush_t.push_back(d_old);
                        write_back_t.erase(d_old);
                        d_old->UnsetDirtyCacheObject(this);
                        d_old->UnsetCacheObject(this);
                        assert(d_old->co_size() == 0);
                    }
                    d_old->UnsetCacheObject(this);
                    diff_t.push_back(d);
                    data_map_t[dreg] = d;
                    d->SetUpCacheObject(this);
                    d->set_cache_type(type);
                }
            }
        }
        for (size_t i = 0; i < write_set_t.size(); ++i) {
            Data *d = write_set_t.at(i);
            GeometricRegion dreg = d->region();
            DMap::iterator it = data_map_t.find(dreg);
            if (it != data_map_t.end()) {
                Data *d_old = it->second;
                if (d_old != d) {
                    if (write_back_t.find(d_old) != write_back_t.end()) {
                        flush_t.push_back(d_old);
                        write_back_t.erase(d_old);
                        d_old->UnsetDirtyCacheObject(this);
                        d_old->UnsetCacheObject(this);
                        assert(d_old->co_size() == 0);
                    }
                    d_old->UnsetCacheObject(this);
                }
            }
            d->InvalidateMappings();
            data_map_t[dreg] = d;
            d->SetUpCacheObject(this);
            write_back_t.insert(d);
            d->SetUpDirtyCacheObject(this);
            d->set_cache_type(t);
        }
        for (size_t i = 0; i < flush_t.size(); ++i) {
          assert(flush_t.at(i)->pending_flag() == 0);
          flush_t.at(i)->set_pending_flag(Data::WRITE);
        }
        for (size_t i = 0; i < sync_t.size(); ++i) {
          /* Data in the sync array may also be present in the flush array. */
          if (sync_t.at(i)->pending_flag() != -1) {
            assert(sync_t.at(i)->pending_flag() == 0);
            sync_t.at(i)->set_pending_flag(Data::WRITE);
          }
        }
        for (size_t i = 0; i < diff_t.size(); ++i) {
          /* Data in the diff array may have also be present in the sync and flush
             arrays. Only set read flag if the data was not marked in write mode. */
          if (diff_t.at(i)->pending_flag() != -1) {
            diff_t.at(i)->set_pending_flag(Data::READ);
          }
        }
    }
}

bool CacheStruct::CheckPendingFlag(
        const std::vector<cache::type_id_t> &var_type,
        const std::vector<DataArray> &read_sets,
        const std::vector<DataArray> &write_sets) {
    size_t num_vars = var_type.size();
    assert(read_sets.size() == num_vars &&
           write_sets.size() == num_vars);
    assert(num_vars <= num_variables_);

    for (size_t t = 0; t < num_vars; ++t) {
      const DataArray &read_set_t = read_sets[t];
      const DataArray &write_set_t = write_sets[t];
      cache::type_id_t type = var_type[t];
      DMap &data_map_t = data_maps_[type];
      DataSet &write_back_t = write_backs_[type];
      for (size_t i = 0; i < read_set_t.size(); ++i) {
        Data *d = read_set_t.at(i);
        if (d->pending_flag() == -1) {
          return false;
        }
        GeometricRegion dreg = d->region();
        DMap::iterator it = data_map_t.find(dreg);
        if (it == data_map_t.end()) {
          if (d->dirty_cache_object()) {
            if (d->pending_flag() != 0) {
              return false;
            }
          }
        } else {
          Data *d_old = it->second;
          if (d_old != d) {
            if (d->dirty_cache_object()) {
              if (d->pending_flag() != 0) {
                return false;
              }
            }
            if (write_back_t.find(d_old) != write_back_t.end()) {
              if (d_old->pending_flag() != 0) {
                return false;
              }
            }
          }
        }
      }
      for (size_t i = 0; i < write_set_t.size(); ++i) {
        Data *d = write_set_t.at(i);
        if (d->pending_flag() != 0) {
          return false;
        }
        GeometricRegion dreg = d->region();
        DMap::iterator it = data_map_t.find(dreg);
        if (it != data_map_t.end()) {
          Data *d_old = it->second;
          if (d_old != d) {
            if (write_back_t.find(d_old) != write_back_t.end()) {
              if (d_old->pending_flag() != 0) {
                return false;
              }
            }
          }
        }
      }
    }
    return true;
}

void CacheStruct::ReleasePendingFlag(
        const std::vector<cache::type_id_t> &var_type,
        std::vector<DataArray> *flush_sets,
        std::vector<DataArray> *diff_sets,
        std::vector<DataArray> *sync_sets) {
    assert(flush_sets != NULL);
    assert(diff_sets != NULL);
    assert(sync_sets != NULL);
    size_t num_vars = var_type.size();
    assert(flush_sets->size() == num_vars &&
           diff_sets->size() == num_vars &&
           sync_sets->size() == num_vars);
    assert(num_vars <= num_variables_);
    for (size_t t = 0; t < num_vars; ++t) {
        DataArray &flush_t = flush_sets->at(t);
        DataArray &diff_t = diff_sets->at(t);
        DataArray &sync_t = sync_sets->at(t);
        for (size_t i = 0; i < flush_t.size(); ++i) {
          flush_t.at(i)->unset_pending_flag(Data::WRITE);
          assert(flush_t.at(i)->pending_flag() == 0);
        }
        for (size_t i = 0; i < sync_t.size(); ++i) {
          /* Check for case where data in sync array is also present in flush array. */
          assert(sync_t.at(i)->pending_flag() <= 0);
          if (sync_t.at(i)->pending_flag() == -1) {
            sync_t.at(i)->unset_pending_flag(Data::WRITE);
          }
          assert(sync_t.at(i)->pending_flag() == 0);
        }
        for (size_t i = 0; i < diff_t.size(); ++i) {
          /* Check for case where data in diff array is also present in flush or sync array */
          if (diff_t.at(i)->pending_flag() > 0) {
            diff_t.at(i)->unset_pending_flag(Data::READ);
          }
        }
    }
}

}  // namespace nimbus

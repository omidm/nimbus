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
 * CacheManager is the interface for application jobs to cache. Application
 * jobs send their request for application objects to the cache manager.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#include <sys/syscall.h>
#include <pthread.h>
#include <cstdio>
#include <ctime>
#include <map>
#include <vector>
#include <string>

#include "data/cache/cache_defs.h"
#include "data/cache/cache_object.h"
#include "data/cache/cache_table.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "worker/cache_manager.h"
#include "worker/data.h"

#define FIRST_UNIQUE_ID 1000
namespace nimbus {

bool CacheManager::print_stat_ = true;

/**
 * \details
 */
CacheManager::CacheManager() {
    unique_id_allocator_ = FIRST_UNIQUE_ID;
    memory_sum_ = 0;
    pool_ = new Pool();
    pthread_mutex_init(&cache_lock, NULL);
    pthread_cond_init(&cache_cond, NULL);
    alloc_log = fopen("cache_objects.txt", "w");
    block_log = fopen("cache_behavior.txt", "w");
    time_log = fopen("cache_time.txt", "w");
}

void CacheManager::WriteImmediately(CacheVar *cache_var,
                                    const DataArray &write_set) {
    PrintTimeStamp("start", "WIV stage");
    DataArray flush_set;
    DataSet &write_back = cache_var->write_back_;
    size_t write_bytes = 0;
    for (size_t i = 0; i < write_set.size(); ++i) {
        Data *d = write_set[i];
        if (write_back.find(d) != write_back.end()) {
            flush_set.push_back(d);
            write_bytes += d->memory_size();
        }
    }
    for (size_t i = 0; i < flush_set.size(); ++i) {
        Data *d = flush_set[i];
        d->UnsetDirtyCacheObject(cache_var);
        write_back.erase(d);
    }
    std::string size_str = "DSW wfcsize " + cache_var->name();
    PrintSizeStamp(size_str.c_str(), write_bytes);
    std::string wfc_str = "GAV wfc " + cache_var->name();
    PrintTimeStamp("start", wfc_str.c_str());
    cache_var->WriteFromCache(flush_set, cache_var->write_region_);
    PrintTimeStamp("end", wfc_str.c_str());
    PrintTimeStamp("end", "WIV stage");
}

/**
 * \detials WriteImmediately(...) checks if data passed to it is in the
 * write_back set for the cache struct instance.
 */
void CacheManager::WriteImmediately(CacheStruct *cache_struct,
                                    const std::vector<cache::type_id_t> &var_type,
                                    const std::vector<DataArray> &write_sets) {
    PrintTimeStamp("start", "WIS stage");
    size_t num_vars = var_type.size();
    if (write_sets.size() != num_vars) {
        dbg(DBG_ERROR, "Mismatch in number of variable types passed to FlushCache\n");
        exit(-1);
    }
    std::vector<DataArray> flush_sets(num_vars);
    std::vector<DataSet> &write_backs = cache_struct->write_backs_;
    size_t write_bytes = 0;
    for (size_t t = 0; t < num_vars; ++t) {
        DataArray &flush_t = flush_sets[t];
        const DataArray &write_set_t = write_sets[t];
        cache::type_id_t type = var_type[t];
        DataSet &write_back_t = write_backs[type];
        for (size_t i = 0; i < write_set_t.size(); ++i) {
            Data *d = write_set_t[i];
            if (write_back_t.find(d) != write_back_t.end()) {
                flush_t.push_back(d);
                write_bytes += d->memory_size();
            }
        }
    }
    for (size_t t = 0; t < num_vars; ++t) {
        const DataArray &flush_t = flush_sets[t];
        cache::type_id_t type = var_type[t];
        DataSet &write_back_t = write_backs[type];
        for (size_t i = 0; i < flush_t.size(); ++i) {
            Data *d = flush_t[i];
            d->UnsetDirtyCacheObject(cache_struct);
            write_back_t.erase(d);
        }
    }
    std::string size_str = "DSW wfcsize " + cache_struct->name();
    PrintSizeStamp(size_str.c_str(), write_bytes);
    std::string wfc_str = "GAV wfc " + cache_struct->name();
    PrintTimeStamp("start", wfc_str.c_str());
    cache_struct->WriteFromCache(var_type, flush_sets,
                                 cache_struct->write_region_);
    PrintTimeStamp("end", wfc_str.c_str());
    PrintTimeStamp("end", "WIS stage");
}

/**
 * \details CacheManager checks if an instance with requested application
 * object region and prototype id is present and available for use. If not, it
 * creates a new instance, and adds the instance to its 2-level map. It then
 * updates the instance as being used in access mode, updates the instance to
 * refelct any unread data, and sets up write set (flush data that may be
 * replaced, set dirty mappings etc.)
 */
CacheVar *CacheManager::GetAppVar(const DataArray &read_set,
                                  const GeometricRegion &read_region,
                                  const DataArray &write_set,
                                  const GeometricRegion &write_region,
                                  const CacheVar &prototype,
                                  const GeometricRegion &region,
                                  cache::CacheAccess access,
                                  void (*aux)(CacheVar*, void*),
                                  void* aux_data) {
    PrintTimeStamp("start", "GAV stage");
    PrintTimeStamp("start", "GAV lock");
    pthread_mutex_lock(&cache_lock);
    PrintTimeStamp("end", "GAV lock");
    CacheVar *cv = NULL;
    // Get a cache object form the cache table.
    if (pool_->find(prototype.id()) == pool_->end()) {
        CacheTable *ct = new CacheTable(cache::VAR);
        (*pool_)[prototype.id()] = ct;
        cv = prototype.CreateNew(region);
        cv->set_unique_id(unique_id_allocator_++);
        assert(cv != NULL);
        ct->AddEntry(region, cv);
    } else {
        CacheTable *ct = (*pool_)[prototype.id()];
        cv = ct->GetClosestAvailable(region, read_set, access);
        if (cv == NULL) {
            cv = prototype.CreateNew(region);
            cv->set_unique_id(unique_id_allocator_++);
            assert(cv != NULL);
            ct->AddEntry(region, cv);
        }
    }
    // cv->AcquireAccess(access);
    DataArray flush, sync, diff;
    CacheObjects sync_co;
    PrintTimeStamp("start", "GAV block");
    while (!cv->CheckPendingFlag(read_set, write_set)) {
      pthread_cond_wait(&cache_cond, &cache_lock);
    }
    PrintTimeStamp("end", "GAV block");
    // Move here.
    cv->AcquireAccess(access);
    cv->SetUpReadWrite(read_set, write_set,
                       &flush, &diff, &sync, &sync_co);
    pthread_mutex_unlock(&cache_lock);
    if (aux != NULL) {
      aux(cv, aux_data);
    }

    GeometricRegion write_region_old = cv->write_region_;
    cv->write_region_ = write_region;
    size_t write_bytes = 0;
    for (size_t i = 0; i < flush.size(); ++i) {
      Data *d = flush[i];
      write_bytes += d->memory_size();
    }
    for (size_t i = 0; i < sync.size(); ++i) {
      Data *d = sync[i];
      write_bytes += d->memory_size();
    }
    std::string size_str = "GAV wfcsize " + cv->name();
    PrintSizeStamp(size_str.c_str(), write_bytes);

    std::string wfc_str = "GAV wfc " + cv->name();
    PrintTimeStamp("start", wfc_str.c_str());
    cv->WriteFromCache(flush, write_region_old);
    for (size_t i = 0; i < sync.size(); ++i) {
        // assert(sync_co[i]->IsAvailable(cache::EXCLUSIVE));
        sync_co[i]->PullData(sync[i]);
    }
    PrintTimeStamp("end", wfc_str.c_str());
    size_t read_bytes = 0;
    for (size_t i = 0; i < diff.size(); ++i) {
      Data *d = diff[i];
      read_bytes += d->memory_size();
    }
    PrintSizeStamp("GAV rfcsize", read_bytes);

    std::string rfc_str = "GAV rfc " + cv->name();
    PrintTimeStamp("start", rfc_str.c_str());
    cv->ReadToCache(diff, read_region);
    PrintTimeStamp("end", rfc_str.c_str());
    PrintTimeStamp("start", "GAV lock");
    pthread_mutex_lock(&cache_lock);
    PrintTimeStamp("end", "GAV lock");
    cv->ReleasePendingFlag(&flush, &diff, &sync, &sync_co);
    cv->unset_pending_flag();
    pthread_cond_broadcast(&cache_cond);
    pthread_mutex_unlock(&cache_lock);
    PrintTimeStamp("end", "GAV stage");
    return cv;
}

/**
 * \details CacheManager checks if an instance with requested application
 * object region and prototype id is present and available for use. If not, it
 * creates a new instance, and adds the instance to its 2-level map. It then
 * updates the instance as being used in access mode, updates the instance to
 * refelct any unread data, and sets up write set (flush data that may be
 * replaced, set dirty mappings etc.)
 */
CacheStruct *CacheManager::GetAppStruct(const std::vector<cache::type_id_t> &var_type,
                                        const std::vector<DataArray> &read_sets,
                                        const GeometricRegion &read_region,
                                        const std::vector<DataArray> &write_sets,
                                        const GeometricRegion &write_region,
                                        const CacheStruct &prototype,
                                        const GeometricRegion &region,
                                        cache::CacheAccess access) {
    PrintTimeStamp("start", "GAS stage");
    PrintTimeStamp("start", "GAS lock");
    pthread_mutex_lock(&cache_lock);
    PrintTimeStamp("end", "GAS lock");
    CacheStruct *cs = NULL;
    if (pool_->find(prototype.id()) == pool_->end()) {
        CacheTable *ct = new CacheTable(cache::STRUCT);
        (*pool_)[prototype.id()] = ct;
        cs = prototype.CreateNew(region);
        cs->set_unique_id(unique_id_allocator_++);
        assert(cs != NULL);
        ct->AddEntry(region, cs);
    } else {
        CacheTable *ct = (*pool_)[prototype.id()];
        cs = ct->GetClosestAvailable(region, var_type, read_sets, access);
        if (cs == NULL) {
            cs = prototype.CreateNew(region);
            cs->set_unique_id(unique_id_allocator_++);
            assert(cs != NULL);
            ct->AddEntry(region, cs);
        }
    }
    // cs->AcquireAccess(access);
    size_t num_var = var_type.size();
    std::vector<DataArray> flush_sets(num_var),
                           sync_sets(num_var),
                           diff_sets(num_var);
    std::vector<CacheObjects> sync_co_sets(num_var);
    PrintTimeStamp("start", "GAS block");
    // Move here.
    while (!cs->CheckPendingFlag(var_type, read_sets, write_sets)) {
      pthread_cond_wait(&cache_cond, &cache_lock);
    }
    PrintTimeStamp("end", "GAS block");
    cs->AcquireAccess(access);
    cs->SetUpReadWrite(var_type, read_sets, write_sets,
                       &flush_sets, &diff_sets, &sync_sets, &sync_co_sets);
    pthread_mutex_unlock(&cache_lock);

    GeometricRegion write_region_old = cs->write_region_;
    cs->write_region_ = write_region;

    size_t write_bytes = 0;
    for (size_t i = 0; i < num_var; ++i) {
      DataArray &flush_t = flush_sets[i];
      for (size_t j = 0; j < flush_t.size(); ++j) {
        Data *d = flush_t[j];
        write_bytes += d->memory_size();
      }
    }
    for (size_t i = 0; i < num_var; ++i) {
      DataArray &sync_t = sync_sets[i];
      for (size_t j = 0; j < sync_t.size(); ++j) {
        Data *d = sync_t[j];
        write_bytes += d->memory_size();
      }
    }
    std::string size_str = "GAS wfcsize " + cs->name();
    PrintSizeStamp(size_str.c_str(), write_bytes);

    std::string wfc_str = "GAV wfc " + cs->name();
    PrintTimeStamp("start", wfc_str.c_str());
    cs->WriteFromCache(var_type, flush_sets, write_region_old);
    for (size_t t = 0; t < num_var; ++t) {
        DataArray &sync_t = sync_sets[t];
        CacheObjects &sync_co_t = sync_co_sets[t];
        for (size_t i = 0; i < sync_t.size(); ++i) {
            // assert(sync_co_t[i]->IsAvailable(cache::EXCLUSIVE));
            sync_co_t[i]->PullData(sync_t[i]);
        }
    }
    PrintTimeStamp("end", wfc_str.c_str());

    size_t read_bytes = 0;
    for (size_t i = 0; i < num_var; ++i) {
      DataArray &diff_t = diff_sets[i];
      for (size_t j = 0; j < diff_t.size(); ++j) {
        Data *d = diff_t[j];
        read_bytes += d->memory_size();
      }
    }
    PrintSizeStamp("GAS rfcsize", read_bytes);

    std::string rfc_str = "GAV rfc " + cs->name();
    PrintTimeStamp("start", rfc_str.c_str());
    cs->ReadToCache(var_type, diff_sets, read_region);
    PrintTimeStamp("end", rfc_str.c_str());
    PrintTimeStamp("start", "GAS lock");
    pthread_mutex_lock(&cache_lock);
    PrintTimeStamp("end", "GAS lock");
    cs->ReleasePendingFlag(var_type,
                           &flush_sets, &diff_sets, &sync_sets, &sync_co_sets);
    cs->unset_pending_flag();
    pthread_cond_broadcast(&cache_cond);
    pthread_mutex_unlock(&cache_lock);
    PrintTimeStamp("end", "GAS stage");
    return cs;
}

/**
 * \details
 */
void CacheManager::SyncData(Data *d) {
    PrintTimeStamp("start", "SD stage");
    CacheObject *co = NULL;
    PrintTimeStamp("start", "SD lock");
    pthread_mutex_lock(&cache_lock);
    PrintTimeStamp("end", "SD lock");
    PrintTimeStamp("start", "SD block");
    while (d->pending_flag() != 0 ||
           (d->dirty_cache_object()
            && d->dirty_cache_object()->pending_flag())) {
       pthread_cond_wait(&cache_cond, &cache_lock);
    }
    PrintTimeStamp("end", "SD block");
    co = d->dirty_cache_object();
    if (!co) {
        pthread_mutex_unlock(&cache_lock);
        PrintTimeStamp("end", "SD stage");
        return;
    }
    // assert(co->IsAvailable(cache::EXCLUSIVE));
    d->set_pending_flag(Data::WRITE);
    // co->set_pending_flag();
    // d->ClearDirtyMappings();
    pthread_mutex_unlock(&cache_lock);

    size_t write_bytes = d->memory_size();
    std::string size_str = "SD wfcsize " + co->name();
    PrintSizeStamp(size_str.c_str(), write_bytes);

    std::string wfc_str = "SD pdata " + co->name();
    PrintTimeStamp("start", wfc_str.c_str());
    co->PullData(d);
    PrintTimeStamp("end", wfc_str.c_str());
    PrintTimeStamp("start", "SD lock");
    pthread_mutex_lock(&cache_lock);
    PrintTimeStamp("end", "SD lock");
    d->ClearDirtyMappings();
    d->unset_pending_flag(Data::WRITE);
    // co->unset_pending_flag();
    pthread_cond_broadcast(&cache_cond);
    pthread_mutex_unlock(&cache_lock);
    PrintTimeStamp("end", "SD stage");
}

/**
 * \details
 */
void CacheManager::InvalidateMappings(Data *d) {
    PrintTimeStamp("start", "IM stage");
    PrintTimeStamp("start", "IM lock");
    pthread_mutex_lock(&cache_lock);
    PrintTimeStamp("end", "IM lock");
    PrintTimeStamp("start", "IM block");
    while (d->pending_flag() != 0) {
        pthread_cond_wait(&cache_cond, &cache_lock);
    }
    PrintTimeStamp("end", "IM block");
    d->InvalidateMappings();
    pthread_mutex_unlock(&cache_lock);
    PrintTimeStamp("end", "IM stage");
}

void CacheManager::ReleaseAccess(CacheObject* cache_object) {
    pthread_mutex_lock(&cache_lock);
    // cache_object->unset_pending_flag();
    // TODO(quhang): use private method and mark friend class.
    cache_object->ReleaseAccessInternal();
    pthread_cond_broadcast(&cache_cond);

    /*
    if (print_stat_) {
      uint64_t data_id =cache_object->unique_id();
      size_t new_size = cache_object->memory_size();
      size_t old_size = 0;
      if (memory_size_map_.find(data_id) == memory_size_map_.end()) {
        memory_size_map_[data_id] = new_size;
      } else {
        old_size = memory_size_map_[data_id];
        memory_size_map_[data_id] = new_size;
      }
      if (new_size != old_size) {
        memory_sum_ = memory_sum_ + new_size - old_size;
        struct timespec t;
        clock_gettime(CLOCK_REALTIME, &t);
        double time_sum = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
        fprintf(alloc_log, "%f %"PRIu64" %s %zu\n",
                time_sum,
                data_id,
                cache_object->name().c_str(),
                new_size);
        fprintf(alloc_log, "%f %zu\n", time_sum, memory_sum_);
        // fflush(alloc_log);
      }
    }
    */

    pthread_mutex_unlock(&cache_lock);
}

/*
void CacheManager::PrintProfile(std::stringstream* output) {
  if (pool_ == NULL) {
    return;
  }
  for (Pool::iterator iter = pool_->begin();
       iter != pool_->end();
       ++iter) {
    *output << "-----------" << std::endl;
    iter->second->PrintProfile(output);
  }
  *output << "-----------" << std::endl;
}
*/

void CacheManager::PrintTimeStamp(const char *status, const char *message) {
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  double time_sum = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  pid_t tid = syscall(SYS_gettid);
  fprintf(time_log, "%d ; %s ; %s ; %f\n", tid, status, message, time_sum);
}

/*
void CacheManager::BlockPrintTimeStamp(const char* message) {
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  double time_sum = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  fprintf(block_log, "%f : %s\n", time_sum, message);
}
*/

void CacheManager::PrintSizeStamp(const char *message, size_t num_bytes) {
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  double time_sum = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  pid_t tid = syscall(SYS_gettid);
  fprintf(time_log, "%d ; %s ; %zu ; %f\n", tid, message, num_bytes, time_sum);
}

}  // namespace nimbus

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
    block_log = fopen("cache_behavior.txt", "w");
    time_log = fopen("cache_time.txt", "w");
}

void CacheManager::DoSetUpWrite(CacheVar* cache_var,
                                const DataArray &write_set,
                                GeometricRegion &write_region) {
    DataArray flush;
    pthread_mutex_lock(&cache_lock);
    BlockPrintTimeStamp("enter");
    while (!cache_var->CheckWritePendingFlag(write_set, write_region)) {
        pthread_cond_wait(&cache_cond, &cache_lock);
    }
    BlockPrintTimeStamp("leave");
    cache_var->SetUpWrite(write_set, write_region, &flush);
    pthread_mutex_unlock(&cache_lock);
    cache_var->PerformSetUpWrite(write_set, write_region, flush);
    pthread_mutex_lock(&cache_lock);
    cache_var->ReleaseWritePendingFlag(write_set, flush);
    pthread_cond_broadcast(&cache_cond);
    pthread_mutex_unlock(&cache_lock);
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
    PrintTimeStamp("start", "GAV");
    PrintTimeStamp("start", "GAV l1");
    pthread_mutex_lock(&cache_lock);
    PrintTimeStamp("end", "GAV l1");
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
    BlockPrintTimeStamp("enter");
    PrintTimeStamp("start", "GAV check");
    while (!cv->CheckPendingFlag(read_set, write_set)) {
      pthread_cond_wait(&cache_cond, &cache_lock);
    }
    PrintTimeStamp("end", "GAV check");
    BlockPrintTimeStamp("leave");
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
    PrintTimeStamp("start", "GAV wfc");
    cv->WriteFromCache(flush, write_region_old);
    for (size_t i = 0; i < sync.size(); ++i) {
        // assert(sync_co[i]->IsAvailable(cache::EXCLUSIVE));
        sync_co[i]->PullData(sync[i]);
    }
    PrintTimeStamp("end", "GAV wfc");
    PrintTimeStamp("start", "GAV rfc");
    cv->ReadToCache(diff, read_region);
    PrintTimeStamp("end", "GAV rfc");
    PrintTimeStamp("start", "GAV l2");
    pthread_mutex_lock(&cache_lock);
    PrintTimeStamp("end", "GAV l2");
    cv->ReleasePendingFlag(&flush, &diff, &sync, &sync_co);
    pthread_cond_broadcast(&cache_cond);
    pthread_mutex_unlock(&cache_lock);
    PrintTimeStamp("end", "GAV");
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
    PrintTimeStamp("start", "GAS");
    PrintTimeStamp("start", "GAS l1");
    pthread_mutex_lock(&cache_lock);
    PrintTimeStamp("end", "GAS l1");
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
    BlockPrintTimeStamp("enter");
    PrintTimeStamp("start", "GAS check");
    // Move here.
    while (!cs->CheckPendingFlag(var_type, read_sets, write_sets)) {
      pthread_cond_wait(&cache_cond, &cache_lock);
    }
    BlockPrintTimeStamp("leave");
    PrintTimeStamp("end", "GAS check");
    cs->AcquireAccess(access);
    cs->SetUpReadWrite(var_type, read_sets, write_sets,
                       &flush_sets, &diff_sets, &sync_sets, &sync_co_sets);
    pthread_mutex_unlock(&cache_lock);

    GeometricRegion write_region_old = cs->write_region_;
    cs->write_region_ = write_region;
    PrintTimeStamp("start", "GAS wfc");
    cs->WriteFromCache(var_type, flush_sets, write_region_old);
    for (size_t t = 0; t < num_var; ++t) {
        DataArray &sync_t = sync_sets[t];
        CacheObjects &sync_co_t = sync_co_sets[t];
        for (size_t i = 0; i < sync_t.size(); ++i) {
            // assert(sync_co_t[i]->IsAvailable(cache::EXCLUSIVE));
            sync_co_t[i]->PullData(sync_t[i]);
        }
    }
    PrintTimeStamp("end", "GAS wfc");
    PrintTimeStamp("start", "GAS rfc");
    cs->ReadToCache(var_type, diff_sets, read_region);
    PrintTimeStamp("end", "GAS rfc");
    PrintTimeStamp("start", "GAS l2");
    pthread_mutex_lock(&cache_lock);
    PrintTimeStamp("end", "GAS l2");
    cs->ReleasePendingFlag(var_type,
                           &flush_sets, &diff_sets, &sync_sets, &sync_co_sets);
    pthread_cond_broadcast(&cache_cond);
    pthread_mutex_unlock(&cache_lock);
    PrintTimeStamp("end", "GAS");
    return cs;
}

/**
 * \details
 */
void CacheManager::SyncData(Data *d) {
    PrintTimeStamp("start", "SD");
    CacheObject *co = NULL;
    PrintTimeStamp("start", "SD l1");
    pthread_mutex_lock(&cache_lock);
    PrintTimeStamp("end", "SD l1");
    BlockPrintTimeStamp("enter");
    PrintTimeStamp("start", "SD check");
    while (d->pending_flag() != 0 ||
           (d->dirty_cache_object()
            && d->dirty_cache_object()->pending_flag())) {
       pthread_cond_wait(&cache_cond, &cache_lock);
    }
    PrintTimeStamp("end", "SD check");
    BlockPrintTimeStamp("leave");
    co = d->dirty_cache_object();
    if (!co) {
        pthread_mutex_unlock(&cache_lock);
        PrintTimeStamp("end", "SD");
        return;
    }
    // assert(co->IsAvailable(cache::EXCLUSIVE));
    d->set_pending_flag(Data::WRITE);
    co->set_pending_flag();
    d->ClearDirtyMappings();
    pthread_mutex_unlock(&cache_lock);

    co->PullData(d);
    PrintTimeStamp("start", "SD l2");
    pthread_mutex_lock(&cache_lock);
    PrintTimeStamp("end", "SD l2");
    d->unset_pending_flag(Data::WRITE);
    co->unset_pending_flag();
    pthread_cond_broadcast(&cache_cond);
    pthread_mutex_unlock(&cache_lock);
    PrintTimeStamp("end", "SD");
}

/**
 * \details
 */
void CacheManager::InvalidateMappings(Data *d) {
    PrintTimeStamp("start", "IM");
    PrintTimeStamp("start", "IM l1");
    pthread_mutex_lock(&cache_lock);
    PrintTimeStamp("end", "IM l1");
    BlockPrintTimeStamp("enter");
    PrintTimeStamp("start", "IM check");
    while (d->pending_flag() != 0) {
        pthread_cond_wait(&cache_cond, &cache_lock);
    }
    PrintTimeStamp("end", "IM check");
    BlockPrintTimeStamp("leave");
    d->InvalidateMappings();
    pthread_mutex_unlock(&cache_lock);
    PrintTimeStamp("end", "IM");
}

void CacheManager::ReleaseAccess(CacheObject* cache_object) {
    pthread_mutex_lock(&cache_lock);
    cache_object->unset_pending_flag();
    // TODO(quhang): use private method and mark friend class.
    cache_object->ReleaseAccessInternal();
    pthread_cond_broadcast(&cache_cond);

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
        fflush(alloc_log);
      }
    }
    pthread_mutex_unlock(&cache_lock);
}

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

void CacheManager::PrintTimeStamp(const char *status, const char *message) {
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  double time_sum = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  pid_t tid = syscall(SYS_gettid);
  fprintf(time_log, "%d ; %s ; %s ; %f\n", tid, status, message, time_sum);
}

void CacheManager::BlockPrintTimeStamp(const char* message) {
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  double time_sum = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  fprintf(block_log, "%f : %s\n", time_sum, message);
}

}  // namespace nimbus

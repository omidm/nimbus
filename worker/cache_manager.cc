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

namespace nimbus {

/**
 * \details
 */
CacheManager::CacheManager() {
    pool_ = new Pool();
    pthread_mutex_init(&cache_lock, NULL);
    pthread_cond_init(&cache_cond, NULL);
    alloc_log = fopen("cache_objects.txt", "w");
    block_log = fopen("cache_behavior.txt", "w");
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
    pthread_mutex_lock(&cache_lock);
    CacheVar *cv = NULL;
    // Get a cache object form the cache table.
    if (pool_->find(prototype.id()) == pool_->end()) {
        CacheTable *ct = new CacheTable(cache::VAR);
        (*pool_)[prototype.id()] = ct;
        cv = prototype.CreateNew(region);
        PrintTimeStamp(sizeof(*cv));
        assert(cv != NULL);
        ct->AddEntry(region, cv);
    } else {
        CacheTable *ct = (*pool_)[prototype.id()];
        cv = ct->GetClosestAvailable(region, read_set, access);
        if (cv == NULL) {
            cv = prototype.CreateNew(region);
            PrintTimeStamp(sizeof(*cv));
            assert(cv != NULL);
            ct->AddEntry(region, cv);
        }
    }
    // cv->AcquireAccess(access);
    DataArray flush, sync, diff;
    CacheObjects sync_co;
    BlockPrintTimeStamp("enter");
    while (!cv->CheckPendingFlag(read_set, write_set)) {
      pthread_cond_wait(&cache_cond, &cache_lock);
    }
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
    cv->WriteFromCache(flush, write_region_old);
    for (size_t i = 0; i < sync.size(); ++i) {
        assert(sync_co[i]->IsAvailable(cache::EXCLUSIVE));
        sync_co[i]->PullData(sync[i]);
    }
    cv->ReadToCache(diff, read_region);
    pthread_mutex_lock(&cache_lock);
    cv->ReleasePendingFlag(&flush, &diff, &sync, &sync_co);
    pthread_cond_broadcast(&cache_cond);
    pthread_mutex_unlock(&cache_lock);
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
    pthread_mutex_lock(&cache_lock);
    CacheStruct *cs = NULL;
    if (pool_->find(prototype.id()) == pool_->end()) {
        CacheTable *ct = new CacheTable(cache::STRUCT);
        (*pool_)[prototype.id()] = ct;
        cs = prototype.CreateNew(region);
        PrintTimeStamp(sizeof(*cs));
        assert(cs != NULL);
        ct->AddEntry(region, cs);
    } else {
        CacheTable *ct = (*pool_)[prototype.id()];
        cs = ct->GetClosestAvailable(region, var_type, read_sets, access);
        if (cs == NULL) {
            cs = prototype.CreateNew(region);
            PrintTimeStamp(sizeof(*cs));
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
    // Move here.
    while (!cs->CheckPendingFlag(var_type, read_sets, write_sets)) {
      pthread_cond_wait(&cache_cond, &cache_lock);
    }
    BlockPrintTimeStamp("leave");
    cs->AcquireAccess(access);
    cs->SetUpReadWrite(var_type, read_sets, write_sets,
                       &flush_sets, &diff_sets, &sync_sets, &sync_co_sets);
    pthread_mutex_unlock(&cache_lock);

    GeometricRegion write_region_old = cs->write_region_;
    cs->write_region_ = write_region;
    cs->WriteFromCache(var_type, flush_sets, write_region_old);
    for (size_t t = 0; t < num_var; ++t) {
        DataArray &sync_t = sync_sets[t];
        CacheObjects &sync_co_t = sync_co_sets[t];
        for (size_t i = 0; i < sync_t.size(); ++i) {
            assert(sync_co_t[i]->IsAvailable(cache::EXCLUSIVE));
            sync_co_t[i]->PullData(sync_t[i]);
        }
    }
    cs->ReadToCache(var_type, diff_sets, read_region);
    pthread_mutex_lock(&cache_lock);
    cs->ReleasePendingFlag(var_type,
                           &flush_sets, &diff_sets, &sync_sets, &sync_co_sets);
    pthread_cond_broadcast(&cache_cond);
    pthread_mutex_unlock(&cache_lock);
    return cs;
}

/**
 * \details
 */
void CacheManager::SyncData(Data *d) {
    CacheObject *co = NULL;
    pthread_mutex_lock(&cache_lock);
    BlockPrintTimeStamp("enter");
    while (d->pending_flag() ||
           (d->dirty_cache_object()
            && d->dirty_cache_object()->pending_flag())) {
       pthread_cond_wait(&cache_cond, &cache_lock);
    }
    BlockPrintTimeStamp("leave");
    co = d->dirty_cache_object();
    if (!co) {
        pthread_mutex_unlock(&cache_lock);
        return;
    }
    assert(co->IsAvailable(cache::EXCLUSIVE));
    d->set_pending_flag();
    co->set_pending_flag();
    d->ClearDirtyMappings();
    pthread_mutex_unlock(&cache_lock);

    co->PullData(d);
    pthread_mutex_lock(&cache_lock);
    d->unset_pending_flag();
    co->unset_pending_flag();
    pthread_cond_broadcast(&cache_cond);
    pthread_mutex_unlock(&cache_lock);
}

/**
 * \details
 */
void CacheManager::InvalidateMappings(Data *d) {
    pthread_mutex_lock(&cache_lock);
    BlockPrintTimeStamp("enter");
    while (d->pending_flag()) {
        pthread_cond_wait(&cache_cond, &cache_lock);
    }
    BlockPrintTimeStamp("leave");
    d->InvalidateMappings();
    pthread_mutex_unlock(&cache_lock);
}

void CacheManager::ReleaseAccess(CacheObject* cache_object) {
    pthread_mutex_lock(&cache_lock);
    cache_object->unset_pending_flag();
    // TODO(quhang): use private method and mark friend class.
    cache_object->ReleaseAccessInternal();
    pthread_cond_broadcast(&cache_cond);
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

void CacheManager::PrintTimeStamp(uint64_t size) {
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  double time_sum = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  fprintf(alloc_log, "%f : %lu\n", time_sum, size);
}

void CacheManager::BlockPrintTimeStamp(const char* message) {
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  double time_sum = t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
  fprintf(block_log, "%f : %s\n", time_sum, message);
}

}  // namespace nimbus

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
 * CacheManager implements application data manager with simple_app_data caching of
 * application data across jobs.
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

#include "data/app_data/app_data_defs.h"
#include "data/app_data/app_object.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "worker/app_data_manager.h"
#include "worker/app_data_managers/cache_manager.h"
#include "worker/app_data_managers/cache_table.h"
#include "worker/data.h"

#define FIRST_UNIQUE_ID 1000

#define RecordTs() timestamps[order++] = GetSizeStamp()
#define OutputTs(x, y) fprintf(time_log, "%d ; %s ; %s ; %f\n", \
                           tid, x, y, timestamps[i++])

namespace nimbus {

FILE* debug_log = NULL;

namespace {
double GetSizeStamp() {
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  return t.tv_sec + .000000001 * static_cast<double>(t.tv_nsec);
}
}  // namespace

/**
 * \details
 */
CacheManager::CacheManager() {
    unique_id_allocator_ = FIRST_UNIQUE_ID;
    pool_ = new Pool();
    pthread_mutex_init(&cache_lock, NULL);
    pthread_cond_init(&cache_cond, NULL);
}

CacheManager::~CacheManager() {}

/**
 * \detials CacheManager writes back all the dirty data from write_sets back to
 * nimbus objects. Note that this is not thread safe, this function does not
 * use locks or pending flags. This is guaranteed to work correctly only under
 * the assumption:
 * The mapping for the cached object and corresponding nimbus objects (dirty)
 * should not be edited simultaneously by other thread. This will hold for a
 * thread if the write_back set is a subset of write set for the job.
 * Locking should probably not add much overhead, but it is not necessary right
 * now.
 */
void CacheManager::WriteImmediately(AppVar *app_var,
                                    const DataArray &write_set) {
    double timestamps[6];
    int order = 0;
    RecordTs();  // start WIV stage
    DataArray flush_set;
    DataSet &write_back = app_var->write_back_;
    // size_t write_bytes = 0;
    RecordTs();  // start WIV mapping
    for (size_t i = 0; i < write_set.size(); ++i) {
        Data *d = write_set[i];
        if (write_back.find(d) != write_back.end()) {
            flush_set.push_back(d);
            // write_bytes += d->memory_size();
        }
    }
    for (size_t i = 0; i < flush_set.size(); ++i) {
        Data *d = flush_set[i];
        d->UnsetDirtyAppObject(app_var);
        write_back.erase(d);
    }
    RecordTs();  // end WIV mapping

    // std::string size_str = "WIV wfcsize " + app_var->name();
    // PrintSizeStamp(size_str.c_str(), write_bytes);

    RecordTs();  // start WIV wfc
    app_var->WriteAppData(flush_set, app_var->write_region_);
    RecordTs();  // end WIV wfc
    RecordTs();  // end WIV stage
    {
      int i = 0;
      pid_t tid = syscall(SYS_gettid);
      OutputTs("start", "WIV stage");
      OutputTs("start", "WIV mapping");
      OutputTs("end", "WIV mapping");
      OutputTs("start", "WIV wfc");
      OutputTs("end", "WIV wfc");
      OutputTs("end", "WIV stage");
      assert(i == order);
      assert(order <= 6);
    }
}

/**
 * \detials CacheManager writes back all the dirty data from write_sets back to
 * nimbus objects. Note that this is not thread safe, this function does not
 * use locks or pending flags. This is guaranteed to work correctly only under
 * the assumption:
 * The mapping for the cache object and corresponding nimbus objects (dirty) is
 * not edited simultaneously by other thread. This will hold for a thread if
 * the write_back set is a subset of write set for the job.
 * Locking should probably not add much overhead, but it is not necessary at
 * right now.
 */
void CacheManager::WriteImmediately(AppStruct *app_struct,
                                    const std::vector<app_data::type_id_t> &var_type,
                                    const std::vector<DataArray> &write_sets) {
    double timestamps[6];
    int order = 0;
    RecordTs();  // start WIS stage
    size_t num_vars = var_type.size();
    if (write_sets.size() != num_vars) {
        dbg(DBG_ERROR, "Mismatch in number of variable types passed to FlushCache\n");
        exit(-1);
    }
    std::vector<DataArray> flush_sets(num_vars);
    std::vector<DataSet> &write_backs = app_struct->write_backs_;
    // size_t write_bytes = 0;
    RecordTs();  // start WIS mapping
    for (size_t t = 0; t < num_vars; ++t) {
        DataArray &flush_t = flush_sets[t];
        const DataArray &write_set_t = write_sets[t];
        app_data::type_id_t type = var_type[t];
        DataSet &write_back_t = write_backs[type];
        for (size_t i = 0; i < write_set_t.size(); ++i) {
            Data *d = write_set_t[i];
            if (write_back_t.find(d) != write_back_t.end()) {
                flush_t.push_back(d);
                // write_bytes += d->memory_size();
            }
        }
    }
    for (size_t t = 0; t < num_vars; ++t) {
        const DataArray &flush_t = flush_sets[t];
        app_data::type_id_t type = var_type[t];
        DataSet &write_back_t = write_backs[type];
        for (size_t i = 0; i < flush_t.size(); ++i) {
            Data *d = flush_t[i];
            d->UnsetDirtyAppObject(app_struct);
            write_back_t.erase(d);
        }
    }
    RecordTs();  // end WIS mapping

    // std::string size_str = "WIS wfcsize " + app_struct->name();
    // PrintSizeStamp(size_str.c_str(), write_bytes);

    RecordTs();  // start WIS wfc
    app_struct->WriteAppData(var_type, flush_sets,
                             app_struct->write_region_);
    RecordTs();  // end WIS wfc
    RecordTs();  // end WIS stage
    {
      int i = 0;
      pid_t tid = syscall(SYS_gettid);
      OutputTs("start", "WIS stage");
      OutputTs("start", "WIS mapping");
      OutputTs("end", "WIS mapping");
      OutputTs("start", "WIS wfc");
      OutputTs("end", "WIS wfc");
      OutputTs("end", "WIS stage");
      assert(i == order);
      assert(order <= 6);
    }
}

/**
 * \details CacheManager checks if an instance with requested application
 * object region and prototype id is present and available for use. If not, it
 * creates a new instance, and adds the instance to its 2-level map. In both
 * cases, it acquires the object, assembles it (mappings, read diff, flush data
 * to be overwritten etc) and then returns.
 * The access mode is overwritten as EXCLUSIVE right now, meaning 2 parallel
 * compute jobs cannot operate on the same cache object.
 */
// TODO(hang/chinmayee): remove *aux, aux_data.
AppVar *CacheManager::GetAppVarV(const DataArray &read_set,
                                 const GeometricRegion &read_region,
                                 const DataArray &write_set,
                                 const GeometricRegion &write_region,
                                 const AppVar &prototype,
                                 const GeometricRegion &region,
                                 app_data::Access access,
                                 void (*aux)(AppVar*, void*),
                                 void* aux_data) {
    // TODO(chinmayee): Remove this when application objects can be shared by compute jobs
    access = app_data::EXCLUSIVE;
    double timestamps[14];
    int order = 0;
    RecordTs();  // start GAV stage
    RecordTs();  // start GAV lock
    pthread_mutex_lock(&cache_lock);
    RecordTs();  // end GAV lock
    AppVar *cv = NULL;
    // Get a cache object form the cache table.
    if (pool_->find(prototype.id()) == pool_->end()) {
        CacheTable *ct = new CacheTable(app_data::VAR);
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
    cv->AcquireAccess(access);
    DataArray flush, sync, diff;
    AppObjects sync_co;
    RecordTs();  // start GAV block
    while (!cv->CheckPendingFlag(read_set, write_set)) {
      pthread_cond_wait(&cache_cond, &cache_lock);
    }
    RecordTs();  // end GAV block
    // Move here.
    RecordTs();  // start GAV mapping
    cv->SetUpReadWrite(read_set, write_set,
                       &flush, &diff, &sync, &sync_co);
    RecordTs();  // end GAV mapping
    pthread_mutex_unlock(&cache_lock);
    if (aux != NULL) {
      aux(cv, aux_data);
    }

    GeometricRegion write_region_old = cv->write_region_;
    cv->write_region_ = write_region;

    // std::string size_str;
    // size_t write_bytes = 0;
    // for (size_t i = 0; i < flush.size(); ++i) {
    //   Data *d = flush[i];
    //   write_bytes += d->memory_size();
    // }
    // for (size_t i = 0; i < sync.size(); ++i) {
    //   Data *d = sync[i];
    //   write_bytes += d->memory_size();
    // }
    // std::string size_str = "GAV wfcsize " + cv->name();
    // PrintSizeStamp(size_str.c_str(), write_bytes);

    RecordTs();  // start GAV wfc
    cv->WriteAppData(flush, write_region_old);
    for (size_t i = 0; i < sync.size(); ++i) {
        sync_co[i]->PullData(sync[i]);
    }
    RecordTs();  // end GAV wfc

    // size_t read_bytes = 0;
    // for (size_t i = 0; i < diff.size(); ++i) {
    //   Data *d = diff[i];
    //   read_bytes += d->memory_size();
    // }
    // size_str = "GAV rtcsize " + cv->name();
    // PrintSizeStamp(size_str.c_str(), read_bytes);

    RecordTs();  // start GAV lock
    pthread_mutex_lock(&cache_lock);
    RecordTs();  // end GAV lock
    cv->ReleasePendingFlag(&flush, &diff, &sync);
    pthread_cond_broadcast(&cache_cond);
    pthread_mutex_unlock(&cache_lock);

    RecordTs();  // start GAV rtc
    cv->ReadAppData(diff, read_region);
    RecordTs();  // end GAV rtc

    RecordTs();  // end GAV stage
    {
      int i = 0;
      pid_t tid = syscall(SYS_gettid);
      OutputTs("start", "GAV stage");
      OutputTs("start", "GAV lock");
      OutputTs("end", "GAV lock");
      OutputTs("start", "GAV block");
      OutputTs("end", "GAV block");
      OutputTs("start", "GAV mapping");
      OutputTs("end", "GAV mapping");
      OutputTs("start", "GAV wfc");
      OutputTs("end", "GAV wfc");
      OutputTs("start", "GAV lock");
      OutputTs("end", "GAV lock");
      OutputTs("start", "GAV rtc");
      OutputTs("end", "GAV rtc");
      OutputTs("end", "GAV stage");
      assert(i == order);
      assert(order <= 14);
    }
    if (diff.size() != 0) {
      size_t max_i = 0;
      for (size_t i = 1; i < diff.size(); ++i)
        if (diff.at(max_i)->region().GetSurfaceArea()
            < diff.at(i)->region().GetSurfaceArea()) {
          max_i = i;
        }
      std::string cv_name = cv->name();
      Data* data = diff.at(max_i);
      std::string region_name = data->region().ToNetworkData();
      fprintf(debug_log,
              "%s %s %f\n",
              cv_name.c_str(),
              region_name.c_str(),
              timestamps[12] - timestamps[11]);
    }
    return cv;
}

/**
 * \details CacheManager checks if an instance with requested application
 * object region and prototype id is present and available for use. If not, it
 * creates a new instance, and adds the instance to its 2-level map. In both
 * cases, it acquires the object, assembles it (mappings, read diff, flush data
 * to be overwritten etc) and then returns.
 * The access mode is overwritten as EXCLUSIVE right now, meaning 2 parallel
 * compute jobs cannot operate on the same cache object.
 */
AppStruct *CacheManager::GetAppStructV(const std::vector<app_data::type_id_t> &var_type,
                                       const std::vector<DataArray> &read_sets,
                                       const GeometricRegion &read_region,
                                       const std::vector<DataArray> &write_sets,
                                       const GeometricRegion &write_region,
                                       const AppStruct &prototype,
                                       const GeometricRegion &region,
                                       app_data::Access access) {
    // TODO(chinmayee): Remove this when application objects can be shared by compute jobs
    access = app_data::EXCLUSIVE;
    double timestamps[14];
    int order = 0;
    RecordTs();  // start GAS stage
    RecordTs();  // start GAS lock
    pthread_mutex_lock(&cache_lock);
    RecordTs();  // end GAS lock
    AppStruct *cs = NULL;
    if (pool_->find(prototype.id()) == pool_->end()) {
        CacheTable *ct = new CacheTable(app_data::STRUCT);
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
    size_t num_var = var_type.size();
    std::vector<DataArray> flush_sets(num_var),
                           sync_sets(num_var),
                           diff_sets(num_var);
    std::vector<AppObjects> sync_co_sets(num_var);
    RecordTs();  // start GAS block
    // Move here.
    cs->AcquireAccess(access);
    while (!cs->CheckPendingFlag(var_type, read_sets, write_sets)) {
      pthread_cond_wait(&cache_cond, &cache_lock);
    }
    RecordTs();  // end GAS block
    RecordTs();  // start GAS mapping
    cs->SetUpReadWrite(var_type, read_sets, write_sets,
                       &flush_sets, &diff_sets, &sync_sets, &sync_co_sets);
    RecordTs();  // end GAS mapping
    pthread_mutex_unlock(&cache_lock);

    GeometricRegion write_region_old = cs->write_region_;
    cs->write_region_ = write_region;

    // std::string size_str = "GAS wfcsize " + cs->name();
    // size_t write_bytes = 0;
    // for (size_t i = 0; i < num_var; ++i) {
    //   DataArray &flush_t = flush_sets[i];
    //   for (size_t j = 0; j < flush_t.size(); ++j) {
    //     Data *d = flush_t[j];
    //     write_bytes += d->memory_size();
    //   }
    // }
    // for (size_t i = 0; i < num_var; ++i) {
    //   DataArray &sync_t = sync_sets[i];
    //   for (size_t j = 0; j < sync_t.size(); ++j) {
    //     Data *d = sync_t[j];
    //     write_bytes += d->memory_size();
    //   }
    // }
    // size_str = "GAS wfcsize " + cs->name();
    // PrintSizeStamp(size_str.c_str(), write_bytes);

    RecordTs();  // start GAS wfc
    cs->WriteAppData(var_type, flush_sets, write_region_old);
    for (size_t t = 0; t < num_var; ++t) {
        DataArray &sync_t = sync_sets[t];
        AppObjects &sync_co_t = sync_co_sets[t];
        for (size_t i = 0; i < sync_t.size(); ++i) {
            sync_co_t[i]->PullData(sync_t[i]);
        }
    }
    RecordTs();  // end GAS wfc

    // size_t read_bytes = 0;
    // for (size_t i = 0; i < num_var; ++i) {
    //   DataArray &diff_t = diff_sets[i];
    //   for (size_t j = 0; j < diff_t.size(); ++j) {
    //     Data *d = diff_t[j];
    //     read_bytes += d->memory_size();
    //   }
    // }
    // size_str = "GAS rtcsize " + cs->name();
    // PrintSizeStamp(size_str.c_str(), read_bytes);

    RecordTs();  // start GAS lock
    pthread_mutex_lock(&cache_lock);
    RecordTs();  // end GAS lock
    cs->ReleasePendingFlag(var_type,
                           &flush_sets, &diff_sets, &sync_sets);
    pthread_cond_broadcast(&cache_cond);
    pthread_mutex_unlock(&cache_lock);

    RecordTs();  // start GAS rtc
    cs->ReadAppData(var_type, diff_sets, read_region);
    RecordTs();  // end GAS rtc

    RecordTs();  // end GAS stage
    {
      int i = 0;
      pid_t tid = syscall(SYS_gettid);
      OutputTs("start", "GAS stage");
      OutputTs("start", "GAS lock");
      OutputTs("end", "GAS lock");
      OutputTs("start", "GAS block");
      OutputTs("end", "GAS block");
      OutputTs("start", "GAS mapping");
      OutputTs("end", "GAS mapping");
      OutputTs("start", "GAS wfc");
      OutputTs("end", "GAS wfc");
      OutputTs("start", "GAS lock");
      OutputTs("end", "GAS lock");
      OutputTs("start", "GAS rtc");
      OutputTs("end", "GAS rtc");
      OutputTs("end", "GAS stage");
      assert(i == order);
      assert(order <= 14);
    }
    return cs;
}

/**
 * \details Pulls data from dirty app object (if there is one), updates
 * mappings and returns.
 */
void CacheManager::SyncData(Data *d) {
    double timestamps[14];
    int order = 0;
    RecordTs();  // start SD stage
    AppObject *co = NULL;
    RecordTs();  // start SD lock
    pthread_mutex_lock(&cache_lock);
    RecordTs();  // end SD lock
    RecordTs();  // start SD block
    while (d->pending_flag() != 0) {
       pthread_cond_wait(&cache_cond, &cache_lock);
    }
    RecordTs();  // end SD block
    co = d->dirty_app_object();
    if (!co) {
        pthread_mutex_unlock(&cache_lock);
        RecordTs();  // end SD stage
        {
          int i = 0;
          pid_t tid = syscall(SYS_gettid);
          OutputTs("start", "SD stage");
          OutputTs("start", "SD lock");
          OutputTs("end", "SD lock");
          OutputTs("start", "SD block");
          OutputTs("end", "SD block");
          OutputTs("end", "SD stage");
          assert(i == order);
          assert(order <= 6);
        }
        return;
    }
    // assert(co->IsAvailable(app_data::EXCLUSIVE));
    RecordTs();  // start SD mapping
    d->set_pending_flag(Data::WRITE);
    RecordTs();  // end SD mapping
    pthread_mutex_unlock(&cache_lock);

    // size_t write_bytes = d->memory_size();
    // std::string size_str = "SD wfcsize " + co->name();
    // PrintSizeStamp(size_str.c_str(), write_bytes);

    RecordTs();  // start SD wfc
    co->PullData(d);
    RecordTs();  // end SD wfc
    RecordTs();  // start SD lock
    pthread_mutex_lock(&cache_lock);
    RecordTs();  // end SD lock
    RecordTs();  // start SD mapping
    d->ClearDirtyMappings();
    d->unset_pending_flag(Data::WRITE);
    RecordTs();  // end SD mapping
    pthread_cond_broadcast(&cache_cond);
    pthread_mutex_unlock(&cache_lock);
    RecordTs();  // end SD stage
    {
      int i = 0;
      pid_t tid = syscall(SYS_gettid);
      OutputTs("start", "SD stage");
      OutputTs("start", "SD lock");
      OutputTs("end", "SD lock");
      OutputTs("start", "SD block");
      OutputTs("end", "SD block");
      OutputTs("start", "SD mapping");
      OutputTs("end", "SD mapping");
      OutputTs("start", "SD wfc");
      OutputTs("end", "SD wfc");
      OutputTs("start", "SD lock");
      OutputTs("end", "SD lock");
      OutputTs("start", "SD mapping");
      OutputTs("end", "SD mapping");
      OutputTs("end", "SD stage");
      assert(i == order);
      assert(order <= 14);
    }
}

/**
 * \details Removes mappings between this data object and any corresponding
 * cached application objects.
 */
void CacheManager::InvalidateMappings(Data *d) {
    double timestamps[8];
    int order = 0;
    RecordTs();  // start IM stage
    RecordTs();  // start IM lock
    pthread_mutex_lock(&cache_lock);
    RecordTs();  // end IM lock
    RecordTs();  // start IM block
    while (d->pending_flag() != 0) {
        pthread_cond_wait(&cache_cond, &cache_lock);
    }
    RecordTs();  // end IM block
    RecordTs();  // start IM mapping
    d->InvalidateMappings();
    RecordTs();  // end IM mapping
    pthread_mutex_unlock(&cache_lock);
    RecordTs();  // end IM stage
    {
      int i = 0;
      pid_t tid = syscall(SYS_gettid);
      OutputTs("start", "IM stage");
      OutputTs("start", "IM lock");
      OutputTs("end", "IM lock");
      OutputTs("start", "IM block");
      OutputTs("end", "IM block");
      OutputTs("start", "IM mapping");
      OutputTs("end", "IM mapping");
      OutputTs("end", "IM stage");
      assert(i == order);
      assert(order <= 8);
    }
}

void CacheManager::ReleaseAccess(AppObject* app_object) {
    pthread_mutex_lock(&cache_lock);
    app_object->ReleaseAccessInternal();
    pthread_cond_broadcast(&cache_cond);
    pthread_mutex_unlock(&cache_lock);
}

void CacheManager::SetLogNames(std::string wid_str) {
    time_log = fopen((wid_str + "_cache_time.txt").c_str(), "w");
    debug_log = fopen((wid_str + "_cache_debug.txt").c_str(), "w");
}


}  // namespace nimbus

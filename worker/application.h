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
  * Nimbus abstraction of an application. Programmers use this base class to
  * write various application served by Nimbus.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  *
  * A job should interact with Nimbus runtime system through interface in
  * application. And the methods offered in applicaiton class are expected to be
  * thread-safe.  --quhang
  *
  */


#ifndef NIMBUS_WORKER_APPLICATION_H_
#define NIMBUS_WORKER_APPLICATION_H_

#include <pthread.h>
#include <map>
#include <string>
#include <vector>
#include "shared/id_maker.h"
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"
#include "shared/scheduler_client.h"
#include "shared/scheduler_command_include.h"
#include "worker/data.h"
#include "worker/job.h"
#include "worker/worker_ldo_map.h"

namespace nimbus {

class Application;
typedef std::map<int, Application*> AppMap;

// forward declaration
class CacheManager;

class Application {
 public:
  Application();
  ~Application();

  virtual void Load();
  virtual void Start(SchedulerClient* client, IDMaker* id_maker, WorkerLdoMap* ldo_map);

  void RegisterJob(std::string name, Job* job);
  void RegisterData(std::string name, Data* data);

  void SpawnComputeJob(const std::string& name,
                       const job_id_t& id,
                       const IDSet<logical_data_id_t>& read,
                       const IDSet<logical_data_id_t>& write,
                       const IDSet<job_id_t>& before,
                       const IDSet<job_id_t>& after,
                       const job_id_t& parent_id,
                       const job_id_t& future_id,
                       const bool& sterile,
                       const GeometricRegion& region,
                       const Parameter& params);


  void SpawnCopyJob(const job_id_t& id,
                    const logical_data_id_t& from_logical_id,
                    const logical_data_id_t& to_logical_id,
                    const IDSet<job_id_t>& before,
                    const IDSet<job_id_t>& after,
                    const job_id_t& parent_id);

  void DefineData(const std::string& name,
                  const logical_data_id_t& logical_data_id,
                  const partition_id_t& partition_id,
                  const IDSet<partition_id_t>& neighbor_partition,
                  const job_id_t& parent_id);

  void DefinePartition(const ID<partition_id_t>& partition_id,
                       const GeometricRegion& r);

  void TerminateApplication(const exit_status_t& exit_status_id);

  Job* CloneJob(std::string name);
  Data* CloneData(std::string name);

  bool GetNewJobID(std::vector<job_id_t>* result, size_t req_num);
  bool GetNewLogicalDataID(std::vector<logical_data_id_t>* result, size_t req_num);

  bool GetPartition(partition_id_t id, GeometricRegion* r);

  const LogicalDataObject* GetLogicalObject(logical_data_id_t id);
  int GetCoveredLogicalObjects(CLdoVector* result,
                               const std::string& variable,
                               const GeometricRegion* r);
  int GetAdjacentLogicalObjects(CLdoVector* result,
                                const std::string& variable,
                                const GeometricRegion* r);
  int GetIntersectingLogicalObjects(CLdoVector* result,
                                    const std::string& variable,
                                    const GeometricRegion* r);

  CacheManager* cache_manager() const;

 private:
  app_id_t id_;
  size_t priority_;
  // Protects job table.
  pthread_mutex_t lock_job_table_;
  JobTable job_table_;
  // Protects data table.
  pthread_mutex_t lock_data_table_;
  DataTable data_table_;
  SchedulerClient* client_;
  IDMaker* id_maker_;
  WorkerLdoMap* ldo_map_;
  CacheManager* cache_manager_;
  class LockGuard {
   public:
    explicit LockGuard(pthread_mutex_t* lock) {
      pthread_mutex_lock(lock);
      lock_ = lock;
    }
    virtual ~LockGuard() {
      pthread_mutex_unlock(lock_);
    }
   private:
    pthread_mutex_t* lock_;
  };
};

}  //  namespace nimbus

#endif  // NIMBUS_WORKER_APPLICATION_H_

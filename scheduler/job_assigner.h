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
  * This is the base job assigner module. It provides methods required for
  * binding a job to a worker and submitting the compute and copy jobs to the
  * worker.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_JOB_ASSIGNER_H_
#define NIMBUS_SCHEDULER_JOB_ASSIGNER_H_

#include <boost/unordered_map.hpp>
#include <boost/thread.hpp>
#include <map>
#include <list>
#include <vector>
#include <string>
#include "shared/nimbus_types.h"
#include "scheduler/job_entry.h"
#include "scheduler/job_profile.h"
#include "scheduler/data_manager.h"
#include "scheduler/job_manager.h"
#include "scheduler/region_map.h"
#include "scheduler/straggler_map.h"
#include "shared/cluster.h"
#include "shared/id_maker.h"
#include "shared/scheduler_server.h"
#include "shared/geometric_region.h"
#include "shared/graph.h"

namespace nimbus {

  class LoadBalancer;

  class JobAssigner {
  public:
    JobAssigner();
    virtual ~JobAssigner();

    void Run();

    virtual void set_id_maker(IDMaker *id_maker);
    virtual void set_server(SchedulerServer* server);
    virtual void set_job_manager(JobManager *job_manager);
    virtual void set_data_manager(DataManager *data_manager);
    virtual void set_load_balancer(LoadBalancer *load_balancer);
    virtual void set_thread_num(size_t thread_num);
    virtual void set_checkpoint_id(checkpoint_id_t checkpoint_id);

    virtual void AssignJobs(const JobEntryList& list);

  protected:
    Log log_;
    Log log_assign_;
    Log log_prepare_;
    Log log_job_manager_;
    Log log_data_manager_;
    Log log_version_manager_;
    Log log_before_set_;
    Log log_assign_stamp_;
    IDMaker *id_maker_;
    SchedulerServer* server_;
    JobManager *job_manager_;
    DataManager *data_manager_;
    LoadBalancer *load_balancer_;
    size_t thread_num_;
    checkpoint_id_t checkpoint_id_;

    JobEntryList job_queue_;
    boost::recursive_mutex job_queue_mutex_;
    boost::condition_variable_any job_queue_cond_;
    size_t pending_assignment_;

    std::list<boost::thread*> job_assigner_threads_;

    virtual void Initialize();

    virtual void JobAssignerThread();

    virtual bool AssignJob(JobEntry* job);

    virtual bool PrepareDataForJobAtWorker(JobEntry* job,
                                           SchedulerWorker* worker,
                                           logical_data_id_t l_id);

    virtual bool AllocateLdoInstanceToJob(JobEntry* job,
                                          LogicalDataObject* ldo,
                                          PhysicalData pd);

    virtual bool SaveJobContextForCheckpoint(JobEntry *job);

    virtual bool CreateDataAtWorker(SchedulerWorker* worker,
                                    LogicalDataObject* ldo,
                                    PhysicalData* created_data);

    virtual bool RemoteCopyData(SchedulerWorker* from_worker,
                                SchedulerWorker* to_worker,
                                LogicalDataObject* ldo,
                                PhysicalData* from_data,
                                PhysicalData* to_data);

    virtual bool LocalCopyData(SchedulerWorker* worker,
                               LogicalDataObject* ldo,
                               PhysicalData* created_data,
                               PhysicalData* to_data);

    virtual bool SaveData(SchedulerWorker* worker,
                          LogicalDataObject* ldo,
                          PhysicalData* from_data,
                          checkpoint_id_t checkpoint_id);

    virtual bool LoadData(SchedulerWorker* worker,
                          LogicalDataObject* ldo,
                          data_version_t version,
                          PhysicalData* to_data,
                          std::string handle);

    virtual bool GetFreeDataAtWorker(SchedulerWorker* worker,
                                     LogicalDataObject* ldo,
                                     PhysicalData* free_data);

    virtual size_t GetObsoleteLdoInstancesAtWorker(SchedulerWorker* worker,
                                                   LogicalDataObject* ldo,
                                                   PhysicalDataList* dest);




    virtual bool SendComputeJobToWorker(SchedulerWorker* worker,
                                        JobEntry* job);

    virtual void PrintLog(JobEntry *job);

  private:
    JobAssigner(const JobAssigner& other) {}
  };

}  // namespace nimbus

#endif  // NIMBUS_SCHEDULER_JOB_ASSIGNER_H_

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
  * This is the base class that serves the scheduler regarding job assignment
  * queries in the cluster. It tries to minimize the completion of simulation
  * by reducing the cost of communication by locality aware data placement and
  * mitigate the effect of stragglers in the system by adapting the job
  * assignment strategies to the dynamic changes of the cloud.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_LOAD_BALANCER_H_
#define NIMBUS_SCHEDULER_LOAD_BALANCER_H_

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
#include "scheduler/job_assigner.h"
#include "shared/cluster.h"
#include "shared/id_maker.h"
#include "shared/scheduler_server.h"
#include "shared/geometric_region.h"
#include "shared/graph.h"

namespace nimbus {

  class LoadBalancer {
  public:
    LoadBalancer();
    virtual ~LoadBalancer();

    ClusterMap* cluster_map();

    virtual bool safe_to_load_balance();
    virtual void set_cluster_map(ClusterMap* cluster_map);
    virtual void set_job_manager(JobManager *job_manager);
    virtual void set_data_manager(DataManager *data_manager);
    virtual void set_job_assigner(JobAssigner *job_assigner);
    virtual void set_assign_batch_size(size_t num);
    virtual void set_split_dimensions(const std::vector<size_t>& split);
    virtual void set_sub_split_dimensions(const std::vector<size_t>& split);
    virtual void set_global_region(const GeometricRegion& region);


    virtual void Run();

    virtual size_t AssignReadyJobs();

    virtual void NotifyJobAssignment(const JobEntry *job);

    virtual void NotifyJobDone(const JobEntry *job);

    virtual void NotifyRegisteredWorker(SchedulerWorker *worker);

    virtual bool NotifyDownWorker(worker_id_t worker_id);

    virtual bool SetWorkerToAssignJob(JobEntry *job);

    virtual bool BalanceLoad(counter_t query_id);

    virtual bool AddWorkerStat(const counter_t& query_id,
                               const worker_id_t& worker_id,
                               const int64_t& run_time,
                               const int64_t& block_time,
                               const int64_t& idle_time);

  protected:
    ClusterMap* cluster_map_;
    JobManager *job_manager_;
    DataManager *data_manager_;
    JobAssigner *job_assigner_;
    size_t assign_batch_size_;
    load_balancing_id_t load_balancing_id_;
    bool safe_to_load_balance_;

    bool enforced_split_;
    bool enforced_sub_split_;
    bool enforced_global_region_;
    std::vector<size_t> split_;
    std::vector<size_t> sub_split_;
    GeometricRegion global_region_;

    class WorkerStat {
      public:
        WorkerStat(const counter_t& query_id,
                   const worker_id_t& worker_id,
                   const int64_t& run_time,
                   const int64_t& block_time,
                   const int64_t& idle_time)
          : query_id_(query_id),
            worker_id_(worker_id),
            run_time_(run_time),
            block_time_(block_time),
            idle_time_(idle_time) {}
        ~WorkerStat() {}

        counter_t query_id_;
        worker_id_t worker_id_;
        int64_t run_time_;
        int64_t block_time_;
        int64_t idle_time_;
    };

    typedef boost::unordered_map<worker_id_t, WorkerStat*> StatQuery;
    typedef std::map<counter_t, StatQuery*> QueryMap;

    QueryMap queries_;

  private:
    LoadBalancer(const LoadBalancer& other) {}

    typedef std::map<worker_id_t, SchedulerWorker*> WorkerMap;
    WorkerMap worker_map_;
  };

}  // namespace nimbus

#endif  // NIMBUS_SCHEDULER_LOAD_BALANCER_H_

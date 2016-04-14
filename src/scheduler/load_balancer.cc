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


#include "src/scheduler/load_balancer.h"

namespace nimbus {

LoadBalancer::LoadBalancer() {
  cluster_map_ = NULL;
  job_manager_ = NULL;
  data_manager_ = NULL;
  job_assigner_ = NULL;
  assign_batch_size_ = 0;
  load_balancing_id_ = NIMBUS_INIT_LOAD_BALANCING_ID;
  // by default it is safe to load balance, unless job manager is memoizing
  // binding, and so we cannot load balance. -omidm
  safe_to_load_balance_ = true;

  enforced_split_ = false;
  enforced_sub_split_ = false;
  enforced_global_region_ = false;
  split_ = std::vector<size_t>(3);
  sub_split_ = std::vector<size_t>(3);
  global_region_ = GeometricRegion();
}

LoadBalancer::~LoadBalancer() {
}

ClusterMap* LoadBalancer::cluster_map() {
  return cluster_map_;
}

bool LoadBalancer::safe_to_load_balance() {
  return safe_to_load_balance_;
}

void LoadBalancer::set_cluster_map(ClusterMap* cluster_map) {
  cluster_map_ = cluster_map;
}

void LoadBalancer::set_job_manager(JobManager *job_manager) {
  job_manager_ = job_manager;
}

void LoadBalancer::set_data_manager(DataManager *data_manager) {
  data_manager_ = data_manager;
}

void LoadBalancer::set_job_assigner(JobAssigner *job_assigner) {
  job_assigner_ = job_assigner;
}

void LoadBalancer::set_assign_batch_size(size_t num) {
  assign_batch_size_ = num;
}

void LoadBalancer::set_split_dimensions(const std::vector<size_t>& split) {
  split_ = split;
  if (split_[0] && split_[1] && split_[2]) {
    enforced_split_ = true;
  }
}

void LoadBalancer::set_sub_split_dimensions(const std::vector<size_t>& sub_split) {
  sub_split_ = sub_split;
  if (sub_split_[0] && sub_split_[1] && sub_split_[2]) {
    enforced_split_ = true;
  }
}

void LoadBalancer::set_global_region(const GeometricRegion& region) {
  global_region_ = region;
  if (global_region_ != GeometricRegion()) {
    enforced_global_region_ = true;
  }
}


void LoadBalancer::Run() {
  dbg(DBG_WARN, "WARNING: Base load balancer is being used for LoadBalancer::Run!\n");
}

size_t LoadBalancer::AssignReadyJobs() {
  JobEntryList list;
  job_manager_->GetJobsReadyToAssign(&list,
                                    assign_batch_size_,
                                    load_balancing_id_,
                                    &safe_to_load_balance_);

  JobEntryList::iterator iter;
  for (iter = list.begin(); iter != list.end(); ++iter) {
    JobEntry* job = *iter;
    if (!SetWorkerToAssignJob(job)) {
      dbg(DBG_ERROR, "ERROR: LoadBalancer: could not get worker to assign job %lu.\n", job->job_id()); // NOLINT
      assert(false);
    }
#ifdef _RUN_MULTI_TENANT_SCENARIO
    if (job->job_type() == JOB_CMPX) {
      static int counter = 10;
      if (--counter == 0) {
        if ((load_balancing_id_ % 2) == 1) {
          ++load_balancing_id_;
        } else {
          --load_balancing_id_;
        }
        counter = 10;
      }
    } else if (job->job_name().find("__MARK_STAT") != std::string::npos) {
      static size_t warmup = 0;
      switch (warmup) {
        case 0:
          ++warmup;
          break;
        case 1:
        case 2:
          ++load_balancing_id_;
          ++warmup;
          break;
        default:
          break;
      }
    }
#endif
  }

  job_assigner_->AssignJobs(list);

  return list.size();
}

bool LoadBalancer::SetWorkerToAssignJob(JobEntry *job) {
  static bool printed_warn = false;
  if (!printed_warn) {
    dbg(DBG_WARN, "WARNING: Base load balancer is being used for SetWorkerToAssignJob, all jobs are assigned to first worker!\n"); // NOLINT
    printed_warn = true;
  }
  assert(worker_map_.size() > 0);
  job->set_assigned_worker(worker_map_.begin()->second);
  return true;
}

bool LoadBalancer::BalanceLoad(counter_t query_id) {
  static bool printed_warn = false;
  if (!printed_warn) {
    dbg(DBG_WARN, "WARNING: Base load balancer is being used for BalanceLoad, no load balancing happens!\n"); // NOLINT
    printed_warn = true;
  }
  QueryMap::iterator iter = queries_.find(query_id);
  assert(iter != queries_.end());
  return true;
}

void LoadBalancer::NotifyJobAssignment(const JobEntry *job) {
  static bool printed_warn = false;
  if (!printed_warn) {
    dbg(DBG_WARN, "WARNING: Base load balancer is being used for NotifyJobAssignment!\n");
    printed_warn = true;
  }
}

void LoadBalancer::NotifyJobDone(const JobEntry *job) {
  static bool printed_warn = false;
  if (!printed_warn) {
    dbg(DBG_WARN, "WARNING: Base load balancer is being used for NotifyJobDone!\n");
    printed_warn = true;
  }
}


void LoadBalancer::NotifyRegisteredWorker(SchedulerWorker *worker) {
  static bool printed_warn = false;
  if (!printed_warn) {
    dbg(DBG_WARN, "WARNING: Base load balancer is being used for NotifyRegisteredWorker!\n");
    printed_warn = true;
  }
  worker_id_t worker_id = worker->worker_id();
  WorkerMap::iterator iter = worker_map_.find(worker_id);
  if (iter == worker_map_.end()) {
    worker_map_[worker_id] = worker;
  } else {
    dbg(DBG_ERROR, "ERROR: LoadBalancer: worker with the same id %lu has already been registered.\n", // NOLINT
        worker_id);
  }
}

bool LoadBalancer::NotifyDownWorker(worker_id_t worker_id) {
  static bool printed_warn = false;
  if (!printed_warn) {
    dbg(DBG_WARN, "WARNING: Base load balancer is being used for NotifyDownWorker!\n");
    printed_warn = true;
  }
  WorkerMap::iterator iter = worker_map_.find(worker_id);
  if (iter != worker_map_.end()) {
    worker_map_.erase(iter);
    return true;
  } else {
    dbg(DBG_ERROR, "ERROR: LoadBalancer: worker id %lu has not been registered.\n", // NOLINT
        worker_id);
    assert(false);
    return false;
  }
}

bool LoadBalancer::AddWorkerStat(const counter_t& query_id,
                                 const worker_id_t& worker_id,
                                 const int64_t& run_time,
                                 const int64_t& block_time,
                                 const int64_t& idle_time) {
  WorkerStat *ws = new WorkerStat(query_id,
                                  worker_id,
                                  run_time,
                                  block_time,
                                  idle_time);
  QueryMap::iterator iter = queries_.find(query_id);
  if (iter != queries_.end()) {
    iter->second->operator[](worker_id) = ws;
  } else {
    StatQuery *sq = new StatQuery();
    queries_[query_id] = sq;
    sq->operator[](worker_id) = ws;
  }

  return true;
}

}  // namespace nimbus

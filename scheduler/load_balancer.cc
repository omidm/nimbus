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


#include "scheduler/load_balancer.h"

namespace nimbus {

LoadBalancer::LoadBalancer() {
  Initialize();
}

void LoadBalancer::Initialize() {
  cluster_map_ = NULL;
  job_manager_ = NULL;
  data_manager_ = NULL;
  job_assigner_ = NULL;
  max_job_to_assign_ = 0;
}

LoadBalancer::~LoadBalancer() {
}

ClusterMap* LoadBalancer::cluster_map() {
  return cluster_map_;
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

void LoadBalancer::set_max_job_to_assign(size_t num) {
  max_job_to_assign_ = num;
}

void LoadBalancer::Run() {
  dbg(DBG_WARN, "WARNING: Base load balancer is being used!!!\n");
}

size_t LoadBalancer::AssignReadyJobs() {
  dbg(DBG_WARN, "WARNING: Base load balancer is being used!!!\n");
  JobEntryList list;
  job_manager_->GetJobsReadyToAssign(&list, max_job_to_assign_);

  JobEntryList::iterator iter;
  for (iter = list.begin(); iter != list.end(); ++iter) {
    JobEntry* job = *iter;
    if (!SetWorkerToAssignJob(job)) {
      dbg(DBG_ERROR, "ERROR: LoadBalancer: could not get worker to assign job %lu.\n", job->job_id()); // NOLINT
      exit(-1);
    }
  }

  job_assigner_->AssignJobs(list);

  return list.size();
}

bool LoadBalancer::SetWorkerToAssignJob(JobEntry *job) {
  dbg(DBG_WARN, "WARNING: Base load balancer is being used, all jobs are assigned to first worker.!!!\n"); // NOLINT
  assert(worker_map_.size() > 0);
  job->set_assigned_worker(worker_map_.begin()->second);
  return true;
}


void LoadBalancer::NotifyJobAssignment(const JobEntry *job) {
  dbg(DBG_WARN, "WARNING: Base load balancer is being used!!!\n");
}

void LoadBalancer::NotifyJobDone(const JobEntry *job) {
  dbg(DBG_WARN, "WARNING: Base load balancer is being used!!!\n");
}


void LoadBalancer::NotifyRegisteredWorker(SchedulerWorker *worker) {
  dbg(DBG_WARN, "WARNING: Base load balancer is being used!!!\n");
  worker_id_t worker_id = worker->worker_id();
  WorkerMap::iterator iter = worker_map_.find(worker_id);
  if (iter == worker_map_.end()) {
    worker_map_[worker_id] = worker;
  } else {
    dbg(DBG_ERROR, "ERROR: LoadBalancer: worker with the same id %lu has already been registered.\n", // NOLINT
        worker_id);
  }
}

}  // namespace nimbus

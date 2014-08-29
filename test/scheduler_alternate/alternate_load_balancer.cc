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
  * This is static version of the load balancer used in scheduler v2.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#include "./alternate_load_balancer.h"

#define SEED_ 123
#define WRITE_FRAME_NAME "write_frame"

namespace nimbus {

AlternateLoadBalancer::AlternateLoadBalancer() {
  seed_ = SEED_;
  log_.set_file_name("load_balancer_log");
}

AlternateLoadBalancer::~AlternateLoadBalancer() {
}

void AlternateLoadBalancer::Run() {
}

size_t AlternateLoadBalancer::AssignReadyJobs() {
  JobEntryList list;
  job_manager_->GetJobsReadyToAssign(&list, max_job_to_assign_);

  JobEntryList::iterator iter;
  for (iter = list.begin(); iter != list.end(); ++iter) {
    JobEntry* job = *iter;
    if (!SetWorkerToAssignJob(job)) {
      dbg(DBG_ERROR, "ERROR: AlternateLoadBalancer: could not get worker to assign job %lu.\n", job->job_id()); // NOLINT
      exit(-1);
    }
  }

  job_assigner_->AssignJobs(list);

  return list.size();
}

bool AlternateLoadBalancer::SetWorkerToAssignJob(JobEntry* job) {
  log_.StartTimer();

  size_t worker_num = worker_map_.size();

  if (worker_num < 1) {
    dbg(DBG_SCHED, "ERROR: there is no worker in scheduler worker list for job assignment");
    return false;
  }

  size_t index;
  if ((job->job_name() == WRITE_FRAME_NAME) ||
      (!job->sterile())) {
    index = 0;
  } else {
    index  = (rand_r(&seed_) % worker_num);
  }
  WorkerMap::iterator it = worker_map_.begin();
  std::advance(it, index);

  log_.StopTimer();
  std::cout
    << "ALTERNATE: Picked worker: " << it->first
    << " for job: " << job->job_name()
    << " took: " << log_.timer()
    << " for union set size of: " << job->union_set_p()->size() << std::endl;

  job->set_assigned_worker(it->second);
  return true;
}

void AlternateLoadBalancer::NotifyJobAssignment(const JobEntry *job) {
}

void AlternateLoadBalancer::NotifyJobDone(const JobEntry *job) {
}


void AlternateLoadBalancer::NotifyRegisteredWorker(SchedulerWorker *worker) {
  worker_id_t worker_id = worker->worker_id();
  WorkerMap::iterator iter = worker_map_.find(worker_id);
  if (iter == worker_map_.end()) {
    worker_map_[worker_id] = worker;
  } else {
    dbg(DBG_ERROR, "ERROR: AlternateLoadBalancer: worker with the same id %lu has already been registered.\n", // NOLINT
        worker_id);
  }
}

}  // namespace nimbus

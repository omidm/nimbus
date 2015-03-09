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
  * This is dynamic load balancer that has load ba;ancer and base class and
  * overrides some of its methods.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#include "./dynamic_load_balancer.h"

#define LB_UPDATE_RATE 2000

namespace nimbus {

DynamicLoadBalancer::DynamicLoadBalancer() {
  worker_num_ = 0;
  update_ = false;
  init_phase_ = true;
  stamp_state_ = -1;
  counter_ = 0;
  log_.set_file_name("load_balancer_log");
  global_region_ = GeometricRegion(0, 0, 0, 0, 0, 0);
}

DynamicLoadBalancer::~DynamicLoadBalancer() {
}

void DynamicLoadBalancer::Run() {
  while (true) {
    boost::unique_lock<boost::recursive_mutex> update_lock(update_mutex_);
    while (!update_) {
      update_cond_.wait(update_lock);
    }
    update_ = false;
    update_cond_.notify_all();

    boost::unique_lock<boost::recursive_mutex> worker_map_lock(worker_map_mutex_);
    boost::unique_lock<boost::recursive_mutex> region_map_lock(region_map_mutex_);

    if (worker_num_ != region_map_.table_size() ||
        global_region_ != data_manager_->global_bounding_region()) {
      if (!data_manager_->initialized_global_bounding_region()) {
        continue;
      } else {
        InitializeRegionMap();
      }
    } else {
      UpdateRegionMap();
    }
  }
}

size_t DynamicLoadBalancer::AssignReadyJobs() {
  JobEntryList list;
  job_manager_->GetJobsReadyToAssign(&list, max_job_to_assign_);

  JobEntryList::iterator iter;
  for (iter = list.begin(); iter != list.end(); ++iter) {
    JobEntry* job = *iter;
    if (!SetWorkerToAssignJob(job)) {
      dbg(DBG_ERROR, "ERROR: DynamicLoadBalancer: could not get worker to assign job %lu.\n", job->job_id()); // NOLINT
      exit(-1);
    }
  }

  job_assigner_->AssignJobs(list);

  return list.size();
}

bool DynamicLoadBalancer::SetWorkerToAssignJob(JobEntry *job) {
  Log log(Log::NO_FILE);
  log.log_StartTimer();

  boost::unique_lock<boost::recursive_mutex> update_lock(update_mutex_);
  while (update_) {
    update_cond_.wait(update_lock);
  }

  boost::unique_lock<boost::recursive_mutex> worker_map_lock(worker_map_mutex_);
  boost::unique_lock<boost::recursive_mutex> region_map_lock(region_map_mutex_);

  if ((worker_num_ != region_map_.table_size()) ||
      (global_region_ != data_manager_->global_bounding_region())) {
    if (data_manager_->initialized_global_bounding_region()) {
      InitializeRegionMap();
    }
  }

  assert(worker_map_.size() > 0);
  assert(worker_num_ > 0);

  GeometricRegion region;

  if (init_phase_ || (job->union_set_p()->size() == 0) || (!job->sterile())) {
    job->set_assigned_worker(worker_map_.begin()->second);
  } else {
    assert(worker_num_ == region_map_.table_size());
    worker_id_t w_id;

    if (!region_map_.QueryWorkerWithMostOverlap(data_manager_, job, &w_id)) {
      dbg(DBG_ERROR, "ERROR: DynamicLoadBalancer: could not find worker for assigning job %lu.\n", job->job_id()); // NOLINT
      return false;
    }

    job->set_assigned_worker(worker_map_[w_id]);

    region_map_.TrackRegionCoverage(data_manager_, job, &w_id);
  }

  log.log_StopTimer();
  std::cout
    << "DYNAMIC: Picked worker: " << job->assigned_worker()->worker_id()
    << " for job: " << job->job_name()
    << " took: " << log.timer()
    << " for union set size of: " << job->union_set_p()->size()
    << std::endl;


  return true;
}


void DynamicLoadBalancer::NotifyJobAssignment(const JobEntry *job) {
  double time = log_.GetTime();

  if (job->job_type() != JOB_COMP) {
    return;
  }

  JobProfile *job_profile =
    new JobProfile(
        job->job_type(),
        job->job_name(),
        job->job_id(),
        job->parent_job_id(),
        job->assigned_worker()->worker_id(),
        job->sterile());

  job_profile->set_assign_time(time);
  job_profile->set_assigned(true);


  Vertex<JobEntry, job_id_t>* vertex;
  job_manager_->GetJobGraphVertex(job->job_id(), &vertex);

  typename Edge<JobEntry, job_id_t>::Iter iter;
  for (iter = vertex->incoming_edges()->begin(); iter != vertex->incoming_edges()->end(); ++iter) {
    JobEntry *j = iter->second->start_vertex()->entry();
    if (!j->done()) {
      job_profile->waiting_set_p()->insert(j->job_id());
    }
  }

  if (job_profile->waiting_set_p()->size() == 0) {
    job_profile->set_ready_time(time);
    job_profile->set_ready(true);
    boost::unique_lock<boost::recursive_mutex> lock_monitor(worker_monitor_mutex_);
    worker_monitor_.AddReadyJob(job->assigned_worker()->worker_id(), job->job_id());
  } else {
    boost::unique_lock<boost::recursive_mutex> lock_monitor(worker_monitor_mutex_);
    worker_monitor_.AddBlockedJob(job->assigned_worker()->worker_id(), job->job_id());
  }

  {
    boost::unique_lock<boost::recursive_mutex> lock_history(job_history_mutex_);
    job_history_[job->job_id()] = job_profile;
  }

  boost::unique_lock<boost::mutex> lock_stamp(stamp_mutex_);
  std::string jname = job->job_name();
  if (jname == "update_ghost_velocities" && (stamp_state_ == 0)) {
    std::cout << "STAMP: FIRST ASSIGNMENT LATENCY: " << log_.timer() << std::endl;
    stamp_state_ = 1;
  }
  if (jname == "projection_main" && (stamp_state_ == 1)) {
    std::cout << "STAMP: ALL ASSIGNMENT LATENCY: " << log_.timer() << std::endl;
    stamp_state_ = -1;
  }
}

void DynamicLoadBalancer::NotifyJobDone(const JobEntry *job) {
  double time = log_.GetTime();

  if (job->job_type() != JOB_COMP) {
    return;
  }

  if (++counter_ >  LB_UPDATE_RATE) {
    log_.log_WriteToFile(worker_monitor_.PrintStats());
    worker_monitor_.ResetWorkerTimers();
    counter_ = 0;
  }

  assert(job->done());
  done_jobs_.push_back(job->job_id());

  JobHistory::iterator it = job_history_.find(job->job_id());
  assert(it != job_history_.end());
  JobProfile *job_profile = it->second;

  job_profile->set_done_time(time);
  job_profile->set_done(true);
  job_profile->set_execute_duration(time - job_profile->ready_time());
  {
    boost::unique_lock<boost::recursive_mutex> lock_monitor(worker_monitor_mutex_);
    worker_monitor_.NotifyJobDone(job_profile->worker_id(), job_profile->job_id());
  }

  Vertex<JobEntry, job_id_t>* vertex;
  job_manager_->GetJobGraphVertex(job->job_id(), &vertex);

  typename Edge<JobEntry, job_id_t>::Iter iter;
  for (iter = vertex->outgoing_edges()->begin(); iter != vertex->outgoing_edges()->end(); ++iter) {
    JobEntry *j = iter->second->end_vertex()->entry();
    it = job_history_.find(j->job_id());
    if (it != job_history_.end()) {
      assert(j->assigned());
      JobProfile *jp = it->second;
      assert(jp->assigned());
      jp->add_log_entry(
          job_profile->worker_id(), job->job_id(), job->job_name(), time);
      jp->waiting_set_p()->remove(job->job_id());
      if (jp->waiting_set_p()->size() == 0) {
        jp->set_ready_time(time);
        jp->set_ready(true);
        {
          boost::unique_lock<boost::recursive_mutex> lock_monitor(worker_monitor_mutex_);
          worker_monitor_.AddReadyJob(jp->worker_id(), jp->job_id());
        }
      }
    }
  }

  // log_.log_WriteToFile(job_profile->Print());

  std::string jname = job->job_name();
  if (jname == "loop_iteration") {
    log_.log_StartTimer();
    stamp_state_ = 0;
  }
}


void DynamicLoadBalancer::NotifyRegisteredWorker(SchedulerWorker *worker) {
  boost::unique_lock<boost::recursive_mutex> update_lock(update_mutex_);
  boost::unique_lock<boost::recursive_mutex> worker_map_lock(worker_map_mutex_);

  {
    boost::unique_lock<boost::recursive_mutex> lock_monitor(worker_monitor_mutex_);
    worker_monitor_.AddWorker(worker->worker_id());
  }

  worker_id_t worker_id = worker->worker_id();
  WorkerMap::iterator iter = worker_map_.find(worker_id);
  if (iter == worker_map_.end()) {
    worker_map_[worker_id] = worker;
    worker_num_ = worker_map_.size();
    update_ = true;
    update_cond_.notify_all();
  } else {
    dbg(DBG_ERROR, "ERROR: DynamicLoadBalancer: worker with the same id %lu has already been registered.\n", // NOLINT
        worker_id);
  }
}

bool DynamicLoadBalancer::NotifyDownWorker(worker_id_t worker_id) {
  boost::unique_lock<boost::recursive_mutex> worker_map_lock(worker_map_mutex_);
  boost::unique_lock<boost::recursive_mutex> region_map_lock(region_map_mutex_);

  {
    boost::unique_lock<boost::recursive_mutex> lock_monitor(worker_monitor_mutex_);
    worker_monitor_.RemoveWorker(worker_id);
  }

  WorkerMap::iterator iter = worker_map_.find(worker_id);
  if (iter != worker_map_.end()) {
    worker_map_.erase(iter);
    worker_num_ = worker_map_.size();
    region_map_.NotifyDownWorker(worker_id);
    return true;
  } else {
    dbg(DBG_ERROR, "ERROR: DynamicLoadBalancer: worker with id %lu has NOT been registered.\n", // NOLINT
        worker_id);
    return false;
  }
}


void DynamicLoadBalancer::InitializeRegionMap() {
  boost::unique_lock<boost::recursive_mutex> worker_map_lock(worker_map_mutex_);
  boost::unique_lock<boost::recursive_mutex> region_map_lock(region_map_mutex_);

  std::vector<worker_id_t> worker_ids;
  WorkerMap::iterator iter = worker_map_.begin();
  for (; iter != worker_map_.end(); ++iter) {
    worker_ids.push_back(iter->first);
  }
  assert(worker_num_ > 0);
  global_region_ = data_manager_->global_bounding_region();

  region_map_.Initialize(worker_ids, global_region_);
  log_.log_WriteToFile(region_map_.Print());

  init_phase_ = false;
}

void DynamicLoadBalancer::UpdateRegionMap() {
  return;
  boost::unique_lock<boost::recursive_mutex> worker_map_lock(worker_map_mutex_);
  boost::unique_lock<boost::recursive_mutex> region_map_lock(region_map_mutex_);

  if (region_map_.BalanceRegions(1, 1)) {
    straggler_map_.ClearRecords();
    log_.log_WriteToFile(region_map_.Print());
  }
}

}  // namespace nimbus

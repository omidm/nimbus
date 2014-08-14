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

#define LB_UPDATE_RATE 50

namespace nimbus {

LoadBalancer::LoadBalancer() {
  Initialize();
}

LoadBalancer::LoadBalancer(ClusterMap* cluster_map) {
  Initialize();
  cluster_map_ = cluster_map;
}

void LoadBalancer::Initialize() {
  worker_num_ = 0;
  global_region_ = GeometricRegion(0, 0, 0, 0, 0, 0);
  update_ = false;
  init_phase_ = true;
  blame_counter_ = 0;
  cluster_map_ = NULL;
  job_manager_ = NULL;
  data_manager_ = NULL;
  log_.set_file_name("load_balancer_log");
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

void LoadBalancer::Run() {
  while (true) {
    boost::adopt_lock_t recursive;
    boost::unique_lock<boost::mutex> update_lock(update_mutex_, recursive);
    while (!update_) {
      update_cond_.wait(update_lock);
    }
    update_ = false;
    update_cond_.notify_all();

    boost::unique_lock<boost::mutex> worker_map_lock(worker_map_mutex_, recursive);
    boost::unique_lock<boost::mutex> region_map_lock(region_map_mutex_, recursive);

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

bool LoadBalancer::GetWorkerToAssignJob(
    JobEntry *job, SchedulerWorker*& worker) {


  Log log;
  log.StartTimer();

  boost::adopt_lock_t recursive;
  boost::unique_lock<boost::mutex> update_lock(update_mutex_, recursive);
  while (update_) {
    update_cond_.wait(update_lock);
  }

  boost::unique_lock<boost::mutex> worker_map_lock(worker_map_mutex_, recursive);
  boost::unique_lock<boost::mutex> region_map_lock(region_map_mutex_, recursive);

  if ((worker_num_ != region_map_.table_size()) ||
      (global_region_ != data_manager_->global_bounding_region())) {
    if (data_manager_->initialized_global_bounding_region()) {
      InitializeRegionMap();
    }
  }

  assert(worker_map_.size() > 0);
  assert(worker_num_ > 0);

  GeometricRegion region;
  bool got_region = job->GetRegion(data_manager_, &region);

  if (init_phase_ || !got_region) {
    worker = worker_map_.begin()->second;
  } else {
    assert(worker_num_ == region_map_.table_size());
    worker_id_t w_id;

    if (!region_map_.QueryWorkerWithMostOverlap(&region, &w_id)) {
      dbg(DBG_ERROR, "ERROR: LoadBalancer: could not find worker for assigning job %lu.\n", job->job_id()); // NOLINT
      return false;
    }

    worker = worker_map_[w_id];
  }

  log.StopTimer();
  std::cout
    << "Picked worker: " << worker->worker_id()
    << " for job: " << job->job_name()
    << " took: " << log.timer()
    << " for union set size of: " << job->union_set_p()->size()
    << std::endl;


  return true;
}


void LoadBalancer::NotifyJobAssignment(
    const JobEntry *job, const SchedulerWorker* worker) {
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
        worker->worker_id(),
        job->sterile());

  job_profile->set_assign_time(time);
  job_profile->set_assigned(true);


  Vertex<JobEntry, job_id_t>* vertex;
  job_manager_->job_graph_p()->GetVertex(job->job_id(), &vertex);

  typename Edge<JobEntry, job_id_t>::Iter iter;
  for (iter = vertex->incoming_edges()->begin(); iter != vertex->incoming_edges()->end(); ++iter) {
    JobEntry *j = iter->second->start_vertex()->entry();
    if (!j->done()) {
      job_profile->waiting_set_p()->insert(j->job_id());
    } else {
      /*
      JobHistory::iterator it = job_history_.find(j->job_id());
      if (it != job_history_.end()) {
        JobProfile *jp = it->second;
        job_profile->add_log_entry(
            jp->worker_id(), jp->job_id(), jp->done_time());
      } else {
        dbg(DBG_WARN, "WARNING: Load balancer, could not find done job in job history.");
        exit(-1);
      }
      */
    }
  }

  if (job_profile->waiting_set_p()->size() == 0) {
    job_profile->set_ready_time(time);
    job_profile->set_ready(true);
  }

  job_history_[job->job_id()] = job_profile;
}

void LoadBalancer::NotifyJobDone(const JobEntry *job) {
  double time = log_.GetTime();

  if (job->job_type() != JOB_COMP) {
    return;
  }

  assert(job->done());
  done_jobs_.push_back(job->job_id());

  JobHistory::iterator it = job_history_.find(job->job_id());
  assert(it != job_history_.end());
  JobProfile *job_profile = it->second;

  job_profile->set_done_time(time);
  job_profile->set_done(true);
  job_profile->set_execute_duration(time - job_profile->ready_time());

  Vertex<JobEntry, job_id_t>* vertex;
  job_manager_->job_graph_p()->GetVertex(job->job_id(), &vertex);

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
      }
    }
  }

  // log_.WriteToFile(job_profile->Print());

  worker_id_t blamed_worker_id;
  if (job_profile->FindBlamedWorker(&blamed_worker_id)) {
    boost::adopt_lock_t recursive;
    boost::unique_lock<boost::mutex> straggler_map_lock(straggler_map_mutex_, recursive);
    straggler_map_.AddRecord(job->assigned_worker(), blamed_worker_id);

    ++blame_map_[blamed_worker_id];

    blame_counter_++;
    if (blame_counter_ > LB_UPDATE_RATE) {
      blame_counter_ = 0;
      update_ = true;
      update_cond_.notify_all();
    }
  }
}


void LoadBalancer::NotifyRegisteredWorker(SchedulerWorker *worker) {
  boost::adopt_lock_t recursive;
  boost::unique_lock<boost::mutex> update_lock(update_mutex_, recursive);
  boost::unique_lock<boost::mutex> worker_map_lock(worker_map_mutex_, recursive);

  worker_id_t worker_id = worker->worker_id();
  WorkerMapIter iter = worker_map_.find(worker_id);
  if (iter == worker_map_.end()) {
    worker_map_[worker_id] = worker;
    worker_num_ = worker_map_.size();
    update_ = true;
    update_cond_.notify_all();
  } else {
    dbg(DBG_ERROR, "ERROR: LoadBalancer: worker with the same id %lu has already been registered.\n", // NOLINT
        worker_id);
  }
}

void LoadBalancer::InitializeRegionMap() {
  boost::adopt_lock_t recursive;
  boost::unique_lock<boost::mutex> worker_map_lock(worker_map_mutex_, recursive);
  boost::unique_lock<boost::mutex> region_map_lock(region_map_mutex_, recursive);

  std::vector<worker_id_t> worker_ids;
  WorkerMapIter iter = worker_map_.begin();
  for (; iter != worker_map_.end(); ++iter) {
    worker_ids.push_back(iter->first);
  }
  assert(worker_num_ > 0);
  global_region_ = data_manager_->global_bounding_region();

  region_map_.Initialize(worker_ids, global_region_);
  log_.WriteToFile(region_map_.Print());

  init_phase_ = false;
}

void LoadBalancer::UpdateRegionMap() {
  boost::adopt_lock_t recursive;
  boost::unique_lock<boost::mutex> worker_map_lock(worker_map_mutex_, recursive);
  boost::unique_lock<boost::mutex> region_map_lock(region_map_mutex_, recursive);
  boost::unique_lock<boost::mutex> straggler_map_lock(straggler_map_mutex_, recursive);

  worker_id_t fast, slow;
  if (straggler_map_.GetMostImbalanceWorkers(&fast, &slow)) {
    std::cout << "LOAD BALANCER: fast worker: " << fast
              << ", slow worker: " << slow << std::endl;
    if (region_map_.BalanceRegions(fast, slow)) {
      straggler_map_.ClearRecords();
      log_.WriteToFile(region_map_.Print());
    }
  }

  worker_id_t worst_worker = 0;
  size_t count = 0;
  std::map<worker_id_t, size_t>::iterator iter = blame_map_.begin();
  for (; iter != blame_map_.end(); ++iter) {
    if (iter->second > count) {
      count = iter->second;
      worst_worker = iter->first;
    }
  }
  std::cout << "LOAD BALANCER: worst worker: " << worst_worker << std::endl;
  blame_map_.clear();
}

}  // namespace nimbus

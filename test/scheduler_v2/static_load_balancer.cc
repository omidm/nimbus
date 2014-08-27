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


#include "./static_load_balancer.h"

#define WEIGHT_NUM 8
#define WEIGHT_X {1, 1, 1, 1, 1, 1, 1, 1}
#define WEIGHT_Y {1, 1, 1, 1, 1, 1, 1, 1}
#define WEIGHT_Z {1, 1, 1, 1, 1, 1, 1, 1}


namespace nimbus {

StaticLoadBalancer::StaticLoadBalancer() {
  stamp_state_ = -1;
  initialized_domains_ = false;
}

StaticLoadBalancer::~StaticLoadBalancer() {
}

void StaticLoadBalancer::Run() {
}

size_t StaticLoadBalancer::AssignReadyJobs() {
  JobEntryList list;
  job_manager_->GetJobsReadyToAssign(&list, max_job_to_assign_);

  JobEntryList::iterator iter;
  for (iter = list.begin(); iter != list.end(); ++iter) {
    JobEntry* job = *iter;
    if (!SetWorkerToAssignJob(job)) {
      dbg(DBG_ERROR, "ERROR: StaticLoadBalancer: could not get worker to assign job %lu.\n", job->job_id()); // NOLINT
      exit(-1);
    }
  }

  job_assigner_->AssignJobs(list);

  return list.size();
}

bool StaticLoadBalancer::SetWorkerToAssignJob(JobEntry* job) {
  log_.StartTimer();

  StaticLoadBalancer::UpdateWorkerDomains();
  assert(worker_num_ == worker_map_.size());
  assert(worker_num_ == worker_domains_.size());

  // initilaize worker ranks.
  WorkerRank worker_rank;
  WorkerMap::iterator wmit = worker_map_.begin();
  for (; wmit != worker_map_.end(); ++wmit) {
    worker_rank[wmit->first] = 0;
  }

  IDSet<logical_data_id_t>::IDSetIter iter;
  for (iter = job->union_set_p()->begin(); iter != job->union_set_p()->end(); ++iter) {
    const LogicalDataObject* ldo;
    ldo = data_manager_->FindLogicalObject(*iter);
    WorkerDomains::iterator wdit = worker_domains_.begin();
    for (; wdit != worker_domains_.end(); ++wdit) {
      if (worker_domains_[wdit->first].Intersects(ldo->region())) {
        ++worker_rank[wdit->first];
      }
    }
  }

  // find the worker that wins the poll.
  WorkerRank::iterator writ = worker_rank.begin();
  assert(writ != worker_rank.end());
  worker_id_t w_id = writ->first;
  size_t count = worker_rank[writ->first];
  for (; writ != worker_rank.end(); ++writ) {
    if (count < worker_rank[writ->first]) {
      count = worker_rank[writ->first];
      w_id = writ->first;
    }
  }

  log_.StopTimer();
  std::cout
    << "OMID: Picked worker: " << w_id
    << " for job: " << job->job_name()
    << " took: " << log_.timer()
    << " for union set size of: " << job->union_set_p()->size() << std::endl;

  wmit = worker_map_.find(w_id);
  if (wmit == worker_map_.end()) {
    return false;
  } else {
    job->set_assigned_worker(wmit->second);
    return true;
  }
}

void StaticLoadBalancer::UpdateWorkerDomains() {
  GeometricRegion global_bounding_region =
    data_manager_->global_bounding_region();

  if ((!initialized_domains_) ||
      (worker_num_ != worker_domains_.size()) ||
      (global_bounding_region_ != global_bounding_region)) {
    global_bounding_region_ = global_bounding_region;
    worker_domains_.clear();

    size_t num_x, num_y, num_z;
    SplitDimensions(worker_num_, &num_x, &num_y, &num_z);

    GenerateDomains(num_x, num_y, num_z,
                    global_bounding_region_);

    initialized_domains_ = true;
  }
}

void StaticLoadBalancer::GenerateDomains(size_t num_x,
                                         size_t num_y,
                                         size_t num_z,
                                         GeometricRegion gbr) {
  size_t weight_x[WEIGHT_NUM] = WEIGHT_X;
  size_t weight_y[WEIGHT_NUM] = WEIGHT_Y;
  size_t weight_z[WEIGHT_NUM] = WEIGHT_Z;

  std::vector<int_dimension_t> width_x;
  size_t weight_sum_x = 0;
  for (size_t i = 0; i < num_x; ++i) {
    weight_sum_x += weight_x[i];
  }
  for (size_t i = 0; i < num_x; ++i) {
    width_x.push_back(gbr.dx() * weight_x[i] / weight_sum_x);
  }
  std::vector<int_dimension_t> marker_x;
  marker_x.push_back(gbr.x());
  for (size_t i = 0; i < num_x; ++i) {
    marker_x.push_back(marker_x[i] + width_x[i]);
  }


  std::vector<int_dimension_t> width_y;
  size_t weight_sum_y = 0;
  for (size_t i = 0; i < num_y; ++i) {
    weight_sum_y += weight_y[i];
  }
  for (size_t i = 0; i < num_y; ++i) {
    width_y.push_back(gbr.dy() * weight_y[i] / weight_sum_y);
  }
  std::vector<int_dimension_t> marker_y;
  marker_y.push_back(gbr.y());
  for (size_t i = 0; i < num_y; ++i) {
    marker_y.push_back(marker_y[i] + width_y[i]);
  }

  std::vector<int_dimension_t> width_z;
  size_t weight_sum_z = 0;
  for (size_t i = 0; i < num_z; ++i) {
    weight_sum_z += weight_z[i];
  }
  for (size_t i = 0; i < num_z; ++i) {
    width_z.push_back(gbr.dz() * weight_z[i] / weight_sum_z);
  }
  std::vector<int_dimension_t> marker_z;
  marker_z.push_back(gbr.z());
  for (size_t i = 0; i < num_z; ++i) {
    marker_z.push_back(marker_z[i] + width_z[i]);
  }

  WorkerMap::iterator it = worker_map_.begin();
  assert((num_x * num_y * num_z) == worker_num_);
  worker_domains_.clear();
  for (size_t i = 0; i < num_x; ++i) {
    for (size_t j = 0; j < num_y; ++j) {
      for (size_t k = 0; k < num_z; ++k) {
        worker_domains_[it->first] =
            GeometricRegion(
              marker_x[i],
              marker_y[j],
              marker_z[k],
              marker_x[i + 1] - marker_x[i],
              marker_y[j + 1] - marker_y[j],
              marker_z[k + 1] - marker_z[k]);
        ++it;
      }
    }
  }
}


void StaticLoadBalancer::SplitDimensions(size_t worker_num,
    size_t *num_x, size_t *num_y, size_t *num_z) {
  switch (worker_num) {
    case 1 :
      *num_x = 1;
      *num_y = 1;
      *num_z = 1;
      break;
    case 2 :
      *num_x = 2;
      *num_y = 1;
      *num_z = 1;
      break;
    case 3 :
      *num_x = 3;
      *num_y = 1;
      *num_z = 1;
      break;
    case 4 :
      *num_x = 4;
      *num_y = 1;
      *num_z = 1;
      break;
    case 5 :
      *num_x = 5;
      *num_y = 1;
      *num_z = 1;
      break;
    case 6 :
      *num_x = 6;
      *num_y = 1;
      *num_z = 1;
      break;
    case 7 :
      *num_x = 7;
      *num_y = 1;
      *num_z = 1;
      break;
    case 8 :
      *num_x = 8;
      *num_y = 1;
      *num_z = 1;
      break;
    default:
      dbg(DBG_ERROR, "ERROR: Do not know how to split!");
      exit(-1);
  }
}

void StaticLoadBalancer::NotifyJobAssignment(const JobEntry *job) {
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

void StaticLoadBalancer::NotifyJobDone(const JobEntry *job) {
  std::string jname = job->job_name();
  if (jname == "loop_iteration") {
    log_.log_StartTimer();
    stamp_state_ = 0;
  }
}


void StaticLoadBalancer::NotifyRegisteredWorker(SchedulerWorker *worker) {
  worker_id_t worker_id = worker->worker_id();
  WorkerMapIter iter = worker_map_.find(worker_id);
  if (iter == worker_map_.end()) {
    worker_map_[worker_id] = worker;
    worker_num_ = worker_map_.size();
  } else {
    dbg(DBG_ERROR, "ERROR: LoadBalancer: worker with the same id %lu has already been registered.\n", // NOLINT
        worker_id);
  }
}

}  // namespace nimbus

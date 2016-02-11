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
#include <stdlib.h>

#define LB_UPDATE_THRESHOLD (double)(0.2) // NOLINT

namespace nimbus {

DynamicLoadBalancer::DynamicLoadBalancer() {
  worker_num_ = 0;
  init_phase_ = true;
  last_global_run_time_ = 0;
  log_.set_file_name("load_balancer_log");
}

DynamicLoadBalancer::~DynamicLoadBalancer() {
}

void DynamicLoadBalancer::Run() {
}


void  DynamicLoadBalancer::SplitDimensions(size_t worker_num) {
  switch (worker_num) {
    case 1 :
      split_[0] = 1;
      split_[1] = 1;
      split_[2] = 1;
      break;
    case 2 :
      split_[0] = 1;
      split_[1] = 2;
      split_[2] = 1;
      break;
    case 3 :
      split_[0] = 3;
      split_[1] = 1;
      split_[2] = 1;
      break;
    case 4 :
      split_[0] = 2;
      split_[1] = 2;
      split_[2] = 1;
      break;
    case 5 :
      split_[0] = 5;
      split_[1] = 1;
      split_[2] = 1;
      break;
    case 6 :
      split_[0] = 2;
      split_[1] = 3;
      split_[2] = 1;
      break;
    case 7 :
      split_[0] = 7;
      split_[1] = 1;
      split_[2] = 1;
      break;
    case 8 :
      split_[0] = 2;
      split_[1] = 2;
      split_[2] = 2;
      break;
    case 64 :
      split_[0] = 4;
      split_[1] = 4;
      split_[2] = 4;
      break;
    case 100 :
      split_[0] = 5;
      split_[1] = 5;
      split_[2] = 4;
      break;
    default:
      dbg(DBG_ERROR, "ERROR: Do not know how to split!");
      assert(false);
  }
}


bool DynamicLoadBalancer::SetWorkerToAssignJob(JobEntry *job) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  if (job->job_name() == "some-random-name-loop_frame") {
    UpdateRegionMap();
  }

  if (job->job_type() == JOB_CMPX) {
    ComplexJobEntry *xj = reinterpret_cast<ComplexJobEntry*>(job);
    char buff[LOG_MAX_BUFF_SIZE];
    snprintf(buff, sizeof(buff), "DYNAMIC: %10.9lf complex job for %s.",
        Log::GetRawTime(), xj->template_entry()->template_name().c_str());
    log_.log_WriteToOutputStream(std::string(buff));
    return true;
  }

  Log log(Log::NO_FILE);
  log.log_StartTimer();


  if (enforced_global_region_) {
    if ((worker_num_ != region_map_.table_size())) {
      BuildRegionMap();
    }
  } else {
    if ((worker_num_ != region_map_.table_size()) ||
        (global_region_ != data_manager_->global_bounding_region())) {
      if (data_manager_->initialized_global_bounding_region()) {
        global_region_ = data_manager_->global_bounding_region();
        BuildRegionMap();
      }
    }
  }

  assert(worker_num_ > 0);
  assert(worker_num_ == worker_map_.size());
  if (!init_phase_) {
    assert(worker_num_ == region_map_.table_size());
  }

  worker_id_t w_id;
  if (init_phase_ || (job->union_set_p()->size() == 0) || (!job->sterile())) {
    job->set_assigned_worker(worker_map_.begin()->second);
    w_id = worker_map_.begin()->second->worker_id();
  } else {
    assert(worker_num_ == region_map_.table_size());

    if (!region_map_.QueryWorkerWithMostOverlap(data_manager_, job, &w_id)) {
      dbg(DBG_ERROR, "ERROR: DynamicLoadBalancer: could not find worker for assigning job %lu.\n", job->job_id()); // NOLINT
      return false;
    }

    job->set_assigned_worker(worker_map_[w_id]);

    region_map_.TrackRegionCoverage(data_manager_, job, &w_id);
  }

  log.log_StopTimer();
  char buff[LOG_MAX_BUFF_SIZE];
  snprintf(buff, sizeof(buff), "DYNAMIC: %10.9lf Picked worker %2.0u for %s.",
      Log::GetRawTime(), w_id, job->job_name().c_str());
  log_.log_WriteToOutputStream(std::string(buff));

  return true;
}

void DynamicLoadBalancer::NotifyRegisteredWorker(SchedulerWorker *worker) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

//  {
//    worker_monitor_.AddWorker(worker->worker_id());
//  }

  worker_id_t worker_id = worker->worker_id();
  WorkerMap::iterator iter = worker_map_.find(worker_id);
  if (iter == worker_map_.end()) {
    worker_map_[worker_id] = worker;
    worker_num_ = worker_map_.size();
  } else {
    dbg(DBG_ERROR, "ERROR: DynamicLoadBalancer: worker with the same id %lu has already been registered.\n", // NOLINT
        worker_id);
    assert(false);
  }
}

bool DynamicLoadBalancer::NotifyDownWorker(worker_id_t worker_id) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

//  {
//    worker_monitor_.RemoveWorker(worker_id);
//  }

  WorkerMap::iterator iter = worker_map_.find(worker_id);
  if (iter != worker_map_.end()) {
    worker_map_.erase(iter);
    worker_num_ = worker_map_.size();

    // Update load balancing id.
    ++load_balancing_id_;

    region_map_.NotifyDownWorker(worker_id);
    return true;
  } else {
    dbg(DBG_ERROR, "ERROR: DynamicLoadBalancer: worker with id %lu has NOT been registered.\n", // NOLINT
        worker_id);
    return false;
  }
}

bool DynamicLoadBalancer::BalanceLoad(counter_t query_id) {
  QueryMap::iterator iter = queries_.find(query_id);
  assert(iter != queries_.end());

  StatQuery *st = iter->second;

  if (worker_map_.size() < 2) {
    std::cout << "There is not at least two workers to load balance.\n";
    return false;
  }

  WorkerMap::iterator w1_iter = worker_map_.begin();
  WorkerMap::iterator w2_iter = worker_map_.begin();
  ++w2_iter;
  for (; w2_iter != worker_map_.end();) {
    worker_id_t w1 = w1_iter->first;
    worker_id_t w2 = w2_iter->first;

    int64_t r1 = 0, r2 = 0;
    {
      StatQuery::iterator it = st->find(w1);
      assert(it != st->end());
      r1 = it->second->run_time_;
    }
    {
      StatQuery::iterator it = st->find(w2);
      assert(it != st->end());
      r2 = it->second->run_time_;
    }
    double load_imbalance = ((double)(std::abs(r1 - r2))) / ((double)(std::max(r1, r2))); // NOLINT
    if (load_imbalance >= LB_UPDATE_THRESHOLD) {
      if (r1 > r2) {
        if (InBlackList(w2, w1)) {
          std::cout << "\n****** LB IN BLACK LIST: "
                    <<  w2 << " CANNOT grow into " << w1
                    << " LIB factor :" << load_imbalance
                    << std::endl;
          ++w1_iter;
          ++w2_iter;
          continue;
        }
        if (CheckPossibleFluctuation(w2, w1, st)) {
          std::cout << "\n****** LB BLACK LISTED: "
                    <<  w2 << " CANNOT grow into " << w1
                    << " LIB factor :" << load_imbalance
                    << std::endl;
          ++w1_iter;
          ++w2_iter;
          continue;
        }
        if (region_map_.BalanceRegions(w2, w1)) {
          log_.log_WriteToFile(region_map_.Print());
          // Update load balancing id.
          ++load_balancing_id_;

          // Remember last exchange stats
          last_exchange_ = std::make_pair(w2, w1);
          last_global_run_time_ = GetGlobalRunTime(st);

          std::cout << "\n****** LB : "
                    <<  w2 << " grows into " << w1
                    << " LIB factor :" << load_imbalance
                    << std::endl;
        }
      } else if (r2 > r1) {
        if (InBlackList(w1, w2)) {
          std::cout << "\n****** LB IN BLACK LIST: "
                    <<  w1 << " CANNOT grow into " << w2
                    << " LIB factor :" << load_imbalance
                    << std::endl;
          ++w1_iter;
          ++w2_iter;
          continue;
        }
        if (CheckPossibleFluctuation(w1, w2, st)) {
          std::cout << "\n****** LB BLACK LISTED: "
                    <<  w1 << " CANNOT grow into " << w2
                    << " LIB factor :" << load_imbalance
                    << std::endl;
          ++w1_iter;
          ++w2_iter;
          continue;
        }
        if (region_map_.BalanceRegions(w1, w2)) {
          log_.log_WriteToFile(region_map_.Print());
          // Update load balancing id.
          ++load_balancing_id_;

          // Remember last exchange stats
          last_exchange_ = std::make_pair(w1, w2);
          last_global_run_time_ = GetGlobalRunTime(st);

          std::cout << "\n****** LB : "
                    <<  w1 << " grows into " << w2
                    << " LIB factor :" << load_imbalance
                    << std::endl;
        }
      }
    }

    ++w1_iter;
    ++w2_iter;
  }

  return true;
}

void DynamicLoadBalancer::BuildRegionMap() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  // Update load balancing id.
  ++load_balancing_id_;

  if (!enforced_split_) {
    SplitDimensions(worker_num_);
  }

  if (!enforced_sub_split_) {
    sub_split_[0] = 1;
    sub_split_[1] = 1;
    sub_split_[2] = 1;
  }

  std::vector<worker_id_t> worker_ids;
  WorkerMap::iterator iter = worker_map_.begin();
  for (; iter != worker_map_.end(); ++iter) {
    worker_ids.push_back(iter->first);
  }
  assert(worker_num_ > 0);

  region_map_.Initialize(worker_ids, split_, sub_split_, global_region_);
  log_.log_WriteToFile(region_map_.Print());

  init_phase_ = false;
}

void DynamicLoadBalancer::UpdateRegionMap() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  if (region_map_.BalanceRegions(1, 2)) {
    log_.log_WriteToFile(region_map_.Print());
    // Update load balancing id.
    ++load_balancing_id_;
    std::cout << "****** LBLBLBLBLB *******\n";
  } else {
    // assert(false);
  }
}

int64_t DynamicLoadBalancer::GetGlobalRunTime(const StatQuery *st) {
  int64_t global_run_time = 0;
  WorkerMap::iterator iter = worker_map_.begin();
  for (; iter != worker_map_.end(); ++iter) {
    StatQuery::const_iterator it = st->find(iter->first);
    assert(it != st->end());
    global_run_time += it->second->run_time_;
  }

  return global_run_time;
}

bool DynamicLoadBalancer::InBlackList(worker_id_t w2, worker_id_t w1) {
  BlackList::iterator iter = black_list_.begin();
  for (; iter != black_list_.end(); ++iter) {
    if ((iter->first == w2) && (iter->second == w1)) {
      return true;
    }
  }

  return false;
}

bool DynamicLoadBalancer::CheckPossibleFluctuation(worker_id_t w2,
                                                   worker_id_t w1,
                                                   const StatQuery *st) {
  if ((last_exchange_.first == w1) && (last_exchange_.second == w2)) {
    if (last_global_run_time_ < GetGlobalRunTime(st)) {
      black_list_.push_back(std::make_pair(w2, w1));
      return true;
    }
  }

  return false;
}

}  // namespace nimbus

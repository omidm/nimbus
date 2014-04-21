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
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "./scheduler_1d.h"

Scheduler1D::Scheduler1D(unsigned int p)
: Scheduler(p) {
  initialized_domains_ = false;
}

bool Scheduler1D::GetWorkerToAssignJob(JobEntry* job, SchedulerWorker*& worker) {
  size_t worker_num = server_->worker_num();
  GeometricRegion global_bounding_region =
    data_manager_->global_bounding_region();

  if (!initialized_domains_ ||
      worker_num_ != worker_num ||
      global_bounding_region_ != global_bounding_region) {
    global_bounding_region_ = global_bounding_region;
    worker_num_ = worker_num;
    worker_domains_.clear();
    int_dimension_t delta = (int_dimension_t)
      ((float) (global_bounding_region_.dy()) / (float) (worker_num_)); // NOLINT
    for (size_t i = 0; i < worker_num_; ++i) {
      worker_domains_.push_back(GeometricRegion(
            global_bounding_region_.x(),
            global_bounding_region_.y() + i * delta,
            global_bounding_region_.z(),
            global_bounding_region_.dx(),
            delta,
            global_bounding_region_.dz()));
    }
    initialized_domains_ = true;
  }

  std::vector<int> workers_rank(worker_num, 0);

  IDSet<logical_data_id_t> union_set = job->union_set();
  IDSet<logical_data_id_t>::IDSetIter iter;
  for (iter = union_set.begin(); iter != union_set.end(); ++iter) {
    const LogicalDataObject* ldo;
    ldo = data_manager_->FindLogicalObject(*iter);
    for (size_t i = 0; i < worker_num; ++i) {
      if (worker_domains_[i].Intersects(ldo->region())) {
        ++workers_rank[i];
      }
    }
  }

  // find the worker that wins the poll.
  worker_id_t w_id = 1;
  int count = workers_rank[0];
  for (size_t i = 1; i < worker_num; ++i) {
    if (count < workers_rank[i]) {
      count = workers_rank[i];
      w_id = i + 1;
    }
  }

  std::cout << "Picked worker: " << w_id << " for job: " << job->job_name() << std::endl;
  return server_->GetSchedulerWorkerById(worker, w_id);
}


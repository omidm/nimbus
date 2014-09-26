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
  * Worker profile used in load balancer.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "scheduler/worker_profile.h"

using namespace nimbus; // NOLINT

WorkerProfile::WorkerProfile(worker_id_t worker_id) {
  worker_id_ = worker_id;
  ready_jobs_ = new IDSet<job_id_t>();
  blocked_jobs_ = new IDSet<job_id_t>();
  ResetTimers();
}

WorkerProfile::~WorkerProfile() {
  delete ready_jobs_;
  delete blocked_jobs_;
}


worker_id_t WorkerProfile::worker_id() {
  return worker_id_;
}

double WorkerProfile::idle_time() {
  boost::unique_lock<boost::mutex> lock(mutex_);
  return idle_timer_.timer();
}

double WorkerProfile::active_time() {
  boost::unique_lock<boost::mutex> lock(mutex_);
  return active_timer_.timer();
}

double WorkerProfile::blocked_time() {
  boost::unique_lock<boost::mutex> lock(mutex_);
  return blocked_timer_.timer();
}

IDSet<job_id_t>* WorkerProfile::ready_jobs() {
  return ready_jobs_;
}

IDSet<job_id_t>* WorkerProfile::blocked_jobs() {
  return blocked_jobs_;
}


void WorkerProfile::ResetTimers() {
  boost::unique_lock<boost::mutex> lock(mutex_);
  if ((ready_jobs_->size() == 0) &&
      (blocked_jobs_->size() == 0)) {
    idle_timer_.Start();
    active_timer_.Reset();
    blocked_timer_.Reset();
  } else if (ready_jobs_->size() == 0) {
    idle_timer_.Reset();
    active_timer_.Reset();
    blocked_timer_.Start();
  } else {
    idle_timer_.Reset();
    active_timer_.Start();
    blocked_timer_.Reset();
  }
}

void WorkerProfile::AddReadyJob(job_id_t job_id) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  ready_jobs_->insert(job_id);
  blocked_jobs_->remove(job_id);
  idle_timer_.Stop();
  blocked_timer_.Stop();
  active_timer_.Resume();
}

void WorkerProfile::AddBlockedJob(job_id_t job_id) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  blocked_jobs_->insert(job_id);
  idle_timer_.Stop();
  if (ready_jobs_->size() == 0) {
    blocked_timer_.Resume();
  }
}

void WorkerProfile::NotifyJobDone(job_id_t job_id) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  assert(ready_jobs_->contains(job_id));
  ready_jobs_->remove(job_id);
  blocked_jobs_->remove(job_id);
  if ((ready_jobs_->size() == 0) &&
      (blocked_jobs_->size() == 0)) {
    idle_timer_.Resume();
    active_timer_.Stop();
    blocked_timer_.Stop();
  } else if (ready_jobs_->size() == 0) {
    active_timer_.Stop();
    blocked_timer_.Resume();
  }
}

std::string WorkerProfile::PrintStatus() {
  boost::unique_lock<boost::mutex> lock(mutex_);
  std::string rval;
  rval += "\n+++++++ Worker Status +++++++\n";

  {
    std::ostringstream ss;
    rval += "worker_id: ";
    ss << worker_id_;
    rval += ss.str();
    ss.str(std::string());
    rval += "\n";
  }

  {
    std::ostringstream ss;
    rval += "idle_time: ";
    ss << idle_time();
    rval += ss.str();
    ss.str(std::string());
    rval += "\n";
  }

  {
    std::ostringstream ss;
    rval += "blocked_time: ";
    ss << blocked_time();
    rval += ss.str();
    ss.str(std::string());
    rval += "\n";
  }

  {
    std::ostringstream ss;
    rval += "active_time: ";
    ss << active_time();
    rval += ss.str();
    ss.str(std::string());
    rval += "\n";
  }

  rval += "\n+++++++++++++++++++++++++++++\n";
  return rval;
}



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
}

WorkerProfile::~WorkerProfile() {
  delete ready_jobs_;
  delete blocked_jobs_;
}


worker_id_t WorkerProfile::worker_id() {
  return worker_id_;
}

double WorkerProfile::idle_time() {
  return idle_timer_.timer();
}

double WorkerProfile::active_time() {
  return active_timer_.timer();
}

double WorkerProfile::blocked_time() {
  return blocked_timer_.timer();
}

IDSet<job_id_t>* WorkerProfile::ready_jobs() {
  return ready_jobs_;
}

IDSet<job_id_t>* WorkerProfile::blocked_jobs() {
  return blocked_jobs_;
}


void WorkerProfile::ResetTimers() {
  idle_timer_.Reset();
  active_timer_.Reset();
  blocked_timer_.Reset();
}

void WorkerProfile::AddReadyJob(job_id_t job_id) {
}

void WorkerProfile::AddBlockedJob(job_id_t job_id) {
}

void WorkerProfile::NotifyJobDone(job_id_t job_id) {
}

std::string WorkerProfile::PrintStatus() {
  std::string rval;
  rval += "\n+++++++ Worker Status +++++++\n";

  std::ostringstream ss;
  rval += "worker_id: ";
  ss << worker_id_;
  rval += ss.str();
  ss.str(std::string());
  rval += "\n";

  rval += "\n+++++++++++++++++++++++++++++\n";
  return rval;
}



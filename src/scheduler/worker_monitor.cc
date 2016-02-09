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
  * Worker monitor keeps track of the workers activity. It captures their idle,
  * blocked, and active time and detects the stragglers and fast workers for
  * load balancing strategy.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#include "src/scheduler/worker_monitor.h"

namespace nimbus {

WorkerMonitor::WorkerMonitor() {
}

WorkerMonitor::~WorkerMonitor() {
  MapIter iter = map_.begin();
  for (; iter != map_.end(); ++iter) {
    delete iter->second;
  }
}

bool WorkerMonitor::AddWorker(worker_id_t worker_id) {
  MapIter iter = map_.find(worker_id);
  if (iter != map_.end()) {
    dbg(DBG_ERROR, "ERROR: WorkerMonitor: worker id %lu already exist!", worker_id);
    assert(false);
    return false;
  }

  map_[worker_id] = new WorkerProfile(worker_id);
  return true;
}

bool WorkerMonitor::RemoveWorker(worker_id_t worker_id) {
  MapIter iter = map_.find(worker_id);
  if (iter == map_.end()) {
    dbg(DBG_ERROR, "ERROR: WorkerMonitor: worker id %lu does not exist!", worker_id);
    assert(false);
    return false;
  }

  map_[worker_id] = new WorkerProfile(worker_id);
  delete iter->second;
  map_.erase(iter);
  return true;
}

bool WorkerMonitor::AddReadyJob(worker_id_t worker_id, job_id_t job_id) {
  MapIter iter = map_.find(worker_id);
  if (iter == map_.end()) {
    dbg(DBG_ERROR, "ERROR: WorkerMonitor: worker id %lu does not exist!", worker_id);
    assert(false);
    return false;
  }

  iter->second->AddReadyJob(job_id);
  return true;
}

bool WorkerMonitor::AddBlockedJob(worker_id_t worker_id, job_id_t job_id) {
  MapIter iter = map_.find(worker_id);
  if (iter == map_.end()) {
    dbg(DBG_ERROR, "ERROR: WorkerMonitor: worker id %lu does not exist!", worker_id);
    assert(false);
    return false;
  }

  iter->second->AddBlockedJob(job_id);
  return true;
}

bool WorkerMonitor::NotifyJobDone(worker_id_t worker_id, job_id_t job_id) {
  MapIter iter = map_.find(worker_id);
  if (iter == map_.end()) {
    dbg(DBG_ERROR, "ERROR: WorkerMonitor: worker id %lu does not exist!", worker_id);
    assert(false);
    return false;
  }

  iter->second->NotifyJobDone(job_id);
  return true;
}

void WorkerMonitor::ResetWorkerTimers() {
  MapIter iter = map_.begin();
  for (; iter != map_.end(); ++iter) {
    iter->second->ResetTimers();
  }
}

size_t WorkerMonitor::GetStragglers(std::list<worker_id_t> *list, size_t percentile) {
  assert(percentile <= 100);
  list->clear();
  return 0;
}

size_t WorkerMonitor::GetFastWorkers(std::list<worker_id_t> *list, size_t percentile) {
  assert(percentile <= 100);
  list->clear();
  return 0;
}

std::string WorkerMonitor::PrintStats() {
  std::string rval;
  rval += "\n+++++++ Worker Monitor Stats +++++++\n";
  MapIter iter = map_.begin();
  for (; iter != map_.end(); ++iter) {
    rval += iter->second->PrintStats();
  }

  rval += "\n++++++++++++++++++++++++++++++++++++\n";
  return rval;
}

}  // namespace nimbus

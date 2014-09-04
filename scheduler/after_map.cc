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
  * AfterMap data structure is meant to keep track of workers that needs to be
  * informes about successful finish of a specific job. In other words, if a
  * worker has to run a job say A, then worker has to be informed about all job
  * dones in the before set of A. So, each job done has to be communicated
  * with all the workers that are running it's after set.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "scheduler/after_map.h"

using namespace nimbus; // NOLINT

AfterMap::AfterMap() {
}

AfterMap::AfterMap(const AfterMap& other) {
  content_ = other.content_;
}

AfterMap::~AfterMap() {
  Iter iter = content_.begin();
  for (; iter != content_.end(); ++iter) {
    delete iter->second;
  }
}

AfterMap::Map AfterMap::content() const {
  return content_;
}

const AfterMap::Map* AfterMap::content_p() const {
  return &content_;
}

void AfterMap::set_content(const AfterMap::Map& content) {
  content_= content;
}

bool AfterMap::AddEntries(JobEntry* job) {
  assert(job->assigned());
  SchedulerWorker *worker = job->assigned_worker();
  assert(worker);

  IDSet<job_id_t>::IDSetIter iter = job->before_set_p()->begin();
  for (; iter != job->before_set_p()->end(); ++iter) {
    AddEntry(*iter, worker);
  }

  return true;
}

bool AfterMap::AddEntry(job_id_t job_id, SchedulerWorker* worker) {
  Iter iter = content_.find(job_id);
  if (iter == content_.end()) {
    Pool *pool = new Pool();
    pool->insert(worker);
    content_[job_id] = pool;
  } else {
    iter->second->insert(worker);
  }

  return true;
}

bool AfterMap::GetWorkersWaitingOn(job_id_t job_id,
                                   std::list<SchedulerWorker*> *list) {
  list->clear();

  Iter iter = content_.find(job_id);
  if (iter != content_.end()) {
    Pool *pool = iter->second;
    Pool::iterator it = pool->begin();
    for (; it != pool->end(); ++it) {
      list->push_back(*it);
    }
  }

  return true;
}

bool AfterMap::RemoveJobRecords(job_id_t job_id) {
  content_.erase(job_id);
  return true;
}

bool AfterMap::RemoveWorkerRecords(SchedulerWorker *worker) {
  Iter iter = content_.begin();
  for (; iter != content_.end();) {
    iter->second->erase(worker);
    if (iter->second->size() == 0) {
      delete iter->second;
      content_.erase(iter++);
    } else {
      ++iter;
    }
  }

  return true;
}

AfterMap& AfterMap::operator=(const AfterMap& right) {
  content_ = right.content_;
  return (*this);
}



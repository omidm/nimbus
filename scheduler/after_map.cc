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
  map_ = new Map();
  late_map_ = new Map();
}

AfterMap::~AfterMap() {
  DestroyMap(map_);
  DestroyMap(late_map_);
}

const AfterMap::Map* AfterMap::map() const {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  return map_;
}

const AfterMap::Map* AfterMap::late_map() const {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  return late_map_;
}

bool AfterMap::AddEntries(JobEntry* job) {
  // Lock will be aciquired.
  assert(job->assigned());
  SchedulerWorker *worker = job->assigned_worker();
  assert(worker);

  IDSet<job_id_t>::IDSetIter iter = job->before_set_p()->begin();
  for (; iter != job->before_set_p()->end(); ++iter) {
    AddEntry(*iter, worker);
  }

  return true;
}

void AfterMap::AddEntryToMap(Map *map, job_id_t job_id, SchedulerWorker *worker) {
  // Already Lock is aciquired.
  Iter iter = map->find(job_id);
  if (iter == map->end()) {
    Pool *pool = new Pool();
    pool->insert(worker);
    map->operator[](job_id) = pool;
  } else {
    iter->second->insert(worker);
  }
}

bool AfterMap::EntryIsInMap(Map *map, job_id_t job_id, SchedulerWorker *worker) {
  // Already Lock is aciquired.
  Iter iter = map->find(job_id);
  if (iter == map->end()) {
    return false;
  } else {
    Pool *pool = iter->second;
    return (pool->find(worker) != pool->end());
  }
}


bool AfterMap::AddEntry(job_id_t job_id, SchedulerWorker* worker) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  if (!JobIsDone(job_id)) {
    AddEntryToMap(map_, job_id, worker);
  } else {
    if (!EntryIsInMap(map_, job_id, worker)) {
      AddEntryToMap(map_, job_id, worker);
      AddEntryToMap(late_map_, job_id, worker);
    }
  }

  return true;
}

bool AfterMap::GetWorkersWaitingOnJob(job_id_t job_id,
                                      std::list<SchedulerWorker*> *list) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  list->clear();

  Iter iter = map_->find(job_id);
  if (iter != map_->end()) {
    Pool *pool = iter->second;
    Pool::iterator it = pool->begin();
    for (; it != pool->end(); ++it) {
      list->push_back(*it);
    }
  }

  return true;
}


bool AfterMap::NotifyJobDone(job_id_t job_id) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  job_done_pool_.insert(job_id);
  return true;
}

bool AfterMap::PullLateMap(Map*& late_map) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  if (late_map_->size() == 0) {
    return false;
  } else {
    late_map = late_map_;
    late_map_ = new Map();
    return true;
  }
}

void AfterMap::DestroyMap(Map *map) {
  Iter iter = map->begin();
  for (; iter != map->end(); ++iter) {
    delete iter->second;
  }
  delete map;
}


void AfterMap::RemoveJobRecordFromMap(Map *map, job_id_t job_id) {
  // Already Lock is aciquired.
  Iter iter = map->find(job_id);
  if (iter != map->end()) {
    delete iter->second;
    map->erase(iter);
  }
}

bool AfterMap::RemoveJobRecords(job_id_t job_id) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  RemoveJobRecordFromMap(map_, job_id);
  RemoveJobRecordFromMap(late_map_, job_id);
  job_done_pool_.erase(job_id);

  return true;
}

void AfterMap::RemoveWorkerRecordFromMap(Map *map, SchedulerWorker *worker) {
  // Already Lock is aciquired.
  Iter iter = map->begin();
  for (; iter != map->end();) {
    iter->second->erase(worker);
    if (iter->second->size() == 0) {
      delete iter->second;
      map->erase(iter++);
    } else {
      ++iter;
    }
  }
}

bool AfterMap::RemoveWorkerRecords(SchedulerWorker *worker) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  RemoveWorkerRecordFromMap(map_, worker);
  RemoveWorkerRecordFromMap(late_map_, worker);

  return true;
}


bool AfterMap::JobIsDone(job_id_t job_id) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  return  (job_done_pool_.find(job_id) != job_done_pool_.end());
}


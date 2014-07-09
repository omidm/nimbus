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
  * Author: Hang Qu <quhang@stanford.edu>
  */

#include <set>

#include "shared/nimbus.h"

#include "worker/job_query.h"

namespace nimbus {

JobQuery::JobQuery(Job* job) {
  job_ = job;
  has_last_barrier_job_ = false;
  query_time_ = 0;
  commit_time_ = 0;
  copy_time_ = 0;
  elimination_time_ = 0;
  spawn_time_ = 0;
}
JobQuery::~JobQuery() {}

bool JobQuery::StageJob(
    const std::string& name, const job_id_t& id,
    const IDSet<logical_data_id_t>& read,
    const IDSet<logical_data_id_t>& write,
    const Parameter& params,
    const bool sterile,
    const bool barrier) {
  struct timespec start_time;
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &start_time);
  IDSet<job_id_t> query_results;
  for (IDSet<logical_data_id_t>::ConstIter index = read.begin();
       index != read.end(); index++) {
    OutstandingAccessors& entry = outstanding_accessors_map_[*index];
    if (entry.has_outstanding_writer) {
        query_results.insert(entry.outstanding_writer);
    }
  }  // RAW
  for (IDSet<logical_data_id_t>::ConstIter index = write.begin();
       index != write.end(); index++) {
    OutstandingAccessors& entry = outstanding_accessors_map_[*index];
    if (entry.has_outstanding_writer) {
        query_results.insert(entry.outstanding_writer);
    }
  }  // WAW
  if (has_last_barrier_job_) {
    query_results.insert(last_barrier_job_id_);
  }
  if (barrier) {
    for (std::list<ShortJobEntry>::iterator index =
         query_results_.begin();
         index != query_results_.end();
         index++) {
       query_results.insert(index->id);
    }
    outstanding_accessors_map_.clear();
    has_last_barrier_job_ = true;
    last_barrier_job_id_ = id;
  }
  clock_gettime(CLOCK_REALTIME, &t);
  query_time_ += difftime(t.tv_sec, start_time.tv_sec)
      + .000000001 * (static_cast<double>(t.tv_nsec - start_time.tv_nsec));
  clock_gettime(CLOCK_REALTIME, &start_time);
  // Put the result in staged areas.
  staged_jobs_.push_back(JobEntry());
  JobEntry& job_entry = staged_jobs_.back();
  job_entry.name = name;
  job_entry.id = id;
  job_entry.read = read;
  job_entry.write = write;
  job_entry.before = query_results;
  job_entry.params = params;
  job_entry.sterile = sterile;
  clock_gettime(CLOCK_REALTIME, &t);
  copy_time_ += difftime(t.tv_sec, start_time.tv_sec)
      + .000000001 * (static_cast<double>(t.tv_nsec - start_time.tv_nsec));
  return true;
}

bool JobQuery::CommitStagedJobs() {
  struct timespec start_time;
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &start_time);
  // Move the job from staged areas to results.
  for (std::list<JobEntry>::iterator iterator = staged_jobs_.begin();
       iterator != staged_jobs_.end();
       iterator++) {
    clock_gettime(CLOCK_REALTIME, &t);
    commit_time_ += difftime(t.tv_sec, start_time.tv_sec)
        + .000000001 * (static_cast<double>(t.tv_nsec - start_time.tv_nsec));
    clock_gettime(CLOCK_REALTIME, &start_time);
    Eliminate(&iterator->before);
    clock_gettime(CLOCK_REALTIME, &t);
    elimination_time_ += difftime(t.tv_sec, start_time.tv_sec)
        + .000000001 * (static_cast<double>(t.tv_nsec - start_time.tv_nsec));
    clock_gettime(CLOCK_REALTIME, &start_time);
    job_->SpawnComputeJob(
        iterator->name, iterator->id, iterator->read, iterator->write,
        iterator->before, IDSet<job_id_t>(), iterator->params,
        iterator->sterile);
    clock_gettime(CLOCK_REALTIME, &t);
    spawn_time_ += difftime(t.tv_sec, start_time.tv_sec)
        + .000000001 * (static_cast<double>(t.tv_nsec - start_time.tv_nsec));
    clock_gettime(CLOCK_REALTIME, &start_time);
    query_results_.push_back(ShortJobEntry());
    ShortJobEntry& temp = query_results_.back();
    temp.name = iterator->name;
    temp.id = iterator->id;
    temp.before = iterator->before;
  }
  // Cleans up outstanding accessor map.
  for (std::list<JobEntry>::iterator iterator = staged_jobs_.begin();
       iterator != staged_jobs_.end();
       iterator++) {
    for (IDSet<logical_data_id_t>::ConstIter index = iterator->write.begin();
         index != iterator->write.end(); index++) {
      OutstandingAccessors& entry = outstanding_accessors_map_[*index];
      entry.has_outstanding_writer = true;
      entry.outstanding_writer = iterator->id;
    }
  }
  staged_jobs_.clear();
  clock_gettime(CLOCK_REALTIME, &t);
  commit_time_ += difftime(t.tv_sec, start_time.tv_sec)
      + .000000001 * (static_cast<double>(t.tv_nsec - start_time.tv_nsec));
  return true;
}

bool JobQuery::CommitJob(const job_id_t& id) {
  // TODO(quhang) not implemented.
  assert(false);
  return true;
}

void JobQuery::Eliminate(IDSet<job_id_t>* before) {
  std::set<job_id_t> temp;
  int64_t scanned = 0;
  for (std::list<ShortJobEntry>::reverse_iterator index =
       query_results_.rbegin();
       index != query_results_.rend();
       index++) {
    if (before->contains(index->id)) {
      if (temp.find(index->id) != temp.end()) {
        before->remove(index->id);
      } else {
        temp.insert(index->id);
        ++scanned;
      }
    }
    if (scanned == before->size()) {
      break;
    }
    if (temp.find(index->id) != temp.end()) {
      for (IDSet<job_id_t>::IDSetIter job_iter = index->before.begin();
           job_iter != index->before.end();
           job_iter++) {
        temp.insert(*job_iter);
      }
    }  // contain statement end
  }
}

void JobQuery::GenerateDotFigure(const std::string& file_name) {
  std::ofstream fout(file_name.c_str(), std::ofstream::out);
  fout << "digraph Workflow {" << std::endl;
  for (std::list<ShortJobEntry>::iterator index = query_results_.begin();
       index != query_results_.end();
       index++) {
    fout << "\t" << "node" << index->id
         << " [label=\"" << index->name << "\"];"
         << std::endl;
    for (IDSet<job_id_t>::IDSetIter prec = index->before.begin();
         prec != index->before.end();
         prec++) {
      fout << "\tnode" << (*prec) << "->node" << index->id << std::endl;
    }
  }
  fout << "}" << std::endl;
  fout.close();
  return;
}

void JobQuery::PrintTimeProfile() {
  printf("\nquery time:%f\ncommit_time:%f\ncopy_time:%f\nelimination_time:%f\n"
         "spawn_time:%f\n", query_time_, commit_time_, copy_time_,
         elimination_time_, spawn_time_);
}

}  // namespace nimbus

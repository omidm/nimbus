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
  * This is the DependencyQuery module to load before set for jobs. It is
  * written based on the code Hang Qu had provided as JobQuery module. It does
  * not have all the optimizations from the initial code and implements a
  * simpler algorithm. The reason is that we need to split different
  * functionalities, to make the code more modular and all optimizations were
  * interleaved in a complex way.
  *
  * Author: Hang Qu <quhang@stanford.edu>
  * Modified by: Omid Mashayekhi <omidm@stanford.edu>
  */

#include <boost/unordered_set.hpp>
#include "shared/nimbus.h"
#include "worker/dependency_query.h"

namespace nimbus {

DependencyQuery::DependencyQuery() {
  has_last_barrier_job_ = false;
  group_id_counter_ = 0;
  total_job_ = 0;
  total_objects_ = 0;
  query_time_ = 0;
  commit_time_ = 0;
  copy_time_ = 0;
  elimination_time_ = 0;
  spawn_time_ = 0;
  e1_time_ = 0;
  e2_time_ = 0;
  e3_time_ = 0;
  e4_time_ = 0;
}
DependencyQuery::~DependencyQuery() {}

bool DependencyQuery::StageJobAndLoadBeforeSet(IDSet<job_id_t> *before_set,
                                               const std::string& name,
                                               const job_id_t& id,
                                               const IDSet<logical_data_id_t>& read,
                                               const IDSet<logical_data_id_t>& write,
                                               const bool barrier) {
  ++total_job_;
  total_objects_ += read.size() + write.size();
  struct timespec start_time;
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &start_time);

  for (IDSet<logical_data_id_t>::ConstIter index = read.begin();
       index != read.end(); index++) {
    OutstandingAccessors& entry = outstanding_accessors_map_[*index];
    if (entry.has_outstanding_writer) {
        before_set->insert(entry.outstanding_writer);
    }
  }  // RAW
  for (IDSet<logical_data_id_t>::ConstIter index = write.begin();
       index != write.end(); index++) {
    if (!read.contains(*index)) {
      OutstandingAccessors& entry = outstanding_accessors_map_[*index];
      if (entry.has_outstanding_writer) {
        before_set->insert(entry.outstanding_writer);
      }
    }
  }  // WAW
  if (has_last_barrier_job_) {
    before_set->insert(last_barrier_job_id_);
  }
  if (barrier) {
    for (std::vector<ShortJobEntry>::iterator index =
         query_log_.begin();
         index != query_log_.end();
         index++) {
       before_set->insert(index->id);
    }
    outstanding_accessors_map_.clear();
    has_last_barrier_job_ = true;
    last_barrier_job_id_ = id;
  }
  clock_gettime(CLOCK_REALTIME, &t);
  query_time_ += difftime(t.tv_sec, start_time.tv_sec)
      + .000000001 * (static_cast<double>(t.tv_nsec - start_time.tv_nsec));

  clock_gettime(CLOCK_REALTIME, &start_time);
  Eliminate(before_set);
  clock_gettime(CLOCK_REALTIME, &t);
  elimination_time_ += difftime(t.tv_sec, start_time.tv_sec)
    + .000000001 * (static_cast<double>(t.tv_nsec - start_time.tv_nsec));

  clock_gettime(CLOCK_REALTIME, &start_time);
  // Put the result in staged areas.
  staged_jobs_.push_back(JobEntry());
  JobEntry& job_entry = staged_jobs_.back();
  job_entry.name = name;
  job_entry.id = id;
  job_entry.read = read;
  job_entry.write = write;
  job_entry.before = *before_set;
  clock_gettime(CLOCK_REALTIME, &t);
  copy_time_ += difftime(t.tv_sec, start_time.tv_sec)
      + .000000001 * (static_cast<double>(t.tv_nsec - start_time.tv_nsec));

  return true;
}

bool DependencyQuery::MarkEndOfStage() {
  struct timespec start_time;
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &start_time);
  // Move the job from staged areas to results.
  for (std::list<JobEntry>::iterator iterator = staged_jobs_.begin();
       iterator != staged_jobs_.end();
       iterator++) {
    clock_gettime(CLOCK_REALTIME, &start_time);
    query_log_.push_back(ShortJobEntry());
    ShortJobEntry& temp = query_log_.back();
    temp.name = iterator->name;
    temp.id = iterator->id;
    temp.before.swap(iterator->before);
    temp.group_id = group_id_counter_;
    job_id_to_rank_[temp.id] = query_log_.size() - 1;
  }
  ++group_id_counter_;
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

  return true;
}

void DependencyQuery::Eliminate(IDSet<job_id_t>* before) {
  struct timespec start_time;
  struct timespec t;

  clock_gettime(CLOCK_REALTIME, &start_time);
  int64_t count = 0;
  if (before->size() <= 1) {
    return;
  }

  int64_t scanned = 0;
  std::vector<boost::unordered_set<RankId> > group_heap;
  boost::unordered_set<job_id_t> group_heap_buffer;
  group_heap.resize(group_id_counter_);

  for (IDSet<job_id_t>::IDSetIter iter = before->begin();
       iter != before->end();
       ++iter) {
    RankId rank_id = job_id_to_rank_[*iter];
    GroupId group_id = query_log_[rank_id].group_id;
    group_heap[group_id].insert(rank_id);
  }

  for (GroupId index = group_id_counter_ - 1; index >= 0; --index) {
    group_heap_buffer.clear();
    for (boost::unordered_set<RankId>::iterator iter =
         group_heap[index].begin();
         iter != group_heap[index].end();
         ++iter) {
      if (before->contains(query_log_[*iter].id)) {
        ++scanned;
        if (scanned == before->size()) {
          return;
        }
      }
      ShortJobEntry& entry = query_log_[*iter];
      group_heap_buffer.insert(entry.before.begin(), entry.before.end());
    }  // Scan the before set.

    for (boost::unordered_set<job_id_t>::iterator iter =
         group_heap_buffer.begin();
         iter != group_heap_buffer.end();
         ++iter) {
      ++count;
      before->remove(*iter);
      RankId temp_rank_id = job_id_to_rank_[*iter];
      GroupId group_id = query_log_[temp_rank_id].group_id;
      group_heap[group_id].insert(temp_rank_id);
    }  // Add to the heap.
  }
}

void DependencyQuery::GenerateDotFigure(const std::string& file_name) {
  return;
  std::ofstream fout(file_name.c_str(), std::ofstream::out);
  fout << "digraph Workflow {" << std::endl;
  for (std::vector<ShortJobEntry>::iterator index = query_log_.begin();
       index != query_log_.end();
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

void DependencyQuery::PrintTimeProfile() {
  return;
  printf("\nquery time:%f\ncommit_time:%f\ncopy_time:%f\nelimination_time:%f\n"
         "spawn_time:%f\ne1_time:%f\ne2_time:%f\ne3_time:%f\ne4_time:%f\n",
         query_time_, commit_time_, copy_time_, elimination_time_, spawn_time_,
         e1_time_, e2_time_, e3_time_, e4_time_);
  printf("\ntotal job:%"PRId64"\ntotal objects:%"PRId64"\n",
         total_job_, total_objects_);
}

}  // namespace nimbus

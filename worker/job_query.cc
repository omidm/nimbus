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

#include <boost/unordered_set.hpp>

#include "shared/nimbus.h"

#include "worker/job_query.h"

namespace nimbus {

JobQuery::JobQuery(Job* job) {
  job_ = job;
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
JobQuery::~JobQuery() {}

bool JobQuery::StageJob(
    const std::string& name, const job_id_t& id,
    IDSet<logical_data_id_t>& read,
    IDSet<logical_data_id_t>& write,
    const Parameter& params,
    const bool sterile,
    const GeometricRegion& region,
    const bool barrier) {
  ++total_job_;
  total_objects_ += read.size() + write.size();
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
    if (!read.contains(*index)) {
      OutstandingAccessors& entry = outstanding_accessors_map_[*index];
      if (entry.has_outstanding_writer) {
        query_results.insert(entry.outstanding_writer);
      }
    }
  }  // WAW
  if (has_last_barrier_job_) {
    query_results.insert(last_barrier_job_id_);
  }
  if (barrier) {
    for (std::vector<ShortJobEntry>::iterator index =
         query_log_.begin();
         index != query_log_.end();
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
  job_entry.read.swap(read);
  job_entry.write.swap(write);
  job_entry.before.swap(query_results);
  job_entry.params = params;
  job_entry.sterile = sterile;
  job_entry.region = region;
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
        iterator->sterile, iterator->region);
    clock_gettime(CLOCK_REALTIME, &t);
    spawn_time_ += difftime(t.tv_sec, start_time.tv_sec)
        + .000000001 * (static_cast<double>(t.tv_nsec - start_time.tv_nsec));

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

void JobQuery::Hint(job_id_t job_id, const GeometricRegion& region,
                    bool bottleneck) {
  hint_map_[job_id] = region;
  if (bottleneck) {
    hint_bottleneck_.insert(job_id);
  }
}

void JobQuery::Eliminate(IDSet<job_id_t>* before) {
  struct timespec start_time;
  struct timespec t;

  clock_gettime(CLOCK_REALTIME, &start_time);
  bool complete_hint = false;
  int64_t count = 0;
  if (before->size() <= 1) {
    return;
  }
  GeometricRegion hint_region(0, 0, 0, 0, 0, 0);
  for (IDSet<job_id_t>::IDSetIter iter = before->begin();
       iter != before->end();
       ++iter) {
    HintMapType::iterator loc = hint_map_.find(*iter);
    if (loc == hint_map_.end()) {
      complete_hint = false;
      break;
    }
    if (hint_bottleneck_.find(*iter) != hint_bottleneck_.end()) {
      continue;
    }
    complete_hint = true;
    hint_region.Union(loc->second);
  }
  clock_gettime(CLOCK_REALTIME, &t);
  e1_time_ += difftime(t.tv_sec, start_time.tv_sec)
      + .000000001 * (static_cast<double>(t.tv_nsec - start_time.tv_nsec));

  clock_gettime(CLOCK_REALTIME, &start_time);
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
  clock_gettime(CLOCK_REALTIME, &t);
  e2_time_ += difftime(t.tv_sec, start_time.tv_sec)
      + .000000001 * (static_cast<double>(t.tv_nsec - start_time.tv_nsec));

  clock_gettime(CLOCK_REALTIME, &start_time);
  for (GroupId index = group_id_counter_ - 1; index >= 0; --index) {
    group_heap_buffer.clear();
    for (boost::unordered_set<RankId>::iterator iter =
         group_heap[index].begin();
         iter != group_heap[index].end();
         ++iter) {
      if (hint_bottleneck_.find(query_log_[*iter].id)
          != hint_bottleneck_.end()) {
        // Clean and finish.
        assert(group_heap[index].size() == 1);
        IDSet<job_id_t> temp_idset;
        for (IDSet<job_id_t>::IDSetIter iter_j = before->begin();
             iter_j != before->end();
             ++iter_j) {
          RankId temp_rank_id = job_id_to_rank_[*iter_j];
          GroupId group_id = query_log_[temp_rank_id].group_id;
          if (group_id >= index) {
            temp_idset.insert(*iter_j);
          }
        }
        before->swap(temp_idset);
        clock_gettime(CLOCK_REALTIME, &t);
        e3_time_ += difftime(t.tv_sec, start_time.tv_sec)
            + .000000001
            * (static_cast<double>(t.tv_nsec - start_time.tv_nsec));

        return;
      }
      if (before->contains(query_log_[*iter].id)) {
        ++scanned;
        if (scanned == before->size()) {
          return;
        }
      }
      ShortJobEntry& entry = query_log_[*iter];
      group_heap_buffer.insert(entry.before.begin(), entry.before.end());
    }  // Scan the before set.
    clock_gettime(CLOCK_REALTIME, &t);
    e3_time_ += difftime(t.tv_sec, start_time.tv_sec)
        + .000000001 * (static_cast<double>(t.tv_nsec - start_time.tv_nsec));

    clock_gettime(CLOCK_REALTIME, &start_time);
    for (boost::unordered_set<job_id_t>::iterator iter =
         group_heap_buffer.begin();
         iter != group_heap_buffer.end();
         ++iter) {
      if (complete_hint &&
          !hint_map_[*iter].Intersects(&hint_region)) {
        continue;
      }
      ++count;
      before->remove(*iter);
      RankId temp_rank_id = job_id_to_rank_[*iter];
      GroupId group_id = query_log_[temp_rank_id].group_id;
      group_heap[group_id].insert(temp_rank_id);
    }  // Add to the heap.
    clock_gettime(CLOCK_REALTIME, &t);
    e4_time_ += difftime(t.tv_sec, start_time.tv_sec)
        + .000000001 * (static_cast<double>(t.tv_nsec - start_time.tv_nsec));

    clock_gettime(CLOCK_REALTIME, &start_time);
  }
}

/*
void JobQuery::Eliminate(IDSet<job_id_t>* before) {
  int64_t count = 0;
  int64_t init_count = 0;
  bool complete_hint = true;
  if (before->size() <= 1) {
    return;
  }
  GeometricRegion hint_region(0, 0, 0, 0, 0, 0);
  for (IDSet<job_id_t>::IDSetIter iter = before->begin();
       iter != before->end();
       ++iter) {
    HintMapType::iterator loc = hint_map_.find(*iter);
    if (loc == hint_map_.end()) {
      complete_hint = false;
      break;
    }
    hint_region.Union(loc->second);
  }
  int64_t scanned = 0;
  std::vector<int64_t> heap;
  boost::unordered_set<int64_t> exist;
  for (IDSet<job_id_t>::IDSetIter iter = before->begin();
       iter != before->end();
       ++iter) {
    heap.push_back(job_id_to_rank_[*iter]);
    exist.insert(job_id_to_rank_[*iter]);
  }
  init_count = heap.size();
  std::make_heap(heap.begin(), heap.end());  // max_heap.
  while (scanned < before->size() && !heap.empty()) {
    pop_heap(heap.begin(), heap.end());
    int64_t top_element = heap.back();
    heap.pop_back();
    exist.erase(top_element);
    ShortJobEntry& entry = query_log_[top_element];
    if (before->contains(entry.id)) {
      ++scanned;
    }
    for (IDSet<job_id_t>::IDSetIter job_iter = entry.before.begin();
         job_iter != entry.before.end();
         job_iter++) {
      if (complete_hint && !hint_map_[*job_iter].Intersects(&hint_region)) {
        continue;
      }
      before->remove(*job_iter);
      int64_t cand_job_id = job_id_to_rank_[*job_iter];
      if (exist.count(cand_job_id) == 0) {
        ++count;
        heap.push_back(cand_job_id);
        push_heap(heap.begin(), heap.end());
        exist.insert(cand_job_id);
      }
    }
  }
  printf("\n%d jobs scanned, %d, %s\n", count, init_count,
         hint_region.ToNetworkData().c_str());
}
*/

void JobQuery::GenerateDotFigure(const std::string& file_name) {
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

void JobQuery::PrintTimeProfile() {
  return;
  printf("\nquery time:%f\ncommit_time:%f\ncopy_time:%f\nelimination_time:%f\n"
         "spawn_time:%f\ne1_time:%f\ne2_time:%f\ne3_time:%f\ne4_time:%f\n",
         query_time_, commit_time_, copy_time_, elimination_time_, spawn_time_,
         e1_time_, e2_time_, e3_time_, e4_time_);
  printf("\ntotal job:%"PRId64"\ntotal objects:%"PRId64"\n",
         total_job_, total_objects_);
}

}  // namespace nimbus

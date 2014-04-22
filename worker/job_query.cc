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
  RESOLVE_WAW = true;
  RESOLVE_WAR = false;
  DISABLE = false;
}
JobQuery::~JobQuery() {}

bool JobQuery::StageJob(
    const std::string& name, const job_id_t& id,
    const IDSet<logical_data_id_t>& read,
    const IDSet<logical_data_id_t>& write,
    const Parameter& params,
    const bool sterile,
    const bool barrier) {
  return StageJob(name, id, read, write,
                  IDSet<job_id_t>(), IDSet<job_id_t>(),
                  params, sterile, barrier);
}
bool JobQuery::StageJob(
    const std::string& name, const job_id_t& id,
    const IDSet<logical_data_id_t>& read,
    const IDSet<logical_data_id_t>& write,
    const IDSet<job_id_t>& before,
    const IDSet<job_id_t>& after,
    const Parameter& params,
    const bool sterile,
    const bool barrier) {
  if (DISABLE) {
    job_->SpawnComputeJob(
        name, id, read, write, before, after, params, sterile);
    return true;
  }
  IDSet<job_id_t> query_results;
  for (IDSet<logical_data_id_t>::ConstIter index = read.begin();
       index != read.end(); index++) {
    OutstandingAccessorsMap::iterator find_entry
        = outstanding_accessors_map_.find(*index);
    if (find_entry != outstanding_accessors_map_.end()) {
      OutstandingAccessorsMap::mapped_type& entry = find_entry->second;
      if (entry.has_outstanding_writer) {
        query_results.insert(entry.outstanding_writer);
      }
    }
  }  // RAW
  for (IDSet<logical_data_id_t>::ConstIter index = write.begin();
       index != write.end(); index++) {
    OutstandingAccessorsMap::iterator find_entry
        = outstanding_accessors_map_.find(*index);
    if (find_entry != outstanding_accessors_map_.end()) {
      OutstandingAccessorsMap::mapped_type& entry = find_entry->second;
      if (RESOLVE_WAW && entry.has_outstanding_writer) {
        query_results.insert(entry.outstanding_writer);
      }  // WAW
      if (RESOLVE_WAR && !entry.outstanding_reader_list.empty()) {
        typedef std::list<job_id_t> ReaderList;
        const ReaderList& reader_list = entry.outstanding_reader_list;
        for (ReaderList::const_iterator read_job_index = reader_list.begin();
             read_job_index != reader_list.end();
             read_job_index++) {
          query_results.insert(*read_job_index);
        }
      }  // WAR
    }
  }
  if (barrier) {
    for (OutstandingAccessorsMap::iterator index =
         outstanding_accessors_map_.begin();
         index != outstanding_accessors_map_.end();
         index++) {
      if (index->second.has_outstanding_writer) {
        query_results.insert(index->second.outstanding_writer);
      }
      typedef std::list<job_id_t> ReaderList;
      for (ReaderList::const_iterator read_job_index =
           index->second.outstanding_reader_list.begin();
           read_job_index != index->second.outstanding_reader_list.end();
           read_job_index++) {
        query_results.insert(*read_job_index);
      }
    }
  }
  // Put the result in staged areas.
  // TODO(quhang) redundent copy happens.
  staged_jobs_.push_back(JobEntry());
  JobEntry& job_entry = staged_jobs_.back();
  job_entry.name = name;
  job_entry.id = id;
  job_entry.read = read;
  job_entry.write = write;
  job_entry.before = query_results;
  job_entry.params = params;
  job_entry.sterile = sterile;

  return true;
}

bool JobQuery::CommitStagedJobs() {
  if (DISABLE) {
    return true;
  }
  // Move the job from staged areas to results.
  for (std::list<JobEntry>::iterator iterator = staged_jobs_.begin();
       iterator != staged_jobs_.end();
       iterator++) {
    // TODO(quhang) eliminate redundent dependencies.
    Eliminate(&iterator->before);
    job_->SpawnComputeJob(
        iterator->name, iterator->id, iterator->read, iterator->write,
        iterator->before, IDSet<job_id_t>(), iterator->params,
        iterator->sterile);
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
    for (IDSet<logical_data_id_t>::ConstIter index = iterator->read.begin();
         index != iterator->read.end(); index++) {
      outstanding_accessors_map_[*index].outstanding_reader_list.push_back(
          iterator->id);
    }
    for (IDSet<logical_data_id_t>::ConstIter index = iterator->write.begin();
         index != iterator->write.end(); index++) {
      outstanding_accessors_map_[*index].has_outstanding_writer = true;
      outstanding_accessors_map_[*index].outstanding_writer = iterator->id;
      outstanding_accessors_map_[*index].outstanding_reader_list.clear();
    }
  }
  staged_jobs_.clear();
  return true;
}

bool JobQuery::CommitJob(const job_id_t& id) {
  if (DISABLE) {
    return true;
  }
  // TODO(quhang) not implemented.
  assert(false);
  return true;
}

// TODO(quhang) too much overhead.
void JobQuery::Eliminate(IDSet<job_id_t>* before) {
  std::set<job_id_t> temp;
  for (std::list<ShortJobEntry>::reverse_iterator index =
       query_results_.rbegin();
       index != query_results_.rend();
       index++) {
    if (before->contains(index->id)) {
      if (temp.find(index->id) != temp.end()) {
        before->remove(index->id);
      }
      temp.insert(index->id);
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
  if (DISABLE) {
    return;
  }
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

}  // namespace nimbus

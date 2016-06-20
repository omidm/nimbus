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
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#include <boost/unordered_set.hpp>
#include "src/shared/nimbus.h"
#include "src/worker/dependency_query.h"

namespace nimbus {

DependencyQuery::DependencyQuery() {
  stage_id_counter_ = 0;
  last_stage_partial_writers_ = new PartialWritersMap();
  current_stage_partial_writers_ = new PartialWritersMap();
  total_job_ = 0;
  pruning_time_ = 0;
  insertion_time_ = 0;
}

DependencyQuery::~DependencyQuery() {
  if ((last_stage_partial_writers_->size() > 0) ||
      (current_stage_partial_writers_->size() > 0)) {
    dbg(DBG_ERROR, "ERROR: there are pending partial writes without reduction.\n"); //NOLINT
    exit(-1);
  }
}

bool DependencyQuery::StageJobAndLoadBeforeSet(IDSet<job_id_t> *before_set,
                                               const std::string& name,
                                               const job_id_t& id,
                                               const IDSet<logical_data_id_t>& read,
                                               const IDSet<logical_data_id_t>& write,
                                               const IDSet<logical_data_id_t>& scratch,
                                               const IDSet<logical_data_id_t>& reduce,
                                               const bool barrier) {
  ++total_job_;
  struct timespec start_time, end_time;
#ifdef __MACH__  // OS X does not have clock_gettime, use clock_get_time
  {
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    start_time.tv_sec = mts.tv_sec;
    start_time.tv_nsec = mts.tv_nsec;
  }
#else
  clock_gettime(CLOCK_REALTIME, &start_time);
#endif

  {
    IDSet<logical_data_id_t>::ConstIter iter = read.begin();
    for (; iter != read.end(); ++iter) {
      LastWriterMap::iterator it = last_writers_.find(*iter);
      if (it != last_writers_.end()) {
        before_set->insert(it->second);
      }
    }
  }

  {
    IDSet<logical_data_id_t>::ConstIter iter = write.begin();
    for (; iter != write.end(); ++iter) {
      if (read.contains(*iter)) {
        continue;
      }
      LastWriterMap::iterator it = last_writers_.find(*iter);
      if (it != last_writers_.end()) {
        before_set->insert(it->second);
      }
    }
  }

  {
    IDSet<logical_data_id_t>::ConstIter iter = reduce.begin();
    for (; iter != reduce.end(); ++iter) {
      PartialWritersMap::iterator it = last_stage_partial_writers_->find(*iter);
      if (it == last_stage_partial_writers_->end()) {
        dbg(DBG_ERROR, "ERROR: no partial write found to reduce!\n");
        exit(-1);
      }
      assert(it->second.size() > 0);
      PartialWriters::iterator i = it->second.begin();
      for (; i != it->second.end(); ++i) {
        before_set->insert(*i);
      }
      last_stage_partial_writers_->erase(it);
    }
  }

  {
    IDSet<logical_data_id_t>::ConstIter iter = scratch.begin();
    for (; iter != scratch.end(); ++iter) {
      current_stage_partial_writers_->operator[](*iter).push_back(id);
    }
  }

  if (last_barrier_.valid) {
    before_set->insert(last_barrier_.id);
  }

  if (barrier) {
    JobList::iterator iter = jobs_.begin();
    for (; iter != jobs_.end(); ++iter) {
       before_set->insert(iter->id);
    }
    if (last_barrier_.valid) {
      if (last_barrier_.stage_id >= stage_id_counter_) {
        dbg(DBG_ERROR, "ERROR: barrier stage id does not make sense, probable mutilple barriers in one stage!\n"); // NOLINT
        exit(-1);
      }
    }
    last_barrier_.valid = true;
    last_barrier_.id = id;
    last_barrier_.stage_id = stage_id_counter_;
  }

#ifdef __MACH__  // OS X does not have clock_gettime, use clock_get_time
  {
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    end_time.tv_sec = mts.tv_sec;
    end_time.tv_nsec = mts.tv_nsec;
  }
#else
  clock_gettime(CLOCK_REALTIME, &end_time);
#endif
  insertion_time_ += difftime(end_time.tv_sec, start_time.tv_sec)
      + .000000001 * (static_cast<double>(end_time.tv_nsec - start_time.tv_nsec));

#ifdef __MACH__  // OS X does not have clock_gettime, use clock_get_time
  {
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    start_time.tv_sec = mts.tv_sec;
    start_time.tv_nsec = mts.tv_nsec;
  }
#else
  clock_gettime(CLOCK_REALTIME, &start_time);
#endif
  Prune(before_set);
#ifdef __MACH__  // OS X does not have clock_gettime, use clock_get_time
  {
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    end_time.tv_sec = mts.tv_sec;
    end_time.tv_nsec = mts.tv_nsec;
  }
#else
  clock_gettime(CLOCK_REALTIME, &end_time);
#endif
  pruning_time_ += difftime(end_time.tv_sec, start_time.tv_sec)
    + .000000001 * (static_cast<double>(end_time.tv_nsec - start_time.tv_nsec));

  staged_jobs_.push_back(Job());
  Job& job = staged_jobs_.back();
  job.id = id;
  job.name = name;
  job.write = write;
  job.before = *before_set;

  return true;
}

bool DependencyQuery::MarkEndOfStage() {
  {
    JobList::iterator iter = staged_jobs_.begin();
    for (; iter != staged_jobs_.end(); ++iter) {
      jobs_.push_back(Job());
      Job& temp = jobs_.back();
      temp.id = iter->id;
      temp.name = iter->name;
      temp.before.swap(iter->before);
      temp.stage_id = stage_id_counter_;
      job_id_to_index_[temp.id] = jobs_.size() - 1;
    }
  }

  ++stage_id_counter_;

  {
    JobList::iterator iter = staged_jobs_.begin();
    for (; iter != staged_jobs_.end(); ++iter) {
      IDSet<logical_data_id_t>::ConstIter it = iter->write.begin();
      for (; it != iter->write.end(); ++it) {
        last_writers_[*it] = iter->id;
      }
    }
  }

  if (last_stage_partial_writers_->size() > 0) {
    dbg(DBG_ERROR, "ERROR: there are pending partial writes from last stage without reduction.\n"); //NOLINT
    exit(-1);
  }
  delete last_stage_partial_writers_;
  last_stage_partial_writers_ = current_stage_partial_writers_;
  current_stage_partial_writers_ = new PartialWritersMap();

  staged_jobs_.clear();

  return true;
}

void DependencyQuery::Prune(IDSet<job_id_t>* before) {
  if (before->size() <= 1) {
    return;
  }
  assert(stage_id_counter_ > 0);

  int64_t scanned = 0;
  std::vector<boost::unordered_set<IndexId> > stage_heap;
  stage_heap.resize(stage_id_counter_);

  {
    IDSet<job_id_t>::IDSetIter iter = before->begin();
    for (; iter != before->end(); ++iter) {
      IndexId index = job_id_to_index_[*iter];
      StageId stage_id = jobs_[index].stage_id;
      stage_heap[stage_id].insert(index);
    }
  }

  for (StageId sid = stage_id_counter_ - 1; sid >= 0; --sid) {
    boost::unordered_set<job_id_t> buffer;
    {
      boost::unordered_set<IndexId>::iterator iter = stage_heap[sid].begin();
      for (; iter != stage_heap[sid].end(); ++iter) {
        if (before->contains(jobs_[*iter].id)) {
          ++scanned;
          if (scanned == before->size()) {
            return;
          }
        }
        Job& entry = jobs_[*iter];
        buffer.insert(entry.before.begin(), entry.before.end());
      }  // Scan the before set.
    }

    {
      boost::unordered_set<job_id_t>::iterator iter = buffer.begin();
      for (; iter != buffer.end(); ++iter) {
        before->remove(*iter);
        if (scanned == before->size()) {
          return;
        }
        IndexId index = job_id_to_index_[*iter];
        StageId stage_id = jobs_[index].stage_id;
        stage_heap[stage_id].insert(index);
      }  // Add to the heap.
    }
  }
}

void DependencyQuery::GenerateDotFigure(const std::string& file_name) {
  return;
  std::ofstream fout(file_name.c_str(), std::ofstream::out);
  fout << "digraph Workflow {" << std::endl;
  JobList::iterator iter = jobs_.begin();
  for (; iter != jobs_.end(); ++iter) {
    fout << "\t" << "node" << iter->id
         << " [label=\"" << iter->name << "\"];"
         << std::endl;
    IDSet<job_id_t>::IDSetIter it = iter->before.begin();
    for (; it != iter->before.end(); ++it) {
      fout << "\tnode" << (*it) << "->node" << iter->id << std::endl;
    }
  }
  fout << "}" << std::endl;
  fout.close();
  return;
}

void DependencyQuery::PrintTimeProfile() {
  return;
  printf("\nDEPENDENCY QUERY: total job:%lu insertion time:%f, pruning time :%f\n",
         total_job_, insertion_time_, pruning_time_);
}

}  // namespace nimbus

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
  * Complex job entry in the job table of the job manager. This job contains a
  * group of compute or explicit copy jobs that are spawned with in a template.
  * The meta data calculated for the template including job dependencies and
  * versioning information is precomputed and is accessible by the complex job.
  * The idea is no avoid expanding the jobs within the template spawning. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "scheduler/complex_job_entry.h"
#include "scheduler/template_entry.h"

using namespace nimbus; // NOLINT

ComplexJobEntry::ComplexJobEntry() {
  job_type_ = JOB_CMPX;
  job_name_ = NIMBUS_COMPLEX_JOB_NAME;
  assign_index_ = 0;
  parent_job_ids_set_ = false;
  job_map_complete_ = false;
}

ComplexJobEntry::ComplexJobEntry(const job_id_t& job_id,
                                 const job_id_t& parent_job_id,
                                 TemplateEntry* template_entry,
                                 const std::vector<job_id_t>& inner_job_ids,
                                 const std::vector<job_id_t>& outer_job_ids,
                                 const std::vector<Parameter>& parameters) {
  job_type_ = JOB_CMPX;
  job_name_ = NIMBUS_COMPLEX_JOB_NAME;
  job_id_ = job_id;
  parent_job_id_ = parent_job_id;
  template_entry_ = template_entry;
  inner_job_ids_ = inner_job_ids;
  outer_job_ids_ = outer_job_ids;
  parameters_ = parameters;
  assign_index_ = 0;
  parent_job_ids_set_ = false;
  job_map_complete_ = false;

  // parent should be explicitally in before set - omidm
  // currentrly before set of complex job is only parent job - omidm
  assert(before_set_.size() == 0);
  before_set_.insert(parent_job_id);

  assignment_dependencies_ = before_set_;
  versioning_dependencies_ = before_set_;
}

ComplexJobEntry::~ComplexJobEntry() {
}

TemplateEntry* ComplexJobEntry::template_entry() const {
  return template_entry_;
}

std::vector<job_id_t> ComplexJobEntry::inner_job_ids() const {
  return inner_job_ids_;
}

const std::vector<job_id_t>* ComplexJobEntry::inner_job_ids_p() const {
  return &inner_job_ids_;
}

std::vector<job_id_t> ComplexJobEntry::outer_job_ids() const {
  return outer_job_ids_;
}

const std::vector<job_id_t>* ComplexJobEntry::outer_job_ids_p() const {
  return &outer_job_ids_;
}

std::vector<Parameter> ComplexJobEntry::parameters() const {
  return parameters_;
}

const std::vector<Parameter>* ComplexJobEntry::parameters_p() const {
  return &parameters_;
}


void ComplexJobEntry::SetParentJobIds() {
  if (parent_job_ids_set_) {
    return;
  }

  std::list<size_t> indices;
  if (!template_entry_->GetParentJobIndices(&indices)) {
    assert(false);
    return;
  }

  parent_job_ids_.clear();
  std::list<size_t>::iterator iter = indices.begin();
  for (; iter != indices.end(); ++iter) {
    assert((*iter) < inner_job_ids_.size());
    parent_job_ids_.push_back(inner_job_ids_[*iter]);
  }

  parent_job_ids_set_ = true;
  return;
}



size_t ComplexJobEntry::GetParentJobIds(std::list<job_id_t>* list) {
  if (!parent_job_ids_set_) {
    SetParentJobIds();
  }

  list->clear();
  *list = parent_job_ids_;
  return parent_job_ids_.size();
}

size_t ComplexJobEntry::GetJobsForAssignment(JobEntryList* list, size_t max_num, bool append) {
  size_t count = 0;
  if (!append) {
    list->clear();
  }

  size_t index = assign_index_;
  for (; (index < inner_job_ids_.size()) && (count < max_num); ++index) {
    ShadowJobEntry* shadow_job;
    ShadowJobEntryMap::iterator it = jobs_.find(inner_job_ids_[index]);
    if (it != jobs_.end()) {
      shadow_job = it->second;
    } else {
      TemplateJobEntry* job = template_entry_->GetJobAtIndex(index);
      shadow_job =
        new ShadowJobEntry(job->job_name(),
                           inner_job_ids_[index],
                           job->read_set_p(),
                           job->write_set_p(),
                           job->union_set_p(),
                           job->vmap_read_diff(),
                           job->vmap_write_diff(),
                           parent_job_id_,
                           0,  // future_job_id, currently not supported - omidm
                           job->sterile(),
                           job->region(),
                           parameters_[index],
                           this);
      jobs_[inner_job_ids_[index]] = shadow_job;
    }

    list->push_back(shadow_job);
    count++;
  }

  assign_index_ += count;
  return count;
}


size_t ComplexJobEntry::GetParentJobs(ShadowJobEntryList* list, bool append) {
  if (!parent_job_ids_set_) {
    SetParentJobIds();
  }

  size_t count = 0;
  if (!append) {
    list->clear();
  }

  std::list<size_t>::iterator iter = parent_job_ids_.begin();
  for (; iter != parent_job_ids_.end(); ++iter) {
    ShadowJobEntry* shadow_job;
    ShadowJobEntryMap::iterator it = jobs_.find(*iter);
    if (it != jobs_.end()) {
      shadow_job = it->second;
    } else {
      size_t index;
      if (!GetJobIndex(*iter, &index)) {
        assert(false);
      }

      TemplateJobEntry* job = template_entry_->GetJobAtIndex(index);
      shadow_job =
        new ShadowJobEntry(job->job_name(),
                           inner_job_ids_[index],
                           job->read_set_p(),
                           job->write_set_p(),
                           job->union_set_p(),
                           job->vmap_read_diff(),
                           job->vmap_write_diff(),
                           parent_job_id_,
                           0,  // future_job_id, currently not supported - omidm
                           job->sterile(),
                           job->region(),
                           parameters_[index],
                           this);
      jobs_[inner_job_ids_[index]] = shadow_job;
    }

    assert(!shadow_job->sterile());
    list->push_back(shadow_job);
    count++;
  }

  return count;
}


bool ComplexJobEntry::GetShadowJobEntry(job_id_t job_id, ShadowJobEntry*& shadow_job) {
  {
    boost::unique_lock<boost::mutex> lock(mutex_);
    if (removed_job_ids_.size() != 0) {
      if (removed_job_ids_.find(job_id) != removed_job_ids_.end()) {
        shadow_job = NULL;
        return false;
      }
    }
  }

  if (job_map_complete_) {
    ShadowJobEntryMap::iterator iter = jobs_.find(job_id);
    if (iter != jobs_.end()) {
      shadow_job = iter->second;
      return true;
    }
  } else {
    size_t index;
    if (GetJobIndex(job_id, &index)) {
      ShadowJobEntryMap::iterator iter = jobs_.find(inner_job_ids_[index]);
      if (iter != jobs_.end()) {
        shadow_job = iter->second;
        return true;
      } else {
        TemplateJobEntry* job = template_entry_->GetJobAtIndex(index);
        ShadowJobEntry* sj =
          new ShadowJobEntry(job->job_name(),
              inner_job_ids_[index],
              job->read_set_p(),
              job->write_set_p(),
              job->union_set_p(),
              job->vmap_read_diff(),
              job->vmap_write_diff(),
              parent_job_id_,
              0,  // future_job_id, currently not supported - omidm
              job->sterile(),
              job->region(),
              parameters_[index],
              this);
        jobs_[inner_job_ids_[index]] = sj;
        shadow_job = sj;
        return true;
      }
    }
  }

  shadow_job = NULL;
  return false;
}


ShadowJobEntryMap ComplexJobEntry::jobs() {
  if (!job_map_complete_) {
    CompleteJobMap();
  }

  return jobs_;
}

const ShadowJobEntryMap* ComplexJobEntry::jobs_p() {
  if (!job_map_complete_) {
    CompleteJobMap();
  }

  return &jobs_;
}

void ComplexJobEntry::CompleteJobMap() {
  if (job_map_complete_) {
    return;
  }

  size_t index = 0;
  for (; index < inner_job_ids_.size(); ++index) {
    ShadowJobEntryMap::iterator it = jobs_.find(inner_job_ids_[index]);
    if (it != jobs_.end()) {
      continue;
    } else {
      TemplateJobEntry* job = template_entry_->GetJobAtIndex(index);
      ShadowJobEntry* shadow_job =
        new ShadowJobEntry(job->job_name(),
                           inner_job_ids_[index],
                           job->read_set_p(),
                           job->write_set_p(),
                           job->union_set_p(),
                           job->vmap_read_diff(),
                           job->vmap_write_diff(),
                           parent_job_id_,
                           0,  // future_job_id, currently not supported - omidm
                           job->sterile(),
                           job->region(),
                           parameters_[index],
                           this);
      jobs_[inner_job_ids_[index]] = shadow_job;
    }
  }

  job_map_complete_ = true;
  return;
}

bool ComplexJobEntry::GetJobIndex(job_id_t job_id, size_t* index) {
  size_t idx = 0;
  std::vector<job_id_t>::iterator iter = inner_job_ids_.begin();
  for (; iter != inner_job_ids_.end(); ++iter) {
    if ((*iter) == job_id) {
      *index = idx;
      return true;
    }

    ++idx;
  }

  return false;
}

bool ComplexJobEntry::DrainedAllJobsForAssignment() {
  return (assign_index_ == inner_job_ids_.size());
}

void ComplexJobEntry::MarkJobAssigned(job_id_t job_id) {
  // TODO(omidm): Implement
}

void ComplexJobEntry::MarkJobDone(job_id_t job_id) {
  done_job_ids_.insert(job_id);
}

bool ComplexJobEntry::AllJobsDone() {
  boost::unique_lock<boost::mutex> lock(mutex_);
  return (done_job_ids_.size() == inner_job_ids_.size());
}

void ComplexJobEntry::MarkJobRemoved(job_id_t job_id) {
  boost::unique_lock<boost::mutex> lock(mutex_);
  removed_job_ids_.insert(job_id);
}

bool ComplexJobEntry::AllJobsRemoved() {
  boost::unique_lock<boost::mutex> lock(mutex_);
  return (removed_job_ids_.size() == inner_job_ids_.size());
}

ComplexJobEntry::Cursor::Cursor() {
}

ComplexJobEntry::Cursor::~Cursor() {
}

ComplexJobEntry::Cursor::State ComplexJobEntry::Cursor::state() {
  return state_;
}

size_t ComplexJobEntry::Cursor::index() {
  return index_;
}

void ComplexJobEntry::Cursor::set_state(State state) {
  state_ = state;
}

void ComplexJobEntry::Cursor::set_index(size_t index) {
  index_ = index;
}








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


#define VMAP_MIN_DENSITY  .2


using namespace nimbus; // NOLINT

ComplexJobEntry::ComplexJobEntry() {
  Initialize();
}

ComplexJobEntry::ComplexJobEntry(const job_id_t& job_id,
                                 const job_id_t& parent_job_id,
                                 TemplateEntry* template_entry,
                                 const std::vector<job_id_t>& inner_job_ids,
                                 const std::vector<job_id_t>& outer_job_ids,
                                 const std::vector<Parameter>& parameters) {
  Initialize();

  job_id_ = job_id;
  parent_job_id_ = parent_job_id;
  template_entry_ = template_entry;
  inner_job_ids_ = inner_job_ids;
  outer_job_ids_ = outer_job_ids;
  parameters_ = parameters;

  double ldid_density =
    static_cast<double>(template_entry_->ldid_count()) /
    (template_entry_->max_ldid() - template_entry_->min_ldid());

  std::cout << "OMID: "
            << template_entry->template_name()
            << " density: "
            << ldid_density << std::endl;

  if (ldid_density > VMAP_MIN_DENSITY) {
  vmap_read_ = boost::shared_ptr<VersionMap>(
      new DenseVersionMap(template_entry_->min_ldid(), template_entry_->max_ldid()));
  } else {
    std::cout << "WARNING: did not use DenseVersionMap for "
              << template_entry->template_name()
              << " logical data id density was "
              << ldid_density << std::endl;
  }

  size_t idx = 0;
  std::vector<job_id_t>::const_iterator iter = inner_job_ids.begin();
  for (; iter != inner_job_ids.end(); ++iter) {
    shadow_job_ids_[*iter] = idx;
    ++idx;
  }

  // parent should be explicitally in before set - omidm
  // currentrly before set of complex job is only parent job - omidm
  assert(before_set_.size() == 0);
  before_set_.insert(parent_job_id);

  assignment_dependencies_ = before_set_;
  versioning_dependencies_ = before_set_;
}

void ComplexJobEntry::Initialize() {
  job_type_ = JOB_CMPX;
  job_name_ = NIMBUS_COMPLEX_JOB_NAME;
  template_entry_ = NULL;
  parent_job_ids_set_ = false;
  parent_job_indices_set_ = false;
  shadow_jobs_complete_ = false;
  initialized_cursor_ = false;
  drained_all_ = false;
}


ComplexJobEntry::~ComplexJobEntry() {
  ShadowJobEntryMap::iterator iter = shadow_jobs_.begin();
  for (; iter != shadow_jobs_.end(); ++iter) {
    delete iter->second;
  }
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
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  if (parent_job_ids_set_) {
    return;
  }

  SetParentJobIndices();

  parent_job_ids_.clear();
  std::list<size_t>::iterator iter = parent_job_indices_.begin();
  for (; iter != parent_job_indices_.end(); ++iter) {
    assert((*iter) < inner_job_ids_.size());
    parent_job_ids_.push_back(inner_job_ids_[*iter]);
  }

  parent_job_ids_set_ = true;
  return;
}

void ComplexJobEntry::SetParentJobIndices() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  if (parent_job_indices_set_) {
    return;
  }

  if (!template_entry_->GetParentJobIndices(&parent_job_indices_)) {
    assert(false);
    return;
  }

  parent_job_indices_set_ = true;
  return;
}

size_t ComplexJobEntry::GetParentJobIds(std::list<job_id_t>* list) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  SetParentJobIds();

  list->clear();
  *list = parent_job_ids_;
  return parent_job_ids_.size();
}

size_t ComplexJobEntry::GetShadowJobsForAssignment(JobEntryList* list,
                                                   size_t max_num,
                                                   bool append) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  if (!initialized_cursor_) {
    template_entry_->InitializeCursor(&cursor_);
    initialized_cursor_ = true;
  }

  size_t count = 0;
  if (!append) {
    list->clear();
  }

  // JobEntryList temp_list;
  while (count < max_num) {
    size_t index = cursor_.index();
    assert(index < inner_job_ids_.size());
    ShadowJobEntry* shadow_job;

    OMIDGetShadowJobEntryByIndex(index, shadow_job);

    list->push_back(shadow_job);
    // temp_list.push_back(shadow_job);
    ++count;

    if (cursor_.state() == Cursor::END_ALL) {
      drained_all_ = true;
      shadow_job->set_to_finalize_binding_template(true);
      break;
    } else if (cursor_.state() == Cursor::END_BATCH) {
      template_entry_->AdvanceCursorForAssignment(&cursor_);
      break;
    } else {
      template_entry_->AdvanceCursorForAssignment(&cursor_);
      // For multi-threaded assignment to work when we Finalize the binding
      // template, the last job needs to be assigned in a batch by it's own and
      // all other jobs needs to be assigned in previous batches. Otherwise,
      // there could be a race - Finalize gets called before all jobs are
      // assigned. -omidm
      if (cursor_.state() == Cursor::END_ALL) {
        break;
      }
    }
  }


//  if (temp_list.size() > 0) {
//    JobEntryList::reverse_iterator iter = temp_list.rbegin();
//    if ((*iter)->job_name() == "--reincorporate_particles" ||
//        (*iter)->job_name() == "--delete_particles" ||
//        (*iter)->job_name() == "--adjust_phi") {
//      for (; iter != temp_list.rend(); ++iter) {
//        list->push_back(*iter);
//      }
//    } else {
//      JobEntryList::iterator it = temp_list.begin();
//      for (; it != temp_list.end(); ++it) {
//        list->push_back(*it);
//      }
//    }
//  }

  return count;
}


size_t ComplexJobEntry::OMIDGetParentShadowJobs(ShadowJobEntryList* list, bool append) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  SetParentJobIndices();

  size_t count = 0;
  if (!append) {
    list->clear();
  }

  std::list<size_t>::iterator iter = parent_job_indices_.begin();
  for (; iter != parent_job_indices_.end(); ++iter) {
    ShadowJobEntry* shadow_job;

    OMIDGetShadowJobEntryByIndex(*iter, shadow_job);

    assert(!shadow_job->sterile());
    list->push_back(shadow_job);
    count++;
  }

  return count;
}

bool ComplexJobEntry::OMIDGetShadowJobEntryById(job_id_t job_id, ShadowJobEntry*& shadow_job) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  IdMap::iterator iter = shadow_job_ids_.find(job_id);
  if (iter == shadow_job_ids_.end()) {
    shadow_job = NULL;
    return false;
  }

  return OMIDGetShadowJobEntryByIndex(iter->second, shadow_job);
}




bool ComplexJobEntry::OMIDGetShadowJobEntryByIndex(size_t index, ShadowJobEntry*& shadow_job) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  assert(index >= 0);
  assert(index < inner_job_ids_.size());
  job_id_t job_id = inner_job_ids_[index];

  ShadowJobEntryMap::iterator iter = shadow_jobs_.find(job_id);
  if (iter != shadow_jobs_.end()) {
    shadow_job = iter->second;
    return true;
  }

  TemplateJobEntry* tj = template_entry_->GetJobAtIndex(index);
  IDSet<job_id_t> before_set;

  // TODO(omidm): This line should always remain commented out. other wise
  // explicit job done from controller is required.
  // HACK
  // if (!tj->sterile())
  //   template_entry_->LoadBeforeSet(&before_set, index, inner_job_ids_, outer_job_ids_);
  // HACK

  ShadowJobEntry* sj =
    new ShadowJobEntry(tj->job_name(),
                       job_id,
                       tj->read_set_p(),
                       tj->write_set_p(),
                       tj->union_set_p(),
                       before_set,
                       tj->vmap_read_diff(),
                       tj->vlist_write_diff(),
                       parent_job_id_,
                       0,  // future_job_id, currently not supported - omidm
                       tj->sterile(),
                       tj->region(),
                       parameters_[index],
                       tj,
                       this);

  shadow_jobs_[job_id] = sj;
  shadow_job = sj;
  return true;
}

ShadowJobEntryMap ComplexJobEntry::shadow_jobs() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  if (!shadow_jobs_complete_) {
    CompleteShadowJobs();
  }

  return shadow_jobs_;
}

const ShadowJobEntryMap* ComplexJobEntry::shadow_jobs_p() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  if (!shadow_jobs_complete_) {
    CompleteShadowJobs();
  }

  return &shadow_jobs_;
}

void ComplexJobEntry::CompleteShadowJobs() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  if (shadow_jobs_complete_) {
    return;
  }

  size_t index = 0;
  for (; index < inner_job_ids_.size(); ++index) {
    ShadowJobEntry* shadow_job;
    OMIDGetShadowJobEntryByIndex(index, shadow_job);
  }

  shadow_jobs_complete_ = true;
  return;
}

bool ComplexJobEntry::DrainedAllShadowJobsForAssignment() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  return drained_all_;
}

bool ComplexJobEntry::ShadowJobContained(job_id_t job_id) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  return (shadow_job_ids_.count(job_id) != 0);
}


bool ComplexJobEntry::ShadowJobSterile(job_id_t job_id) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  assert(shadow_job_ids_.count(job_id) == 1);

  SetParentJobIds();

  std::list<job_id_t>::iterator iter = parent_job_ids_.begin();
  for (; iter != parent_job_ids_.end(); ++iter) {
    if (*iter == job_id) {
      return false;
    }
  }

  return true;
}

void ComplexJobEntry::MarkShadowJobAssigned(job_id_t job_id) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  assigned_shadow_job_ids_.insert(job_id);
}

void ComplexJobEntry::MarkShadowJobDone(job_id_t job_id) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  done_shadow_job_ids_.insert(job_id);
}

bool ComplexJobEntry::ShadowJobAssigned(job_id_t job_id) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  assert(shadow_job_ids_.count(job_id) == 1);

  return (assigned_shadow_job_ids_.count(job_id) != 0);
}

bool ComplexJobEntry::ShadowJobDone(job_id_t job_id) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  assert(shadow_job_ids_.count(job_id) == 1);

  return (done_shadow_job_ids_.count(job_id) != 0);
}

bool ComplexJobEntry::AllShadowJobsAssigned() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  return (assigned_shadow_job_ids_.size() == shadow_job_ids_.size());
}

bool ComplexJobEntry::AllShadowJobsDone() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  return (done_shadow_job_ids_.size() == shadow_job_ids_.size());
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

size_t ComplexJobEntry::Cursor::pivot() {
  return pivot_;
}

void ComplexJobEntry::Cursor::set_state(State state) {
  state_ = state;
}

void ComplexJobEntry::Cursor::set_index(size_t index) {
  index_ = index;
}

void ComplexJobEntry::Cursor::set_pivot(size_t pivot) {
  pivot_ = pivot;
}







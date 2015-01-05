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


size_t ComplexJobEntry::GetParentJobIds(std::list<job_id_t>* list) {
  std::list<size_t> indices;
  if (!template_entry_->GetParentJobIndices(&indices)) {
    assert(false);
    return 0;
  }

  size_t count = 0;
  list->clear();
  std::list<size_t>::iterator iter = indices.begin();
  for (; iter != indices.end(); ++iter) {
    list->push_back(inner_job_ids_[*iter]);
    ++count;
  }

  return count;
}

size_t ComplexJobEntry::GetJobsForAssignment(JobEntryList* list, size_t max_num, bool append) {
  size_t count = 0;
  if (!append) {
    list->clear();
  }

  size_t index = assign_index_;
  for (; (index < inner_job_ids_.size()) && (count < max_num); ++index) {
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
    list->push_back(shadow_job);
    count++;
  }

  assign_index_ += count;
  return count;
}


bool ComplexJobEntry::DrainedAllJobsForAssignment() {
  return (assign_index_ == inner_job_ids_.size());
}

void ComplexJobEntry::MarkJobAssigned(job_id_t job_id) {
  // TODO(omidm): Implement
}

void ComplexJobEntry::MarkJobDone(job_id_t job_id) {
  // TODO(omidm): Implement
}



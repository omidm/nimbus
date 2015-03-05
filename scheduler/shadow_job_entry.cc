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
  * Shadow job entry is a version of job entry that contanis only required
  * information to assign the job properly. It includes version maps and
  * read/write set. This data id=s read from the main job entry and shared
  * among multiple shadoes. This way the creation of shadow jobs is very cost
  * effective, cause only a couple od pointers are passed to each instance.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "scheduler/shadow_job_entry.h"
#include "scheduler/complex_job_entry.h"

using namespace nimbus; // NOLINT


ShadowJobEntry::ShadowJobEntry() {
  job_type_ = JOB_SHDW;
  complex_job_ = NULL;
}

ShadowJobEntry::ShadowJobEntry(const job_id_t& job_id) {
  job_type_ = JOB_SHDW;
  job_id_ = job_id;
  complex_job_ = NULL;
}

ShadowJobEntry::ShadowJobEntry(const std::string& job_name,
                               const job_id_t& job_id,
                               const IDSet<logical_data_id_t>* read_set_p,
                               const IDSet<logical_data_id_t>* write_set_p,
                               const IDSet<logical_data_id_t>* union_set_p,
                               const IDSet<job_id_t> before_set,
                               boost::shared_ptr<VersionMap> vmap_read_diff,
                               boost::shared_ptr<VersionMap> vmap_write_diff,
                               const job_id_t& parent_job_id,
                               const job_id_t& future_job_id,
                               const bool& sterile,
                               const GeometricRegion& region,
                               const Parameter& params,
                               TemplateJobEntry* template_job,
                               ComplexJobEntry* complex_job) {
  job_type_ = JOB_SHDW;
  job_name_ = job_name;
  job_id_ = job_id;
  read_set_p_ = read_set_p;
  write_set_p_ = write_set_p;
  union_set_p_ = union_set_p;
  before_set_ = before_set;
  vmap_read_diff_ = vmap_read_diff;
  vmap_write_diff_ = vmap_write_diff;
  parent_job_id_ = parent_job_id;
  future_job_id_ = future_job_id;
  sterile_ = sterile;
  region_ = region;
  params_ = params;
  template_job_ = template_job;
  complex_job_ = complex_job;
}

ShadowJobEntry::~ShadowJobEntry() {
}


const IDSet<logical_data_id_t>* ShadowJobEntry::read_set_p() const {
  return read_set_p_;
}

const IDSet<logical_data_id_t>* ShadowJobEntry::write_set_p() const {
  return write_set_p_;
}

const IDSet<logical_data_id_t>* ShadowJobEntry::union_set_p() const {
  return union_set_p_;
}

boost::shared_ptr<VersionMap> ShadowJobEntry::vmap_read_diff() const {
  return vmap_read_diff_;
}

boost::shared_ptr<VersionMap> ShadowJobEntry::vmap_write_diff() const {
  return vmap_write_diff_;
}

ComplexJobEntry* ShadowJobEntry::complex_job() {
  return complex_job_;
}

TemplateJobEntry* ShadowJobEntry::template_job() {
  return template_job_;
}

void ShadowJobEntry::set_read_set_p(const IDSet<logical_data_id_t>* read_set_p) {
  read_set_p_ = read_set_p;
}

void ShadowJobEntry::set_write_set_p(const IDSet<logical_data_id_t>* write_set_p) {
  write_set_p_ = write_set_p;
}


void ShadowJobEntry::set_union_set_p(const IDSet<logical_data_id_t>* union_set_p) {
  union_set_p_ = union_set_p;
}

void ShadowJobEntry::set_vmap_read_diff(boost::shared_ptr<VersionMap> vmap_read_diff) {
  vmap_read_diff_ = vmap_read_diff;
}

void ShadowJobEntry::set_vmap_write_diff(boost::shared_ptr<VersionMap> vmap_write_diff) {
  vmap_write_diff_ = vmap_write_diff;
}

void ShadowJobEntry::set_complex_job(ComplexJobEntry* complex_job) {
  complex_job_ = complex_job;
}

void ShadowJobEntry::set_template_job(TemplateJobEntry* template_job) {
  template_job_ = template_job;
}

bool ShadowJobEntry::IsReadyForCompleteVersioning() {
  assert(complex_job_);
  return complex_job_->IsReadyForCompleteVersioning();
}

job_depth_t ShadowJobEntry::job_depth() const {
  return complex_job_->job_depth();
}

bool ShadowJobEntry::GetPhysicalReadSet(IDSet<physical_data_id_t>* set) {
  set->clear();
  IDSet<logical_data_id_t>::IDSetIter it;
  for (it = read_set_p_->begin(); it != read_set_p_->end(); ++it) {
    if (physical_table_.count(*it) == 0)
      return false;
    set->insert(physical_table_[*it]);
  }
  return true;
}

bool ShadowJobEntry::GetPhysicalWriteSet(IDSet<physical_data_id_t>* set) {
  set->clear();
  IDSet<logical_data_id_t>::IDSetIter it;
  for (it = write_set_p_->begin(); it != write_set_p_->end(); ++it) {
    if (physical_table_.count(*it) == 0)
      return false;
    set->insert(physical_table_[*it]);
  }
  return true;
}

bool ShadowJobEntry::LookUpMetaBeforeSet(JobEntry* job) {
  // if query job is not ShadowJob then we assume it is in before set path -omidm
  if (job->job_type() != JOB_SHDW) {
    return true;
  }
  ShadowJobEntry* sj = reinterpret_cast<ShadowJobEntry*>(job);

  // if query job is not in the same complex job then we assume it is in before set path -omidm
  if (sj->complex_job() != complex_job_) {
    return true;
  }

  TemplateJobEntry* tj = sj->template_job();

  return template_job_->meta_before_set()->LookUpBeforeSetChain(tj->job_id(), tj->job_depth());
}



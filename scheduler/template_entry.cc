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
  * This is TemplateEntry module to hold and instantiate the templates.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "scheduler/template_entry.h"

using namespace nimbus; // NOLINT


TemplateEntry::TemplateEntry() {
  finalized_ = false;
  // TODO(omidm): currently we do not support future job id in templates!
  future_job_id_ptr_ = boost::shared_ptr<job_id_t>(new job_id_t(0));
}

TemplateEntry::~TemplateEntry() {
  // Used shared ptr for allocated pointers. -omidm
}

bool TemplateEntry::finalized() {
  return finalized_;
}

boost::shared_ptr<VersionMap> TemplateEntry::vmap_base() const {
  return vmap_base_;
}

void TemplateEntry::set_vmap_base(boost::shared_ptr<VersionMap> vmap_base) {
  vmap_base_ = vmap_base;
}

bool TemplateEntry::Finalize() {
  if (finalized_) {
    dbg(DBG_WARN, "WARNING: template has been already finalized!\n");
    return true;
  }

  if (job_id_ptrs_map_.size() != job_id_ptrs_.size()) {
    dbg(DBG_ERROR, "ERROR: referenced jobs are not equal to defined jobs!\n");
    return false;
  }

  size_t index = 0;
  parent_job_indices_.clear();
  EntryList::iterator iter = entry_list_.begin();
  for (; iter != entry_list_.end(); ++iter) {
    if (!iter->sterile_) {
      parent_job_indices_.push_back(index);
    }
    ++index;
  }

  finalized_ = true;
  return true;
}

bool TemplateEntry::CleanPartiallyFilledTemplate() {
  if (finalized_) {
    dbg(DBG_ERROR, "ERROR: template has been finalized and cannot get cleaned!\n");
    return false;
  }

  entry_list_.clear();
  compute_jobs_.clear();
  job_id_ptrs_.clear();
  job_id_ptrs_map_.clear();
  return true;
}

bool TemplateEntry::Instantiate(JobManager *job_manager,
                                const std::vector<job_id_t>& inner_job_ids,
                                const std::vector<job_id_t>& outer_job_ids,
                                const std::vector<Parameter>& parameters,
                                const job_id_t& parent_job_id) {
  if (!finalized_) {
    dbg(DBG_ERROR, "ERROR: template has NOT been finalized and cannot get instantiated!\n");
    return false;
  }

  assert(entry_list_.size() == job_id_ptrs_.size());

  if (inner_job_ids.size() != job_id_ptrs_.size()) {
    dbg(DBG_ERROR, "ERROR: number of provided ids does not match the required ids!\n");
    return false;
  }

  if (parameters.size() != job_id_ptrs_.size()) {
    dbg(DBG_ERROR, "ERROR: number of provided parameters does not match the required ids!\n");
    return false;
  }


  // Set the job_id pointers to the new values.
  size_t index = 0;
  PtrList::iterator piter = job_id_ptrs_.begin();
  for (; piter != job_id_ptrs_.end(); ++piter) {
    *(*piter) = inner_job_ids[index];
    ++index;
  }

  index = 0;
  EntryList::iterator iter = entry_list_.begin();
  for (; iter != entry_list_.end(); ++iter) {
    IDSet<job_id_t> before_set;
    {
      // TODO(omidm) Does accesing a field in class make a copy?
      PtrSet::iterator it = iter->before_set_ptrs_.begin();
      for (; it != iter->before_set_ptrs_.end(); ++it) {
        before_set.insert(*(*it));
      }
    }

    IDSet<job_id_t> after_set;
    {
      // TODO(omidm) Does accesing a field in class make a copy?
      PtrSet::iterator it = iter->after_set_ptrs_.begin();
      for (; it != iter->after_set_ptrs_.end(); ++it) {
        after_set.insert(*(*it));
      }
    }

    job_manager->AddComputeJobEntry(iter->job_name_,
                                    *(job_id_ptrs_[index]),
                                    iter->read_set_,
                                    iter->write_set_,
                                    before_set,
                                    after_set,
                                    parent_job_id,
                                    *(iter->future_job_id_ptr_),
                                    iter->sterile_,
                                    iter->region_,
                                    parameters[index]);
    ++index;
  }

  return true;
}

size_t TemplateEntry::GetParentJobIndices(std::list<size_t>* list) {
  assert(finalized_);
  list->clear();
  *list = parent_job_indices_;
  return list->size();
}


TemplateJobEntry* TemplateEntry::GetJobAtIndex(size_t index) {
  if (index >= compute_jobs_.size()) {
    return NULL;
  }

  return compute_jobs_[index];
}


bool TemplateEntry::GetComplexJobEntry(ComplexJobEntry*& complex_job,
                                       const job_id_t& job_id,
                                       const job_id_t& parent_job_id,
                                       const std::vector<job_id_t>& inner_job_ids,
                                       const std::vector<job_id_t>& outer_job_ids,
                                       const std::vector<Parameter>& parameters) {
  if (!finalized_) {
    dbg(DBG_ERROR, "ERROR: template has NOT been finalized and cannot get instantiated!\n");
    return false;
  }

  assert(entry_list_.size() == job_id_ptrs_.size());

  if (inner_job_ids.size() != job_id_ptrs_.size()) {
    dbg(DBG_ERROR, "ERROR: number of provided ids does not match the required ids!\n");
    return false;
  }

  if (parameters.size() != job_id_ptrs_.size()) {
    dbg(DBG_ERROR, "ERROR: number of provided parameters does not match the required ids!\n");
    return false;
  }

  complex_job = new ComplexJobEntry(job_id,
                                    parent_job_id,
                                    this,
                                    inner_job_ids,
                                    outer_job_ids,
                                    parameters);

  return true;
}

TemplateJobEntry* TemplateEntry::AddComputeJob(const std::string& job_name,
                                               const job_id_t& job_id,
                                               const IDSet<logical_data_id_t>& read_set,
                                               const IDSet<logical_data_id_t>& write_set,
                                               const IDSet<job_id_t>& before_set,
                                               const IDSet<job_id_t>& after_set,
                                               const job_id_t& parent_job_id,
                                               const job_id_t& future_job_id,
                                               const bool& sterile,
                                               const GeometricRegion& region) {
  if (finalized_) {
    dbg(DBG_ERROR, "ERROR: template has been finalized and cannot add compute job!\n");
    return NULL;
  }

  boost::shared_ptr<job_id_t> job_id_ptr;
  {
    PtrMap::iterator iter = job_id_ptrs_map_.find(job_id);
    if (iter == job_id_ptrs_map_.end()) {
      job_id_ptr = boost::shared_ptr<job_id_t>(new job_id_t(job_id));
      job_id_ptrs_map_[job_id] = job_id_ptr;
    } else {
      job_id_ptr = iter->second;
    }
  }

  PtrSet before_set_ptrs;
  {
    IDSet<job_id_t>::IDSetIter it = before_set.begin();
    for (; it != before_set.end(); ++it) {
      boost::shared_ptr<job_id_t> ptr;
      PtrMap::iterator iter = job_id_ptrs_map_.find(*it);
      if (iter == job_id_ptrs_map_.end()) {
        ptr = boost::shared_ptr<job_id_t>(new job_id_t(*it));
        job_id_ptrs_map_[*it] = ptr;
      } else {
        ptr = iter->second;
      }
      before_set_ptrs.insert(ptr);
    }
  }

  PtrSet after_set_ptrs;
  {
    IDSet<job_id_t>::IDSetIter it = after_set.begin();
    for (; it != after_set.end(); ++it) {
      boost::shared_ptr<job_id_t> ptr;
      PtrMap::iterator iter = job_id_ptrs_map_.find(*it);
      if (iter == job_id_ptrs_map_.end()) {
        ptr = boost::shared_ptr<job_id_t>(new job_id_t(*it));
        job_id_ptrs_map_[*it] = ptr;
      } else {
        ptr = iter->second;
      }
      after_set_ptrs.insert(ptr);
    }
  }

  TemplateJobEntry *job =
    new TemplateJobEntry(job_name,
                         read_set,
                         write_set,
                         sterile,
                         region);

  compute_jobs_.push_back(job);

  TemplateComputeJobEntry entry(job_name,
                                job_id_ptr,
                                read_set,
                                write_set,
                                before_set_ptrs,
                                after_set_ptrs,
                                future_job_id_ptr_,
                                sterile,
                                region);

  entry_list_.push_back(entry);
  job_id_ptrs_.push_back(job_id_ptr);

  return job;
}

bool TemplateEntry::AddExplicitCopyJob() {
  dbg(DBG_ERROR, "ERROR: explicit copy jobs from application are not supported yet!.\n");
  exit(-1);
  return false;
}





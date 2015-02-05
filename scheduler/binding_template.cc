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
  * This is BindingTemplate module to hold and instantiate copy and compute
  * commands sent to workers to bind physical data to a batch of jobs in a
  * template complex job.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "scheduler/binding_template.h"
#include "scheduler/template_entry.h"
#include "scheduler/job_manager.h"

using namespace nimbus; // NOLINT

BindingTemplate::BindingTemplate(TemplateEntry *template_entry) {
  finalized_ = false;
  template_entry_ = template_entry;
  future_job_id_ptr_ = JobIdPtr(new job_id_t(0));
}

BindingTemplate::~BindingTemplate() {
}

bool BindingTemplate::finalized() {
  return finalized_;
}

size_t BindingTemplate::copy_job_num() {
  assert(finalized_);
  return copy_job_id_list_.size();
}

size_t BindingTemplate::compute_job_num() {
  assert(finalized_);
  return compute_job_id_list_.size();
}

bool BindingTemplate::Finalize() {
  assert(!finalized_);

  assert(phy_id_map_.size() == phy_id_list_.size());

  assert(copy_job_id_map_.size() == copy_job_id_list_.size());
  assert(compute_job_id_map_.size() == compute_job_id_list_.size());
  assert(compute_job_id_map_.size() == template_entry_->compute_jobs_num());

  assert(end_pattern_map_.size() == end_pattern_list_.size());
  assert(entry_pattern_map_.size() == entry_pattern_list_.size());
  assert(entry_pattern_map_.size() == end_pattern_map_.size());

  PatternList::iterator iter = end_pattern_list_.begin();
  for (; iter != end_pattern_list_.end(); ++iter) {
    assert((*iter)->version_type_ == REGULAR);
  }

  finalized_ = true;
  return true;
}


bool BindingTemplate::Instantiate(const std::vector<job_id_t>& compute_job_ids,
                                  const std::vector<job_id_t>& copy_job_ids,
                                  const std::vector<physical_data_id_t> physical_ids,
                                  SchedulerServer *server) {
  assert(compute_job_ids.size() == compute_job_id_list_.size());
  assert(copy_job_ids.size() == copy_job_id_list_.size());
  assert(physical_ids.size() == phy_id_list_.size());

  assert(false);
  return true;
}

bool BindingTemplate::TrackDataObject(const worker_id_t& worker_id,
                                      const logical_data_id_t& ldid,
                                      const physical_data_id_t& pdid,
                                      VERSION_TYPE version_type,
                                      data_version_t version_diff_from_base) {
  PhyIdPtrMap::iterator iter =  phy_id_map_.find(pdid);
  if (iter != phy_id_map_.end()) {
    return true;
  }

  PhyIdPtr pdid_ptr = PhyIdPtr(new physical_data_id_t(pdid));
  phy_id_map_[pdid] = pdid_ptr;
  phy_id_list_.push_back(pdid_ptr);

  {
    PatternEntry *pattern =
      new PatternEntry(worker_id, ldid, version_type, version_diff_from_base);
    entry_pattern_map_[pdid] = pattern;
    entry_pattern_list_.push_back(pattern);
  }

  {
    PatternEntry *pattern =
      new PatternEntry(worker_id, ldid, version_type, version_diff_from_base);
    end_pattern_map_[pdid] = pattern;
    end_pattern_list_.push_back(pattern);
  }

  return true;
}

bool BindingTemplate::UpdateDataObject(const physical_data_id_t& pdid,
                                       data_version_t version_diff_from_base) {
  PatternMap::iterator iter =  end_pattern_map_.find(pdid);
  if (iter == end_pattern_map_.end()) {
    assert(false);
    return false;
  }

  PatternEntry *pe = iter->second;
  pe->version_type_ = REGULAR;
  pe->version_diff_from_base_ = version_diff_from_base;

  return true;
}

bool BindingTemplate::AddComputeJobCommand(ComputeJobCommand* command,
                                           worker_id_t w_id) {
  JobIdPtr job_id_ptr = GetComputeJobIdPtr(command->job_id().elem());

  PhyIdPtrSet read_set;
  {
    IDSet<physical_data_id_t>::IDSetIter iter = command->read_set_p()->begin();
    for (; iter != command->read_set_p()->end(); ++iter) {
      PhyIdPtr phy_id_ptr = GetExistingPhyIdPtr(*iter);
      read_set.insert(phy_id_ptr);
    }
  }

  PhyIdPtrSet write_set;
  {
    IDSet<physical_data_id_t>::IDSetIter iter = command->write_set_p()->begin();
    for (; iter != command->write_set_p()->end(); ++iter) {
      PhyIdPtr phy_id_ptr = GetExistingPhyIdPtr(*iter);
      write_set.insert(phy_id_ptr);
    }
  }

  JobIdPtrSet before_set;
  {
    IDSet<job_id_t>::IDSetIter iter = command->before_set_p()->begin();
    for (; iter != command->before_set_p()->end(); ++iter) {
      JobIdPtr job_id_ptr;
      if (IDMaker::SchedulerProducedJobID(*iter)) {
        if (!GetCopyJobIdPtrIfExisted(*iter, &job_id_ptr)) {
          continue;
        }
      } else {
        job_id_ptr = GetComputeJobIdPtr(*iter);
      }
      before_set.insert(job_id_ptr);
    }
  }

  JobIdPtrSet after_set;
  {
    IDSet<job_id_t>::IDSetIter iter = command->after_set_p()->begin();
    for (; iter != command->after_set_p()->end(); ++iter) {
      JobIdPtr job_id_ptr;
      if (IDMaker::SchedulerProducedJobID(*iter)) {
        if (!GetCopyJobIdPtrIfExisted(*iter, &job_id_ptr)) {
          continue;
        }
      } else {
        job_id_ptr = GetComputeJobIdPtr(*iter);
      }
      before_set.insert(job_id_ptr);
    }
  }

  ComputeJobCommandTemplate *cm =
    new ComputeJobCommandTemplate(command->job_name(),
                                  job_id_ptr,
                                  read_set,
                                  write_set,
                                  before_set,
                                  after_set,
                                  future_job_id_ptr_,
                                  command->sterile(),
                                  command->region(),
                                  w_id);

  command_templates_.push_back(cm);

  return true;
}

bool BindingTemplate::AddLocalCopyCommand(LocalCopyCommand* command,
                                          worker_id_t w_id) {
  JobIdPtr job_id_ptr = GetCopyJobIdPtr(command->job_id().elem());

  PhyIdPtr from_physical_data_id_ptr =
    GetExistingPhyIdPtr(command->from_physical_data_id().elem());

  PhyIdPtr to_physical_data_id_ptr =
    GetExistingPhyIdPtr(command->to_physical_data_id().elem());

  JobIdPtrSet before_set;
  {
    IDSet<job_id_t>::IDSetIter iter = command->before_set_p()->begin();
    for (; iter != command->before_set_p()->end(); ++iter) {
      JobIdPtr job_id_ptr;
      if (IDMaker::SchedulerProducedJobID(*iter)) {
        if (!GetCopyJobIdPtrIfExisted(*iter, &job_id_ptr)) {
          continue;
        }
      } else {
        job_id_ptr = GetComputeJobIdPtr(*iter);
      }
      before_set.insert(job_id_ptr);
    }
  }

  LocalCopyCommandTemplate *cm =
    new LocalCopyCommandTemplate(job_id_ptr,
                                 from_physical_data_id_ptr,
                                 to_physical_data_id_ptr,
                                 before_set,
                                 w_id);

  command_templates_.push_back(cm);

  return true;
}

bool BindingTemplate::AddRemoteCopySendCommand(RemoteCopySendCommand* command,
                                               worker_id_t w_id) {
  JobIdPtr job_id_ptr = GetCopyJobIdPtr(command->job_id().elem());

  JobIdPtr receive_job_id_ptr = GetCopyJobIdPtr(command->receive_job_id().elem());

  PhyIdPtr from_physical_data_id_ptr =
    GetExistingPhyIdPtr(command->from_physical_data_id().elem());

  JobIdPtrSet before_set;
  {
    IDSet<job_id_t>::IDSetIter iter = command->before_set_p()->begin();
    for (; iter != command->before_set_p()->end(); ++iter) {
      JobIdPtr job_id_ptr;
      if (IDMaker::SchedulerProducedJobID(*iter)) {
        if (!GetCopyJobIdPtrIfExisted(*iter, &job_id_ptr)) {
          continue;
        }
      } else {
        job_id_ptr = GetComputeJobIdPtr(*iter);
      }
      before_set.insert(job_id_ptr);
    }
  }

  RemoteCopySendCommandTemplate *cm =
    new RemoteCopySendCommandTemplate(job_id_ptr,
                                      receive_job_id_ptr,
                                      from_physical_data_id_ptr,
                                      command->to_worker_id(),
                                      command->to_ip(),
                                      command->to_port(),
                                      before_set,
                                      w_id);

  command_templates_.push_back(cm);

  return true;
}

bool BindingTemplate::AddRemoteCopyReceiveCommand(RemoteCopyReceiveCommand* command,
                                                  worker_id_t w_id) {
  JobIdPtr job_id_ptr = GetCopyJobIdPtr(command->job_id().elem());

  PhyIdPtr to_physical_data_id_ptr =
    GetExistingPhyIdPtr(command->to_physical_data_id().elem());

  JobIdPtrSet before_set;
  {
    IDSet<job_id_t>::IDSetIter iter = command->before_set_p()->begin();
    for (; iter != command->before_set_p()->end(); ++iter) {
      JobIdPtr job_id_ptr;
      if (IDMaker::SchedulerProducedJobID(*iter)) {
        if (!GetCopyJobIdPtrIfExisted(*iter, &job_id_ptr)) {
          continue;
        }
      } else {
        job_id_ptr = GetComputeJobIdPtr(*iter);
      }
      before_set.insert(job_id_ptr);
    }
  }

  RemoteCopyReceiveCommandTemplate *cm =
    new RemoteCopyReceiveCommandTemplate(job_id_ptr,
                                         to_physical_data_id_ptr,
                                         before_set,
                                         w_id);

  command_templates_.push_back(cm);

  return true;
}

BindingTemplate::JobIdPtr BindingTemplate::GetCopyJobIdPtr(job_id_t job_id) {
  JobIdPtr job_id_ptr;

  JobIdPtrMap::iterator iter = copy_job_id_map_.find(job_id);
  if (iter != copy_job_id_map_.end()) {
    job_id_ptr = iter->second;
  } else {
    job_id_ptr = JobIdPtr(new job_id_t(job_id));
    copy_job_id_map_[job_id] = job_id_ptr;
    copy_job_id_list_.push_back(job_id_ptr);
  }

  return job_id_ptr;
}

bool BindingTemplate::GetCopyJobIdPtrIfExisted(job_id_t job_id, JobIdPtr *job_id_ptr) {
  JobIdPtrMap::iterator iter = copy_job_id_map_.find(job_id);
  if (iter != copy_job_id_map_.end()) {
    *job_id_ptr = iter->second;
    return true;
  }

  return false;
}

BindingTemplate::JobIdPtr BindingTemplate::GetComputeJobIdPtr(job_id_t job_id) {
  JobIdPtr job_id_ptr;

  JobIdPtrMap::iterator iter = compute_job_id_map_.find(job_id);
  if (iter != compute_job_id_map_.end()) {
    job_id_ptr = iter->second;
  } else {
    job_id_ptr = JobIdPtr(new job_id_t(job_id));
    compute_job_id_map_[job_id] = job_id_ptr;
    compute_job_id_list_.push_back(job_id_ptr);
  }

  return job_id_ptr;
}

BindingTemplate::PhyIdPtr BindingTemplate::GetExistingPhyIdPtr(physical_data_id_t pdid) {
  PhyIdPtr phy_id_ptr;

  PhyIdPtrMap::iterator iter = phy_id_map_.find(pdid);
  if (iter != phy_id_map_.end()) {
    phy_id_ptr = iter->second;
  } else {
    assert(false);
  }

  return phy_id_ptr;
}



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

const BindingTemplate::PatternMetaData* BindingTemplate::patterns_meta_data_p() const {
  return &patterns_meta_data_;
}

const BindingTemplate::PatternList* BindingTemplate::entry_pattern_list_p() const {
  return &entry_pattern_list_;
}

const BindingTemplate::PatternList* BindingTemplate::end_pattern_list_p() const {
  return &end_pattern_list_;
}

bool BindingTemplate::Finalize(const std::vector<job_id_t>& compute_job_ids) {
  assert(!finalized_);
  assert(compute_job_id_map_.size() == template_entry_->compute_jobs_num());
  assert(compute_job_id_map_.size() == compute_job_ids.size());

  compute_job_id_list_.clear();
  {
    std::vector<job_id_t>::const_iterator iter = compute_job_ids.begin();
    for (; iter != compute_job_ids.end(); ++iter) {
      JobIdPtrMap::iterator it = compute_job_id_map_.find(*iter);
      assert(it != compute_job_id_map_.end());
      compute_job_id_list_.push_back(it->second);
    }
  }

  copy_job_id_list_.clear();
  {
    JobIdPtrMap::iterator iter = copy_job_id_map_.begin();
    for (; iter != copy_job_id_map_.end(); ++iter) {
      copy_job_id_list_.push_back(iter->second);
    }
  }

  phy_id_list_.clear();
  end_pattern_list_.clear();
  entry_pattern_list_.clear();
  {
    PatternSorted::iterator iter = ordered_entry_patterns_.begin();
    for (; iter != ordered_entry_patterns_.end(); ++iter) {
      size_t regular_num = 0;
      size_t wild_card_num = 0;
      {
        PatternList::iterator it = iter->second.first->begin();
        for (; it != iter->second.first->end(); ++it) {
          ++regular_num;
          physical_data_id_t pdid = (*it)->pdid_;
          {
            PhyIdPtrMap::iterator i = phy_id_map_.find(pdid);
            assert(i != phy_id_map_.end());
            phy_id_list_.push_back(i->second);
          }
          {
            PatternMap::iterator i = end_pattern_map_.find(pdid);
            assert(i != end_pattern_map_.end());
            end_pattern_list_.push_back(i->second);
          }
          {
            PatternMap::iterator i = entry_pattern_map_.find(pdid);
            assert(i != entry_pattern_map_.end());
            entry_pattern_list_.push_back(i->second);
            assert(i->second == (*it));
          }
        }
      }
      {
        PatternList::iterator it = iter->second.second->begin();
        for (; it != iter->second.second->end(); ++it) {
          ++wild_card_num;
          physical_data_id_t pdid = (*it)->pdid_;
          {
            PhyIdPtrMap::iterator i = phy_id_map_.find(pdid);
            assert(i != phy_id_map_.end());
            phy_id_list_.push_back(i->second);
          }
          {
            PatternMap::iterator i = end_pattern_map_.find(pdid);
            assert(i != end_pattern_map_.end());
            end_pattern_list_.push_back(i->second);
          }
          {
            PatternMap::iterator i = entry_pattern_map_.find(pdid);
            assert(i != entry_pattern_map_.end());
            entry_pattern_list_.push_back(i->second);
            assert(i->second == (*it));
          }
        }
      }

      patterns_meta_data_.push_back(std::make_pair(regular_num, wild_card_num));
    }
  }

  {
    PatternList::iterator iter = end_pattern_list_.begin();
    for (; iter != end_pattern_list_.end(); ++iter) {
      assert((*iter)->version_type_ == REGULAR);
    }
  }

  assert(copy_job_id_map_.size() == copy_job_id_list_.size());
  assert(compute_job_id_map_.size() == compute_job_id_list_.size());
  assert(phy_id_map_.size() == phy_id_list_.size());
  assert(phy_id_map_.size() == end_pattern_map_.size());
  assert(phy_id_map_.size() == entry_pattern_map_.size());
  assert(end_pattern_map_.size() == end_pattern_list_.size());
  assert(entry_pattern_map_.size() == entry_pattern_list_.size());
  assert(ordered_entry_patterns_.size() == patterns_meta_data_.size());

  finalized_ = true;
  return true;
}

bool BindingTemplate::Instantiate(const std::vector<job_id_t>& compute_job_ids,
                                  const std::vector<job_id_t>& copy_job_ids,
                                  const std::vector<physical_data_id_t> physical_ids,
                                  SchedulerServer *server) {
  assert(finalized_);
  assert(compute_job_ids.size() == compute_job_id_list_.size());
  assert(copy_job_ids.size() == copy_job_id_list_.size());
  assert(physical_ids.size() == phy_id_list_.size());

  {
    size_t idx = 0;
    JobIdPtrList::iterator iter = compute_job_id_list_.begin();
    for (; iter != compute_job_id_list_.end(); ++iter) {
      *(*iter) = compute_job_ids[idx];
      ++idx;
    }
  }

  {
    size_t idx = 0;
    JobIdPtrList::iterator iter = copy_job_id_list_.begin();
    for (; iter != copy_job_id_list_.end(); ++iter) {
      *(*iter) = copy_job_ids[idx];
      ++idx;
    }
  }

  {
    size_t idx = 0;
    PhyIdPtrList::iterator iter = phy_id_list_.begin();
    for (; iter != phy_id_list_.end(); ++iter) {
      *(*iter) = physical_ids[idx];
      ++idx;
    }
  }







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

  {
    PatternEntry *pattern =
      new PatternEntry(worker_id, ldid, pdid, version_type, version_diff_from_base);
    entry_pattern_map_[pdid] = pattern;

    PatternSorted::iterator iter = ordered_entry_patterns_.find(ldid);
    if (iter != ordered_entry_patterns_.end()) {
      if (version_type == REGULAR) {
        iter->second.first->push_back(pattern);
      } else {
        iter->second.second->push_back(pattern);
      }
    } else {
      PatternList *regular = new PatternList();
      PatternList *wild_card = new PatternList();
      if (version_type == REGULAR) {
        regular->push_back(pattern);
      } else {
        wild_card->push_back(pattern);
      }
      ordered_entry_patterns_[ldid] = std::make_pair(regular, wild_card);
    }
  }

  {
    PatternEntry *pattern =
      new PatternEntry(worker_id, ldid, pdid, version_type, version_diff_from_base);
    end_pattern_map_[pdid] = pattern;
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



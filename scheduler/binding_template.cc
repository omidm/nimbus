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
#include "scheduler/job_manager.h"

using namespace nimbus; // NOLINT

BindingTemplate::BindingTemplate() {
  finalized_ = false;
  future_job_id_ptr_ = JobIdPtr(new job_id_t(0));
}

BindingTemplate::~BindingTemplate() {
}

bool BindingTemplate::TrackDataObject(const logical_data_id_t& ldid,
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

  PatternEntry *pattern =
    new PatternEntry(ldid, version_type, version_diff_from_base);
  entry_pattern_.push_back(pattern);

  return true;
}

bool BindingTemplate::AddComputeJobCommand(ComputeJobCommand* command,
                                           worker_id_t w_id) {
  JobIdPtr job_id_ptr = GetJobIdPtr(command->job_id().elem());

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
      JobIdPtr job_id_ptr = GetJobIdPtr(*iter);
      before_set.insert(job_id_ptr);
    }
  }

  JobIdPtrSet after_set;
  {
    IDSet<job_id_t>::IDSetIter iter = command->after_set_p()->begin();
    for (; iter != command->after_set_p()->end(); ++iter) {
      JobIdPtr job_id_ptr = GetJobIdPtr(*iter);
      after_set.insert(job_id_ptr);
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

  compute_job_commands_.push_back(cm);

  return false;
}

bool BindingTemplate::AddLocalCopyCommand(LocalCopyCommand* command,
                                          worker_id_t w_id) {
  return false;
}

bool BindingTemplate::AddRemoteCopySendCommand(RemoteCopySendCommand* command,
                                               worker_id_t w_id) {
  return false;
}

bool BindingTemplate::AddRemoteCopyReceiveCommand(RemoteCopyReceiveCommand* command,
                                                  worker_id_t w_id) {
  return false;
}

BindingTemplate::JobIdPtr BindingTemplate::GetJobIdPtr(job_id_t job_id) {
  JobIdPtr job_id_ptr;

  JobIdPtrMap::iterator iter = job_id_map_.find(job_id);
  if (iter != job_id_map_.end()) {
    job_id_ptr = iter->second;
  } else {
    job_id_ptr = JobIdPtr(new job_id_t(job_id));
    job_id_map_[job_id] = job_id_ptr;
    job_id_list_.push_back(job_id_ptr);
  }

  return job_id_ptr;
}

BindingTemplate::JobIdPtr BindingTemplate::GetExistingJobIdPtr(job_id_t job_id) {
  JobIdPtr job_id_ptr;

  JobIdPtrMap::iterator iter = job_id_map_.find(job_id);
  if (iter != job_id_map_.end()) {
    job_id_ptr = iter->second;
  } else {
    assert(false);
  }

  return job_id_ptr;
}

BindingTemplate::PhyIdPtr BindingTemplate::GetPhyIdPtr(physical_data_id_t pdid) {
  PhyIdPtr phy_id_ptr;

  PhyIdPtrMap::iterator iter = phy_id_map_.find(pdid);
  if (iter != phy_id_map_.end()) {
    phy_id_ptr = iter->second;
  } else {
    phy_id_ptr = JobIdPtr(new physical_data_id_t(pdid));
    phy_id_map_[pdid] = phy_id_ptr;
    phy_id_list_.push_back(phy_id_ptr);
  }

  return phy_id_ptr;
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






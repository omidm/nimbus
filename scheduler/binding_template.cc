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

  PatternEntry *pattern =
    new PatternEntry(ldid, version_type, version_diff_from_base);

  entry_pattern_.push_back(pattern);

  PhyIdPtr pdid_ptr = PhyIdPtr(new physical_data_id_t(pdid));
  phy_id_map_[pdid] = pdid_ptr;
  phy_id_list_.push_back(pdid_ptr);

  return true;
}

bool BindingTemplate::AddComputeJobCommand(const ComputeJobCommand& command,
                                           worker_id_t w_id) {
  return false;
}

bool BindingTemplate::AddLocalCopyCommand(const LocalCopyCommand& command,
                                          worker_id_t w_id) {
  return false;
}

bool BindingTemplate::AddRemoteCopySendCommand(const RemoteCopySendCommand& command,
                                               worker_id_t w_id) {
  return false;
}

bool BindingTemplate::AddRemoteCopyReceiveCommand(const RemoteCopyReceiveCommand& command,
                                                  worker_id_t w_id) {
  return false;
}




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
  * Scheduler Version Manager. This module serves the job manager by keeping
  * track of the version numbers of each logical data object that are needed by
  * the jobs in the job graph.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "scheduler/version_manager.h"

using namespace nimbus; // NOLINT

VersionManager::VersionManager() {
  ldo_map_p_ = NULL;
}

VersionManager::~VersionManager() {
}

bool VersionManager::AddJobEntry(JobEntry *job) {
  if (job->sterile()) {
    IDSet<logical_data_id_t>::ConstIter it;
    for (it = job->read_set_p()->begin(); it != job->read_set_p()->end(); ++it) {
      IndexIter iter = index_.find(*it);
      if (iter == index_.end()) {
        VersionEntry ve(*it);
        ve.AddJobEntry(job);
        index_.insert(std::make_pair(*it, ve));
      } else {
        iter->second.AddJobEntry(job);
      }
    }
  } else {
    std::map<logical_data_id_t, LogicalDataObject*>::const_iterator it;
    for (it = ldo_map_p_->begin(); it != ldo_map_p_->end(); ++it) {
      IndexIter iter = index_.find(it->first);
      if (iter == index_.end()) {
        VersionEntry ve(it->first);
        ve.AddJobEntry(job);
        index_.insert(std::make_pair(it->first, ve));
      } else {
        iter->second.AddJobEntry(job);
      }
    }
  }

  return true;
}

size_t VersionManager::GetJobsNeedDataVersion(
    JobEntryList* list, VLD vld) {
  IndexIter iter = index_.find(vld.first);
  if (iter == index_.end()) {
    list->clear();
    return 0;
  } else {
    return iter->second.GetJobsNeedVersion(list, vld.second);
  }
}

bool VersionManager::RemoveJobEntry(JobEntry* job) {
  assert(job->versioned());

  VersionMap::ConstIter it = job->vmap_read()->content_p()->begin();
  for (; it != job->vmap_read()->content_p()->end(); ++it) {
    IndexIter iter = index_.find(it->first);
    if (iter == index_.end()) {
      dbg(DBG_ERROR, "Version manager: ldid %lu is not in the index.\n", *it);
      exit(-1);
      return false;
    } else {
      iter->second.RemoveJobEntry(job);
    }
  }

  return true;
}

void VersionManager::set_ldo_map_p(
    const std::map<logical_data_id_t, LogicalDataObject*>* ldo_map_p) {
  ldo_map_p_ = ldo_map_p;
}



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
/*
  if (job->sterile()) {
    log_sterile_.ResumeTimer();
    IDSet<logical_data_id_t>::ConstIter it;
    for (it = job->read_set_p()->begin(); it != job->read_set_p()->end(); ++it) {
      data_version_t version;
      log_lookup_.ResumeTimer();
      lookup_count_++;
      bool found = LookUpVersion(job, *it, &version);
      log_lookup_.StopTimer();
      if (found) {
        job->vmap_read()->set_entry(*it, version);
      }
    }
    log_sterile_.StopTimer();
  } else {
    log_nonsterile_.ResumeTimer();
    boost::shared_ptr<VersionMap> vmap = boost::shared_ptr<VersionMap>(new VersionMap());
    std::map<logical_data_id_t, LogicalDataObject*>::const_iterator it;
    for (it = ldo_map_p_->begin(); it != ldo_map_p_->end(); ++it) {
      data_version_t version;
      log_lookup_.ResumeTimer();
      lookup_count_++;
      bool found = LookUpVersion(job, it->first, &version);
      log_lookup_.StopTimer();
      if (found) {
        job->vmap_read()->set_entry(it->first, version);
      }
    }
    log_nonsterile_.StopTimer();
  }
*/





  return false;
}

size_t VersionManager::GetJobsNeedDataVersion(
    JobEntryList* list, VLD vld) {
  return 0;
}

bool VersionManager::RemoveJobEntry(JobEntry* job) {
  return false;
}


void VersionManager::set_ldo_map_p(
    const std::map<logical_data_id_t, LogicalDataObject*>* ldo_map_p) {
  ldo_map_p_ = ldo_map_p;
}



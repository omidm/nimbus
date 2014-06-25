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
  * Scheduler Version Entry. It holds the meta data for each version of the
  * data containing logical data id, version number, and a pointer to a job
  * related to the version, and the type of relation meaning whether job needs
  * the version as input or dumps it as output.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "scheduler/version_entry.h"

using namespace nimbus; // NOLINT

VersionEntry::VersionEntry(logical_data_id_t ldid) {
  ldid_ = ldid;
}


VersionEntry::VersionEntry(const VersionEntry& other) {
  ldid_ = other.ldid_;
  pending_jobs_ = other.pending_jobs_;
  index_ = other.index_;
}

VersionEntry::~VersionEntry() {
}

VersionEntry& VersionEntry::operator= (const VersionEntry& right) {
  ldid_ = right.ldid_;
  pending_jobs_ = right.pending_jobs_;
  index_ = right.index_;
  return *this;
}


bool VersionEntry::AddJobEntry(JobEntry *job) {
  pending_jobs_.insert(job);
  return true;
}

size_t VersionEntry::GetJobsNeedVersion(
    JobEntryList* list, data_version_t version) {
  size_t count = 0;
  list->clear();

  /*
   * TODO(omidm): for now just assume that when called all jobs are versioned.
   *  1. Partial versioning (future jobs)
   *  2. Versioning right before asignment.
   * causes the following to fail.
   */

  BucketIter iter = pending_jobs_.begin();
  for (; iter != pending_jobs_.end();) {
    assert((*iter)->versioned());
    assert((*iter)->job_name() != "localcopy");
    data_version_t ver;
    if ((*iter)->vmap_read()->query_entry(ldid_, &ver)) {
      index_[ver].insert(*iter);
    } else {
      dbg(DBG_ERROR, "Version Entry: ldid %lu in read set of job %lu is not versioned.\n",
          ldid_, (*iter)->job_id());
      exit(-1);
      return 0;
    }
    pending_jobs_.erase(iter++);
  }

  assert(pending_jobs_.size() == 0);

  IndexIter iiter = index_.find(version);
  if (iiter == index_.end()) {
    return count;
  } else {
    BucketIter it = iiter->second.begin();
    for (; it != iiter->second.end(); ++it) {
      if (!(*it)->assigned()) {
        list->push_back(*it);
        ++count;
      }
    }
  }

  return count;
}

bool VersionEntry::RemoveJobEntry(JobEntry *job) {
  assert(job->versioned());

  pending_jobs_.erase(job);

  data_version_t ver;
  if (job->vmap_read()->query_entry(ldid_, &ver)) {
    index_[ver].erase(job);
    if (index_[ver].size() == 0) {
      index_.erase(ver);
    }
  } else {
    dbg(DBG_ERROR, "Version Entry: ldid %lu is not versioned for job %lu.\n",
        ldid_, job->job_id());
    exit(-1);
    return false;
  }

  return true;
}



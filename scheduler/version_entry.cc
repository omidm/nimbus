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
  ldl_.set_ldid(ldid);
}

VersionEntry::VersionEntry(const VersionEntry& other) {
  ldid_ = other.ldid_;
  pending_reader_jobs_ = other.pending_reader_jobs_;
  pending_writer_jobs_ = other.pending_writer_jobs_;
  index_ = other.index_;
  ldl_ = other.ldl_;
}

VersionEntry::~VersionEntry() {
  IndexIter it = index_.begin();
  for (; it != index_.end(); ++it) {
    delete it->second;
  }
}

VersionEntry& VersionEntry::operator= (const VersionEntry& right) {
  ldid_ = right.ldid_;
  pending_reader_jobs_ = right.pending_reader_jobs_;
  pending_writer_jobs_ = right.pending_writer_jobs_;
  index_ = right.index_;
  ldl_ = right.ldl_;
  return *this;
}

void VersionEntry::InitializeLdl(
    const job_id_t& job_id,
    const job_depth_t& job_depth) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  ldl_.AppendLdlEntry(job_id, NIMBUS_INIT_DATA_VERSION, job_depth, false);
}



bool VersionEntry::AddJobEntryReader(JobEntry *job) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  if (job->sterile()) {
    pending_reader_jobs_.insert(job);
  } else {
    BucketIter iter = seen_parent_jobs_.find(job);
    if (iter == seen_parent_jobs_.end()) {
      seen_parent_jobs_.insert(job);
      pending_reader_jobs_.insert(job);
    }
  }

  return true;
}

bool VersionEntry::AddJobEntryWriter(JobEntry *job) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  pending_writer_jobs_.insert(job);
  return true;
}

size_t VersionEntry::GetJobsNeedVersion(
    JobEntryList* list, data_version_t version) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  if (!UpdateLdl()) {
    dbg(DBG_ERROR, "ERROR: Could not update the ldl for ldid %lu.\n", ldid_);
  }

  size_t count = 0;
  list->clear();

  BucketIter iter = pending_reader_jobs_.begin();
  for (; iter != pending_reader_jobs_.end();) {
    if ((*iter)->assigned()) {
      pending_reader_jobs_.erase(iter++);
      continue;
    }

    data_version_t partial_version;
    if ((*iter)->vmap_partial()->query_entry(ldid_, &partial_version)) {
      IndexIter it = index_.find(partial_version);
      assert(it != index_.end());
      it->second->erase(*iter);
    }

    data_version_t new_version;
    if ((*iter)->vmap_read()->query_entry(ldid_, &new_version)) {
    } else if (LookUpVersion(*iter, &new_version)) {
      if ((*iter)->IsReadyForCompleteVersioning()) {
        (*iter)->vmap_read()->set_entry(ldid_, new_version);
      } else {
        (*iter)->vmap_partial()->set_entry(ldid_, new_version);
      }
    } else {
      dbg(DBG_ERROR, "ERROR: Version Entry: ldid %lu in read set of job %lu could not be not versioned.\n", // NOLINT
          ldid_, (*iter)->job_id());
      exit(-1);
    }

    IndexIter it = index_.find(new_version);
    if (it != index_.end()) {
      it->second->insert(*iter);
    } else {
      Bucket *b = new Bucket();
      b->insert(*iter);
      index_.insert(std::make_pair(new_version, b));
    }


    if ((*iter)->IsReadyForCompleteVersioning()) {
      pending_reader_jobs_.erase(iter++);
    } else {
      ++iter;
    }
  }



  IndexIter iiter = index_.find(version);
  if (iiter == index_.end()) {
    return count;
  } else {
    BucketIter it = iiter->second->begin();
    for (; it != iiter->second->end(); ++it) {
      if ((!(*it)->assigned()) ||
          ((!(*it)->sterile()) && (!(*it)->done()))) {
        list->push_back(*it);
        ++count;
      }
    }
  }

  return count;
}

bool VersionEntry::RemoveJobEntry(JobEntry *job) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  assert(job->versioned());

  seen_parent_jobs_.erase(job);
  pending_reader_jobs_.erase(job);
  pending_writer_jobs_.erase(job);

  {
    data_version_t ver;
    if (job->vmap_read()->query_entry(ldid_, &ver)) {
      IndexIter it = index_.find(ver);
      if (it != index_.end()) {
        it->second->erase(job);
        if (it->second->size() == 0) {
          delete it->second;
          index_.erase(it);
        }
      }
    }
  }

  if (!job->sterile()) {
    data_version_t ver;
    if (job->vmap_partial()->query_entry(ldid_, &ver)) {
      IndexIter it = index_.find(ver);
      if (it != index_.end()) {
        it->second->erase(job);
        if (it->second->size() == 0) {
          delete it->second;
          index_.erase(it);
        }
      }
    }
  }

  return true;
}

bool VersionEntry::LookUpVersion(
    JobEntry *job,
    data_version_t *version) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  if (!UpdateLdl()) {
    dbg(DBG_ERROR, "Could not update the ldl for ldid %lu.\n", ldid_);
  }

  if (job->job_type() == JOB_SHDW) {
    ShadowJobEntry* sj = reinterpret_cast<ShadowJobEntry*>(job);
    ComplexJobEntry* xj = sj->complex_job();

    data_version_t base_version;
    if (xj->vmap_read()->query_entry(ldid_, &base_version)) {
    } else if (LookUpVersion(xj, &base_version)) {
      xj->vmap_read()->set_entry(ldid_, base_version);
    } else {
      dbg(DBG_ERROR, "ERROR: could not version the base complex job %lu.\n", xj->job_id());
      return false;
    }

    data_version_t diff_version;
    if (!sj->vmap_read_diff()->query_entry(ldid_, &diff_version)) {
      dbg(DBG_ERROR, "ERROR: could not get diff version for shadow job %lu.\n", sj->job_id());
      return false;
    }

    *version = base_version + diff_version;
    return true;
  }

  return ldl_.LookUpVersion(job->meta_before_set(), version);
}

bool VersionEntry::is_empty() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  return ((index_.size() == 0) &&
          (seen_parent_jobs_.size() == 0) &&
          (pending_reader_jobs_.size() == 0) &&
          (pending_writer_jobs_.size() == 0));
}


bool CompareJobDepths(JobEntry* i, JobEntry* j) {
  return ((i->job_depth()) < (j->job_depth()));
}

bool VersionEntry::UpdateLdl() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  if (pending_writer_jobs_.size() == 0) {
    return true;
  }

  std::vector<JobEntry*> depth_sorted;

  BucketIter iter = pending_writer_jobs_.begin();
  for (; iter != pending_writer_jobs_.end();) {
    if ((*iter)->IsReadyForCompleteVersioning()) {
      depth_sorted.push_back(*iter);
      pending_writer_jobs_.erase(iter++);
    } else {
      ++iter;
    }
  }

  if (depth_sorted.size() == 0) {
    return true;
  }

  std::sort(depth_sorted.begin(), depth_sorted.end(), CompareJobDepths);

  data_version_t version = ldl_.last_version();
  std::vector<JobEntry*>::iterator it = depth_sorted.begin();
  for (; it != depth_sorted.end(); ++it) {
    if ((*it)->read_set_p()->contains(ldid_)) {
      (*it)->vmap_read()->set_entry(ldid_, version);
    }
    ++version;
    (*it)->vmap_write()->set_entry(ldid_, version);

    if (!ldl_.AppendLdlEntry(
          (*it)->job_id(),
          version,
          (*it)->job_depth(),
          (*it)->sterile())) {
      return false;
    }
  }

  return true;
}



bool VersionEntry::InsertCheckPointLdlEntry(
    const job_id_t& job_id,
    const data_version_t& version,
    const job_depth_t& job_depth) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  return ldl_.InsertCheckpointLdlEntry(job_id, version, job_depth);
}


bool VersionEntry::CleanLdl(const IDSet<job_id_t>& snap_shot) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  return ldl_.CleanChain(snap_shot);
}

void VersionEntry::Reinitialize() {
  seen_parent_jobs_.clear();
  pending_reader_jobs_.clear();
  pending_writer_jobs_.clear();
  ClearIndex();
  ReinitializeLdl();
}

void VersionEntry::ReinitializeLdl() {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  return ldl_.Reinitialize();
}

void VersionEntry::ClearIndex() {
  IndexIter it = index_.begin();
  for (; it != index_.end(); ++it) {
    delete it->second;
  }
  index_.clear();
}


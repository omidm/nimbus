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
  parent_removed_ = false;
}

VersionManager::~VersionManager() {
  IndexIter it = index_.begin();
  for (; it != index_.end(); ++it) {
    delete it->second;
  }
}

bool VersionManager::AddJobEntry(JobEntry *job) {
  if (job->sterile()) {
    IDSet<logical_data_id_t>::ConstIter it;
    for (it = job->read_set_p()->begin(); it != job->read_set_p()->end(); ++it) {
      IndexIter iter = index_.find(*it);
      if (iter == index_.end()) {
        dbg(DBG_ERROR, "ERROR: ldid %lu appeared in read set of %lu is not defined yet.\n",
            *it, job->job_id());
      } else {
        iter->second->AddJobEntryReader(job);
      }
    }
  } else {
    LdoMap::const_iterator it;
    for (it = ldo_map_p_->begin(); it != ldo_map_p_->end(); ++it) {
      IndexIter iter = index_.find(it->first);
      if (iter == index_.end()) {
        dbg(DBG_ERROR, "ERROR: ldid %lu appeared in ldo_map set of %lu is not defined yet.\n",
            *it, job->job_id());
      } else {
        iter->second->AddJobEntryReader(job);
      }
    }
  }

  {
    IDSet<logical_data_id_t>::ConstIter it;
    for (it = job->write_set_p()->begin(); it != job->write_set_p()->end(); ++it) {
      IndexIter iter = index_.find(*it);
      if (iter == index_.end()) {
        dbg(DBG_ERROR, "ERROR: ldid %lu appeared in ldo_map set of %lu is not defined yet.\n",
            *it, job->job_id());
      } else {
        iter->second->AddJobEntryWriter(job);
      }
    }
  }

  ChildCounterIter it = child_counter_.find(job->parent_job_id());
  if (it != child_counter_.end()) {
    ++it->second;
  }

  return true;
}

bool VersionManager::ResolveJobDataVersions(JobEntry *job) {
  assert(job->IsReadyForCompleteVersioning());

  if (job->sterile()) {
    IDSet<logical_data_id_t>::ConstIter it;
    for (it = job->read_set_p()->begin(); it != job->read_set_p()->end(); ++it) {
      data_version_t version;
      if (job->vmap_read()->query_entry(*it, &version)) {
        continue;
      }
      if (LookUpVersion(job, *it, &version)) {
        job->vmap_read()->set_entry(*it, version);
      } else {
        return false;
      }
    }
  } else {
    LdoMap::const_iterator it;
    for (it = ldo_map_p_->begin(); it != ldo_map_p_->end(); ++it) {
      data_version_t version;
      if (job->vmap_read()->query_entry(it->first, &version)) {
        continue;
      }
      if (LookUpVersion(job, it->first, &version)) {
        job->vmap_read()->set_entry(it->first, version);
      } else {
        return false;
      }
    }
  }

  {
    IDSet<logical_data_id_t>::ConstIter it;
    for (it = job->write_set_p()->begin(); it != job->write_set_p()->end(); ++it) {
      data_version_t version;
      if (job->vmap_write()->query_entry(*it, &version)) {
        continue;
      }
      if (job->vmap_read()->query_entry(*it, &version)) {
        job->vmap_write()->set_entry(*it, version + 1);
      } else {
        if (LookUpVersion(job, *it, &version)) {
          job->vmap_write()->set_entry(*it, version + 1);
        } else {
          return false;
        }
      }
    }
  }

  if (!job->sterile()) {
    // Clear meta before set and update the ldl.

    job->meta_before_set()->Clear();

    LdoMap::const_iterator it;
    for (it = ldo_map_p_->begin(); it != ldo_map_p_->end(); ++it) {
      data_version_t version;
      if (job->vmap_write()->query_entry(it->first, &version)) {
      } else if (job->vmap_read()->query_entry(it->first, &version)) {
      } else {
        dbg(DBG_ERROR, "ERROR: Version Manager: ldid %lu is not versioned for parent job %lu.\n",
            it->first, job->job_id());
        exit(-1);
      }

      InsertParentLdlEntry(
          it->first, job->job_id(), version, job->job_depth());
    }

    boost::unique_lock<boost::mutex> lock(child_counter_mutex_);
    ChildCounterIter cit = child_counter_.find(job->job_id());
    if (cit == child_counter_.end()) {
      live_parents_.insert(job->job_id());
      child_counter_[job->job_id()] = (counter_t) (1);
    } else {
      assert(false);
    }
  }

  boost::unique_lock<boost::mutex> lock(child_counter_mutex_);
  ChildCounterIter it = child_counter_.find(job->parent_job_id());
  if (it != child_counter_.end()) {
    --it->second;
    if (it->second == 0) {
      live_parents_.remove(it->first);
      child_counter_.erase(it);
      parent_removed_ = true;
    }
  }

  return true;
}


bool VersionManager::InsertParentLdlEntry(
    const logical_data_id_t ldid,
    const job_id_t& job_id,
    const data_version_t& version,
    const job_depth_t& job_depth) {
  IndexIter iter = index_.find(ldid);
  if (iter == index_.end()) {
    return false;
  } else {
    return iter->second->InsertParentLdlEntry(job_id, version, job_depth);
  }
}


bool VersionManager::LookUpVersion(
    JobEntry *job,
    logical_data_id_t ldid,
    data_version_t *version) {
  IndexIter iter = index_.find(ldid);
  if (iter == index_.end()) {
    return false;
  } else {
    return iter->second->LookUpVersion(job, version);
  }
}



size_t VersionManager::GetJobsNeedDataVersion(
    JobEntryList* list, VLD vld) {
  IndexIter iter = index_.find(vld.first);
  if (iter == index_.end()) {
    list->clear();
    return 0;
  } else {
    return iter->second->GetJobsNeedVersion(list, vld.second);
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
      iter->second->RemoveJobEntry(job);
      if (iter->second->is_empty()) {
        delete iter->second;
        index_.erase(iter);
        std::cout << "version entry got empty!!\n";
      }
    }
  }

  if (!(job->sterile())) {
    ChildCounterIter it = child_counter_.find(job->job_id());
    if (it != child_counter_.end()) {
      --it->second;
      if (it->second == 0) {
        live_parents_.remove(it->first);
        child_counter_.erase(it);
        parent_removed_ = true;
      }
    } else {
      assert(false);
    }
  }

  return true;
}

bool VersionManager::DefineData(
    const logical_data_id_t ldid,
    const job_id_t& job_id,
    const job_depth_t& job_depth) {
  IndexIter iter = index_.find(ldid);
  if (iter == index_.end()) {
    VersionEntry *ve = new VersionEntry(ldid);
    ve->InitializeLdl(job_id, job_depth);
    index_.insert(std::make_pair(ldid, ve));
    return true;
  } else {
    dbg(DBG_ERROR, "ERROR: defining logical data id %lu, which already exist.\n", ldid);
    exit(-1);
    return false;
  }
}

bool VersionManager::CleanUp() {
  if (parent_removed_) {
    Log log;
    log.StartTimer();

    IndexIter iter = index_.begin();
    for (; iter != index_.end(); ++iter) {
      iter->second->CleanLdl(live_parents_);
    }
    parent_removed_ = false;

    log.StopTimer();
    std::cout << "Version Manager Clean up: " << log.timer() << std::endl;
  }

  return true;
}


void VersionManager::set_ldo_map_p(const LdoMap* ldo_map_p) {
  ldo_map_p_ = ldo_map_p;
}



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

#define DEFAULT_SNAP_SHOT_RATE 20;

using namespace nimbus; // NOLINT

VersionManager::VersionManager() {
  ldo_map_p_ = NULL;
  snap_shot_pending_ = false;
  non_sterile_counter_ = 0;
  snap_shot_rate_ = DEFAULT_SNAP_SHOT_RATE;
}

VersionManager::~VersionManager() {
  Index::iterator it = index_.begin();
  for (; it != index_.end(); ++it) {
    delete it->second;
  }
}

void VersionManager::set_snap_shot_rate(size_t rate) {
  snap_shot_rate_ = rate;
}

bool VersionManager::AddJobEntry(JobEntry *job) {
  if (job->sterile()) {
    IDSet<logical_data_id_t>::ConstIter it;
    for (it = job->read_set_p()->begin(); it != job->read_set_p()->end(); ++it) {
      Index::iterator iter = index_.find(*it);
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
      Index::iterator iter = index_.find(it->first);
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
      Index::iterator iter = index_.find(*it);
      if (iter == index_.end()) {
        dbg(DBG_ERROR, "ERROR: ldid %lu appeared in ldo_map set of %lu is not defined yet.\n",
            *it, job->job_id());
      } else {
        iter->second->AddJobEntryWriter(job);
      }
    }
  }

  DetectNewJob(job);

  return true;
}

bool VersionManager::ResolveJobDataVersions(JobEntry *job) {
  assert(job->IsReadyForCompleteVersioning());

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

  DetectVersionedJob(job);

  if (!job->sterile()) {
    boost::unique_lock<boost::recursive_mutex> lock(snap_shot_mutex_);
    if (++non_sterile_counter_ == snap_shot_rate_) {
      GetSnapShot();
      non_sterile_counter_ = 0;
    }
  }

  return true;
}

bool VersionManager::ResolveEntireContextForJob(JobEntry *job) {
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
  return true;
}

bool VersionManager::CreateCheckPoint(JobEntry *job) {
  assert(!job->sterile());

  // Resolve entire ldo_map for job.
  // Note: it only resolves the read versions, and vmap_write may be empty.
  ResolveEntireContextForJob(job);

  // Clear meta before set.
  job->meta_before_set()->Clear();

  // Add the checkpoint to each ldl.
  LdoMap::const_iterator it;
  for (it = ldo_map_p_->begin(); it != ldo_map_p_->end(); ++it) {
    data_version_t version;
    if (job->vmap_read()->query_entry(it->first, &version)) {
      if (job->write_set_p()->contains(it->first)) {
        version = version + 1;
      }
    } else {
      dbg(DBG_ERROR, "ERROR: Version Manager: ldid %lu is not versioned for checkpoint job %lu.\n",
          it->first, job->job_id());
      exit(-1);
    }

    InsertCheckPointLdlEntry(
        it->first, job->job_id(), version, job->job_depth());
  }

  return true;
}




bool VersionManager::DetectNewJob(JobEntry *job) {
  boost::unique_lock<boost::recursive_mutex> lock(snap_shot_mutex_);
  if (job->job_name() != NIMBUS_MAIN_JOB_NAME) {
    ChildCounter::iterator it = child_counter_.find(job->parent_job_id());
    assert(it != child_counter_.end());
    ++it->second;
  }

  if (!job->sterile()) {
    assert(child_counter_.find(job->job_id()) == child_counter_.end());
    child_counter_[job->job_id()] = (counter_t) (1);
    assert(parent_map_.find(job->job_id()) == parent_map_.end());
    parent_map_[job->job_id()] = job;
  }
  return true;
}

bool VersionManager::DetectVersionedJob(JobEntry *job) {
  boost::unique_lock<boost::recursive_mutex> lock(snap_shot_mutex_);
  if (job->job_name() != NIMBUS_MAIN_JOB_NAME) {
    ChildCounter::iterator it = child_counter_.find(job->parent_job_id());
    assert(it != child_counter_.end());
    --it->second;
    if (it->second == 0) {
      child_counter_.erase(it);
    }
  }
  return true;
}

bool VersionManager::DetectDoneJob(JobEntry *job) {
  boost::unique_lock<boost::recursive_mutex> lock(snap_shot_mutex_);
  if (!job->sterile()) {
    ChildCounter::iterator it = child_counter_.find(job->job_id());
    assert(it != child_counter_.end());
    --it->second;
    if (it->second == 0) {
      child_counter_.erase(it);
    }
  }
  return true;
}

bool VersionManager::InsertCheckPointLdlEntry(
    const logical_data_id_t ldid,
    const job_id_t& job_id,
    const data_version_t& version,
    const job_depth_t& job_depth) {
  Index::iterator iter = index_.find(ldid);
  if (iter == index_.end()) {
    return false;
  } else {
    return iter->second->InsertCheckPointLdlEntry(job_id, version, job_depth);
  }
  return true;
}

bool VersionManager::GetSnapShot() {
  boost::unique_lock<boost::recursive_mutex> lock(snap_shot_mutex_);
  if (snap_shot_pending_) {
    CleanUp();
  }

  snap_shot_.clear();
  ChildCounter::iterator iter = child_counter_.begin();
  for (; iter != child_counter_.end(); ++iter) {
    assert(iter->second > 0);
    snap_shot_.insert(iter->first);
  }
  snap_shot_pending_ = true;

  return true;
}

bool VersionManager::LookUpVersion(
    JobEntry *job,
    logical_data_id_t ldid,
    data_version_t *version) {
  Index::iterator iter = index_.find(ldid);
  if (iter == index_.end()) {
    return false;
  } else {
    return iter->second->LookUpVersion(job, version);
  }
}



size_t VersionManager::GetJobsNeedDataVersion(
    JobEntryList* list, VLD vld) {
  Index::iterator iter = index_.find(vld.first);
  if (iter == index_.end()) {
    list->clear();
    return 0;
  } else {
    return iter->second->GetJobsNeedVersion(list, vld.second);
  }
}

bool VersionManager::NotifyJobDone(JobEntry* job) {
  assert(job->versioned());
  DetectDoneJob(job);
  return true;
}

bool VersionManager::RemoveJobEntry(JobEntry* job) {
  assert(job->versioned());

  if (job->sterile()) {
    IDSet<logical_data_id_t>::ConstIter it;
    for (it = job->read_set_p()->begin(); it != job->read_set_p()->end(); ++it) {
      Index::iterator iter = index_.find(*it);
      if (iter != index_.end()) {
        iter->second->RemoveJobEntry(job);
        if (iter->second->is_empty()) {
          delete iter->second;
          index_.erase(iter);
          std::cout << "version entry got empty!!\n";
        }
      }
    }
  } else {
    LdoMap::const_iterator it;
    for (it = ldo_map_p_->begin(); it != ldo_map_p_->end(); ++it) {
      Index::iterator iter = index_.find(it->first);
      if (iter != index_.end()) {
        iter->second->RemoveJobEntry(job);
        if (iter->second->is_empty()) {
          delete iter->second;
          index_.erase(iter);
          std::cout << "version entry got empty!!\n";
        }
      }
    }
  }


  {
    boost::unique_lock<boost::recursive_mutex> lock(snap_shot_mutex_);
    if (snap_shot_.contains(job->job_id())) {
      CleanUp();
    }

    if (!job->sterile()) {
      parent_map_.erase(job->job_id());
    }
  }

  return true;
}

bool VersionManager::DefineData(
    const logical_data_id_t ldid,
    const job_id_t& job_id,
    const job_depth_t& job_depth) {
  Index::iterator iter = index_.find(ldid);
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
  boost::unique_lock<boost::recursive_mutex> lock(snap_shot_mutex_);
  if (snap_shot_pending_) {
    Log log(Log::NO_FILE);
    log.StartTimer();

    IDSet<job_id_t>::IDSetIter iter = snap_shot_.begin();
    for (; iter != snap_shot_.end(); ++iter) {
      ParentMap::iterator it = parent_map_.find(*iter);
      assert(it != parent_map_.end());
      CreateCheckPoint(it->second);
    }

    Index::iterator i = index_.begin();
    for (; i != index_.end(); ++i) {
      i->second->CleanLdl(snap_shot_);
    }

    snap_shot_.clear();
    snap_shot_pending_ = false;

    log.StopTimer();
    std::cout << "Version Manager Clean up: " << log.timer() << std::endl;
  }

  return true;
}


void VersionManager::set_ldo_map_p(const LdoMap* ldo_map_p) {
  ldo_map_p_ = ldo_map_p;
}



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

  batch_size_ = 0;
  batch_state_ = INIT;
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
  FlushPendingBatch();
  {
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

void VersionManager::FlushPendingBatch() {
  boost::unique_lock<boost::mutex> lock(batch_update_mutex_);

  switch (batch_state_) {
    case INIT:
      // std::cout << "LAZY: FLUSH on INIT\n";
      assert(batch_size_ == 0);
      break;
    case DETECTING:
      std::cout << "LAZY: FLUSH on DETECTING\n";
      assert(batch_size_ == 0);
      batch_state_ = INIT;
      break;
    case ACTIVE:
      std::cout << "LAZY: FLUSH on ACTIVE\n";
      assert(batch_size_ > 0);
      InsertComplexJobBatchInLdl();
      batch_size_ = 0;
      batch_state_ = INIT;
      break;
    default:
      assert(false);
      break;
  }
}

bool VersionManager::AddComplexJobEntry(ComplexJobEntry *complex_job) {
  std::string name = complex_job->template_entry()->template_name();
  switch (batch_state_) {
    case INIT:
      std::cout << "LAZY: INIT\n";
      assert(batch_size_ == 0);
      batch_state_ = DETECTING;
      last_complex_job_ = complex_job;
      InsertComplexJobInLdl(complex_job);
      break;
    case DETECTING:
      std::cout << "LAZY: DETECTING\n";
      assert(batch_size_ == 0);
      if (name == last_complex_job_->template_entry()->template_name()) {
        batch_size_ = 1;
        batch_state_ = ACTIVE;
        last_complex_job_ = complex_job;
        base_batch_mbs_ = complex_job->meta_before_set();
        std::cout << "LAZY: TURNED TO ACTIVE\n";
      } else {
        batch_state_ = INIT;
        InsertComplexJobInLdl(complex_job);
      }
      break;
    case ACTIVE:
      std::cout << "LAZY: ACTIVE\n";
      assert(batch_size_ > 0);
      if (name == last_complex_job_->template_entry()->template_name()) {
        ++batch_size_;
        last_complex_job_ = complex_job;
      } else {
        InsertComplexJobBatchInLdl();
        batch_size_ = 0;
        batch_state_ = INIT;
        InsertComplexJobInLdl(complex_job);
      }
      break;
    default:
      assert(false);
      break;
  }

  complex_jobs_[complex_job->job_id()] = complex_job;
  DetectNewComplexJob(complex_job);

//  Log log(Log::NO_FILE);
//  log.StartTimer();
//
//  const ShadowJobEntryMap* jobs =  complex_job->jobs_p();
//
//  ShadowJobEntryMap::const_iterator iter = jobs->begin();
//  for (; iter != jobs->end(); ++iter) {
//    ShadowJobEntry* job = iter->second;
//    {
//      IDSet<logical_data_id_t>::ConstIter it;
//      for (it = job->read_set_p()->begin(); it != job->read_set_p()->end(); ++it) {
//        Index::iterator iter = index_.find(*it);
//        if (iter == index_.end()) {
//          dbg(DBG_ERROR, "ERROR: ldid %lu appeared in read set of %lu is not defined yet.\n",
//              *it, job->job_id());
//        } else {
//          iter->second->AddJobEntryReader(job);
//        }
//      }
//    }
//
//    DetectNewJob(job);
//  }
//
//  log.StopTimer();
//  std::cout << "SPAWN ADD COMPLEX JOB: " << log.timer() << std::endl;

  return true;
}

bool VersionManager::ResolveJobDataVersions(JobEntry *job) {
  FlushPendingBatch();
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
  FlushPendingBatch();
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

bool VersionManager::ResolveJobDataVersionsForPattern(JobEntry *job,
                          const BindingTemplate::PatternList* patterns) {
  FlushPendingBatch();
  BindingTemplate::PatternList::const_iterator it;
  for (it = patterns->begin(); it != patterns->end(); ++it) {
    logical_data_id_t ldid = (*it)->ldid_;
    data_version_t version;
    if (job->vmap_read()->query_entry(ldid, &version)) {
      continue;
    }
    if (LookUpVersion(job, ldid, &version)) {
      job->vmap_read()->set_entry(ldid, version);
    } else {
      return false;
    }
  }
  return true;
}

bool VersionManager::MemoizeVersionsForTemplate(JobEntry *job) {
  FlushPendingBatch();
  TemplateJobEntry* tj = job->template_job();
  assert(job->memoize());
  assert(tj);
  TemplateEntry* te = tj->template_entry();

  assert(job->IsReadyForCompleteVersioning());

  // Set the meta before set for later before set relation lookups.
  tj->set_meta_before_set(job->meta_before_set());
  tj->set_job_depth(job->job_depth());

  if (job->sterile()) {
    IDSet<logical_data_id_t>::ConstIter it;
    for (it = job->union_set_p()->begin(); it != job->union_set_p()->end(); ++it) {
      data_version_t version;
      if (job->vmap_read()->query_entry(*it, &version)) {
      } else if (LookUpVersion(job, *it, &version)) {
        job->vmap_read()->set_entry(*it, version);
      } else {
        return false;
      }

      data_version_t base_version;
      if (!te->vmap_base()->query_entry(*it, &base_version)) {
        return false;
      }

      data_version_t diff_version = version - base_version;
      assert(diff_version >= 0);
      tj->vmap_read_diff()->set_entry(*it, diff_version);

      // Memoize the acceess pattern
      if (tj->read_set_p()->contains(*it)) {
        te->AddToAccessPattern(*it, diff_version, tj->index());
      }

      if (tj->write_set_p()->contains(*it)) {
        diff_version = diff_version + 1;
      }
      if (diff_version > 0) {
        tj->vlist_write_diff()->push_back(std::make_pair(*it, diff_version));
      }
    }
  } else {
    LdoMap::const_iterator it;
    for (it = ldo_map_p_->begin(); it != ldo_map_p_->end(); ++it) {
      data_version_t version;
      if (job->vmap_read()->query_entry(it->first, &version)) {
      } else if (LookUpVersion(job, it->first, &version)) {
        job->vmap_read()->set_entry(it->first, version);
      } else {
        return false;
      }

      data_version_t base_version;
      if (!te->vmap_base()->query_entry(it->first, &base_version)) {
        return false;
      }

      data_version_t diff_version = version - base_version;
      assert(diff_version >= 0);
      tj->vmap_read_diff()->set_entry(it->first, diff_version);

      // Memoize the acceess pattern
      if (tj->read_set_p()->contains(it->first)) {
        te->AddToAccessPattern(it->first, diff_version, tj->index());
      }

      if (tj->write_set_p()->contains(it->first)) {
        diff_version = diff_version + 1;
      }
      if (diff_version > 0) {
        tj->vlist_write_diff()->push_back(std::make_pair(it->first, diff_version));
      }
    }
  }

  return true;
}

bool VersionManager::InsertComplexJobBatchInLdl() {
  assert(batch_size_ > 0);
  std::cout << "LAZY: Flush a Batch of size: " << batch_size_ << std::endl;
  assert(batch_state_ == ACTIVE);
  ShadowJobEntryList list;
  last_complex_job_->OMIDGetParentShadowJobs(&list);
  // For now only one parent job per complex job is allowd!
  assert(list.size() == 1);
  ShadowJobEntry* sj = *(list.begin());
  ComplexJobEntry* xj = sj->complex_job();
  assert(xj == last_complex_job_);

  VersionList::const_iterator iter = sj->vlist_write_diff()->begin();
  for (; iter != sj->vlist_write_diff()->end(); ++iter) {
    data_version_t diff_version = iter->second;
    data_version_t base_version;
    if (xj->vmap_read()->query_entry(iter->first, &base_version)) {
    } else if (LookUpVersionByMetaBeforeSet(base_batch_mbs_, iter->first, &base_version)) {
      // inceremnt the base based on the batch size
      xj->vmap_read()->set_entry(iter->first, base_version + (batch_size_ - 1) * diff_version);
    } else {
      dbg(DBG_ERROR, "ERROR: complex job %lu as batch could not be versioned for ldid %lu.", xj->job_id(), iter->first); //NOLINT
      assert(false);
      return false;
    }

    // inceremnt the base based on the batch size
    InsertCheckPointLdlEntry(
        iter->first, xj->job_id(), base_version + batch_size_ * diff_version, xj->job_depth());
  }

  return true;
}

bool VersionManager::InsertComplexJobInLdl(ComplexJobEntry *job) {
  ShadowJobEntryList list;
  job->OMIDGetParentShadowJobs(&list);
  // For now only one parent job per complex job is allowd!
  assert(list.size() == 1);
  ShadowJobEntry* sj = *(list.begin());
  ComplexJobEntry* xj = sj->complex_job();
  assert(xj == job);

  VersionList::const_iterator iter = sj->vlist_write_diff()->begin();
  for (; iter != sj->vlist_write_diff()->end(); ++iter) {
    data_version_t base_version;
    if (xj->vmap_read()->query_entry(iter->first, &base_version)) {
    } else if (LookUpVersion(xj, iter->first, &base_version)) {
      xj->vmap_read()->set_entry(iter->first, base_version);
    } else {
      dbg(DBG_ERROR, "ERROR: complex job %lu could not be versioned for ldid %lu.", xj->job_id(), iter->first); //NOLINT
      assert(false);
      return false;
    }

    data_version_t diff_version = iter->second;

    InsertCheckPointLdlEntry(
        iter->first, xj->job_id(), base_version + diff_version, xj->job_depth());
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
      assert(false);
    }

    InsertCheckPointLdlEntry(
        it->first, job->job_id(), version, job->job_depth());
  }

  return true;
}


bool VersionManager::DetectNewComplexJob(ComplexJobEntry *xj) {
  boost::unique_lock<boost::recursive_mutex> lock(snap_shot_mutex_);
  {
    ChildCounter::iterator it = child_counter_.find(xj->parent_job_id());
    assert(it != child_counter_.end());
    it->second += xj->inner_job_ids_p()->size();
  }

  ShadowJobEntryList list;
  xj->OMIDGetParentShadowJobs(&list);
  // For now complex job can have only one parent job - omidm
  assert(list.size() == 1);
  JobEntry *j = *(list.begin());
  assert(!j->sterile());
  {
    assert(child_counter_.find(j->job_id()) == child_counter_.end());
    child_counter_[j->job_id()] = (counter_t) (1);
    assert(parent_map_.find(j->job_id()) == parent_map_.end());
    parent_map_[j->job_id()] = j;
  }

  return true;
}

bool VersionManager::DetectNewJob(JobEntry *job) {
  boost::unique_lock<boost::recursive_mutex> lock(snap_shot_mutex_);
  // if (job->job_name() != NIMBUS_MAIN_JOB_NAME) {
  if (job->parent_job_id() != NIMBUS_KERNEL_JOB_ID) {
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
  // if (job->job_name() != NIMBUS_MAIN_JOB_NAME) {
  if (job->parent_job_id() != NIMBUS_KERNEL_JOB_ID) {
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
  // TODO(omidm): fix snap shot logic when there are complex jobs!
  return false;

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

bool VersionManager::LookUpVersionByMetaBeforeSet(
    boost::shared_ptr<MetaBeforeSet> mbs,
    logical_data_id_t ldid,
    data_version_t *version) {
  Index::iterator iter = index_.find(ldid);
  if (iter == index_.end()) {
    return false;
  } else {
    return iter->second->LookUpVersionByMetaBeforeSet(mbs, version);
  }
}

size_t VersionManager::GetJobsNeedDataVersion(
    JobEntryList* list, VLD vld) {
  FlushPendingBatch();
  size_t count = 0;
  list->clear();

  // TODO(omidm): This is not necessarily true, if the non-sterile job is not
  // bottleneck and runs and spawns the new jobs before it siblings finish.
  assert(complex_jobs_.size() <= 1);

  {
    ComplexJobMap::iterator iter = complex_jobs_.begin();
    for (; iter != complex_jobs_.end(); ++iter) {
      ComplexJobEntry *xj = iter->second;
      assert(!xj->done());

      TemplateEntry *te = xj->template_entry();

      data_version_t base_version;
      if (xj->vmap_read()->query_entry(vld.first, &base_version)) {
      } else if (LookUpVersion(xj, vld.first, &base_version)) {
        xj->vmap_read()->set_entry(vld.first, base_version);
      } else {
        dbg(DBG_ERROR, "ERROR: could not version the base complex job %lu.\n", xj->job_id());
        assert(false);
      }

      if (vld.second < base_version) {
        continue;
      }

      data_version_t diff_version = vld.second - base_version;

      std::list<size_t> indices;
      te->QueryAccessPattern(vld.first, diff_version, &indices);
      std::list<size_t>::iterator it = indices.begin();
      for (; it != indices.end(); ++it) {
        ShadowJobEntry *sj;
        if (xj->OMIDGetShadowJobEntryByIndex(*it, sj)) {
          bool done = xj->ShadowJobDone(sj->job_id());
          if ((!sj->assigned()) ||
              ((!sj->sterile()) && (!done))) {
            list->push_back(sj);
            ++count;
          }
        }
      }
    }
  }

  {
    Index::iterator iter = index_.find(vld.first);
    if (iter == index_.end()) {
      return count;
    } else {
      {
        boost::unique_lock<boost::recursive_mutex> lock(snap_shot_mutex_);
        ParentMap::iterator it = parent_map_.begin();
        for (; it != parent_map_.end(); ++it) {
          if (!(it->second->done())) {
            iter->second->AddJobEntryReader(it->second);
          }
        }
      }

      count += iter->second->GetJobsNeedVersion(list, vld.second, true);
    }
  }

  return count;
}

bool VersionManager::NotifyJobDone(JobEntry* job) {
  if (job->job_type() == JOB_CMPX) {
    // TODO(omidm): does child_counter_ gets cleaned for the parent shadow job
    // with in the complex job (DetectDoneJob ...).
    complex_jobs_.erase(job->job_id());
    return true;
  }

  assert(job->versioned());
  DetectDoneJob(job);
  return true;
}

bool VersionManager::RemoveJobEntry(JobEntry* job) {
  // TODO(omidm): may need better clean up for complex and """*SHADOW*""" jobs.
  if (job->job_type() == JOB_CMPX) {
    assert(complex_jobs_.find(job->job_id()) == complex_jobs_.end());

    ComplexJobEntry *xj = reinterpret_cast<ComplexJobEntry*>(job);
    ShadowJobEntryList list;
    xj->OMIDGetParentShadowJobs(&list);
    // For now complex job can have only one parent job - omidm
    assert(list.size() == 1);
    JobEntry *j = *(list.begin());
    assert(!j->sterile());

    {
      boost::unique_lock<boost::recursive_mutex> lock(snap_shot_mutex_);
      if (snap_shot_.contains(j->job_id())) {
        CleanUp();
      }

      assert(parent_map_.find(j->job_id()) != parent_map_.end());
      parent_map_.erase(j->job_id());
    }

    LdoMap::const_iterator it;
    for (it = ldo_map_p_->begin(); it != ldo_map_p_->end(); ++it) {
      Index::iterator iter = index_.find(it->first);
      if (iter != index_.end()) {
        // Set versioned to comply with version entry assumption.
        // Note that with active binding template job is never versioned -omidm
        j->set_versioned(true);
        iter->second->RemoveJobEntry(j);
        // Even empty you cannot remove, otherwise the ldl would show as undefined!
        // if (iter->second->is_empty()) {
        //   delete iter->second;
        //   index_.erase(iter);
        //   std::cout << "version entry got empty!!\n";
        // }
      }
    }

    return true;
  }



  assert(job->versioned());

  {
    boost::unique_lock<boost::recursive_mutex> lock(snap_shot_mutex_);
    if (snap_shot_.contains(job->job_id())) {
      CleanUp();
    }

    if (!job->sterile()) {
      parent_map_.erase(job->job_id());
    }
  }

  if (job->sterile()) {
    IDSet<logical_data_id_t>::ConstIter it;
    for (it = job->read_set_p()->begin(); it != job->read_set_p()->end(); ++it) {
      Index::iterator iter = index_.find(*it);
      if (iter != index_.end()) {
        iter->second->RemoveJobEntry(job);
        // Even empty you cannot remove, otherwise the ldl would show as undefined!
        // if (iter->second->is_empty()) {
        //   delete iter->second;
        //   index_.erase(iter);
        //   std::cout << "version entry got empty!!\n";
        // }
      }
    }
  } else {
    LdoMap::const_iterator it;
    for (it = ldo_map_p_->begin(); it != ldo_map_p_->end(); ++it) {
      Index::iterator iter = index_.find(it->first);
      if (iter != index_.end()) {
        iter->second->RemoveJobEntry(job);
        // Even empty you cannot remove, otherwise the ldl would show as undefined!
        // if (iter->second->is_empty()) {
        //   delete iter->second;
        //   index_.erase(iter);
        //   std::cout << "version entry got empty!!\n";
        // }
      }
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
    assert(false);
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

void VersionManager::Reinitialize(const JobEntryList *list) {
  FlushPendingBatch();
  snap_shot_pending_ = false;
  non_sterile_counter_ = 0;
  snap_shot_.clear();
  parent_map_.clear();
  child_counter_.clear();
  complex_jobs_.clear();

  {
    Index::iterator iter = index_.begin();
    for (; iter != index_.end(); ++iter) {
      iter->second->Reinitialize();
      assert(iter->second->is_empty());
    }
  }

  JobEntryList::const_iterator iter = list->begin();
  for (; iter != list->end(); ++iter) {
    // Add the checkpoint to each ldl.
    LdoMap::const_iterator it;
    for (it = ldo_map_p_->begin(); it != ldo_map_p_->end(); ++it) {
      data_version_t version;
      if ((*iter)->vmap_read()->query_entry(it->first, &version)) {
        if ((*iter)->write_set_p()->contains(it->first)) {
          version = version + 1;
        }
      } else {
        dbg(DBG_ERROR, "ERROR: Version Manager: ldid %lu is not versioned for checkpoint job %lu.\n", // NOLINT
            it->first, (*iter)->job_id());
        assert(false);
      }

      InsertCheckPointLdlEntry(
          it->first, (*iter)->job_id(), version, (*iter)->job_depth());
    }
  }
}

void VersionManager::set_ldo_map_p(const LdoMap* ldo_map_p) {
  ldo_map_p_ = ldo_map_p;
}



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
}

VersionManager::~VersionManager() {
  VersionIndex::iterator iter = version_index_.begin();
  for (; iter != version_index_.end(); ++iter) {
    VersionEntryList* list = (*iter).second;
    delete list;
  }
}

bool VersionManager::AddVersionEntry(
    logical_data_id_t logical_id, data_version_t version,
    JobEntry* job_entry, VersionEntry::Relation relation) {
  VersionEntryList* list;
  if (version_index_.find(logical_id) == version_index_.end()) {
    list = new VersionEntryList();
    version_index_[logical_id] = list;
  } else {
    list = version_index_[logical_id];
  }
  VersionEntry ve(logical_id, version, job_entry, relation);
  list->push_back(ve);
  return true;
}

bool VersionManager::AddVersionEntry(const VersionEntry& ve) {
  logical_data_id_t logical_id = ve.logical_id();
  VersionEntryList* list;
  if (version_index_.find(logical_id) == version_index_.end()) {
    list = new VersionEntryList();
    version_index_[logical_id] = list;
  } else {
    list = version_index_[logical_id];
  }
  list->push_back(ve);
  return true;
}

bool VersionManager::AddJobVersionTables(JobEntry* job_entry) {
  if (!job_entry->versioned()) {
    return false;
  }

  JobEntry::VersionTable::iterator iter;

  JobEntry::VersionTable  vt_in = job_entry->version_table_in();
  for (iter = vt_in.begin(); iter != vt_in.end(); ++iter) {
    AddVersionEntry(iter->first, iter->second, job_entry, VersionEntry::IN);
  }

  JobEntry::VersionTable  vt_out = job_entry->version_table_out();
  for (iter = vt_out.begin(); iter != vt_out.end(); ++iter) {
    AddVersionEntry(iter->first, iter->second, job_entry, VersionEntry::OUT);
  }
  return true;
}

bool VersionManager::RemoveVersionEntry(
    logical_data_id_t logical_id, data_version_t version,
    JobEntry* job_entry, VersionEntry::Relation relation) {
  if (version_index_.find(logical_id) == version_index_.end()) {
    dbg(DBG_ERROR,"ERROR: the logical id %lu is not in the version index.\n", logical_id); // NOLINT 
    return false;
  } else {
    VersionEntryList* list = version_index_[logical_id];
    VersionEntryList::iterator iter = list->begin();
    for (; iter != list->end(); ++iter) {
      if ((iter->version() == version) &&
          (iter->job_entry() == job_entry) &&
          (iter->relation() == relation)) {
        list->erase(iter);
        return true;
      }
    }
    dbg(DBG_ERROR,"ERROR: could not match any entry in version index.\n"); // NOLINT
    return false;
  }
}

bool VersionManager::RemoveVersionEntry(const VersionEntry& ve) {
  logical_data_id_t logical_id = ve.logical_id();
  if (version_index_.find(logical_id) == version_index_.end()) {
    dbg(DBG_ERROR,"ERROR: the logical id %lu is not in the version index.\n", logical_id); // NOLINT 
    return false;
  } else {
    VersionEntryList* list = version_index_[logical_id];
    VersionEntryList::iterator iter = list->begin();
    for (; iter != list->end(); ++iter) {
      if ((iter->version() == ve.version()) &&
          (iter->job_entry() == ve.job_entry()) &&
          (iter->relation() == ve.relation())) {
        list->erase(iter);
        return true;
      }
    }
    dbg(DBG_ERROR,"ERROR: could not match any entry in version index.\n"); // NOLINT
    return false;
  }
}

bool VersionManager::RemoveJobVersionTables(JobEntry* job_entry) {
  if (!job_entry->versioned()) {
    return false;
  }

  JobEntry::VersionTable::iterator iter;

  JobEntry::VersionTable  vt_in = job_entry->version_table_in();
  for (iter = vt_in.begin(); iter != vt_in.end(); ++iter) {
    RemoveVersionEntry(iter->first, iter->second, job_entry, VersionEntry::IN);
  }

  JobEntry::VersionTable  vt_out = job_entry->version_table_out();
  for (iter = vt_out.begin(); iter != vt_out.end(); ++iter) {
    RemoveVersionEntry(iter->first, iter->second, job_entry, VersionEntry::OUT);
  }
  return true;
}

size_t VersionManager::GetJobsNeedDataVersion(
    JobEntryList* result, VersionedLogicalData vld) {
  size_t num = 0;
  result->clear();
  if (version_index_.find(vld.first) == version_index_.end()) {
    return num;
  } else {
    VersionEntryList* list = version_index_[vld.first];
    VersionEntryList::iterator iter = list->begin();
    for (; iter != list->end(); ++iter) {
      JobEntry* j = iter->job_entry();
      if ((iter->version() == vld.second) &&
          (iter->relation() == VersionEntry::IN) &&
          !(j->assigned()) &&
          ((j->read_set_p()->contains(vld.first)) || !(j->sterile()))) {
        result->push_back(iter->job_entry());
        ++num;
      }
    }
  }

  return num;
}

size_t VersionManager::GetJobsOutputDataVersion(
    JobEntryList* result, VersionedLogicalData vld) {
  size_t num = 0;
  result->clear();
  if (version_index_.find(vld.first) == version_index_.end()) {
    return num;
  } else {
    VersionEntryList* list = version_index_[vld.first];
    VersionEntryList::iterator iter = list->begin();
    for (; iter != list->end(); ++iter) {
      if ((iter->version() == vld.second) &&
          (iter->relation() == VersionEntry::OUT)) {
        result->push_back(iter->job_entry());
        ++num;
      }
    }
  }

  return num;
}

size_t VersionManager::RemoveObsoleteVersionEntriesOfLdo(
    logical_data_id_t logical_id) {
  size_t num = 0;
  if (version_index_.find(logical_id) == version_index_.end()) {
    return num;
  } else {
    VersionEntryList* list = version_index_[logical_id];
    VersionEntryList::iterator iter = list->begin();
    for (; iter != list->end();) {
      // TODO(omidm): figure out when it is obsolete.
      if (false) {
        list->erase(iter++);
        ++num;
      } else {
        ++iter;
      }
    }
  }
  return num;
}

size_t VersionManager::RemoveAllObsoleteVersionEntries() {
  VersionIndex::iterator it = version_index_.begin();
  size_t num = 0;
  for (; it != version_index_.end(); ++it) {
    VersionEntryList* list = (*it).second;
    VersionEntryList::iterator iter = list->begin();
    for (; iter != list->end();) {
      // TODO(omidm): figure out when it is obsolete.
      if (false) {
        list->erase(iter++);
        ++num;
      } else {
        ++iter;
      }
    }
  }
  return num;
}




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
  * Scheduler Job Manager object. This module serves the scheduler by providing
  * facilities about jobs ready to be maped, and their dependencies.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "scheduler/job_manager.h"

using namespace nimbus; // NOLINT

JobManager::JobManager() {
  // Add the SCHED job which is the parent of main, create and copy jobs that
  // are spawned by the scheduler.
  IDSet<job_id_t> job_id_set;
  IDSet<logical_data_id_t> logical_data_id_set;
  Parameter params;
  JobEntry* job = new JobEntry(JOB_SCHED, "kernel", (job_id_t)(0), (job_id_t)(0));
  if (!job_graph_.AddJobEntry(job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add scheduler kernel job in job manager constructor.\n");
  } else {
    job->set_versioned(true);
    job->set_assigned(true);
  }
}

JobManager::~JobManager() {
  JobEntry* job;
  if (JobManager::GetJobEntry((job_id_t)(0), job)) {
    delete job;
  }
}

bool JobManager::AddJobEntry(const JobType& job_type,
    const std::string& job_name,
    const job_id_t& job_id,
    const IDSet<logical_data_id_t>& read_set,
    const IDSet<logical_data_id_t>& write_set,
    const IDSet<job_id_t>& before_set,
    const IDSet<job_id_t>& after_set,
    const job_id_t& parent_job_id,
    const Parameter& params) {
  JobEntry* job = new JobEntry(job_type, job_name, job_id, read_set, write_set,
      before_set, after_set, parent_job_id, params);
  if (job_graph_.AddJobEntry(job)) {
    return true;
  } else {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    return false;
  }
}

bool JobManager::AddJobEntry(const JobType& job_type,
    const std::string& job_name,
    const job_id_t& job_id,
    const job_id_t& parent_job_id,
    const bool& versioned,
    const bool& assigned) {
  JobEntry* job = new JobEntry(job_type, job_name, job_id, parent_job_id);
  if (job_graph_.AddJobEntry(job)) {
    job->set_versioned(versioned);
    job->set_assigned(assigned);
    return true;
  } else {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    return false;
  }
}

bool JobManager::GetJobEntry(job_id_t job_id, JobEntry*& job) {
  return job_graph_.GetJobEntry(job_id, job);
}

bool JobManager::RemoveJobEntry(JobEntry* job) {
  if (job_graph_.RemoveJobEntry(job)) {
    delete job;
    return true;
  } else {
    return false;
  }
}

bool JobManager::RemoveJobEntry(job_id_t job_id) {
  JobEntry* job;
  if (GetJobEntry(job_id, job)) {
    job_graph_.RemoveJobEntry(job);
    delete job;
    return true;
  } else {
    return false;
  }
}

size_t JobManager::GetJobsReadyToAssign(JobEntryList* list, size_t max_num) {
  while (ResolveVersions() > 0) {
    continue;
  }

  size_t num = 0;
  list->clear();
  JobGraph::Iter iter = job_graph_.Begin();
  for (; (iter != job_graph_.End()) && (num < max_num); ++iter) {
    JobEntry* job = iter->second;
    if (job->versioned() && !job->assigned()) {
      bool before_set_done = true;
      IDSet<job_id_t>::IDSetIter it;
      IDSet<job_id_t> before_set = job->before_set();
      for (it = before_set.begin(); it != before_set.end(); ++it) {
        JobEntry* j;
        job_id_t id = *it;
        if (GetJobEntry(id, j)) {
          if (!(j->done())) {
            dbg(DBG_SCHED, "Job in befor set (id: %lu) is not done yet.\n", id);
            before_set_done = false;
            break;
          }
        } else {
          dbg(DBG_ERROR, "ERROR: Job in befor set (id: %lu) is not in the graph.\n", id);
          before_set_done = false;
          break;
        }
      }
      if (before_set_done) {
        // job->set_assigned(true); No, we are not sure yet thet it will be assignd!
        list->push_back(job);
        ++num;
      }
    }
  }
  return num;
}

size_t JobManager::RemoveObsoleteJobEntries() {
  return 0;
}

void JobManager::JobDone(job_id_t job_id) {
  JobEntry* job;
  if (GetJobEntry(job_id, job)) {
    job->set_done(true);
  } else {
    dbg(DBG_WARN, "WARNING: done job with id %lu is not in the graph.\n", job_id);
  }
}

void JobManager::DefineData(job_id_t job_id, logical_data_id_t ldid) {
  JobEntry* job;
  if (GetJobEntry(job_id, job)) {
    JobEntry::VersionTable vt;
    vt = job->version_table();
    JobEntry::VTIter it = vt.begin();
    bool new_logical_data = true;
    for (; it != vt.end(); ++it) {
      if (it->first == ldid) {
        new_logical_data = false;
        dbg(DBG_ERROR, "ERROR: defining logical data id %lu, which already exist.\n", ldid);
        break;
      }
    }
    if (new_logical_data) {
      vt[ldid] = (data_version_t)(0);
      job->set_version_table(vt);
    }
  } else {
    dbg(DBG_WARN, "WARNING: parent of define data with job id %lu is not in the graph.\n", job_id);
  }
}

bool JobManager::ResolveJobDataVersions(JobEntry* job) {
  if (job->versioned())
    return true;

  JobEntry* j;
  JobEntry::VersionTable version_table, vt;

  job_id_t parent_id = job->parent_job_id();
  if (GetJobEntry(parent_id, j)) {
    if (j->versioned()) {
      vt = j->version_table();
      JobEntry::VTIter it = vt.begin();
      for (; it != vt.end(); ++it) {
        if (j->write_set().contains(it->first)) {
          version_table[it->first] =  ++(it->second);
        } else {
          version_table[it->first] =  (it->second);
        }
      }
    } else {
      dbg(DBG_ERROR, "ERROR: parent job (id: %lu) is not versioned yet.\n", parent_id);
      return false;
    }
  } else {
    dbg(DBG_ERROR, "ERROR: parent job (id: %lu) is not in job graph.\n", parent_id);
    return false;
  }

  IDSet<job_id_t>::IDSetIter iter;

  IDSet<job_id_t> before_set = job->before_set();
  for (iter = before_set.begin(); iter != before_set.end(); ++iter) {
    job_id_t id = (*iter);
    if (GetJobEntry(id, j)) {
      if (j->versioned()) {
        vt = j->version_table();
        JobEntry::VTIter it = vt.begin();
        for (; it != vt.end(); ++it) {
          if (version_table.count(it->first) == 0) {
            if (j->write_set().contains(it->first)) {
              version_table[it->first] =  ++(it->second);
            } else {
              version_table[it->first] =  (it->second);
            }
          } else {
            if (j->write_set().contains(it->first)) {
              version_table[it->first] =
                std::max(++(it->second), version_table[it->first]);
            } else {
              version_table[it->first] =
                std::max((it->second), version_table[it->first]);
            }
          }
        }
      } else {
        dbg(DBG_SCHED, "Job in befor set (id: %lu) is not versioned yet.\n", id);
        return false;
      }
    } else {
      dbg(DBG_SCHED, "Job in befor set (id: %lu) is not in the graph.\n", id);
      return false;
    }
  }

  IDSet<logical_data_id_t> read_set = job->read_set();
  for (iter = read_set.begin(); iter != read_set.end(); ++iter) {
    if (version_table.count(*iter) == 0) {
      dbg(DBG_ERROR, "ERROR: parent and before set could not resolve read id %lu.\n", *iter);
      return false;
    }
  }
  IDSet<logical_data_id_t> write_set = job->write_set();
  for (iter = write_set.begin(); iter != write_set.end(); ++iter) {
    if (version_table.count(*iter) == 0) {
      dbg(DBG_ERROR, "ERROR: parent and before set could not resolve write id %lu.\n", *iter);
      return false;
    }
  }

  job->set_versioned(true);
  job->set_version_table(version_table);
  return true;
}

size_t JobManager::ResolveVersions() {
  size_t num_new_versioned = 0;
  JobGraph::Iter iter = job_graph_.Begin();
  for (; iter != job_graph_.End(); ++iter) {
    if (!iter->second->versioned()) {
      if (ResolveJobDataVersions(iter->second))
        ++num_new_versioned;
    }
  }
  return num_new_versioned;
}


size_t JobManager::GetJobsNeedDataVersion(JobEntryList* list,
    JobEntry::VersionedLogicalData vld) {
  size_t num = 0;
  list->clear();
  JobGraph::Iter iter = job_graph_.Begin();
  for (; iter != job_graph_.End(); ++iter) {
    JobEntry* job = iter->second;
    if (job->versioned() && !job->assigned()) {
      JobEntry::VersionTable version_table = job->version_table();
      if (version_table.count(vld.first) != 0) {
        if ((version_table[vld.first] == vld.second) &&
            (job->union_set().contains(vld.first))) {
          list->push_back(job);
          ++num;
        }
      }
    }
  }
  return num;
}

void JobManager::UpdateJobBeforeSet(JobEntry* job) {
  IDSet<job_id_t> before_set = job->before_set();
  UpdateBeforeSet(&before_set);
  job->set_before_set(before_set);
}

void JobManager::UpdateBeforeSet(IDSet<job_id_t>* before_set) {
  IDSet<job_id_t>::IDSetIter it;
  for (it = before_set->begin(); it != before_set->end();) {
    JobEntry* j;
    job_id_t id = *it;
    if (GetJobEntry(id, j)) {
      if (j->done()) {
        before_set->remove(it++);
      } else {
        ++it;
      }
    } else {
      ++it;
    }
  }
}




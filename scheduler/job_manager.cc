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
  // Add the KERNEL job.
  if (AddKernelJobEntry() == NULL) {
    dbg(DBG_ERROR, "ERROR: could not add scheduler kernel job in job manager constructor.\n");
    exit(-1);
  }

  ldo_map_p_ = NULL;
  after_map_ = NULL;
  log_.set_file_name("log_job_manager");
}

JobManager::~JobManager() {
  // TODO(omidm): do you need to call remove obsolete?
  JobEntry* job;
  if (JobManager::GetJobEntry(NIMBUS_KERNEL_JOB_ID, job)) {
    delete job;
  }
}

void JobManager::set_after_map(AfterMap* after_map) {
  after_map_ = after_map;
}

void JobManager::set_ldo_map_p(const LdoMap* ldo_map_p) {
  ldo_map_p_ = ldo_map_p;
  version_manager_.set_ldo_map_p(ldo_map_p);
}

Graph<JobEntry, job_id_t>* JobManager::job_graph_p() {
  return &job_graph_;
}

JobEntry* JobManager::AddComputeJobEntry(const std::string& job_name,
                                         const job_id_t& job_id,
                                         const IDSet<logical_data_id_t>& read_set,
                                         const IDSet<logical_data_id_t>& write_set,
                                         const IDSet<job_id_t>& before_set,
                                         const IDSet<job_id_t>& after_set,
                                         const job_id_t& parent_job_id,
                                         const job_id_t& future_job_id,
                                         const bool& sterile,
                                         const Parameter& params) {
  JobEntry* job =
    new ComputeJobEntry(job_name,
                        job_id,
                        read_set,
                        write_set,
                        before_set,
                        after_set,
                        parent_job_id,
                        future_job_id,
                        sterile,
                        params);

  if (!job_graph_.AddVertex(job_id, job)) {
    dbg(DBG_SCHED, "Filling possible future job (id: %lu) in job manager.\n", job_id);
    delete job;
    GetJobEntry(job_id, job);
    if (job->future()) {
      job->set_job_type(JOB_COMP);
      job->set_job_name(job_name);
      job->set_read_set(read_set);
      job->set_write_set(write_set);
      job->set_before_set(before_set, true);
      job->set_after_set(after_set);
      job->set_parent_job_id(parent_job_id, true);
      job->set_params(params);
      job->set_sterile(sterile);
      job->set_future(false);
      dbg(DBG_SCHED, "Filled the information for future job (id: %lu).\n", job_id);
    } else {
      dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
      exit(-1);
      return NULL;
    }
  }

  if (!AddJobEntryIncomingEdges(job)) {
    job_graph_.RemoveVertex(job_id);
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    exit(-1);
    return NULL;
  }
  ReceiveMetaBeforeSetDepthVersioningDependency(job);
  PassMetaBeforeSetDepthVersioningDependency(job);

  version_manager_.AddJobEntry(job);

  return job;
}

JobEntry* JobManager::AddExplicitCopyJobEntry() {
  dbg(DBG_ERROR, "ERROR: explicit copy jobs from application are not supported yet!.\n");
  exit(-1);
  return NULL;
}

JobEntry* JobManager::AddKernelJobEntry() {
  JobEntry* job = new KernelJobEntry();

  if (!job_graph_.AddVertex(NIMBUS_KERNEL_JOB_ID, job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add kernel job in job manager.\n");
    exit(-1);
    return NULL;
  }
  job->set_done(true);

  return job;
}

JobEntry* JobManager::AddMainJobEntry(const job_id_t& job_id) {
  JobEntry* job = new MainJobEntry(job_id);

  if (!job_graph_.AddVertex(job_id, job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    exit(-1);
    return NULL;
  }

  if (!AddJobEntryIncomingEdges(job)) {
    job_graph_.RemoveVertex(job_id);
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    exit(-1);
    return NULL;
  }
  ReceiveMetaBeforeSetDepthVersioningDependency(job);
  PassMetaBeforeSetDepthVersioningDependency(job);

  version_manager_.AddJobEntry(job);

  jobs_ready_to_assign_[job_id] = job;

  return job;
}

JobEntry* JobManager::AddCreateDataJobEntry(const job_id_t& job_id) {
  JobEntry* job = new CreateDataJobEntry(job_id);

  boost::unique_lock<boost::mutex> lock(job_graph_mutex_);
  if (!job_graph_.AddVertex(job_id, job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    exit(-1);
    return NULL;
  }

  return job;
}

JobEntry* JobManager::AddLocalCopyJobEntry(const job_id_t& job_id) {
  JobEntry* job = new LocalCopyJobEntry(job_id);

  boost::unique_lock<boost::mutex> lock(job_graph_mutex_);
  if (!job_graph_.AddVertex(job_id, job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    exit(-1);
    return NULL;
  }

  return job;
}

JobEntry* JobManager::AddRemoteCopySendJobEntry(const job_id_t& job_id) {
  JobEntry* job = new RemoteCopySendJobEntry(job_id);

  boost::unique_lock<boost::mutex> lock(job_graph_mutex_);
  if (!job_graph_.AddVertex(job_id, job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    exit(-1);
    return NULL;
  }

  return job;
}

JobEntry* JobManager::AddRemoteCopyReceiveJobEntry(const job_id_t& job_id) {
  JobEntry* job = new RemoteCopyReceiveJobEntry(job_id);

  boost::unique_lock<boost::mutex> lock(job_graph_mutex_);
  if (!job_graph_.AddVertex(job_id, job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    exit(-1);
    return NULL;
  }

  return job;
}

JobEntry* JobManager::AddFutureJobEntry(const job_id_t& job_id) {
  JobEntry* job = new FutureJobEntry(job_id);

  boost::unique_lock<boost::mutex> lock(job_graph_mutex_);
  if (!job_graph_.AddVertex(job_id, job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager as future job.\n", job_id);
    exit(-1);
    return NULL;
  }

  return job;
}

bool JobManager::AddJobEntryIncomingEdges(JobEntry *job) {
  Edge<JobEntry, job_id_t> *edge;
  JobEntry *j;

  if (job_graph_.AddEdge(job->parent_job_id(), job->job_id(), &edge)) {
    j = edge->start_vertex()->entry();
    assert(j->versioned());
  } else {
    dbg(DBG_ERROR, "ERROR: could not add edge from parent (id: %lu) for job (id: %lu) in job manager.\n", // NOLINT
        job->parent_job_id(), job->job_id());
    return false;
  }

  IDSet<job_id_t>::ConstIter it;
  for (it = job->before_set_p()->begin(); it != job->before_set_p()->end(); ++it) {
    if (*it == job->parent_job_id()) {
      continue;
    }

    if (job_graph_.AddEdge(*it, job->job_id(), &edge)) {
      j = edge->start_vertex()->entry();
    } else {
      dbg(DBG_SCHED, "Adding possible future job (id: %lu) in job manager.\n", *it);
      AddFutureJobEntry(*it);
      if (!job_graph_.AddEdge(*it, job->job_id())) {
        return false;
      }
    }
  }

  return true;
}


void JobManager::ReceiveMetaBeforeSetDepthVersioningDependency(JobEntry* job) {
  job_depth_t depth = NIMBUS_INIT_JOB_DEPTH;
  job_id_t job_id = job->job_id();
  Vertex<JobEntry, job_id_t>* vertex;
  job_graph_.GetVertex(job_id, &vertex);
  typename Edge<JobEntry, job_id_t>::Iter it;
  for (it = vertex->incoming_edges()->begin(); it != vertex->incoming_edges()->end(); ++it) {
    JobEntry *j = it->second->start_vertex()->entry();
    if (j->IsReadyForCompleteVersioning()) {
      job->meta_before_set()->table_p()->insert(
          std::pair<job_id_t, boost::shared_ptr<MetaBeforeSet> > (j->job_id(), j->meta_before_set())); // NOLINT

      if (depth < j->job_depth()) {
        depth = j->job_depth();
      }

      job->remove_versioning_dependency(j->job_id());
    }
  }
  job->set_job_depth(depth + 1);
}


void JobManager::PassMetaBeforeSetDepthVersioningDependency(JobEntry* job) {
  if (job->IsReadyForCompleteVersioning()) {
    job_id_t job_id = job->job_id();
    Vertex<JobEntry, job_id_t>* vertex;
    job_graph_.GetVertex(job_id, &vertex);
    typename Edge<JobEntry, job_id_t>::Iter it;
    for (it = vertex->outgoing_edges()->begin(); it != vertex->outgoing_edges()->end(); ++it) {
      JobEntry *j = it->second->end_vertex()->entry();

      j->meta_before_set()->table_p()->insert(
          std::pair<job_id_t, boost::shared_ptr<MetaBeforeSet> > (job->job_id(), job->meta_before_set())); // NOLINT

      if (j->job_depth() <= job->job_depth()) {
        j->set_job_depth(job->job_depth() + 1);
      }

      j->remove_versioning_dependency(job_id);
      if (j->IsReadyForCompleteVersioning()) {
        PassMetaBeforeSetDepthVersioningDependency(j);
      }
    }
  }
}


bool JobManager::GetJobEntry(job_id_t job_id, JobEntry*& job) {
  Vertex<JobEntry, job_id_t>* vertex;
  boost::unique_lock<boost::mutex> lock(job_graph_mutex_);
  if (job_graph_.GetVertex(job_id, &vertex)) {
    job = vertex->entry();
    return true;
  } else {
    job = NULL;
    return false;
  }
}

bool JobManager::RemoveJobEntry(JobEntry* job) {
  if (job_graph_.RemoveVertex(job->job_id())) {
    if (job->parent_job_id() != NIMBUS_KERNEL_JOB_ID ||
        job->job_name() == NIMBUS_MAIN_JOB_NAME) {
      version_manager_.RemoveJobEntry(job);
    }
    delete job;
    return true;
  } else {
    return false;
  }
}

bool JobManager::RemoveJobEntry(job_id_t job_id) {
  JobEntry* job;
  if (GetJobEntry(job_id, job)) {
    assert(job_id == job->job_id());
    job_graph_.RemoveVertex(job_id);
    if (job->parent_job_id() != NIMBUS_KERNEL_JOB_ID ||
        job->job_name() == NIMBUS_MAIN_JOB_NAME) {
      version_manager_.RemoveJobEntry(job);
    }
    delete job;
    return true;
  } else {
    return false;
  }
}

bool JobManager::ResolveJobDataVersions(JobEntry *job) {
  if (!job->IsReadyForCompleteVersioning()) {
    dbg(DBG_ERROR, "ERROR: job %lu is not reaqdy for complete versioing.\n", job->job_id());
    exit(-1);
    return false;
  }

  if (version_manager_.ResolveJobDataVersions(job)) {
    job->set_versioned(true);
    return true;
  } else {
    dbg(DBG_ERROR, "ERROR: could not version job %lu.\n", job->job_id());
    exit(-1);
    return false;
  }
}

size_t JobManager::NumJobsReadyToAssign() {
  return jobs_ready_to_assign_.size();
}

size_t JobManager::GetJobsReadyToAssign(JobEntryList* list, size_t max_num) {
  size_t num = 0;
  list->clear();

  JobEntryMap::iterator iter = jobs_ready_to_assign_.begin();
  for (; (iter != jobs_ready_to_assign_.end()) && (num < max_num);) {
    JobEntry *job = iter->second;
    assert(!job->assigned());

    jobs_pending_to_assign_[iter->first] = job;
    jobs_ready_to_assign_.erase(iter++);

    list->push_back(job);
    ++num;
  }

  return num;
}

size_t JobManager::RemoveObsoleteJobEntries(size_t max_to_remove) {
  size_t num = 0;

  JobEntryList jobs_to_remove;
  {
    boost::unique_lock<boost::mutex> lock(job_graph_mutex_);
    for (size_t i = 0; (i < max_to_remove) && (i < jobs_done_.size()); ++i) {
      jobs_to_remove.push_back(jobs_done_.front());
      jobs_done_.pop_front();
    }
  }

  JobEntryList::iterator iter = jobs_to_remove.begin();
  for (; iter != jobs_to_remove.end(); ++iter) {
    assert((*iter)->done());
    RemoveJobEntry(*iter);
    dbg(DBG_SCHED, "removed job with id %lu from job manager.\n", (*iter)->job_id());
    ++num;
  }

  version_manager_.CleanUp();

  return num;
}

void JobManager::NotifyJobAssignment(JobEntry *job) {
  job->set_assigned(true);
  job->set_assigned_worker_id(job->assigned_worker()->worker_id());

  // AfterMap has internal locking.
  after_map_->AddEntries(job);

  if (job->job_type() != JOB_COMP) {
    return;
  }

  {
    boost::unique_lock<boost::recursive_mutex> job_queue_lock(job_queue_mutex_);
    jobs_pending_to_assign_.erase(job->job_id());
  }

  if (job->sterile()) {
    job_id_t job_id = job->job_id();
    Vertex<JobEntry, job_id_t>* vertex;
    {
      boost::unique_lock<boost::mutex> job_graph_lock(job_graph_mutex_);
      job_graph_.GetVertex(job_id, &vertex);
    }
    typename Edge<JobEntry, job_id_t>::Iter it;
    for (it = vertex->outgoing_edges()->begin(); it != vertex->outgoing_edges()->end(); ++it) {
      JobEntry *j = it->second->end_vertex()->entry();
      boost::unique_lock<boost::recursive_mutex> job_queue_lock(job_queue_mutex_);
      j->remove_assignment_dependency(job_id);
      if (j->IsReadyToAssign()) {
        jobs_ready_to_assign_[j->job_id()] = j;
      }
    }
  }
}

void JobManager::NotifyJobDone(JobEntry *job) {
  job->set_done(true);
  job_id_t job_id = job->job_id();

  // Initialy let the version manager know that the job is done.
  version_manager_.NotifyJobDone(job);

  // Put the job in the list for removal.
  {
    boost::unique_lock<boost::mutex> lock(job_graph_mutex_);
    jobs_done_.push_back(job);
  }

  if (!job->sterile()) {
    Vertex<JobEntry, job_id_t>* vertex;
    job_graph_.GetVertex(job_id, &vertex);
    typename Edge<JobEntry, job_id_t>::Iter it;
    for (it = vertex->outgoing_edges()->begin(); it != vertex->outgoing_edges()->end(); ++it) {
      JobEntry *j = it->second->end_vertex()->entry();
      j->remove_assignment_dependency(job_id);
      if (j->IsReadyToAssign()) {
        jobs_ready_to_assign_[j->job_id()] = j;
      }
    }
  }
}

void JobManager::DefineData(job_id_t job_id, logical_data_id_t ldid) {
  JobEntry* job;
  if (GetJobEntry(job_id, job)) {
    if (job->sterile()) {
      dbg(DBG_ERROR, "ERROR: sterile job cannot define data.\n");
    }
    if (!version_manager_.DefineData(ldid, job_id, job->job_depth())) {
      dbg(DBG_ERROR, "ERROR: could not define data in ldl_map for ldid %lu.\n", ldid);
    }
  } else {
    dbg(DBG_ERROR, "ERROR: parent of define data with job id %lu is not in the graph.\n", job_id);
    exit(-1);
  }
}

size_t JobManager::GetJobsNeedDataVersion(JobEntryList* list,
    VersionedLogicalData vld) {
  return version_manager_.GetJobsNeedDataVersion(list, vld);
}

bool JobManager::AllJobsAreDone() {
  bool all_done = true;

  typename Vertex<JobEntry, job_id_t>::Iter iter = job_graph_.begin();
  for (; iter != job_graph_.end(); ++iter) {
    JobEntry* job = iter->second->entry();
    if (!job->done()) {
      all_done = false;
      break;
    }
  }
  return all_done;
}

void JobManager::UpdateJobBeforeSet(JobEntry* job) {
  UpdateBeforeSet(job->before_set_p());
}

void JobManager::UpdateBeforeSet(IDSet<job_id_t>* before_set) {
  IDSet<job_id_t>::IDSetIter it;
  for (it = before_set->begin(); it != before_set->end();) {
    JobEntry* j;
    job_id_t id = *it;
    if (GetJobEntry(id, j)) {
      if ((j->done()) || (id == NIMBUS_KERNEL_JOB_ID)) {
        before_set->remove(it++);
      } else {
        ++it;
      }
    } else {
      // if the job is not in the table it is already done and removed.
      before_set->remove(it++);
    }
  }
}

bool JobManager::CausingUnwantedSerialization(JobEntry* job,
    const logical_data_id_t& l_id, const PhysicalData& pd) {
  bool result = false;

  if (!job->write_set_p()->contains(l_id)) {
    return result;
  }

  IDSet<job_id_t>::IDSetIter iter;
  for (iter = pd.list_job_read_p()->begin(); iter != pd.list_job_read_p()->end(); iter++) {
    JobEntry *j;
    if (GetJobEntry(*iter, j)) {
      if ((!j->done()) &&
          (j->job_type() == JOB_COMP) &&
          (!job->before_set_p()->contains(*iter))) {
        result = true;
        break;
      }
    }
  }

  return result;
}


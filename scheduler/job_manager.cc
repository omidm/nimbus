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
  if (!AddKernelJobEntry()) {
    dbg(DBG_ERROR, "ERROR: could not add scheduler kernel job in job manager constructor.\n");
    exit(-1);
  }

  ldo_map_p_ = NULL;

  pass_version_in_progress_ = 0;
}

JobManager::~JobManager() {
  // TODO(omidm): do you need to call remove obsolete?
  JobEntry* job;
  if (JobManager::GetJobEntry(NIMBUS_KERNEL_JOB_ID, job)) {
    delete job;
  }
}

bool JobManager::AddComputeJobEntry(
    const std::string& job_name,
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
      return false;
    }
  }

  if (!AddJobEntryIncomingEdges(job)) {
    job_graph_.RemoveVertex(job_id);
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    exit(-1);
    return false;
  }
  ReceiveMetaBeforeSetDepthVersioningDependency(job);
  PassMetaBeforeSetDepthVersioningDependency(job);

  version_manager_.AddJobEntry(job);

  return true;
}

bool JobManager::AddExplicitCopyJobEntry() {
  dbg(DBG_ERROR, "ERROR: explicit copy jobs from application are not supported yet!.\n");
  exit(-1);
  return false;
}

bool JobManager::AddKernelJobEntry() {
  JobEntry* job = new KernelJobEntry();

  if (!job_graph_.AddVertex(NIMBUS_KERNEL_JOB_ID, job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add kernel job in job manager.\n");
    exit(-1);
    return false;
  }
  job->set_done(true);

  return true;
}

bool JobManager::AddMainJobEntry(const job_id_t& job_id) {
  JobEntry* job = new MainJobEntry(job_id);

  if (!job_graph_.AddVertex(job_id, job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    exit(-1);
    return false;
  }

  if (!AddJobEntryIncomingEdges(job)) {
    job_graph_.RemoveVertex(job_id);
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    exit(-1);
    return false;
  }
  ReceiveMetaBeforeSetDepthVersioningDependency(job);
  PassMetaBeforeSetDepthVersioningDependency(job);

  version_manager_.AddJobEntry(job);

  jobs_ready_to_assign_[job_id] = job;

  return true;
}

bool JobManager::AddCreateDataJobEntry(const job_id_t& job_id) {
  JobEntry* job = new CreateDataJobEntry(job_id);

  boost::unique_lock<boost::mutex> lock(job_graph_mutex_);
  if (!job_graph_.AddVertex(job_id, job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    exit(-1);
    return false;
  }

  return true;
}

bool JobManager::AddLocalCopyJobEntry(const job_id_t& job_id) {
  JobEntry* job = new LocalCopyJobEntry(job_id);

  boost::unique_lock<boost::mutex> lock(job_graph_mutex_);
  if (!job_graph_.AddVertex(job_id, job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    exit(-1);
    return false;
  }

  return true;
}

bool JobManager::AddRemoteCopySendJobEntry(const job_id_t& job_id) {
  JobEntry* job = new RemoteCopySendJobEntry(job_id);

  boost::unique_lock<boost::mutex> lock(job_graph_mutex_);
  if (!job_graph_.AddVertex(job_id, job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    exit(-1);
    return false;
  }

  return true;
}

bool JobManager::AddRemoteCopyReceiveJobEntry(const job_id_t& job_id) {
  JobEntry* job = new RemoteCopyReceiveJobEntry(job_id);

  boost::unique_lock<boost::mutex> lock(job_graph_mutex_);
  if (!job_graph_.AddVertex(job_id, job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    exit(-1);
    return false;
  }

  return true;
}

bool JobManager::AddFutureJobEntry(const job_id_t& job_id) {
  JobEntry* job = new FutureJobEntry(job_id);

  boost::unique_lock<boost::mutex> lock(job_graph_mutex_);
  if (!job_graph_.AddVertex(job_id, job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager as future job.\n", job_id);
    exit(-1);
    return false;
  }
  return true;
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
  log_version_.log_ResetTimer();
  log_merge_.log_ResetTimer();
  log_lookup_.log_ResetTimer();
  log_sterile_.log_ResetTimer();
  log_nonsterile_.log_ResetTimer();
  lookup_count_ = 0;
  // while (ResolveDataVersions() > 0) {
  //   continue;
  // }
  if (lookup_count_ > 0) {
    std::cout << "Versioning in get ready jobs: versioning: " <<
      log_version_.timer() << " merge: " << log_merge_.timer() <<
      " lookup: " << log_lookup_.timer() <<
      " lookup_count: " << lookup_count_ <<
      " sterile: " << log_sterile_.timer() <<
      " nonsterile: " << log_nonsterile_.timer() << std::endl;
  }


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

//  size_t num = 0;
//  list->clear();
//  typename Vertex<JobEntry, job_id_t>::Iter iter = job_graph_.begin();
//  for (; (iter != job_graph_.end()) && (num < max_num); ++iter) {
//    Vertex<JobEntry, job_id_t>* vertex = iter->second;
//    JobEntry* job = vertex->entry();
//    if (job->versioned() && !job->assigned()) {
//      // Job is already versioned so it has the information from parent and
//      // beforeset already, they may not be in the graph at this point though. -omidm
//      JobEntry* j;
//
//      if (GetJobEntry(job->parent_job_id(), j)) {
//        if (!(j->done())) {
//          continue;
//        }
//      }
//
//      bool before_set_done_or_sterile = true;
//
//      Edge<JobEntry, job_id_t>::Iter it;
//      typename Edge<JobEntry, job_id_t>::Map* incoming_edges = vertex->incoming_edges();
//      for (it = incoming_edges->begin(); it != incoming_edges->end(); ++it) {
//        j = it->second->start_vertex()->entry();
//        /*
//         * Due to the current graph traversal we should only assign the job if
//         * before set is done otherwise by leveraging the sterile flag we could
//         * end up flooding worker with lots of jobs that depend on a job that
//         * has not been assigned yet and so it causes a lot of latency. the
//         * graph traversal right now is based on iterator of map (so job ids)
//         * we need to traverse based on graph shape. In addition to efficiency
//         * issues it could cause problem since the before set may not be
//         * assigned yet and so we may not find the data version for the job in
//         * the system. -omidm
//         */
//         // if (!(j->done())) {
//        /*
//         * For now and sake of current water multiple since the number of jobs
//         * are not too large and we have huge number of data partitions the
//         * application would benefit if scheduler can use the sterile flag and
//         * scheduler jobs in advance. Note that since still the job ids may not
//         * be in order we have to make sure that the jobs in before set are
//         * already assigned, otherwise we may not find the data version we want
//         * for the job in the system. -omidm
//         */
//        if (!(j->done()) && !(j->sterile() && j->assigned())) {
//        /*
//         * The ultimate goal is to turn it to this after we have built the
//         * graph traversal, since the job in before set is assigned for sure
//         * before the job in after set if we travers properly. -omidm
//         */
//        // if (!(j->done()) && !(j->sterile())) {
//          before_set_done_or_sterile = false;
//          break;
//        }
//      }
//
//      if (before_set_done_or_sterile) {
//        // job->set_assigned(true); No, we are not sure yet thet it will be assignd!
//        list->push_back(job);
//        ++num;
//      }
//    }
//  }
//  return num;
}

size_t JobManager::RemoveObsoleteJobEntries() {
  log_version_.log_ResetTimer();
  log_merge_.log_ResetTimer();
  log_lookup_.log_ResetTimer();
  log_sterile_.log_ResetTimer();
  log_nonsterile_.log_ResetTimer();
  lookup_count_ = 0;
  // while (ResolveDataVersions() > 0) {
  //   continue;
  // }
  if (lookup_count_ > 0) {
    std::cout << "Versioning in get ready jobs: versioning: " <<
      log_version_.timer() << " merge: " << log_merge_.timer() <<
      " lookup: " << log_lookup_.timer() <<
      " lookup_count: " << lookup_count_ <<
      " sterile: " << log_sterile_.timer() <<
      " nonsterile: " << log_nonsterile_.timer() << std::endl;
  }

  size_t num = 0;

  if (NumJobsReadyToAssign() > 0) {
    return num;
  }

  JobEntryMap::iterator iter;
  for (iter = jobs_done_.begin(); iter != jobs_done_.end();) {
    assert(iter->second->done());
    RemoveJobEntry(iter->second);
    dbg(DBG_SCHED, "removed job with id %lu from job manager.\n", iter->first);
    ++num;
    jobs_done_.erase(iter++);
  }

  version_manager_.CleanUp();

  return num;
}

void JobManager::NotifyJobAssignment(JobEntry *job) {
  job->set_assigned(true);
  job->set_assigned_worker_id(job->assigned_worker()->worker_id());

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
  jobs_done_[job_id] = job;
  jobs_need_version_.erase(job_id);

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
//    if (!ldl_map_.DefineData(ldid, job_id, job->job_depth(), job->sterile())) {
//      dbg(DBG_ERROR, "ERROR: could not define data in ldl_map for ldid %lu.\n", ldid);
//    }
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


size_t JobManager::ResolveDataVersions() {
  log_version_.log_ResumeTimer();

  size_t num = 0;

  size_t passed_num = 0;

  std::map<job_id_t, JobEntryList> new_pass_version;
  std::map<job_id_t, JobEntryList>::iterator iter;
  for (iter = pass_version_.begin(); iter != pass_version_.end(); ++iter) {
    JobEntry *job;
    job_id_t job_id = iter->first;
    if (GetJobEntry(job_id, job)) {
      JobEntry *parent_job;
      GetJobEntry(job->parent_job_id(), parent_job);
      if (!parent_job->done()) {
        continue;
      }
      ++passed_num;

      PassDataVersionToJob(job, iter->second);
      jobs_need_version_[job_id] = job;
      if (JobVersionIsComplete(job)) {
        /*
         * Logical data lineage approach with meta before set.
         */
        if (job->sterile()) {
          log_sterile_.log_ResumeTimer();
          boost::shared_ptr<VersionMap> vmap = boost::shared_ptr<VersionMap>(new VersionMap());
          IDSet<logical_data_id_t>::ConstIter itw;
          for (itw = job->write_set_p()->begin(); itw != job->write_set_p()->end(); ++itw) {
            data_version_t version;
            if (job->vmap_read()->query_entry(*itw, &version)) {
              ldl_map_.AppendLdlEntry(
                  *itw, job->job_id(), version + 1, job->job_depth(), job->sterile());
              vmap->set_entry(*itw, version + 1);
            } else {
              log_lookup_.log_ResumeTimer();
              lookup_count_++;
              bool found = LookUpVersion(job, *itw, &version);
              log_lookup_.log_StopTimer();
              if (found) {
                ldl_map_.AppendLdlEntry(
                    *itw, job->job_id(), version + 1, job->job_depth(), job->sterile());
                vmap->set_entry(*itw, version + 1);
              } else {
                dbg(DBG_ERROR, "ERROR: could not resolve data id: %lu.\n", *itw); // NOLINT
                exit(-1);
              }
            }
          }
          job->set_vmap_write(vmap);

          log_sterile_.log_StopTimer();
        } else {
          log_nonsterile_.log_ResumeTimer();
          boost::shared_ptr<VersionMap> vmap = boost::shared_ptr<VersionMap>(new VersionMap());
          vmap->set_content(job->vmap_read()->content());
          IDSet<logical_data_id_t>::ConstIter itw;
          for (itw = job->write_set_p()->begin(); itw != job->write_set_p()->end(); ++itw) {
            data_version_t version;
            if (job->vmap_read()->query_entry(*itw, &version)) {
              ldl_map_.AppendLdlEntry(
                  *itw, job->job_id(), version + 1, job->job_depth(), job->sterile());
              vmap->set_entry(*itw, version + 1);
            } else {
              log_lookup_.log_ResumeTimer();
              lookup_count_++;
              bool found = LookUpVersion(job, *itw, &version);
              log_lookup_.log_StopTimer();
              if (found) {
                ldl_map_.AppendLdlEntry(
                    *itw, job->job_id(), version + 1, job->job_depth(), job->sterile());
                vmap->set_entry(*itw, version + 1);
              } else {
                dbg(DBG_ERROR, "ERROR: could not resolve data id: %lu.\n", *itw); // NOLINT
                exit(-1);
              }
            }
          }
          job->set_vmap_write(vmap);

          // Clear meta before set and update the ldl.
          log_merge_.log_ResumeTimer();

          job->meta_before_set()->Clear();

          LdoMap::const_iterator it;
          for (it = ldo_map_p_->begin(); it != ldo_map_p_->end(); ++it) {
            data_version_t v_out;
            vmap->query_entry(it->first, &v_out);
            ldl_map_.InsertParentLdlEntry(
                it->first, job->job_id(), v_out, job->job_depth(), job->sterile());
          }
          log_merge_.log_StopTimer();

          log_nonsterile_.log_StopTimer();
        }

        job->set_versioned(true);
        Vertex<JobEntry, job_id_t>* vertex;
        job_graph_.GetVertex(job_id, &vertex);
        typename Edge<JobEntry, job_id_t>::Iter it;
        for (it = vertex->outgoing_edges()->begin(); it != vertex->outgoing_edges()->end(); ++it) {
          new_pass_version[it->first].push_back(job);
        }
        ++num;
      }
    } else {
      dbg(DBG_ERROR, "ERROR: Job (id: %lu) is not in the graph to receive the versions.\n", iter->first); // NOLINT
      exit(-1);
    }
  }


  assert((passed_num == 0) || (passed_num == pass_version_.size()));
  if (passed_num == pass_version_.size()) {
    pass_version_.clear();
    pass_version_ = new_pass_version;
  }


  log_version_.log_StopTimer();
  return num;
}

void JobManager::PassDataVersionToJob(
    JobEntry *job, const JobEntryList& source_jobs) {
  assert(!job->future() && !job->versioned());
  assert(source_jobs.size() > 0);

  if (job->partial_versioned()) {
    job->meta_before_set()->InvalidateNegativeQueryCache();
  } else {
    job->set_vmap_read(
        boost::shared_ptr<VersionMap>(new VersionMap()));
    job->set_meta_before_set(
        boost::shared_ptr<MetaBeforeSet>(new MetaBeforeSet()));
  }

  job_depth_t depth = (*source_jobs.begin())->job_depth();
  JobEntryList::const_iterator iter;
  for (iter = source_jobs.begin(); iter != source_jobs.end(); ++iter) {
    JobEntry* j = (*iter);
    assert(j->versioned());
    job->meta_before_set()->table_p()->insert(
        std::pair<job_id_t, boost::shared_ptr<MetaBeforeSet> > (j->job_id(), j->meta_before_set())); // NOLINT
    if (depth < j->job_depth()) {
      depth = j->job_depth();
    }
    job->add_job_passed_versions(j->job_id());
  }
  job->set_job_depth(depth + 1);

  if (job->sterile()) {
    log_sterile_.log_ResumeTimer();
    IDSet<logical_data_id_t>::ConstIter it;
    for (it = job->read_set_p()->begin(); it != job->read_set_p()->end(); ++it) {
      data_version_t version;
      log_lookup_.log_ResumeTimer();
      lookup_count_++;
      bool found = LookUpVersion(job, *it, &version);
      log_lookup_.log_StopTimer();
      if (found) {
        job->vmap_read()->set_entry(*it, version);
      }
    }
    log_sterile_.log_StopTimer();
  } else {
    log_nonsterile_.log_ResumeTimer();
    boost::shared_ptr<VersionMap> vmap = boost::shared_ptr<VersionMap>(new VersionMap());
    LdoMap::const_iterator it;
    for (it = ldo_map_p_->begin(); it != ldo_map_p_->end(); ++it) {
      data_version_t version;
      log_lookup_.log_ResumeTimer();
      lookup_count_++;
      bool found = LookUpVersion(job, it->first, &version);
      log_lookup_.log_StopTimer();
      if (found) {
        job->vmap_read()->set_entry(it->first, version);
      }
    }
    log_nonsterile_.log_StopTimer();
  }

  job->set_partial_versioned(true);
}

bool JobManager::LookUpVersion(JobEntry *job,
    logical_data_id_t ldid, data_version_t *version) {
  return ldl_map_.LookUpVersion(
      ldid, job->meta_before_set(), version);
}

bool JobManager::JobVersionIsComplete(JobEntry *job) {
  IDSet<job_id_t> need = job->need_set();
  return (need.size() == 0);
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

void JobManager::set_ldo_map_p(const LdoMap* ldo_map_p) {
  ldo_map_p_ = ldo_map_p;
  version_manager_.set_ldo_map_p(ldo_map_p);
}

Graph<JobEntry, job_id_t>* JobManager::job_graph_p() {
  return &job_graph_;
}





void JobManager::WaitToPassAllVersions() {
  // TODO(omidm): Implement!
  // boost::lock_guard<boost::mutex> l(mtx);
  boost::unique_lock<boost::mutex> lock(pass_version_mutex_);
  while (pass_version_in_progress_ > 0 ||
         pass_version_.size() > 0) {
    pass_version_draw_cond_.wait(lock);
  }
}






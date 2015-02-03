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
#include "scheduler/binding_template.h"

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
  checkpoint_creation_rate_ = DEFAULT_CHECKPOINT_CREATION_RATE;
  non_sterile_counter_ = 0;
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

bool JobManager::AddJobEntryToJobGraph(job_id_t job_id, JobEntry *job) {
  boost::unique_lock<boost::mutex> lock(job_graph_mutex_);
  return job_graph_.AddVertex(job_id, job);
}

bool JobManager::RemoveJobEntryFromJobGraph(job_id_t job_id) {
  boost::unique_lock<boost::mutex> lock(job_graph_mutex_);
  return job_graph_.RemoveVertex(job_id);
}

bool JobManager::ClearJobGraph() {
  boost::unique_lock<boost::mutex> lock(job_graph_mutex_);
  std::list<job_id_t> list_to_remove;
  typename Vertex<JobEntry, job_id_t>::Iter iter = job_graph_.begin();
  for (; iter != job_graph_.end(); ++iter) {
    JobEntry* job = iter->second->entry();
    if (job->job_id() != NIMBUS_KERNEL_JOB_ID) {
      list_to_remove.push_back(job->job_id());
      delete job;
    }
  }

  std::list<job_id_t>::iterator it = list_to_remove.begin();
  for (; it != list_to_remove.end(); ++it) {
    job_graph_.RemoveVertex(*it);
  }

  return true;
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
                                         const GeometricRegion& region,
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
                        region,
                        params);

  if (!AddJobEntryToJobGraph(job_id, job)) {
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
      job->set_region(region);
      job->set_future(false);
      dbg(DBG_SCHED, "Filled the information for future job (id: %lu).\n", job_id);
    } else {
      dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
      exit(-1);
      return NULL;
    }
  }

  if (!AddJobEntryIncomingEdges(job)) {
    RemoveJobEntryFromJobGraph(job_id);
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    exit(-1);
    return NULL;
  }

  ReceiveMetaBeforeSetDepthVersioningDependency(job);
  PassMetaBeforeSetDepthVersioningDependency(job);

  version_manager_.AddJobEntry(job);

  if (!job->sterile()) {
    non_sterile_jobs_[job_id] = job;
  }

  return job;
}

bool JobManager::AddComplexJobEntry(ComplexJobEntry* complex_job) {
  job_id_t job_id = complex_job->job_id();

  if (!AddJobEntryToJobGraph(job_id, complex_job)) {
    dbg(DBG_ERROR, "ERROR: could not add complex job (id: %lu) in job manager.\n", job_id);
    delete complex_job;
    return false;
  }

  if (!AddJobEntryIncomingEdges(complex_job)) {
    RemoveJobEntryFromJobGraph(job_id);
    delete complex_job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    return false;
  }

  ReceiveMetaBeforeSetDepthVersioningDependency(complex_job);
  PassMetaBeforeSetDepthVersioningDependency(complex_job);

  version_manager_.AddComplexJobEntry(complex_job);

  std::list<job_id_t> list;
  complex_job->GetParentJobIds(&list);
  // For now complex job could have only one parent job, otherwise versioning goes wrong!
  assert(list.size() == 1);
  complex_jobs_[job_id] = complex_job;

  // TODO(omidm): what does it mean for checkpointing logic?!!
  if (!complex_job->sterile()) {
    non_sterile_jobs_[job_id] = complex_job;
  }

  return true;
}

JobEntry* JobManager::AddExplicitCopyJobEntry() {
  dbg(DBG_ERROR, "ERROR: explicit copy jobs from application are not supported yet!.\n");
  exit(-1);
  return NULL;
}

JobEntry* JobManager::AddKernelJobEntry() {
  JobEntry* job = new KernelJobEntry();

  if (!AddJobEntryToJobGraph(NIMBUS_KERNEL_JOB_ID, job)) {
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

  if (!AddJobEntryToJobGraph(job_id, job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    exit(-1);
    return NULL;
  }

  if (!AddJobEntryIncomingEdges(job)) {
    RemoveJobEntryFromJobGraph(job_id);
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

  if (!AddJobEntryToJobGraph(job_id, job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    exit(-1);
    return NULL;
  }

  return job;
}

JobEntry* JobManager::AddLocalCopyJobEntry(const job_id_t& job_id) {
  JobEntry* job = new LocalCopyJobEntry(job_id);

  if (!AddJobEntryToJobGraph(job_id, job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    exit(-1);
    return NULL;
  }

  return job;
}

JobEntry* JobManager::AddRemoteCopySendJobEntry(const job_id_t& job_id) {
  JobEntry* job = new RemoteCopySendJobEntry(job_id);

  if (!AddJobEntryToJobGraph(job_id, job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    exit(-1);
    return NULL;
  }

  return job;
}

JobEntry* JobManager::AddRemoteCopyReceiveJobEntry(const job_id_t& job_id) {
  JobEntry* job = new RemoteCopyReceiveJobEntry(job_id);

  if (!AddJobEntryToJobGraph(job_id, job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager.\n", job_id);
    exit(-1);
    return NULL;
  }

  return job;
}

JobEntry* JobManager::AddFutureJobEntry(const job_id_t& job_id) {
  JobEntry* job = new FutureJobEntry(job_id);

  if (!AddJobEntryToJobGraph(job_id, job)) {
    delete job;
    dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) in job manager as future job.\n", job_id);
    exit(-1);
    return NULL;
  }

  return job;
}

Edge<JobEntry, job_id_t>* JobManager::AddEdgeToJobGraph(job_id_t from, job_id_t to) {
  boost::unique_lock<boost::mutex> lock(job_graph_mutex_);
  Edge<JobEntry, job_id_t> *edge;
  if (job_graph_.AddEdge(from, to, &edge)) {
    return edge;
  } else {
    return NULL;
  }
}

bool JobManager::GetComplexJobContainer(const job_id_t& job_id,
                                        ComplexJobEntry*& complex_job) {
  ComplexJobEntryMap::iterator iter = complex_jobs_.begin();
  for (; iter != complex_jobs_.end(); ++iter) {
    ShadowJobEntry* shadow_job;
    if (iter->second->GetShadowJobEntry(job_id, shadow_job)) {
      assert(job_id == shadow_job->job_id());
      complex_job = iter->second;
      return true;
    }
  }

  complex_job = NULL;
  return false;
}

bool JobManager::AddJobEntryIncomingEdges(JobEntry *job) {
  JobEntry *j;
  Edge<JobEntry, job_id_t> *edge;

  ComplexJobEntry* complex_job;
  edge = AddEdgeToJobGraph(job->parent_job_id(), job->job_id());
  if (edge != NULL) {
    j = edge->start_vertex()->entry();
    assert(j->versioned());
  } else if (GetComplexJobContainer(job->parent_job_id(), complex_job)) {
    edge = AddEdgeToJobGraph(complex_job->job_id(), job->job_id());
    assert(edge);
    j = edge->start_vertex()->entry();
    // TODO(omidm): uncomment this check!
    // assert(j->versioned());
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

    edge = AddEdgeToJobGraph(*it, job->job_id());
    if (edge != NULL) {
      j = edge->start_vertex()->entry();
    } else if (GetComplexJobContainer(*it, complex_job)) {
      // cannot have before set in a complex job - omidm
      assert(false);
    } else {
      dbg(DBG_SCHED, "Adding possible future job (id: %lu) in job manager.\n", *it);
      AddFutureJobEntry(*it);
      edge = AddEdgeToJobGraph(*it, job->job_id());
      if (edge == NULL) {
        return false;
      }
    }
  }

  return true;
}

bool JobManager::GetJobGraphVertex(job_id_t job_id, Vertex<JobEntry, job_id_t> **vertex) {
  boost::unique_lock<boost::mutex> lock(job_graph_mutex_);
  return job_graph_.GetVertex(job_id, vertex);
}

void JobManager::ReceiveMetaBeforeSetDepthVersioningDependency(JobEntry* job) {
  job_depth_t depth = NIMBUS_INIT_JOB_DEPTH;
  job_id_t job_id = job->job_id();
  Vertex<JobEntry, job_id_t>* vertex;
  bool got_vertex = GetJobGraphVertex(job_id, &vertex);
  assert(got_vertex);
  typename Edge<JobEntry, job_id_t>::Iter it;
  for (it = vertex->incoming_edges()->begin(); it != vertex->incoming_edges()->end(); ++it) {
    JobEntry *j = it->second->start_vertex()->entry();
    if (j->IsReadyForCompleteVersioning()) {
      job->meta_before_set()->table_p()->insert(
          std::pair<job_id_t, boost::shared_ptr<MetaBeforeSet> > (j->job_id(), j->meta_before_set())); // NOLINT

      if (depth < j->job_depth()) {
        depth = j->job_depth();
      }

      if (j->job_type() == JOB_CMPX) {
        ComplexJobEntry* complex_job = reinterpret_cast<ComplexJobEntry*>(j);
        std::list<job_id_t> list;
        complex_job->GetParentJobIds(&list);
        // For now complex job could have only one parent job, otherwise versioning goes wrong!
        assert(list.size() == 1);
        job->remove_versioning_dependency(*(list.begin()));
      } else {
        job->remove_versioning_dependency(j->job_id());
      }
    }
  }
  job->set_job_depth(depth + 1);
}


void JobManager::PassMetaBeforeSetDepthVersioningDependency(JobEntry* job) {
  if (job->IsReadyForCompleteVersioning()) {
    job_id_t job_id = job->job_id();
    Vertex<JobEntry, job_id_t>* vertex;
    bool got_vertex = GetJobGraphVertex(job_id, &vertex);
    assert(got_vertex);
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
    ComplexJobEntryMap::iterator iter = complex_jobs_.begin();
    for (; iter != complex_jobs_.end(); ++iter) {
      ShadowJobEntry* shadow_job;
      if (iter->second->GetShadowJobEntry(job_id, shadow_job)) {
        job = shadow_job;
        return true;
      }
    }

    job = NULL;
    return false;
  }
}

bool JobManager::RemoveJobEntry(JobEntry* job) {
  if (RemoveJobEntryFromJobGraph(job->job_id())) {
    if (job->parent_job_id() != NIMBUS_KERNEL_JOB_ID ||
        // job->job_name() == NIMBUS_MAIN_JOB_NAME) {
        !job->sterile()) {
      version_manager_.RemoveJobEntry(job);
    }
    if (job->job_type() == JOB_CMPX) {
      ComplexJobEntry* xj = reinterpret_cast<ComplexJobEntry*>(job);
      assert(xj->AllJobsRemoved());
    }
    delete job;
    return true;
  } else if (job->job_type() == JOB_SHDW) {
    ShadowJobEntry* sj = reinterpret_cast<ShadowJobEntry*>(job);
    version_manager_.RemoveJobEntry(job);
    sj->complex_job()->MarkJobRemoved(job->job_id());
    delete job;
    return true;
  } else {
    return false;
  }
}

bool JobManager::ResolveJobDataVersions(JobEntry *job) {
  if (!job->IsReadyForCompleteVersioning()) {
    dbg(DBG_ERROR, "ERROR: job %lu is not reaqdy for complete versioing.\n", job->job_id());
    assert(false);
    return false;
  }

  if (version_manager_.ResolveJobDataVersions(job)) {
    job->set_versioned(true);
    return true;
  } else {
    dbg(DBG_ERROR, "ERROR: could not version job %lu.\n", job->job_id());
    assert(false);
    return false;
  }
}


bool JobManager::ResolveEntireContextForJob(JobEntry *job) {
  if (!job->IsReadyForCompleteVersioning()) {
    dbg(DBG_ERROR, "ERROR: job %lu is not reaqdy for complete versioing.\n", job->job_id());
    exit(-1);
    return false;
  }

  if (version_manager_.ResolveEntireContextForJob(job)) {
    job->set_versioned_entire_context(true);
    return true;
  } else {
    dbg(DBG_ERROR, "ERROR: could not version job %lu for entire context.\n", job->job_id());
    exit(-1);
    return false;
  }
}

bool JobManager::MemoizeVersionsForTemplate(JobEntry *job) {
  if (!job->IsReadyForCompleteVersioning()) {
    dbg(DBG_ERROR, "ERROR: job %lu is not reaqdy for complete versioing.\n", job->job_id());
    exit(-1);
    return false;
  }

  if (version_manager_.MemoizeVersionsForTemplate(job)) {
    return true;
  } else {
    dbg(DBG_ERROR, "ERROR: could not memoize versions for job %lu.\n", job->job_id());
    exit(-1);
    return false;
  }
}





bool JobManager::CompleteJobForCheckpoint(checkpoint_id_t checkpoint_id,
                                          const JobEntry *job) {
  return checkpoint_manager_.CompleteJobForCheckpoint(checkpoint_id, job);
}

bool JobManager::AddSaveDataJobToCheckpoint(checkpoint_id_t checkpoint_id,
                                            job_id_t job_id,
                                            logical_data_id_t ldid,
                                            data_version_t version,
                                            worker_id_t worker_id) {
  return checkpoint_manager_.AddSaveDataJobToCheckpoint(checkpoint_id,
                                                        job_id,
                                                        ldid, version,
                                                        worker_id);
}

bool JobManager::GetHandleToLoadData(checkpoint_id_t checkpoint_id,
                                     logical_data_id_t ldid,
                                     data_version_t version,
                                     WorkerHandleList *handles) {
    return checkpoint_manager_.GetHandleToLoadData(checkpoint_id,
                                                   ldid,
                                                   version,
                                                   handles);
}


bool JobManager::NotifySaveDataJobDoneForCheckpoint(checkpoint_id_t checkpoint_id,
                                                    job_id_t job_id,
                                                    std::string handle) {
  return checkpoint_manager_.NotifySaveDataJobDoneForCheckpoint(checkpoint_id,
                                                                job_id,
                                                                handle);
}



bool JobManager::RewindFromLastCheckpoint(checkpoint_id_t *checkpoint_id) {
  if (!checkpoint_manager_.GetCheckpointToRewind(checkpoint_id)) {
    dbg(DBG_ERROR, "ERROR: could not get any checkpoint to rewind!\n");
    exit(-1);
    return false;
  }

  JobEntryList list;
  checkpoint_manager_.GetJobListFromCheckpoint(*checkpoint_id, &list);

  after_map_->Clear();

  version_manager_.Reinitialize(&list);

  ClearJobGraph();

  {
    boost::unique_lock<boost::mutex> lock(jobs_done_mutex_);
    jobs_done_.clear();
    non_sterile_jobs_.clear();
  }

  {
    boost::unique_lock<boost::recursive_mutex> lock(job_queue_mutex_);
    jobs_ready_to_assign_.clear();
    jobs_pending_to_assign_.clear();
    JobEntryList::iterator iter = list.begin();
    for (; iter != list.end(); ++iter) {
      IDSet<job_id_t> empty_set;


      JobEntry* job =
        new ComputeJobEntry((*iter)->job_name(),
                            (*iter)->job_id(),
                            (*iter)->read_set(),
                            (*iter)->write_set(),
                            empty_set,
                            empty_set,
                            NIMBUS_KERNEL_JOB_ID,
                            (*iter)->future_job_id(),
                            (*iter)->sterile(),
                            (*iter)->region(),
                            (*iter)->params());

      if (!AddJobEntryToJobGraph((*iter)->job_id(), job)) {
          dbg(DBG_ERROR, "ERROR: could not add job (id: %lu) to job graph.\n", (*iter)->job_id());
          exit(-1);
          return false;
      }

      assert((*iter)->versioned_entire_context());
      job->set_vmap_read((*iter)->vmap_read());
      job->set_vmap_write((*iter)->vmap_write());
      job->set_job_depth((*iter)->job_depth());
      job->MarkJobAsCompletelyResolved();

      jobs_ready_to_assign_[(*iter)->job_id()] = job;

      version_manager_.AddJobEntry(job);

      assert(!job->sterile());
      non_sterile_jobs_[(*iter)->job_id()] = job;
    }
  }

  return true;
}

size_t JobManager::NumJobsReadyToAssign() {
  return jobs_ready_to_assign_.size();
}

size_t JobManager::GetJobsReadyToAssign(JobEntryList* list, size_t max_num) {
  size_t num = 0;
  list->clear();

  // JobEntryList temp_list;
  JobEntryMap::iterator iter = jobs_ready_to_assign_.begin();
  for (; (iter != jobs_ready_to_assign_.end()) && (num < max_num);) {
    JobEntry* job = iter->second;
    assert(!job->assigned());

    if (job->job_type() == JOB_CMPX) {
      ComplexJobEntry* complex_job = reinterpret_cast<ComplexJobEntry*>(job);
      TemplateEntry *te = complex_job->template_entry();
      BindingTemplate *bt;
      bool create_bt = true;
      if (te->QueryBindingRecord(STATIC_BINDING_RECORD, bt)) {
        create_bt = false;
        if (bt->finalized()) {
          complex_job->set_binding_template(bt);
          jobs_pending_to_assign_[complex_job->job_id()] = complex_job;
          jobs_ready_to_assign_.erase(iter++);
          num += bt->compute_job_num();
          continue;
        }
      }

      if (create_bt) {
        bt = new BindingTemplate();
        if (!te->AddBindingRecord(STATIC_BINDING_RECORD, bt)) {
          assert(false);
        }
      }

      JobEntryList l;
      size_t c_num = complex_job->GetJobsForAssignment(&l, max_num - num);
      assert(c_num > 0);
      num += c_num;
      JobEntryList::iterator it = l.begin();
      for (; it != l.end(); ++it) {
        (*it)->set_memoize_binding(true);
        (*it)->set_binding_template(bt);
        list->push_back(*it);
        jobs_pending_to_assign_[(*it)->job_id()] = *it;
      }
      if (complex_job->DrainedAllJobsForAssignment()) {
        jobs_ready_to_assign_.erase(iter++);
      } else {
        ++iter;
      }

    } else {
      list->push_back(job);
      // temp_list.push_back(job);
      jobs_pending_to_assign_[iter->first] = job;
      jobs_ready_to_assign_.erase(iter++);
      ++num;
    }
  }

//  if (temp_list.size() > 0) {
//    JobEntryList::reverse_iterator it = temp_list.rbegin();
//    if ((*it)->job_name() == "--reincorporate_particles" ||
//        (*it)->job_name() == "--delete_particles" ||
//        (*it)->job_name() == "--adjust_phi") {
//      for (; it != temp_list.rend(); ++it) {
//        list->push_back(*it);
//      }
//    } else {
//      JobEntryList::iterator i = temp_list.begin();
//      for (; i != temp_list.end(); ++i) {
//        list->push_back(*i);
//      }
//    }
//  }

  return num;
}

size_t JobManager::RemoveObsoleteJobEntries(size_t max_to_remove) {
  size_t num = 0;

  JobEntryList jobs_to_remove;
  {
    boost::unique_lock<boost::mutex> lock(jobs_done_mutex_);
    for (size_t i = 0; (i < max_to_remove) && (i < jobs_done_.size()); ++i) {
      jobs_to_remove.push_back(jobs_done_.front());
      jobs_done_.pop_front();
    }
  }

  JobEntryList::iterator iter = jobs_to_remove.begin();
  for (; iter != jobs_to_remove.end(); ++iter) {
    assert((*iter)->done());
    dbg(DBG_SCHED, "removing job with id %lu from job manager.\n", (*iter)->job_id());
    RemoveJobEntry(*iter);
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

  if ((job->job_type() != JOB_COMP) &&
      (job->job_type() != JOB_CMPX)) {
    return;
  }

  {
    boost::unique_lock<boost::recursive_mutex> job_queue_lock(job_queue_mutex_);
    jobs_pending_to_assign_.erase(job->job_id());
  }

  if (job->sterile()) {
    job_id_t job_id = job->job_id();
    Vertex<JobEntry, job_id_t>* vertex;
    bool got_vertex = GetJobGraphVertex(job_id, &vertex);
    assert(got_vertex);
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

  if (!job->sterile()) {
    ComplexJobEntry* complex_job;
    Vertex<JobEntry, job_id_t>* vertex;
    if (!GetJobGraphVertex(job_id, &vertex)) {
      if (!GetComplexJobContainer(job_id, complex_job)) {
        assert(false);
      }
      bool got_vertex = GetJobGraphVertex(complex_job->job_id(), &vertex);
      assert(got_vertex);
    }
    typename Edge<JobEntry, job_id_t>::Iter it;
    for (it = vertex->outgoing_edges()->begin(); it != vertex->outgoing_edges()->end(); ++it) {
      JobEntry *j = it->second->end_vertex()->entry();
      j->remove_assignment_dependency(job_id);
      if (j->IsReadyToAssign()) {
        jobs_ready_to_assign_[j->job_id()] = j;
      }
    }
  }


  if (NIMBUS_FAULT_TOLERANCE_ACTIVE) {
    if (!job->sterile()) {
      non_sterile_jobs_.erase(job_id);
      if ((++non_sterile_counter_) == checkpoint_creation_rate_) {
        // Initiate checkpoint creation.
        non_sterile_counter_ = 0;
        checkpoint_id_t checkpoint_id;
        checkpoint_manager_.CreateNewCheckpoint(&checkpoint_id);
        JobEntryMap::iterator it = non_sterile_jobs_.begin();
        for (; it != non_sterile_jobs_.end(); ++it) {
          // TODO(omidm): these checks may not necessarily pass, need general way
          // to enforce these constraints. leave it like this for now, since the
          // job graphs are not complicated.
          assert(!it->second->assigned());
          it->second->set_checkpoint_id(checkpoint_id);
          checkpoint_manager_.AddJobToCheckpoint(checkpoint_id, it->second);
        }
      }
    }
  }

  // After traversing the out going endges then put in the dueue for removal.
  {
    boost::unique_lock<boost::mutex> lock(jobs_done_mutex_);
    jobs_done_.push_back(job);

    // TODO(omidm): right now because the jobs_done is treated as a FIFO, the
    // coplex jobs will be removed after all shadow jobs are removed such that
    // the complex job pointer in shadow jobs is always valid. But you are
    // gonna need a cleaner logic for that (using AllJobsRemoved method). - omidm
    if (job->job_type() == JOB_SHDW) {
      ShadowJobEntry* sj = reinterpret_cast<ShadowJobEntry*>(job);
      ComplexJobEntry* xj = sj->complex_job();
      xj->MarkJobDone(job_id);
      if (xj->AllJobsDone()) {
        xj->set_done(true);
        complex_jobs_.erase(xj->job_id());
        // Let the version manager clean the complex job from its track.
        version_manager_.NotifyJobDone(xj);
        jobs_done_.push_back(xj);
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

  boost::unique_lock<boost::mutex> lock(job_graph_mutex_);
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
    } else if (IDMaker::SchedulerProducedJobID(id)) {
      ++it;
    } else {
      // if the job is not in the table or a copy/create job it is already done and removed.
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
          (j->job_type() == JOB_COMP || j->job_type() == JOB_SHDW) &&
          // (!job->before_set_p()->contains(*iter))) {
          // the job may not be in the immediate before set but still there is
          // indirect dependency. so we need to look through meta before set.
          (!job->LookUpMetaBeforeSet(j))) {
        result = true;
        break;
      }
    }
  }

  return result;
}

